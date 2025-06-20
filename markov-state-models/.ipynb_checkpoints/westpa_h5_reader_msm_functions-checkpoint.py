import numpy as np
from collections import Counter

import westpa
from westpa.analysis import Run

import matplotlib.pyplot as plt
from decimal import Decimal

#Jonathan Borowsky
#4/22/25
#Grabe lab

#---------------------------------------------------------------------------------------------------------
#load paired progress coordinates of every walker and its parent
#---------------------------------------------------------------------------------------------------------

#see https://westpa.readthedocs.io/en/latest/documentation/analysis/index.html
#and https://westpa.readthedocs.io/en/latest/users_guide/hdf5.html

def h5_2_pcs(h5path, miniter, maxiter):
    
    run = Run.open(h5path)

    #set maximum iteration automatically
    if maxiter == -1:
        maxiter = run.num_iterations

    #collect matching arrays of the progress coordinate of each walker and its parent in the previous round
    #i.e. parent_pcs[i] contains the pc of the parent of the walker whose pc is recorded at pcs[i]
    pcs = []
    parent_pcs = []
    weights = []     #for parent walkers

    #initial progress coordinate (only collected if miniter <= 1)
    pc_init = []

    #loop over all WESTPA rounds
    for iteration in run:

        #check if current iteration is in the correct range
        if iteration.number < miniter:
            if iteration.number == 1:
                #initial pc
                pc_init = [[wa.parent.pcoord] for wa in iteration.walkers][0][-1]            
            continue
        elif iteration.number >= maxiter:
            break

        #collect PCs for current round
        pcs.append(iteration.pcoords)

        #collect PCs for parents of walkers in current round
        #note that these are collected as lists of 2d arrays rather than 3d arrays
        #    they are stacked rather than concatenated at the end to achieve the same output shape
        if iteration.number == 1:
            #for the first round there is no parent walker in the regular sense 
            #     but the pcoord is recorded under a slightly different variable name
            #     weight is similarly recorded differently
            parent_pcs.append(np.stack([[wa.parent.pcoord, wa.parent.pcoord] for wa in iteration.walkers]))
            weights.append([wa.parent.basis_state.probability for wa in iteration.walkers])
            
            #initial pc
            pc_init = parent_pcs[-1][-1]
            
        else:
            #not everything appears in .__dict__
            parent_pcs.append(np.stack([wa.parent.pcoords for wa in iteration.walkers]))
            weights.append([wa.parent.weight for wa in iteration.walkers])


    return pcs, parent_pcs, weights, pc_init, maxiter


#---------------------------------------------------------------------------------------------------------
#digitize and flatten a list of walker PCs
#---------------------------------------------------------------------------------------------------------

#returns a list of the bin assignments of the coordinates in pcs
def digitize_flatten_trj(pcs, binbounds, n_discrete_pc_vals):

    #put walkers in discrete bins
    connectivity_bins = np.digitize(pcs[:,-1,0], binbounds)
    #path_bins = [int(round(i)) for i in pcs[:,-1,1]]

    #maximum number of bins along the continuous PC axis; these may not all be occupied
    #n_connectivity_bins = len(binbounds)+1

    #flatten data in 2d PC space into a 1d PC space consisting of blocks along the discrete PC
    #bins_flattened = [pb*n_connectivity_bins + cb for pb, cb in zip(path_bins, connectivity_bins)]

    return connectivity_bins #bins_flattened


#---------------------------------------------------------------------------------------------------------
#digitize and flatten parent and child walker PCs and return a list of discrete transitions
#---------------------------------------------------------------------------------------------------------

def pcs_2_transitions(pcs, parent_pcs, binbounds, n_discrete_pc_vals):

    walker_bins_flattened = digitize_flatten_trj(pcs, binbounds, n_discrete_pc_vals)
    parent_bins_flattened = digitize_flatten_trj(parent_pcs, binbounds, n_discrete_pc_vals)
    
    return np.stack((parent_bins_flattened, walker_bins_flattened)).transpose()


#---------------------------------------------------------------------------------------------------------
#get a list of transitions and their weights using the methods above
#---------------------------------------------------------------------------------------------------------

def h5_2_transitions(h5path, miniter, maxiter, binbounds, n_discrete_pc_vals):

    pcs, parent_pcs, weights, pc_init, maxiter = h5_2_pcs(h5path, miniter, maxiter)
    transitions = pcs_2_transitions(np.concatenate(pcs), np.concatenate(parent_pcs), binbounds, n_discrete_pc_vals)

    return transitions, weights, pc_init, maxiter


#---------------------------------------------------------------------------------------------------------
#weight walkers based on the MSM weights of the bins they occupy and the number of walkers in each bin 
#so that the total weight of all walkers in each bin equals that bin's MSM probability
#---------------------------------------------------------------------------------------------------------

def get_walker_mi_weights(pcs, binbounds, n_discrete_pc_vals, pyem):

    #put walkers in discrete bins
    connectivity_bins = np.digitize(pcs[:,0,0], binbounds)
    path_bins = [int(round(i)) for i in pcs[:,0,1]]

    #maximum number of bins along the continuous PC axis; these may not all be occupied
    n_connectivity_bins = len(binbounds)+1

    #flatten data in 2d PC space into a 1d PC space consisting of blocks along the discrete PC
    bins_flattened = [pb*n_connectivity_bins + cb for pb, cb in zip(path_bins, connectivity_bins)]

    bin_counts = Counter(bins_flattened)

    # print(np.where(pyem.active_set == 231)[0])
    # print(np.where(pyem.active_set == 231)[0][0])
    try:
        walker_mi_weights = [pyem.stationary_distribution[np.where(pyem.active_set == walker_bin)[0][0]]/bin_counts[walker_bin] for walker_bin in bins_flattened]
        
    except IndexError as e:
        print("The following error probably arises from insufficient sampling. This leads to some states for which data are present being ergodically trimmed from the MSM, causing an index mismatch. All walkers will be given equal weights instead of MSM-derived ones:\n")
        print(e)
        print("")
        walker_mi_weights = [1/(len(bins_flattened)*bin_counts[walker_bin]) for walker_bin in bins_flattened]
        
    return walker_mi_weights


#--------------------------------------------------------------
# Plot PyEMMA markov state model energies
#--------------------------------------------------------------
#TODO write method spec.

# Map PyEMMA MSM state indices back to progress coordinate values, convert them to energies,
# and plot a line for the energy as a function of the continuous PC at each discrete PC value

# set std=True for bayesian MSM only
def plot_2d_pc_webins(pyem, nbins, binset, pclims, discrete_pc_vals, plottitle="", zeromin=False):

    bins = binset

    # -------------------------------------------------------
    # convert populations to energies
    energy = [-np.log(p) for p in pyem.stationary_distribution]
    if zeromin:
        min_energy = min(energy)
        energy = [ei - min_energy for ei in energy]

    # -------------------------------------------------------
    # align PC values for wires at each discrete PC value and plot them against each other
    # 'align' here means to account for the index shift
    # due to PyEMMA discarding disconnected energy landscape components,
    # causing it to output fewer equilibrium probabilities than the original number of bins

    # empty list of bin indices for each PC 2 value
    # note that these are the bins used by h5_2_transitions_forpyemma() to discretize the PC for MSM construction
    # which need not match the original westpa bins
    bin_vals = [[] for x in range(discrete_pc_vals)]
    # empty lists of MSM state energies and standard deviations thereof for each PC 2 value
    energies = [[] for x in range(discrete_pc_vals)]
    # empty lists of MSM state occupancies and standard deviations thereof for each PC 2 value
    occupancies = [[] for x in range(discrete_pc_vals)]

    # pyem.active_set is a list of the input indices of the bins which are part of the largest connected component
    # state i in the MSM returned by PyEMMA therefore corresponds to bin pyem.active_set[i]
    # note that these are bin indices in the flattened data sent to PyEMMA, not just PC 1 bin indices
    # also note that they are 1-indexed
    for i, b in enumerate(pyem.active_set):
        # the value of PC 2
        #  Though naturally distributed in 2D PC1-PC2 space, bins were flattened
        #  for transition matrix (and MSM) construction, so the transition matrix provided to PyEMMA
        #  has dimensions of (discrete_pc_vals*nbins) x (discrete_pc_vals*nbins),
        #  with all the bins sharing each PC2 value grouped together in order.
        #  The first nbins PyEMMA input states therefore have the lowest PC 2 value,
        #  the next nbins PyEMMA input states the second lowest, etc.
        #  I think 1 is subtracted because pyem.active_set treats the input states as 1-indexed.
        #  This assumes a 0-indexed integer PC 2 value.
        pc2_val = (b - 1) // nbins

        # separate MSM data by PC 2 value

        # determine the bin index of the current bin along PC 1 and add it to the list
        # note that (b-1)%nbins = b-1-nbins*pc2_val
        bin_vals[pc2_val].append((b - 1) % nbins)

        # add MSM state energies and occupancies to lists for subsequent plotting
        energies[pc2_val].append(energy[i])

        occupancies[pc2_val].append(pyem.stationary_distribution[i])

    #get x coordinates of bins in active set for each discrete PC value
    bin_x_all = [[bins[b] for b in bin_vals[i]] for i in range(discrete_pc_vals)]

    return bin_x_all, energies
