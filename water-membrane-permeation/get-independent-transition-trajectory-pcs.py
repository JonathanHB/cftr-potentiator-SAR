import numpy as np
import h5py
import os
import sys
sys.path.insert(0,'../independent-transitions/')
#put westpa_msm_functions.py in the same folder I believe
from get_independent_transitions import get_independent_transitions
from walker_ancestors import walker_ancestors
#Jonathan Borowsky
#Grabe lab
#062025

#-----------------------------------------------methods-----------------------------------------------

#------------------------
# Parameters: 
# pc: a progress coordinate (currently a 1d array [or list] of floats, 
# but could be changed to a 2d array if the PC had multiple time points 
# and the indexing in get_unique_transitions() was adjusted accordingly)
#------------------------
# Returns: 
# macrostate: An integer. The possible values should be consecutive integers counting up from 0, and -1.
#             takes values in the range [-1,4] in this implementation
# note that this method must return a number other than -1 for the macrostate of the initial structure

def pc_2_macrostate(pc):
    pc0_min = 2.5
    pc0_max = 20 
    
    if pc[0] < pc0_min: 
        macrostate = 0
    elif pc[0] > pc0_max:
        macrostate = 1
    else:
        return -1
        
    return macrostate


def get_transition_representatives(h5path, macrostate_classifier, n_macrostates, maxround):

    transitions, wbf, itt, pcs_by_tss = get_independent_transitions(h5path, macrostate_classifier, n_macrostates, minround=0, maxround=maxround)
    tcounts = np.array([[len(j) for j in i] for i in transitions])
    print(tcounts)
    
    #most_advanced_descendant_of_each_independent_transition
    madoeit = []

    for k in pcs_by_tss.keys():
        #print(k)
        pcs = [i[2] for i in pcs_by_tss[k][1]]

        #print(len(pcs))
        #print(max(pcs))
        #print(min(pcs))

        #print(pcs_by_tss[k][1][np.argmax(pcs)])
        madoeit.append(pcs_by_tss[k][1][np.argmax(pcs)])
        #print()
    #transition_rep_sets = [tr_0_i for tr_0_i in transitions[0][1:] if len(tr_0_i) > 0]

    return [madoeit]

n_macrostates=2


#---------------------------------paths and target walker information---------------------------------

pathdict = {"nonlip_glpg_1": ["/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_cftr_1_degrabo", "../../jhb-simulation-data/wstp_nonlip_glpg_1_continued", "west-050225.h5", 2000],
            "nonlip_glpg_2": ["/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_cftr_2_wynton",  "../../jhb-simulation-data/wstp_nonlip_glpg_2_continued", "west-050225.h5", 1000],
            "lip_glpg_1":    ["/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_lip_glpg_1",     "../../jhb-simulation-data/wstp_lip_glpg_1_continued",    "west-050125.h5", 2000],
            "lip_glpg_2":    ["../../chloe-cftr-project/wstp_lip_glpg_2",                                  "../../jhb-simulation-data/wstp_lip_glpg_2_continued",    "west-050125.h5", 2000]}

abspath = os.getcwd()

#get PCs for reactive trajectories with each walker
for k in pathdict.keys():

    we_data_paths = pathdict[k]

    we_data_path_1 = we_data_paths[0]
    we_data_path_2 = we_data_paths[1]
    west_fn = we_data_paths[2]
    maxround = we_data_paths[3]

    h5path = f'{we_data_path_2}/{west_fn}'


    #-----------------collect trajectory ancestors for each transition representative---------------------
    print(f"reading {h5path}")

    #get transition representatives
    transition_representatives = get_transition_representatives(h5path, pc_2_macrostate, n_macrostates, maxround)
    print(transition_representatives)

    #extract trajectory files of the transition representatives
    for ms_ind, tr_set in enumerate(transition_representatives):
        for tr in tr_set:
            print(tr)

            (ids, pcs) = walker_ancestors(h5path, tr[0]+1, tr[1])

            np.save(f"{k}_round_{tr[0]+1}_walker_{tr[1]}.npy", pcs)
            #print(pcs)


