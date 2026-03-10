#Jonathan Borowsky
#Grabe lab
#062025

import h5py

def walker_ancestors(h5path, walker_round, walker_num):
    
    #load h5 file
    with h5py.File(h5path, 'r') as f:
        
        pcoords = []
        walker_ids = []
        
        #determine names of westpa rounds
        #note that round 1's name is at index 0
        iterations = [iter for iter in f["iterations"]]
        
        for iter_ind in range(walker_round-1, -1, -1):
        
            #this extracts only the iteration name, not the iteration data
            iter_name = iterations[iter_ind] 
            #using the iteration name to extract the data
            
            #log walker ID and progress coordinate
            #zeros are for trimming excess nested array layers
            pcoords.append(f["iterations"][iter_name]["pcoord"][walker_num][0][0])
            walker_ids.append(walker_num)
            
            #update round and walker numbers
            #current_round -= 1
            if iter_ind == 0:
                break
            walker_num = f["iterations"][iter_name]["seg_index"][walker_num][1]

    walker_ids.reverse()
    pcoords.reverse()
    return [walker_ids, pcoords]