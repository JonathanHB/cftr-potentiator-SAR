#Jonathan Borowsky
#Grabe lab
#062025

import os
import sys
from walker_ancestors import walker_ancestors

#---------------------------------paths and target walker information---------------------------------

#specify input file and walker
we_data_path_1 = sys.argv[1] #data from Chloe's initial simulation
we_data_path_2 = sys.argv[2] #data from my extension
west_path = sys.argv[3]

h5path = f'{we_data_path_2}/{west_path}'

walker_round = int(sys.argv[4])
walker_num = int(sys.argv[5])

#-----------------------------------------------------------------------------------------------------

#make directory to store relevant traj segs
trj_seg_dir = f"{str(walker_round).zfill(6)}-{str(walker_num).zfill(6)}-ancestors-2.5A-20A"
if not os.path.exists(trj_seg_dir):
    os.mkdir(trj_seg_dir)

#get ancestors of target walkers
walker_ancestry = walker_ancestors(h5path, walker_round, walker_num)
walker_ids = walker_ancestry[0]

print(walker_ancestry)

#extract trajectory files of the ancestors of the target trajectory
for round in range(1, walker_round+1):
    archive = f"round-{str(round).zfill(6)}-segs"

    use_path_2 = False
    
    if os.path.exists(f"{we_data_path_1}/traj_segs/{archive}.tar.gz"):
        os.system(f"tar -zxf {we_data_path_1}/traj_segs/{archive}.tar.gz traj_segs/{str(round).zfill(6)}/{str(walker_ids[round-1]).zfill(6)}/traj_comp.xtc")
        
        if not os.path.exists(f"traj_segs/{str(round).zfill(6)}/{str(walker_ids[round-1]).zfill(6)}/traj_comp.xtc"):
            use_path_2 = True
            
    else:
        use_path_2 = True

    if use_path_2:
        os.system(f"tar -zxf {we_data_path_2}/traj_segs/{archive}.tar.gz traj_segs/{str(round).zfill(6)}/{str(walker_ids[round-1]).zfill(6)}/traj_comp.xtc")

    #the untarred file comes out named 'traj_segs' and I can't figure out how to untar only one file directly from a .tar.gz to a nonstandard new location
    os.system(f"mv traj_segs/{str(round).zfill(6)}/{str(walker_ids[round-1]).zfill(6)}/traj_comp.xtc {trj_seg_dir}/{str(round).zfill(6)}-{str(walker_ids[round-1]).zfill(6)}-traj_comp.xtc") 

