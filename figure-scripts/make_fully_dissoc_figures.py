import pymol
from pymol import cmd

import os

color_by_run = False
outpath = "/home/jonathan/Documents/grabelab/cftr/figures/new"
upperpath = "/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation"
lip = False

colors = ["g", "c", "m"]
if lip:
    folders = [["lip_glpg_1", "001798-000087-ancestors-2.5A-20A"], ["lip_glpg_1", "001998-000130-ancestors-2.5A-20A"], ["lip_glpg_2", "001986-000211-ancestors-2.5A-20A"]]
    prefix = "CFTRi-C10"
else:
    folders = [["nonlip_glpg_1", "001655-000156-ancestors-2.5A-20A"], ["nonlip_glpg_1", "001913-000187-ancestors-2.5A-20A"], ["nonlip_glpg_2", "000691-000198-ancestors-2.5A-20A"]]
    prefix = "ABBV-974"


frame_data = [["unbound", [-2, -2, -2]]]

state_ind = 0

#------------------------------------load data--------------------------------------------

cmd.delete("all")

cmd.load(f"{upperpath}/nonlip_glpg_1/topology/input.gro", "viewref") #this is used to make sure the protein is aligned the same in both C10 and pyrazole derivatives

cmd.load(f"{upperpath}/{folders[0][0]}/topology/input.gro", "ref")

cmd.align("ref", "object viewref and (resi 77-149 or resi 192-245 or resi 298-362 or resi 988-1034 or resi 857-889 or resi 900-942 or resi 1094-1154)")

#load frames of interest from manageably sized .xtc files from individual WE rounds,
#  converting from subsampled frame indices from VMD to WE rounds and frame indices in .xtc files
for folder, frame, color in zip(folders, frame_data[state_ind][1], colors):
    if frame == -1:
        continue

    cmd.load(f"{upperpath}/{folder[0]}/topology/input.gro", f"{folder[0]}-{folder[1]}")

    we_round = int(round(frame*10/3))
    subframe = int(round(3*((frame*10/3) % 1)))

    if frame == -2:
        we_round = int(folder[1].split("-")[0])-1
        print(we_round)
        subframe = 2

    #print([f for f in os.listdir(f"{upperpath}/{folder[0]}/2.5A-20A/{folder[1]}/") if f[-14:] == "-traj_comp.xtc"])
    trj_files = [f for f in os.listdir(f"{upperpath}/{folder[0]}/2.5A-20A/{folder[1]}/") if (f.endswith("-traj_comp.xtc") and int(f[0:6]) == we_round+1)]

    cmd.load_traj(f"{upperpath}/{folder[0]}/2.5A-20A/{folder[1]}/{trj_files[0]}", f"{folder[0]}-{folder[1]}", start=subframe+1, stop=subframe+1)

    #align on TMD
    cmd.align(f"{folder[0]}-{folder[1]}", "ref and (resi 77-149 or resi 192-245 or resi 298-362 or resi 988-1034 or resi 857-889 or resi 900-942 or resi 1094-1154)")

    #the linter is wrong about 'util' being syntactically incorrect here
    if color == "g":
        util.cbag(f"{folder[0]}-{folder[1]}")
    elif color == "c":
        util.cbac(f"{folder[0]}-{folder[1]}")
    elif color == "m":
        util.cbam(f"{folder[0]}-{folder[1]}")
        
#------------------------------------graphics--------------------------------------------

#secondary structure assignment
cmd.dss("ref")

#hide/show and rescale
cmd.hide("everything")

cmd.show("cart", "ref and poly")

cmd.show("sticks", "ref and resi 229+233+236+304+305+308+309+312+313+316+928+930+931+932 and not name C+N+O")
cmd.show("spheres", "ref and resi 229+233+236+304+305+308+309+312+313+316+928+930+931+932 and name CA")
cmd.set("sphere_scale", 0.3, "ref and resi 229+233+236+304+305+308+309+312+313+316+928+930+931+932 and name CA")

cmd.show("spheres", "ref and resi 873+933 and not name C+N+O") #926, 931, and 932 form interactions in some cases, but not in others

cmd.show("spheres", "resname LJP and not ref and not viewref")

cmd.hide("sticks", "elem H") # and not (resname LJP and name H12+H13+H14)")
cmd.hide("spheres", "elem H")

cmd.show("spheres", "ref and name P31")
cmd.set("sphere_scale", 0.6, "ref and name P31")

#coloring
util.cbaw("poly")
util.cbao("resn LJP and not ref and not viewref")

if color_by_run:
    util.cbag("ref and resi 229+233+236+304+305+308+309+312+313+316+928+930+931+932+873+933")
    util.cbac(f"resn LJP and {folders[2][0]}-{folders[2][1]} and not ref and not viewref")


cmd.color("grey40", "ref and name P31")

#view settings

cmd.set("orthoscopic", "on")

cmd.set_view((\
    -0.889825165,    0.001684323,    0.456256568,\
    -0.456212103,   -0.018101078,   -0.889667928,\
     0.006769278,   -0.999811947,    0.016872900,\
     0.002988905,    0.005761471, -217.890579224,\
    39.281814575,   48.769935608,  116.061782837,\
   126.229881287,  309.521667480,   20.000000000 ))

# cmd.set_view((\
#     -0.889825165,    0.001684323,    0.456256568,\
#     -0.456212103,   -0.018101078,   -0.889667928,\
#      0.006769278,   -0.999811947,    0.016872900,\
#      0.003174372,    0.005853221, -180.625808716,\
#     43.932701111,   51.183319092,  116.882499695,\
#     88.953842163,  272.245513916,   20.000000000 ))

cmd.png(f"{outpath}/{prefix}_end_states_coloredbyrun.png", width=2400, height=1800, ray=True)
