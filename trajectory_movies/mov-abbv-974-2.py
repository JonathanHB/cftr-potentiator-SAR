import pymol
from pymol import cmd
import os

def tmd_query_pymol():
    segment_resis = [[77, 149], [192, 245], [298, 362], [988, 1034], [857, 889], [900, 942], [1094, 1154]]
    return " or ".join([f"resi {sr[0]}-{sr[1]}" for sr in segment_resis])

cmd.delete("all")

cmd.load("/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/nonlip_glpg_2/topology/input.gro")
cmd.load_traj("/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/nonlip_glpg_2/2.5A-20A/000691-000198-trj-pbcmol-centered-tmd-rot-s10.xtc")

util.cbag()
cmd.hide("sticks")
#water hydrogen spheres look bad in smoothed trajectories because the water rotation time and smoothing window length are similar
cmd.hide("spheres", "elem H") 

cmd.show("spheres", "resn Cl- or resn Na+")
cmd.color("palegreen", "resn Cl-")
#cmd.color("grey", "resn Na+")

cmd.show("spheres", "resn LJP and not elem H")
cmd.show("spheres", "resn LJP and name H14")
util.cbac("resn LJP")

resis = "873+933+229+233+236+304+305+308+309+312+313+316+928+930+931+932+869"

cmd.show("spheres", f"input and resi {resis} and not name C+N+O and not elem H")
cmd.show("spheres", "input and resi 931 and name C+N+O")
cmd.show("spheres", "resn Arg and resi 933 and name HE+HH11+HH12+HH21+HH22")


util.cbay(f"input and resi {resis}")

cmd.color("firebrick", "elem O and not resn TP3")


cmd.smooth()

cmd.set("orthoscopic", "1")

views = {"membrane_plane": (\
     0.929763377,   -0.094949186,   -0.355696768,\
    -0.360433549,   -0.037993312,   -0.932006478,\
     0.074980721,    0.994754553,   -0.069548033,\
     0.000168606,   -0.000308543, -134.784576416,\
    50.596939087,   79.852500916,   85.034301758,\
    27.431268692,  160.785385132,   20.000000000 ),

    "outside_end": (\
    -0.948640406,    0.303268284,   -0.090073213,\
    -0.303331971,   -0.952790082,   -0.013304352,\
    -0.089855865,    0.014701462,    0.995845735,\
    -0.000085017,    0.000250384, -181.889022827,\
    34.538902283,   51.821651459,   88.789581299,\
   178.132415771,  228.513900757,  -20.000000000 )}

trjlen = 206

vk = "outside_end" #"membrane_plane" #
cmd.set_view(views[vk])
cmd.viewport(1000,720) #will not work if in fullscreen or maximized window
print("view set")

cmd.set("movie_fps", 15)

#turn off for view debugging
make_movie = True
if make_movie:
    cmd.do(f"mset 1-{trjlen}")
    cmd.do(f"movie.produce /home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/trajectory_movies/movies/abbv-974_we2_{vk}.mpg, quality=70")


#waters 39575 and 41931 hitchhike with LJP from the protein (wherein they begin; one begins bonded to pyrazole) into bulk water

#these are in the tutorial but caused all but 1 frame not to render
# cmd.do("mview store, 1, state=1, object=input")
# cmd.do(f"mview store, {trjlen}, state={trjlen}, object=input")