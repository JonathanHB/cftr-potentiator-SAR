import pymol
from pymol import cmd
import os

def tmd_query_pymol():
    segment_resis = [[77, 149], [192, 245], [298, 362], [988, 1034], [857, 889], [900, 942], [1094, 1154]]
    return " or ".join([f"resi {sr[0]}-{sr[1]}" for sr in segment_resis])

cmd.delete("all")

cmd.load("/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/nonlip_glpg_1/topology/input.gro")
cmd.load_traj("/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/nonlip_glpg_1/001913-000187-trj-pbcmol-centered-tmd-rot-s10.xtc")

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
cmd.bg_color("black")

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
   178.132415771,  228.513900757,  -20.000000000 ),
    "membrane_plane_small": (\
     0.929763377,   -0.094949186,   -0.355696768,\
    -0.360433549,   -0.037993312,   -0.932006478,\
     0.074980721,    0.994754553,   -0.069548033,\
     0.000162765,   -0.000322100, -112.460617065,\
    51.537811279,   79.235702515,   87.557777405,\
     5.171977997,  138.526092529,   20.000000000 ),
    "outside_end_small":(\
    -0.934112728,    0.345433086,   -0.090073213,\
    -0.345682442,   -0.938255370,   -0.013304352,\
    -0.089107640,    0.018709399,    0.995845735,\
    -0.000092924,    0.000237629, -156.497390747,\
    35.693954468,   47.645637512,   88.835281372,\
   152.737869263,  203.119354248,  -20.000000000 )
    }


trjlen = 573

#for vk in views.keys(): #this does not work with moviemaker

vk = "outside_end_small" #"membrane_plane_small" #
cmd.set_view(views[vk])
width = 720
height = 512
cmd.viewport(width, height) #will not work if in fullscreen or maximized window
print("view set")

#number of movie frames per second; should not affect file size
fps = 10 #when exporting .mpg files, pymol says the lowest legal frame rate is 23.976 fps and uses that instead
#make the movie using every trj_frame_increment-th trajectory frame
#i.e. a value of 2 takes every other frame
trj_frame_increment = 2
#seems to work differently for .mov and .mpg files and only sometimes affects the size of the former???
quality = 20

#turn off for view debugging
make_movie = True
if make_movie:
    cmd.set("movie_fps", fps)

    nstates = cmd.count_states()
    mset_arg = " ".join([str(1+i) for i in range(1, nstates, trj_frame_increment)])
    cmd.mset(mset_arg)
    cmd.do(f"movie.produce /home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/trajectory_movies/movies/abbv-974_we1_{vk}_q{quality}_{fps}fps_{width}x{height}px.mov, quality={quality}")


#waters 39575 and 41931 hitchhike with LJP from the protein (wherein they begin; one begins bonded to pyrazole) into bulk water

#these are in the tutorial but caused all but 1 frame not to render
# cmd.do("mview store, 1, state=1, object=input")
# cmd.do(f"mview store, {trjlen}, state={trjlen}, object=input")