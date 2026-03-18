import pymol
from pymol import cmd
import os

def tmd_query_pymol():
    segment_resis = [[77, 149], [192, 245], [298, 362], [988, 1034], [857, 889], [900, 942], [1094, 1154]]
    return " or ".join([f"resi {sr[0]}-{sr[1]}" for sr in segment_resis])

cmd.delete("all")

#make the movie using every trj_frame_increment-th trajectory frame
#i.e. a value of 2 takes every other frame
trj_frame_increment = 2

cmd.load("/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/lip_glpg_2/topology/input.gro")
cmd.load_traj("/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/lip_glpg_2/001986-000211-trj-pbcmol-centered-tmd-rot-s10.xtc", interval = trj_frame_increment)

util.cbag()
cmd.hide("sticks")
#water hydrogen spheres look bad in smoothed trajectories because the water rotation time and smoothing window length are similar
cmd.hide("spheres", "elem H") 

cmd.show("spheres", "resn Cl- or resn Na+")
cmd.color("palegreen", "resn Cl-")
#cmd.color("grey", "resn Na+")

cmd.show("spheres", "resn LJP and not elem H")
cmd.show("spheres", "resn LJP and name H13")
util.cbac("resn LJP")

resis = "873+933+229+233+236+304+305+308+309+312+313+316+928+930+931+932+869"

cmd.show("spheres", f"input and resi {resis} and not name C+N+O and not elem H")
cmd.show("spheres", "input and resi 931 and name C+N+O")
cmd.show("spheres", "resn Arg and resi 933 and name HE+HH11+HH12+HH21+HH22")


util.cbay(f"input and resi {resis}")

cmd.color("firebrick", "elem O and not resn TP3")

#view from extracellular end looking towards cytoplasm
# cmd.set_view((\
#     -0.882278562,   -0.464601129,   -0.075704679,\
#      0.466403872,   -0.884541273,   -0.007132889,\
#     -0.063650310,   -0.041601930,    0.997103751,\
#     -0.000160344,    0.000113741, -192.223220825,\
#     54.672969818,   56.323814392,   80.269195557,\
#    179.534149170,  204.912185669,  -20.000000000 ))

cmd.smooth()

cmd.set("orthoscopic", "1")
cmd.bg_color("white")

views = {"membrane_plane": (\
    -0.037405413,    0.026684867,    0.998940945,\
     0.998739898,   -0.032339491,    0.038263123,\
     0.033326045,    0.999118030,   -0.025440669,\
    -0.000245899,   -0.000372943, -165.307769775,\
    63.780288696,   34.571521759,   81.286361694,\
    10.415435791,  238.533187866,  -20.000000000 ),
    "outside_end": (\
     0.297796905,    0.950370848,   -0.090073213,\
    -0.954521239,    0.297836810,   -0.013304352,\
     0.014183471,    0.089939050,    0.995845735,\
     0.000267893,   -0.000010717, -181.879959106,\
    54.426780701,   44.406246185,   90.498329163,\
   178.872360229,  227.773956299,  -20.000000000 ),
    "membrane_plane_small": (\
    -0.037405413,    0.026684867,    0.998940945,\
     0.998739898,   -0.032339491,    0.038263123,\
     0.033326045,    0.999118030,   -0.025440669,\
    -0.000286520,   -0.000380380, -136.630310059,\
    63.619785309,   38.534919739,   81.208564758,\
   -18.268718719,  209.849044800,  -20.000000000 ),
    "outside_end_small": (\
     0.470223397,    0.877939403,   -0.090073213,\
    -0.882009268,    0.471038520,   -0.013304352,\
     0.030747997,    0.085701607,    0.995845735,\
     0.000237361,   -0.000040471, -149.483673096,\
    61.426898956,   40.454261780,   91.073120117,\
   146.470474243,  195.372085571,  -20.000000000 )
   }

vk = "outside_end_small" #"membrane_plane_small" #"membrane_plane" #
cmd.set_view(views[vk])
width = 720
height = 512
cmd.viewport(width, height) #will not work if in fullscreen or maximized window
print("view set")

#number of movie frames per second; should not affect file size
fps = 10 #when exporting .mpg files, pymol says the lowest legal frame rate is 23.976 fps and uses that instead
#seems to work differently for .mov and .mpg files but affects size of both
quality = 30

#turn off for view debugging
make_movie = True
if make_movie:
    cmd.set("movie_fps", fps)
    cmd.mset(f"1-{cmd.count_states()}")
    cmd.do(f"movie.produce /home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/trajectory_movies/movies/cftri-c10_we2_{vk}_q{quality}_{fps}fps_{width}x{height}px_bgwhite.mov, quality={quality}")

























#trjlen = 596

#for vk in views.keys(): #this does not work with moviemaker

# vk = "membrane_plane_small" #"outside_end" #
# cmd.set_view(views[vk])
# cmd.viewport(1000,720) #will not work if in fullscreen or maximized window
# print("view set")

# cmd.set("movie_fps", 15)

# #turn off for view debugging
# make_movie = True
# if make_movie:
#     cmd.do(f"mset 1-{trjlen}")
#     cmd.do(f"movie.produce /home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/trajectory_movies/movies/abbv-cftri-c10-we2_{vk}.mpg, quality=70")


#waters 39575 and 41931 hitchhike with LJP from the protein (wherein they begin; one begins bonded to pyrazole) into bulk water

#these are in the tutorial but caused all but 1 frame not to render
# cmd.do("mview store, 1, state=1, object=input")
# cmd.do(f"mview store, {trjlen}, state={trjlen}, object=input")