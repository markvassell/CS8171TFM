# from rosetta import *
from pyrosetta.rosetta import *
from pyrosetta import *
# from Bio.PDB import PDBParser
# import nglview as nv
from math import exp
from random import randint
init()
#pose = pose_from_pdb("../Data/T0880/T0880TS023_1.pdb")
# pose = pose_from_pdb("~/home/markvassell/Project_2/CS8171TFM/Data/T0880/T0880TS023_1.pdb")

# print (pose.is_fullatom())
# print(pose.residue(5))
scorefxn = get_fa_scorefxn()
#print(scorefxn.show(pose))

def prob_function(current_energy, pass_energy, current_temp):
	return exp(-(current_energy - pass_energy)/current_temp) 



# print(prob_function(100,200,100))

current_temp = 100
structure = pose_from_pdb("../Data/T0880/T0880TS023_1.pdb")
print(scorefxn.show(structure))
sbest = "initial structure" # same as structure to start off with
ebest = "Energy of initial structure" # energy of the initial structure

while current_temp > 0:
	# s_new = 'new s' # Randomly pick a fragment and replace the angle to get new s
	# if prob_function('Energy of old structure', 'Energy of new structure', current_temp) > .85:
	# 	s = new_s
	# if ebest > energy(s_new):
	# 	sbest = s_new
	print(prob_function(randint(100,500),randint(100,500), current_temp))
	current_temp -= 1
	 

# Write best structure to pdb
# pose.dump_pdb("../Results/output_file.pdb") 






