from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from rosetta import *
from pyrosetta import *
from Bio import SeqIO
import random
from math import exp

# parser = PDBParser(PERMISSIVE=1)
# structure_id = "T0779"
# filename = "../Data/T0779.pdb"
# structure = parser.get_structure(structure_id, filename)

pyrosetta.init()

for seq_record in SeqIO.parse("../Data/T0866/T0866.fasta", "fasta"):
    seq_T0866 = seq_record.seq

for seq_record in SeqIO.parse("../Data/T0868/T0868.fasta", "fasta"):
    seq_T0868 = seq_record.seq

for seq_record in SeqIO.parse("../Data/T0880/T0880.fasta", "fasta"):
    seq_T0880 = seq_record.seq

def prob_function(new_energy, current_energy, current_temp):
   pro =  exp(-(new_energy - current_energy)/current_temp)
   num = random.randint(1, 100)*0.01
   if(pro > num):
       return True
   else:
       return False

def function(filename, fragmentFile):
    for seq_record in SeqIO.parse(filename, "fasta"):
        seq = seq_record.seq

    seqString = ''.join(str(e) for e in seq)
    init_pose = pose_from_sequence(seqString)

    fragments = open(fragmentFile, "r")
    current_count = 1
    while current_count <= len(seq):
        positionLine = fragments.readline()
        temp = positionLine.split(" ")
        temp = filter(None, temp)
        if temp[0] == "position:" and int(temp[1]) == current_count:
            fragments.readline()
            for j in range(0, 9):
                temp = fragments.readline().split(" ");
                temp = filter(None, temp)
                if current_count <= len(seq):
                    init_pose.set_phi(current_count, float(temp[5]))
                    init_pose.set_psi(current_count, float(temp[6]))
                    init_pose.set_omega(current_count, float(temp[7]))
                    current_count = current_count + 1
                else:
                    break

    init_pose.dump_pdb("initial_structure.pdb")

    temperature = 100
    full_scorefxn = create_score_function('ref2015')

    current_energy = full_scorefxn(init_pose)
    current_s = init_pose
    ebest = current_energy
    sbest = init_pose

    while temperature > 0:

        new_s = current_s
        # random replace one fragment
        position_random = random.randint(1, len(seq) - 9)
        fragments = open(fragmentFile, "r")
        while 1:
            positionLine = fragments.readline()
            temp = positionLine.split(" ")
            temp = filter(None, temp)
            if temp[0] == "position:" and int(temp[1]) == position_random:
                fragment_random = random.randint(1, 200)
                for i in range(0, fragment_random * 10):
                    fragments.readline()

                fragments.readline()
                for j in range(0, 9):
                    temp = fragments.readline().split(" ");
                    temp = filter(None, temp)
                    print(temp)
                    new_s.set_phi(position_random + j, float(temp[5]))
                    new_s.set_psi(position_random + j, float(temp[6]))
                    new_s.set_omega(position_random + j, float(temp[7]))

                break

        new_score = full_scorefxn(new_s)

        if prob_function(new_score, current_energy, temperature):
            current_s = new_s
            current_energy = new_score
            current_s.dump_pdb("temperature:" + str(temperature) + ".pdb")

        if ebest > new_score:
            sbest = new_s
            ebest = new_score


        temperature -= 0.1

    sbest.dump_pdb("best_structure.pdb")


function("../Data/T0866/T0866.fasta","../Data/T0880/t000_.200.9mers")




