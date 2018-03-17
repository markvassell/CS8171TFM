from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from rosetta import *
from pyrosetta import *
from Bio import SeqIO
import random

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


def function(filename, fragmentFile):
    for seq_record in SeqIO.parse(filename, "fasta"):
        seq = seq_record.seq

    seqString = ''.join(str(e) for e in seq)
    pose = pose_from_sequence(seqString)

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
                    pose.set_phi(current_count, float(temp[5]))
                    pose.set_psi(current_count, float(temp[6]))
                    pose.set_omega(current_count, float(temp[7]))
                    current_count = current_count + 1
                else:
                    break

    pose.dump_pdb(str(current_count) + ".pdb")

    # position_random =  random.randint(1, len(seq) - 9)
    # fragments = open(fragmentFile, "r")
    #
    # while 1:
    #     positionLine = fragments.readline()
    #     temp = positionLine.split(" ")
    #     temp = filter(None, temp)
    #     if temp[0] == "position:" and int(temp[1]) == position_random:
    #         fragment_random = random.randint(1, 200)
    #         for i in range(0,fragment_random*10):
    #             fragments.readline()
    #
    #         fragments.readline()
    #         for j in range(0, 9):
    #             temp = fragments.readline().split(" ");
    #             temp = filter(None, temp)
    #             print(temp)
    #             pose.set_phi(position_random + j, float(temp[5]))
    #             pose.set_psi(position_random + j, float(temp[6]))
    #             pose.set_omega(position_random + j, float(temp[7]))
    #
    #         break


    # full_scorefxn = create_score_function('ref2015')
    # pose_score = full_scorefxn(pose)


function("../Data/T0866/T0866.fasta","../Data/T0880/t000_.200.9mers")




