from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from rosetta import *
from pyrosetta import *
from Bio import SeqIO
import random
import math
import sys
import matplotlib.pyplot as plt

# parser = PDBParser(PERMISSIVE=1)
# structure_id = "T0779"
# filename = "../Data/T0779.pdb"
# structure = parser.get_structure(structure_id, filename)

pyrosetta.init()

def prob_function(new_energy, current_energy, current_temp):
    if new_energy < current_energy:
        return True
    else:
        prob = math.exp((current_energy - new_energy)/current_temp*0.01)
        random_num = random.randint(1, 10)*0.1
        # print('prosibility = ', prob, 'random number = ', random_num)
        if prob > random_num:
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
        try:
            if temp[0] == "position:":
                if int(temp[1]) == current_count and int(temp[1]) < len(seq)-2:
                    fragments.readline()
                    for j in range(0, 3):
                        temp = fragments.readline().split(" ");
                        temp = filter(None, temp)
                        init_pose.set_phi(current_count, float(temp[5]))
                        init_pose.set_psi(current_count, float(temp[6]))
                        init_pose.set_omega(current_count, float(temp[7]))
                        current_count = current_count + 1

                elif(int(temp[1]) == len(seq) - 2):
                    print('secondary')
                    print(temp)
                    fragments.readline()
                    for j in range(0, 3):
                        temp = fragments.readline().split(" ");
                        temp = filter(None, temp)
                        print(temp)
                        init_pose.set_phi(int(temp[1])+j, float(temp[5]))
                        init_pose.set_psi(int(temp[1])+j, float(temp[6]))
                        init_pose.set_omega(int(temp[1])+j, float(temp[7]))

                    break

        except:
            print('Caught error: ', sys.exc_info()[0])
            break


    init_pose.dump_pdb("initial_structure.pdb")

    temperature = 100
    full_scorefxn = create_score_function('ref2015')

    current_energy = full_scorefxn(init_pose)
    current_s = init_pose
    ebest = current_energy
    sbest = init_pose

    energy_history = []

    while temperature > 0:
        new_s = current_s
        # random replace one fragment
        position_random = random.randint(1, len(seq) - 3)
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
                for j in range(0, 3):
                    temp = fragments.readline().split(" ");
                    temp = filter(None, temp)
                    try:
                        new_s.set_phi(position_random + j, float(temp[5]))
                        new_s.set_psi(position_random + j, float(temp[6]))
                        new_s.set_omega(position_random + j, float(temp[7]))
                    except:
                        continue

                break

        new_score = full_scorefxn(new_s)
        # print("new score = ",new_score,"current_energy", current_energy,temperature)
        if prob_function(new_score, current_energy, temperature):
            current_s = new_s
            current_energy = new_score
            energy_history.append(current_energy)
            current_s.dump_pdb("temperature" + str(temperature) + ".pdb")

        print('new score', new_score, 'e best',ebest)

        if ebest > new_score:
            sbest = new_s
            ebest = new_score

        temperature -= 0.1

    sbest.dump_pdb("best_structure.pdb")
    print(energy_history)



function("../Data/T0866/T0866.fasta","../Data/T0866/t000_.200.3mers")




