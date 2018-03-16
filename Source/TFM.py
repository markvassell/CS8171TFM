from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from rosetta import *
from pyrosetta import *
from Bio import SeqIO


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

    fragments = open(fragmentFile, "r")

    seqString = ''.join(str(e) for e in seq)
    pose = pose_from_sequence(seqString)

    loopNum = (len(seq) - 1) // 9 + 1

    for i in range(0, loopNum):
        positionLine = fragments.readline()
        tp2 = positionLine.split(" ");
        while tp2[0] != 'position:':
            tp = pos.readline()
            tp2 = tp.split(" ");
        pos.readline()
        for j in range(0, 9):
            tp = pos.readline()
            tp2 = tp.split(" ");
            while '' in tp2:
                tp2.remove('')
            # print((tp2))
            # break
            L2.append(tp2[5])
            L3.append(tp2[6])
    for i in range(1, pose.total_residue() + 1):
        pose.set_phi(i, -180)
        pose.set_psi(i, 180)
        pose.set_omega(i, 180)

    full_scorefxn = create_score_function('ref2015')
    # pose_score = full_scorefxn(pose)
    print('ScoreFunction', full_scorefxn(pose))

    # for i in range(1, pose.total_residue() + 1):
    #     pose.set_phi(i, -180)
    #     pose.set_psi(i, 180)
    #     pose.set_omega(i, 180)
    #
    # pose.dump_pdb('1.pdb')





function("../Data/T0866/T0866.fasta")




