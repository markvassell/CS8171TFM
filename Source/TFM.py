from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser


for seq_record in SeqIO.parse("../Data/T0779.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))



parser = PDBParser(PERMISSIVE=1)
structure_id = "T0779"
filename = "../Data/T0779.pdb"

structure = parser.get_structure(structure_id, filename)



