from Bio import Entrez, SeqIO
from Bio.SeqUtils import molecular_weight as mw
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from Bio.Blast import NCBIWWW
from Bio import SearchIO
from Bio.PDB import PDBParser
import nglview as nv



# Analyse SARS-CoV-2

"""
The GenBank accession number for the complete genome sequence of
the severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2), 
which is the virus responsible for the COVID-19 pandemic.

            **MN908947**
"""
Entrez.email = "samaniloqman91@gmail.com"  # access email







path1 = "/home/sam/Documents/projects/biopy/files/6YYT.pdb"
parser = PDBParser()

structure = parser.get_structure(id="6YYT", file=path1)

print(structure)
""" <Structure id=6YYT> """


for chain in structure[0]:
    print(chain.id)
"""
A
B
C
D
P
Q
T
U
"""


nv.show_biopython(structure, gui=True)


