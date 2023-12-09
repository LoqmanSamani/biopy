from Bio import pairwise2
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices as mat
import random
import math



align1 = pairwise2.align.globalxx("ACGTAGCTACGAAG", "ACGTACGTGG")  # globalxx - matches score 1, mismatches 0 and no gap penalty.

for align in align1:
    print(pairwise2.format_alignment(*align))

"""
ACGTAGC-TACGAAG
||||| | |  |  |
ACGTA-CGT--G--G
  Score=9

ACGTA-GCTACGAAG
||||| | |  |  |
ACGTACG-T--G--G
  Score=9
"""



# another way to do the same

aligner = PairwiseAligner()

align2 = aligner.align(seqA="ACGTAGCTAG", seqB="ACGTGG")

for align in align2:
    print(align)

"""
target            0 ACGTAGCTAG 10
                  0 ||||-|---| 10
query             0 ACGT-G---G  6
"""


# globalmx - # matches score 2, mismatches -1. No gap penalty.
align3 = pairwise2.align.globalmx("AGCGATGGCTGACTGAC", "ACTGCTGACGTCGACGG", match=2, mismatch=-1)

for align in align3:
    print(pairwise2.format_alignment(*align))

"""
AGCGATGGCTGAC-T-GAC--
| |  | |||||| | |||  
A-C--T-GCTGACGTCGACGG
  Score=26

AGCGATGGCTGAC-T-GAC--
| |  || ||||| | |||  
A-C--TG-CTGACGTCGACGG
  Score=26
"""


# second method

aligner2 = PairwiseAligner()

aligner.gap_score = -1
aligner.match = 2
aligner.mismatch = -1
aligner.query_end_gap_score = 0

align4 = aligner2.align(seqA="AGCGATGGCTGACTGAC", seqB="ACTGCTGACGTCGACGG")

for align in align4:
    print(align)


"""
target            0 AGCGATGGCTGAC-T-GAC-- 17
                  0 |-|--||-|||||-|-|||-- 21
query             0 A-C--TG-CTGACGTCGACGG 17

target            0 AGCGATGGCTGAC-T-GAC-- 17
                  0 |-|--|-||||||-|-|||-- 21
query             0 A-C--T-GCTGACGTCGACGG 17
"""



# globalxs - matches score 1, mismatches 0, opening gap -2, extended gap -1
align5 = pairwise2.align.globalxs("AGTGATGCTGACGGTAATC", "ATCATCATGCACGT", open=-2, extend=-1)

for align in align5:
    print(pairwise2.format_alignment(*align))

"""
AGTGATGCTGACGGTAATC
| |.||..||...||    
A-TCATCATGCACGT----
  Score=1
"""


random.seed(42)

seq1 = "".join([random.choice(["A", "C", "G", "G", "T"]) for _ in range(35)])
seq2 = "".join([random.choice(["A", "C", "G", "G", "T"]) for _ in range(25)])

print(seq1)
print(seq2)
"""
AAGCCCATATGAAACCTTATCTGCGTGACGGGCCG
AAGAGGTGAGTAGATGTGTCAACGA
"""




bl62 = mat.load("BLOSUM62")
print(bl62)
"""
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
     A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    B    Z    X    *
A  4.0 -1.0 -2.0 -2.0  0.0 -1.0 -1.0  0.0 -2.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0  1.0  0.0 -3.0 -2.0  0.0 -2.0 -1.0  0.0 -4.0
R -1.0  5.0  0.0 -2.0 -3.0  1.0  0.0 -2.0  0.0 -3.0 -2.0  2.0 -1.0 -3.0 -2.0 -1.0 -1.0 -3.0 -2.0 -3.0 -1.0  0.0 -1.0 -4.0
N -2.0  0.0  6.0  1.0 -3.0  0.0  0.0  0.0  1.0 -3.0 -3.0  0.0 -2.0 -3.0 -2.0  1.0  0.0 -4.0 -2.0 -3.0  3.0  0.0 -1.0 -4.0
D -2.0 -2.0  1.0  6.0 -3.0  0.0  2.0 -1.0 -1.0 -3.0 -4.0 -1.0 -3.0 -3.0 -1.0  0.0 -1.0 -4.0 -3.0 -3.0  4.0  1.0 -1.0 -4.0
C  0.0 -3.0 -3.0 -3.0  9.0 -3.0 -4.0 -3.0 -3.0 -1.0 -1.0 -3.0 -1.0 -2.0 -3.0 -1.0 -1.0 -2.0 -2.0 -1.0 -3.0 -3.0 -2.0 -4.0
Q -1.0  1.0  0.0  0.0 -3.0  5.0  2.0 -2.0  0.0 -3.0 -2.0  1.0  0.0 -3.0 -1.0  0.0 -1.0 -2.0 -1.0 -2.0  0.0  3.0 -1.0 -4.0
E -1.0  0.0  0.0  2.0 -4.0  2.0  5.0 -2.0  0.0 -3.0 -3.0  1.0 -2.0 -3.0 -1.0  0.0 -1.0 -3.0 -2.0 -2.0  1.0  4.0 -1.0 -4.0
G  0.0 -2.0  0.0 -1.0 -3.0 -2.0 -2.0  6.0 -2.0 -4.0 -4.0 -2.0 -3.0 -3.0 -2.0  0.0 -2.0 -2.0 -3.0 -3.0 -1.0 -2.0 -1.0 -4.0
H -2.0  0.0  1.0 -1.0 -3.0  0.0  0.0 -2.0  8.0 -3.0 -3.0 -1.0 -2.0 -1.0 -2.0 -1.0 -2.0 -2.0  2.0 -3.0  0.0  0.0 -1.0 -4.0
I -1.0 -3.0 -3.0 -3.0 -1.0 -3.0 -3.0 -4.0 -3.0  4.0  2.0 -3.0  1.0  0.0 -3.0 -2.0 -1.0 -3.0 -1.0  3.0 -3.0 -3.0 -1.0 -4.0
L -1.0 -2.0 -3.0 -4.0 -1.0 -2.0 -3.0 -4.0 -3.0  2.0  4.0 -2.0  2.0  0.0 -3.0 -2.0 -1.0 -2.0 -1.0  1.0 -4.0 -3.0 -1.0 -4.0
K -1.0  2.0  0.0 -1.0 -3.0  1.0  1.0 -2.0 -1.0 -3.0 -2.0  5.0 -1.0 -3.0 -1.0  0.0 -1.0 -3.0 -2.0 -2.0  0.0  1.0 -1.0 -4.0
M -1.0 -1.0 -2.0 -3.0 -1.0  0.0 -2.0 -3.0 -2.0  1.0  2.0 -1.0  5.0  0.0 -2.0 -1.0 -1.0 -1.0 -1.0  1.0 -3.0 -1.0 -1.0 -4.0
F -2.0 -3.0 -3.0 -3.0 -2.0 -3.0 -3.0 -3.0 -1.0  0.0  0.0 -3.0  0.0  6.0 -4.0 -2.0 -2.0  1.0  3.0 -1.0 -3.0 -3.0 -1.0 -4.0
P -1.0 -2.0 -2.0 -1.0 -3.0 -1.0 -1.0 -2.0 -2.0 -3.0 -3.0 -1.0 -2.0 -4.0  7.0 -1.0 -1.0 -4.0 -3.0 -2.0 -2.0 -1.0 -2.0 -4.0
S  1.0 -1.0  1.0  0.0 -1.0  0.0  0.0  0.0 -1.0 -2.0 -2.0  0.0 -1.0 -2.0 -1.0  4.0  1.0 -3.0 -2.0 -2.0  0.0  0.0  0.0 -4.0
T  0.0 -1.0  0.0 -1.0 -1.0 -1.0 -1.0 -2.0 -2.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0  1.0  5.0 -2.0 -2.0  0.0 -1.0 -1.0  0.0 -4.0
W -3.0 -3.0 -4.0 -4.0 -2.0 -2.0 -3.0 -2.0 -2.0 -3.0 -2.0 -3.0 -1.0  1.0 -4.0 -3.0 -2.0 11.0  2.0 -3.0 -4.0 -3.0 -2.0 -4.0
Y -2.0 -2.0 -2.0 -3.0 -2.0 -1.0 -2.0 -3.0  2.0 -1.0 -1.0 -2.0 -1.0  3.0 -3.0 -2.0 -2.0  2.0  7.0 -1.0 -3.0 -2.0 -1.0 -4.0
V  0.0 -3.0 -3.0 -3.0 -1.0 -2.0 -2.0 -3.0 -3.0  3.0  1.0 -2.0  1.0 -1.0 -2.0 -2.0  0.0 -3.0 -1.0  4.0 -3.0 -2.0 -1.0 -4.0
B -2.0 -1.0  3.0  4.0 -3.0  0.0  1.0 -1.0  0.0 -3.0 -4.0  0.0 -3.0 -3.0 -2.0  0.0 -1.0 -4.0 -3.0 -3.0  4.0  1.0 -1.0 -4.0
Z -1.0  0.0  0.0  1.0 -3.0  3.0  4.0 -2.0  0.0 -3.0 -3.0  1.0 -1.0 -3.0 -1.0  0.0 -1.0 -3.0 -2.0 -2.0  1.0  4.0 -1.0 -4.0
X  0.0 -1.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0  0.0  0.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -4.0
* -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0  1.0
"""

# globaldx - matching/mismatching scores read from blosum62 matrix, no gap penalty

align6 = pairwise2.align.globaldx(seq1, seq2, match_dict=bl62)

for align in align6:
    print(pairwise2.format_alignment(*align))

"""
AAGCCCATATGAAACCTTATC-TGC-GTGACG--G-G-C--CG-
|||     | |           ||  || | |  | | |  || 
AAG-----A-G----------GTG-AGT-A-GATGTGTCAACGA
  Score=92

AAGCCCATATGAAACCTTATC-TGC-GTGACG--G-G-C--CG-
|||   |   |           ||  || | |  | | |  || 
AAG---A---G----------GTG-AGT-A-GATGTGTCAACGA
  Score=92
"""


# globalmc - matches score 5, mismatches -4, gap penalty defined through function gap_function

def gap_function(x, y):  # x is gap position in seq, y is gap length
     if y == 0:  # No gap
        return 0
     elif y == 1:  # Gap open penalty
        return -2

     return - (2 + y/4.0 + math.log(y)/2.0)


align7 = pairwise2.align.globalmc(seq1, seq2, match=5, mismatch=-4, gap_A_fn=gap_function, gap_B_fn=gap_function)

for align in align7:
    print(pairwise2.format_alignment(*align))

"""
AAGCCCATA--TGA--A-ACCTTATCTGCGTG---ACGGGCCG-
|||     |  |||  | |  |       |||   |     || 
AAG-----AGGTGAGTAGA--T-------GTGTCAA-----CGA
  Score=51.3286

AAGCCCATA--TGA--A-ACCTTATCTGCGTG---ACGGGCCG-
|||     |  |||  | |       |  |||   |     || 
AAG-----AGGTGAGTAGA-------T--GTGTCAA-----CGA
  Score=51.3286

AAGCCCATA--TGA--A-ACCTTATCTGCGTG---ACGGGCCG-
|||     |  |||  | |       ||  ||   |     || 
AAG-----AGGTGAGTAGA-------TG--TGTCAA-----CGA
  Score=51.3286

AAGCCCATA--TGA--A-ACCTTATCTGCGTG---ACGGGCCG-
|||     |  |||  | |  |       |||   ||     | 
AAG-----AGGTGAGTAGA--T-------GTGTCAAC-----GA
  Score=51.3286

AAGCCCATA--TGA--A-ACCTTATCTGCGTG---ACGGGCCG-
|||     |  |||  | |       |  |||   ||     | 
AAG-----AGGTGAGTAGA-------T--GTGTCAAC-----GA
  Score=51.3286

AAGCCCATA--TGA--A-ACCTTATCTGCGTG---ACGGGCCG-
|||     |  |||  | |       ||  ||   ||     | 
AAG-----AGGTGAGTAGA-------TG--TGTCAAC-----GA
  Score=51.3286

AAGCCCATA--TGA--A-ACCTTATCTGCGTG---ACGGGCCG-
|||     |  |||  | |  |       |||   |||      
AAG-----AGGTGAGTAGA--T-------GTGTCAACG-----A
  Score=51.3286

AAGCCCATA--TGA--A-ACCTTATCTGCGTG---ACGGGCCG-
|||     |  |||  | |       |  |||   |||      
AAG-----AGGTGAGTAGA-------T--GTGTCAACG-----A
  Score=51.3286

AAGCCCATA--TGA--A-ACCTTATCTGCGTG---ACGGGCCG-
|||     |  |||  | |       ||  ||   |||      
AAG-----AGGTGAGTAGA-------TG--TGTCAACG-----A
  Score=51.3286

"""

