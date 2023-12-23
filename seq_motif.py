from Bio import motifs
from Bio.Seq import Seq


# Create a simple DNA motif

samples = [Seq("AATGCGAT"),
           Seq("TTGCCAGG"),
           Seq("GGTACGTA"),
           Seq("TTCATACT"),
           Seq("GGGCCCAT"),
           Seq("TTTCGATG")]

mot = motifs.create(samples)

print(mot)
"""
AATGCGAT
TTGCCAGG
GGTACGTA
TTCATACT
GGGCCCAT
TTTCGATG
"""

print(mot.instances)

print(mot.counts)
"""
      0      1      2      3      4      5      6      7
A:   1.00   1.00   0.00   2.00   0.00   3.00   2.00   1.00
C:   0.00   0.00   1.00   3.00   4.00   1.00   1.00   0.00
G:   2.00   2.00   2.00   1.00   1.00   2.00   1.00   2.00
T:   3.00   3.00   3.00   0.00   1.00   0.00   2.00   3.00
"""

print(mot.counts["A"])
""" [1, 1, 0, 2, 0, 3, 2, 1] """

print(mot.counts[:, 2])
""" {'A': 0, 'C': 1, 'G': 2, 'T': 3} """

print(mot.consensus)
""" TTTCCAAT """

print(mot.anticonsensus)
""" CCATATCC """

print(mot.degenerate_consensus)
""" KKKMCRNK """

