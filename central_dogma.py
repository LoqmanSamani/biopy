from Bio.Seq import Seq
from Bio import SeqUtils
import random


alphabet = ["A", "C", "G", "T"]

random.seed(42)  # ensure reproducibility

seq1 = Seq("".join(random.choice(alphabet) for _ in range(100)))  # generate a random sequence with the length of 100

print(seq1)

"""
AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTAT
"""

print(seq1.count("AA"))  # Count the number of non-overlapping occurrences of 'AA'

print(seq1.count_overlap("AA"))  # Count the number of overlapping occurrences of 'AA'

print(SeqUtils.GC(seq1))  # GC percentage

print(SeqUtils.gc_fraction(seq1))  # GC percentile




#  Transcription:  DNA --> RNA


forward = Seq("".join(random.choice(alphabet) for _ in range(100)))
reverse = forward.reverse_complement()
complement = forward.complement()


print(f"5_{forward}_3")
print(f"5_{reverse}_3")

"""
5_TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATA_3
5_TATGACGAGAGGGTCATTTAAACATTGAAGAGGACGTCTTGTTTGGTCTGGTAACACGGACGAGGTATTGTGCCTTACTGCGGCTGCTAGACTCCTTCAA_3
"""

print(f"5_{forward}_3")
print(f"3_{complement}_5")

"""
5_TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATA_3
3_AACTTCCTCAGATCGTCGGCGTCATTCCGTGTTATGGAGCAGGCACAATGGTCTGGTTTGTTCTGCAGGAGAAGTTACAAATTTACTGGGAGAGCAGTAT_5
"""


rna1 = forward.transcribe()

print(forward)
print(rna1)

"""
TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATA
UUGAAGGAGUCUAGCAGCCGCAGUAAGGCACAAUACCUCGUCCGUGUUACCAGACCAAACAAGACGUCCUCUUCAAUGUUUAAAUGACCCUCUCGUCAUA
"""

#  back transcription

forward1 = rna1.back_transcribe()

print(f"Original template:    {forward}")
print(f"Transcribed template: {rna1}")
print(f"Regenerated template: {forward1}")

"""
Original template:    TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATA
Transcribed template: UUGAAGGAGUCUAGCAGCCGCAGUAAGGCACAAUACCUCGUCCGUGUUACCAGACCAAACAAGACGUCCUCUUCAAUGUUUAAAUGACCCUCUCGUCAUA
Regenerated template: TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATA
"""


#  Translation:  DNA --> polypeptide chain  or RNA --> polypeptide chain

poly1 = forward.translate()
poly2 = rna1.translate()

print(f"DNA-Template: {forward}")
print(f"Polypeptide from DNA: {poly1}")
print("__________________________________________________________________________________")
print(f"RNA: {rna1}")
print(f"Polypeptide from RNA: {poly2}")


"""
DNA-Template: TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAATGACCCTCTCGTCATA
Polypeptide from DNA: LKESSSRSKAQYLVRVTRPNKTSSSMFK*PSRH
__________________________________________________________________________________

RNA: UUGAAGGAGUCUAGCAGCCGCAGUAAGGCACAAUACCUCGUCCGUGUUACCAGACCAAACAAGACGUCCUCUUCAAUGUUUAAAUGACCCUCUCGUCAUA
Polypeptide from RNA: LKESSSRSKAQYLVRVTRPNKTSSSMFK*PSRH
"""





