from Bio.SeqIO import parse  # import parse function from SeqIO module
from Bio.SeqIO.FastaIO import SimpleFastaParser  # this will be used in case of big input data
import time





# Method 1

seqs = {}  # a dictionary to store sequences from a fasta file

path = "/home/sam/Documents/projects/biopy/files/sample.txt"  # path to the fasta file



with open(path, "r") as file:  # open the file

    for sequence in parse(handle=file, format="fasta"):  # iterate through the file

        # store each sequence as an item, key(sequence-id) and values (length of the sequence and the sequence), in the dictionary
        seqs[sequence.id] = [len(sequence.seq), str(sequence.seq)]



#  iterate through the dictionary (seqs) and print each key, value pair
for item in seqs:

    print(item)
    print(seqs[item])



"""
Output:

gi|2765658|emb|Z78533.1|CIZ78533
[740, 'CGTAACAAGGTTTCCGTGAGGTGGACGGATGCTGGCAGCAGCTGCCGTGCGAATCCCCCATGCT...']
gi|2765657|emb|Z78532.1|CCZ78532
[753, 'CCCAAGTCTTTGAACGCAAGGTCAGGCGGGGCCATGCTGGCAGCAGCTCCCGCTGAGTTGAGGC...']

...

"""






# Method 2

path1 = "/home/sam/Documents/projects/biopy/files/samp4.fasta"  # use another file

seqs1 = []  # use a list to store the sequences

with open(file=path, mode="r") as file1:

    for seq1 in parse(handle=file1, format="fasta"):

        seqs1.append(str(seq1.seq))

        print(seq1.id)
        print(seq1.seq)
        print(len(seq1.seq))





"""
output:

gi|2765565|emb|Z78440.1|PPZ78440
CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCAC...TTTTTCCCCCAATC
744
gi|2765564|emb|Z78439.1|PBZ78439
CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTT...GCTTTAGTTGGGCC
592
...

"""



print(len(seqs1))




"""
output:
94
"""


print(seqs1)






path2 = "/home/sam/Documents/projects/xenoMrna.fa.gz"  # the file-size = 6.7G (very big)


start = time.time()

count = 0  # count the number of sequences in the file

with open(file=path2, mode="r") as file2:

    for seq2_id, seq2 in SimpleFastaParser(file2):

        count += 1


end = time.time()


process_time = end - start  # calculate the process time


print("Process-time: ", process_time)
print("Number of Sequences: ", count)



"""
output:

Process-time: 544.55364549243774
Number of Sequences: 26345479

"""







# Another relative big file

path3 = "/home/sam/Documents/projects/biopy/files/samp6.fasta"


start3 = time.time()

count3 = 0

with open(file=path3, mode="r") as file3:

    for seq3_id, seq3 in SimpleFastaParser(file3):

        count3 += 1

end3 = time.time()

pro_time = end3 - start3

print(f"process-time: {pro_time}")
print(f"Number of Sequences: {count3}")








# Exercises

# 1) Create a function that takes the name of a FASTA file as input and returns its content (sequences) as list

path4 = "/home/sam/Documents/projects/biopy/files/samp1.fasta"


def read_fasta(path, ids=None, length=None, seqs=True):

    seqs2 = []

    with open(file=path, mode="r") as exercise1:

        for seq4 in parse(handle=exercise1, format="fasta"):

            if seqs:

                if ids and length:

                    seqs2.append([seq4.id, len(str(seq4.seq)), str(seq4.seq)])

                elif ids:

                    seqs2.append([seq4.id, str(seq4.seq)])

                elif length:

                    seqs2.append([len(seq4.seq), str(seq4.seq)])

            elif not seqs:

                if ids and length:

                    seqs2.append([seq4.id, len(str(seq4.seq))])

                elif ids:

                    seqs2.append(seq4.id)

                elif length:

                    seqs2.append(len(seq4.seq))

    return seqs2



result1 = read_fasta(path4, ids=None, length=None)

print(result1)



# 2) Create a function that takes the name of a FASTA file as input and returns a list of FASTA identifiers


result2 = read_fasta(path4, ids=True, seqs=False)

print(result2)

