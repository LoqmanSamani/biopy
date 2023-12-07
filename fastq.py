from Bio import SeqIO
import pandas as pd
from collections import Counter
import gzip
import matplotlib.pyplot as plt
import seaborn as sns


path = "/home/sam/Documents/projects/biopy/files/samp3.fastq"

sequences = {}

with open(file=path, mode="r") as file:

    for sequence in SeqIO.parse(handle=file, format="fastq"):

        sequences[sequence.id] = [len(sequence.seq), str(sequence.seq)]


print(sequences)

"""

{
  'SRR26411560.1.1': [301, 'CCTACGGGTGGCAGCAGTGAGGAATATTGGTCAATGGTCGGAAGTCTGAACCAGCCATGCCGCGTGCAGGATGAATGCCTTATGGGTTGTAAACTGCTTTTATATGGGAAGAATAAGCAGTACGCGTACTTTGATGACGGTACCATATGAATAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTCATACGGAGGATGCGAGCGTTATTCGGAATCATTGGGTTTAAAGGGTCTGTAGGCGGGCTATTAAGTCCGGGGTGAAAGGTTTCAGCTTAACTGAGAAATTGCCTTTGATTCTGG'],
  'SRR26411560.1.2': [300, 'GACTACCAGGGTCTCTAATCCTGTTCGCTCCCCACGCTTTCGTCCCTCAGCGTCAATTGTCTCCCAGTCCCCTGCCTTCGCAATCGGTGTTCTGCGTCATATCTAAGCATTTCACCGCTACACTACACATTCCAAGAACCTCACAGACACTCAAGTCTACCAGTATCAAAGGCAATTTCTCAGTTAAGCCGAAACCTTTCACCCCTGACTTAATACCCCGCCTACCGACCCTTTCAACCCATTGATTCCGAATACCGCCCGCACCCCCCGGATTACCCCGGCTGCTGGCACTGTATTCAC'], 
  'SRR26411560.2.1': [301, 'CCTACGGGTGGCAGCAGTGGGGAATATTGGACAATGGGGGGAACCCTGATCCAGCAATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTTAAGCTCTTTTAGTAGGGAGGAAAGGTCCGTCGTTTATACCTGTGCAAGTGTCGCTCCCCACCGATGAAGCACCGGCTACCTCCGCGCCAGCAGCCGCGTTATTACGGTGGGTTCGTGCGCTATTCGGAATTTTTTGCCGTAAAGCGCACGCAGGCGGTTGCCCAAGTCAGTTGTTATCCCCCCGGCCCTAACCCGTGACCTTTCTTTGAAA'], 
  'SRR26411560.2.2': [300, 'GACTACTCGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACATCAGCGTCAGTTTCGACCCAGTAAGCTGCCTTCGCCATCGGTGTTCCTCCTAATATCTACGCATTTCACCGCTCCCCCTCAAATTCCGCCTTCCTCTCCCACACTCTCGTCGCCCAGTTTCAACTGCAGTTCCCAGTTTAAGCCCGGGGCTTTCCCCCCTGACTTGCTCACCCGCCTGCGTCCCCTTTCAGCACCTTAATTCCGACTCACCCCCGCACCCCCCGTATTACCGTGCCTGCTGCCACGCATTTTCC'], 
  'SRR26411560.3.1': [301, 'CCTACGGGTGGCAGCAGTGGGGAATATTGGACAATGGGCGCAACCCTGATCCAGCAATGCCGCGTGTGTGATGAAGGCCTTCGGGTTGTAAACCTCTTTTTGTAGGGAGGAAAGGTCCTTTGTTCTTACCTTTTTATCTGACGCTACCTCTAGACGATGCACCGCCTCCCTCCGTTCCAGCAGCCCCGTTCTTACGCTCCGTGCGTGCTTTTTTCGGATTTACTGTGCGTTAAGCGCACGCAGTTGGTTGCCCAAGTCACTTGTGAACGCCCCGTCCTCCTCCTGCGTACGTCCTTTGAAA'], 
  'SRR26411560.3.2': [300, 'GACTACTCGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACATCAGCGTCAGTTCCCACCCAGTAACCTGCCTTCGCCATCGGTGTTCCTCCTCATATCTACGCATTTCACCGCTACACCTGAACTTCCGCTCTCCCCTCCCACACTCTCGTCCCCCAGTTTCAAATGCACTTCCCACTTTTACCCCGCGGCTTTCACACCTTCCTTGCCCACCCCCCTGCCTCCCCTTTCCGCCCCTTAATTCCAAACCACCCCACCCCCCCCACGTTTTCCCGCGCTTCTTGACGTGCATTTAC'], 
  ...
}

"""

# Download a fastq file and analyse it

"""
(base) sam@sam-IdeaPad-L340-15IRH-Gaming:~/Documents/projects$ rm -f SRR003265.filt.fastq.gz 2>/dev/null
(base) sam@sam-IdeaPad-L340-15IRH-Gaming:~/Documents/projects$ wget -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz
--2023-12-07 20:05:09--  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz
           => ‘SRR003265.filt.fastq.gz’
Resolving ftp.1000genomes.ebi.ac.uk (ftp.1000genomes.ebi.ac.uk)... 193.62.193.167
Connecting to ftp.1000genomes.ebi.ac.uk (ftp.1000genomes.ebi.ac.uk)|193.62.193.167|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /vol1/ftp/phase3/data/NA18489/sequence_read ... done.
==> SIZE SRR003265.filt.fastq.gz ... 28919712
==> PASV ... done.    ==> RETR SRR003265.filt.fastq.gz ... done.
Length: 28919712 (28M) (unauthoritative)

SRR003265.filt.fastq.gz                                100%[===========================================================================================================================>]  27.58M  10.0MB/s    in 2.7s    

2023-12-07 20:05:13 (10.0 MB/s) - ‘SRR003265.filt.fastq.gz’ saved [28919712]

"""

path2 = "/home/sam/Documents/projects/biopy/files/samp7.gz"

records = SeqIO.parse(gzip.open(filename=path2, mode='rt', encoding='utf-8'), format='fastq')

record = next(records)

print(record)
print(record.id, record.description, record.seq)
print(record.letter_annotations)

"""
ID: SRR003265.31
Name: SRR003265.31
Description: SRR003265.31 3042NAAXX:3:1:1252:1819 length=51
Number of features: 0
Per letter annotation for: phred_quality
Seq('GGGAAAAGAAAAACAAACAAACAAAAACAAAACACAGAAACAAAAAAACCA')
SRR003265.31 SRR003265.31 3042NAAXX:3:1:1252:1819 length=51 GGGAAAAGAAAAACAAACAAACAAAAACAAAACACAGAAACAAAAAAACCA
{'phred_quality': [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 30, 23, 40, 32, 35,
 29, 40, 16, 40, 40, 32, 35, 31, 40, 40, 39, 22, 40, 24, 20, 28, 31, 12, 31, 10, 22, 28, 13, 26, 20, 23, 23]}
"""


# records1 = list(records)

# print(records1)

"""
id='SRR003265.9654212', name='SRR003265.9654212', description='SRR003265.9654212 3042NAAXX:3:99:436:1041 length=51', dbxrefs=[]),SeqRecord(seq=Seq('ATGTCTGCGACCACACGATGCAGGTGTAAACGTATTGTTGTACACTGGCTC'), 
id='SRR003265.9654213', name='SRR003265.9654213', description='SRR003265.9654213 3042NAAXX:3:99:219:762 length=51', dbxrefs=[]), SeqRecord(seq=Seq('CATCCCAGTTACTTGGGGGGCCAAGGCACGATAAATGTTTGCACCACGGAG'), 
id='SRR003265.9654214', name='SRR003265.9654214', description='SRR003265.9654214 3042NAAXX:3:99:39:750 length=51', dbxrefs=[]), SeqRecord(seq=Seq('TGAGAACTAGATTAGAATGAAGTTTAACTGTGAACTAGTGAGTTGGTATAT'),
... 
"""



# The following code calculates the total percentage for each nucleotide.

counts = {}

for record in records:
    count = Counter(record.seq)
    for key in count.keys():
        if key not in counts:
            counts[key] = count[key]
        else:
            counts[key] += count[key]

num_nuc = sum(counts.values())

new_counts = {}

for key in counts.keys():
    new_counts[key] = counts[key] / num_nuc


print(new_counts)

"""

{
   'G': 0.20676842412896096,
   'A': 0.2859598043556052,
   'T': 0.295796307602681,
   'C': 0.21003681594631063,
   'N': 0.0014386479664422215
}

"""





# Distribution of nucleotide reads


def pos_count(data, nucleotide):

    count_nuc = {}

    for rec in data:
        for index, letter in enumerate(str(rec.seq)):
            if letter == nucleotide:
                nuc = 1
            else:
                nuc = 0

            if index not in count_nuc:
                count_nuc[index] = nuc
            else:
                count_nuc[index] += nuc

    return count_nuc





path3 = "/home/sam/Documents/projects/biopy/files/samp7.gz"

sequences = SeqIO.parse(gzip.open(filename=path3, mode='rt', encoding='utf-8'), format='fastq')

count_a = pos_count(data=sequences, nucleotide="A")




path4 = "/home/sam/Documents/projects/biopy/files/samp7.gz"

sequences1 = SeqIO.parse(gzip.open(filename=path4, mode='rt', encoding='utf-8'), format='fastq')

count_t = pos_count(data=sequences1, nucleotide="T")



path5 = "/home/sam/Documents/projects/biopy/files/samp7.gz"

sequences2 = SeqIO.parse(gzip.open(filename=path5, mode='rt', encoding='utf-8'), format='fastq')

count_c = pos_count(data=sequences2, nucleotide="C")




path6 = "/home/sam/Documents/projects/biopy/files/samp7.gz"

sequences3 = SeqIO.parse(gzip.open(filename=path6, mode='rt', encoding='utf-8'), format='fastq')

count_g = pos_count(data=sequences3, nucleotide="G")





path7 = "/home/sam/Documents/projects/biopy/files/samp7.gz"

sequences4 = SeqIO.parse(gzip.open(filename=path7, mode='rt', encoding='utf-8'), format='fastq')

count_n = pos_count(data=sequences4, nucleotide="N")


print(count_a)
print(count_t)
print(count_c)
print(count_g)
print(count_n)




"""

count_a:

{0: 152145, 1: 155972, 2: 151872, 3: 149034, 4: 154439, 5: 155548, 6: 152621, 7: 151544, 8: 149383, 9: 153301, 10: 151083,
 11: 150111, 12: 149784, 13: 149006, 14: 148620, 15: 148471, 16: 149472, 17: 147655, 18: 148528, 19: 147385, 20: 148280,
 21: 147784, 22: 146868, 23: 147502, 24: 147539, 25: 146249, 26: 146350, 27: 145667, 28: 144766, 29: 143747, 30: 143644, 
 31: 142783, 32: 144022, 33: 143806, 34: 143159, 35: 142746, 36: 142523, 37: 142289, 38: 141624, 39: 141675, 40: 141227, 
 41: 140103, 42: 140042, 43: 139296, 44: 138261, 45: 137566, 46: 136796, 47: 135910, 48: 134705, 49: 131936, 50: 127126}
 
 

count_t:
 
{0: 160065, 1: 138195, 2: 145048, 3: 138774, 4: 139001, 5: 131352, 6: 141270, 7: 141815, 8: 143626, 9: 142690, 10: 144314,
 11: 147136, 12: 147666, 13: 147431, 14: 146957, 15: 147821, 16: 148077, 17: 146938, 18: 146885, 19: 146252, 20: 146330, 
 21: 146601, 22: 146904, 23: 145839, 24: 146217, 25: 146873, 26: 146971, 27: 147483, 28: 148028, 29: 148709, 30: 147059, 
 31: 149010, 32: 147067, 33: 149650, 34: 149611, 35: 152633, 36: 152734, 37: 155440, 38: 155024, 39: 155130, 40: 155378, 
 41: 156777, 42: 157750, 43: 158810, 44: 160130, 45: 160753, 46: 163989, 47: 164178, 48: 166372, 49: 166748, 50: 181374}
 
 
 
count_c:
 
{0: 98497, 1: 105842, 2: 104313, 3: 104964, 4: 102753, 5: 108969, 6: 106073, 7: 107464, 8: 109547, 9: 107253, 10: 106438, 
 11: 106461, 12: 106919, 13: 108889, 14: 110100, 15: 108710, 16: 108138, 17: 110523, 18: 109600, 19: 111468, 20: 108558, 
 21: 110179, 22: 110609, 23: 110738, 24: 110131, 25: 110105, 26: 110343, 27: 109449, 28: 109229, 29: 108307, 30: 109508, 
 31: 107461, 32: 108220, 33: 107374, 34: 109436, 35: 107385, 36: 107221, 37: 105414, 38: 105863, 39: 105351, 40: 105370, 
 41: 105231, 42: 104781, 43: 105183, 44: 103259, 45: 104517, 46: 101036, 47: 102118, 48: 101659, 49: 104107, 50: 92990}



count_g:

{0: 97519, 1: 108217, 2: 106993, 3: 115454, 4: 112033, 5: 112357, 6: 108262, 7: 107403, 8: 105670, 9: 104982, 10: 106391,
 11: 104518, 12: 103857, 13: 102900, 14: 102549, 15: 103224, 16: 102539, 17: 103110, 18: 103213, 19: 103121, 20: 105058,
 21: 103662, 22: 103845, 23: 104147, 24: 104339, 25: 104938, 26: 104514, 27: 103835, 28: 103731, 29: 104355, 30: 103823,
 31: 104994, 32: 104255, 33: 103408, 34: 104352, 35: 104456, 36: 104087, 37: 104524, 38: 105261, 39: 105491, 40: 105483,
 41: 105895, 42: 105270, 43: 104436, 44: 105431, 45: 104976, 46: 105539, 47: 105556, 48: 105067, 49: 104784, 50: 105510}
    
 
 
count_n:
   
{0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 
 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 61, 26: 48, 27: 1792, 28: 2472, 29: 3108, 30: 4192, 31: 3978, 
 32: 4662, 33: 3988, 34: 1668, 35: 1006, 36: 1661, 37: 559, 38: 454, 39: 579, 40: 768, 41: 220, 42: 383, 43: 501, 
 44: 1145, 45: 414, 46: 866, 47: 464, 48: 423, 49: 651, 50: 1226}

"""


a = {
     0: 152145, 1: 155972, 2: 151872, 3: 149034, 4: 154439, 5: 155548, 6: 152621, 7: 151544, 8: 149383, 9: 153301, 10: 151083,
     11: 150111, 12: 149784, 13: 149006, 14: 148620, 15: 148471, 16: 149472, 17: 147655, 18: 148528, 19: 147385, 20: 148280,
     21: 147784, 22: 146868, 23: 147502, 24: 147539, 25: 146249, 26: 146350, 27: 145667, 28: 144766, 29: 143747, 30: 143644,
     31: 142783, 32: 144022, 33: 143806, 34: 143159, 35: 142746, 36: 142523, 37: 142289, 38: 141624, 39: 141675, 40: 141227,
     41: 140103, 42: 140042, 43: 139296, 44: 138261, 45: 137566, 46: 136796, 47: 135910, 48: 134705, 49: 131936, 50: 127126
}


keys1 = list(a.keys())
values1 = list(a.values())


plt.bar(keys1, values1)
plt.xlabel('Index')
plt.ylabel('Count')
plt.title('Nucleotide (A) Counts at Each Position')
plt.show()


t = {
     0: 160065, 1: 138195, 2: 145048, 3: 138774, 4: 139001, 5: 131352, 6: 141270, 7: 141815, 8: 143626, 9: 142690, 10: 144314,
     11: 147136, 12: 147666, 13: 147431, 14: 146957, 15: 147821, 16: 148077, 17: 146938, 18: 146885, 19: 146252, 20: 146330,
     21: 146601, 22: 146904, 23: 145839, 24: 146217, 25: 146873, 26: 146971, 27: 147483, 28: 148028, 29: 148709, 30: 147059,
     31: 149010, 32: 147067, 33: 149650, 34: 149611, 35: 152633, 36: 152734, 37: 155440, 38: 155024, 39: 155130, 40: 155378,
     41: 156777, 42: 157750, 43: 158810, 44: 160130, 45: 160753, 46: 163989, 47: 164178, 48: 166372, 49: 166748, 50: 181374
}


keys2 = list(t.keys())
values2 = list(t.values())


plt.bar(keys2, values2)
plt.xlabel('Index')
plt.ylabel('Count')
plt.title('Nucleotide (T) Counts at Each Position')
plt.show()



c = {
     0: 98497, 1: 105842, 2: 104313, 3: 104964, 4: 102753, 5: 108969, 6: 106073, 7: 107464, 8: 109547, 9: 107253, 10: 106438,
     11: 106461, 12: 106919, 13: 108889, 14: 110100, 15: 108710, 16: 108138, 17: 110523, 18: 109600, 19: 111468, 20: 108558,
     21: 110179, 22: 110609, 23: 110738, 24: 110131, 25: 110105, 26: 110343, 27: 109449, 28: 109229, 29: 108307, 30: 109508,
     31: 107461, 32: 108220, 33: 107374, 34: 109436, 35: 107385, 36: 107221, 37: 105414, 38: 105863, 39: 105351, 40: 105370,
     41: 105231, 42: 104781, 43: 105183, 44: 103259, 45: 104517, 46: 101036, 47: 102118, 48: 101659, 49: 104107, 50: 92990
}

keys3 = list(c.keys())
values3 = list(c.values())


plt.bar(keys3, values3)
plt.xlabel('Index')
plt.ylabel('Count')
plt.title('Nucleotide (C) Counts at Each Position')
plt.show()



g = {
     0: 97519, 1: 108217, 2: 106993, 3: 115454, 4: 112033, 5: 112357, 6: 108262, 7: 107403, 8: 105670, 9: 104982, 10: 106391,
     11: 104518, 12: 103857, 13: 102900, 14: 102549, 15: 103224, 16: 102539, 17: 103110, 18: 103213, 19: 103121, 20: 105058,
     21: 103662, 22: 103845, 23: 104147, 24: 104339, 25: 104938, 26: 104514, 27: 103835, 28: 103731, 29: 104355, 30: 103823,
     31: 104994, 32: 104255, 33: 103408, 34: 104352, 35: 104456, 36: 104087, 37: 104524, 38: 105261, 39: 105491, 40: 105483,
     41: 105895, 42: 105270, 43: 104436, 44: 105431, 45: 104976, 46: 105539, 47: 105556, 48: 105067, 49: 104784, 50: 105510
}


keys4 = list(g.keys())
values4 = list(g.values())


plt.bar(keys4, values4)
plt.xlabel('Index')
plt.ylabel('Count')
plt.title('Nucleotide (G) Counts at Each Position')
plt.show()



n = {
     0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0,
     18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 61, 26: 48, 27: 1792, 28: 2472, 29: 3108, 30: 4192, 31: 3978,
     32: 4662, 33: 3988, 34: 1668, 35: 1006, 36: 1661, 37: 559, 38: 454, 39: 579, 40: 768, 41: 220, 42: 383, 43: 501,
     44: 1145, 45: 414, 46: 866, 47: 464, 48: 423, 49: 651, 50: 1226
}


keys5 = list(n.keys())
values5 = list(n.values())


plt.bar(keys5, values5)
plt.xlabel('Index')
plt.ylabel('Count')
plt.title('Nucleotide (N) Counts at Each Position')
plt.show()




path8 = "/home/sam/Documents/projects/biopy/files/samp7.gz"

sequences5 = SeqIO.parse(gzip.open(filename=path8, mode='rt', encoding='utf-8'), format='fastq')

phred_quality = []
length = []
# max_length = max(length)  # 51

for sequence in sequences5:

    phred_quality.append(sequence.letter_annotations["phred_quality"])
    length.append(len(sequence.letter_annotations["phred_quality"]))


# Create a DataFrame from the data
data1 = pd.DataFrame(phred_quality, columns=range(1, 51 + 1))


# Plot boxplot using Seaborn
plt.figure(figsize=(15, 8))
sns.boxplot(data=data1, width=0.5)
plt.xlabel('Position')
plt.ylabel('Phred Quality')
plt.title('Boxplot of Phred Quality at Each Position')
plt.show()


