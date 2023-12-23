from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.SeqUtils import GC


"""
 !wget https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.fasta
wget -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz https://files.rcsb.org/download/6YYT.pdb https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.fasta
--2023-12-23 22:15:58--  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz
           => ‘SRR003265.filt.fastq.gz.1’
Resolving ftp.1000genomes.ebi.ac.uk (ftp.1000genomes.ebi.ac.uk)... 193.62.193.167
Connecting to ftp.1000genomes.ebi.ac.uk (ftp.1000genomes.ebi.ac.uk)|193.62.193.167|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /vol1/ftp/phase3/data/NA18489/sequence_read ... done.
==> SIZE SRR003265.filt.fastq.gz ... 28919712
==> PASV ... done.    ==> RETR SRR003265.filt.fastq.gz ... done.
Length: 28919712 (28M) (unauthoritative)

SRR003265.filt.fastq.gz.1                              100%[===========================================================================================================================>]  27.58M  12.5MB/s    in 2.2s

2023-12-23 22:16:03 (12.5 MB/s) - ‘SRR003265.filt.fastq.gz.1’ saved [28919712]

--2023-12-23 22:16:03--  https://files.rcsb.org/download/6YYT.pdb
Resolving files.rcsb.org (files.rcsb.org)... 128.6.159.245
Connecting to files.rcsb.org (files.rcsb.org)|128.6.159.245|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: unspecified [application/octet-stream]
Saving to: ‘6YYT.pdb’

6YYT.pdb                                                   [   <=>                                                                                                                      ] 954.91K  1.55MB/s    in 0.6s

2023-12-23 22:16:04 (1.55 MB/s) - ‘6YYT.pdb’ saved [977832]

--2023-12-23 22:16:04--  https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.fasta
Resolving raw.githubusercontent.com (raw.githubusercontent.com)... 2606:50c0:8003::154, 2606:50c0:8001::154, 2606:50c0:8000::154, ...
Connecting to raw.githubusercontent.com (raw.githubusercontent.com)|2606:50c0:8003::154|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 76480 (75K) [text/plain]
Saving to: ‘ls_orchid.fasta’

ls_orchid.fasta                                        100%[===========================================================================================================================>]  74.69K  --.-KB/s    in 0.02s

2023-12-23 22:16:05 (3.67 MB/s) - ‘ls_orchid.fasta’ saved [76480/76480]

FINISHED --2023-12-23 22:16:05--
Total wall clock time: 6.3s
Downloaded: 3 files, 29M in 2.8s (10.1 MB/s)

"""

path = "/home/sam/Documents/projects/biopy/files/ls_orchid.fasta"

seqs = []

for sequence in SeqIO.parse(handle=path, format="fasta"):
    seqs.append(sequence.seq)


print(seqs[:10])
"""
[
Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), 
Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC'), 
Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA'), 
Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT'), 
Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA'), 
Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC'), 
Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT'), 
Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA'), 
Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC'), 
Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
]
"""


lengths = []
for seq in seqs:
    lengths.append(len(list(seq)))


plt.hist(lengths, bins=100, color="red")
plt.xlabel("Sequence length (bp)")
plt.ylabel("Frequency")
plt.title("Sequence Length vs. Frequency")
plt.show()


print(len(seqs))
""" 94 """

gc_contents = sorted([GC(seq) for seq in seqs])

plt.plot(gc_contents)
plt.xlabel("Sequence")
plt.ylabel("GC Content")
plt.title("Sequence GC %")
plt.show()

