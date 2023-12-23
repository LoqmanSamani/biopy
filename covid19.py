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


handle = Entrez.efetch(
    db="nucleotide",
    id="MN908947",
    rettype="gb",
    retmode="text"
)

records = list(SeqIO.parse(handle=handle, format="gb"))

handle.close()  # close the connection

print(records)

"""
[
  SeqRecord(seq=Seq('ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGT...AAA'), 
  id='MN908947.3', 
  name='MN908947', 
  description='Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, 
  complete genome', 
  dbxrefs=[])
]
"""

covid_dna = records[0].seq  # extract the DNA sequence from the file

print(covid_dna)
"""
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGA...GACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

print(len(covid_dna))
""" 29903 """

covid_weight = mw(covid_dna)  # Molecular weight refers to the sum of the atomic weights of all atoms in a molecule
print(covid_weight)
""" 9241219.214400413  daltons (Da) """

gc_content = GC(covid_dna)
print(gc_content)
""" 37.97277865097148 """

nuc_count = {
    "A": covid_dna.count("A"),
    "G": covid_dna.count("G"),
    "C": covid_dna.count("C"),
    "T": covid_dna.count("T")
}
print(nuc_count)
"""
'A': 8954 
'G': 5863 
'C': 5492 
'T': 9594
"""

custom = sns.color_palette(palette="husl", n_colors=len(nuc_count.keys()))
sns.barplot(data=nuc_count, palette=custom, width=0.6)
plt.xlabel("Nucleotide")
plt.ylabel("Frequency")
plt.title("Nucleotide Frequency of COVID 19")
plt.show()

covid_mrna = covid_dna.transcribe()
print(covid_mrna)
"""
AUUAAAGGUUUAUACCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUC...AUGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

covid_aa = covid_dna.translate()
print(covid_aa)
print(len(covid_aa))
"""
The asterisk (*) represents a stop codon
IKGLYLPR*QTNQLSISCRSVL*TNFKICVAVTRLHA*CTHAV*LI...PM*F**LLRRMTKKKKKKKKKK

9967
"""

aa_frequency = Counter(covid_aa)
print(aa_frequency)
"""
Counter({
   'L': 886, 
   'S': 810, 
   '*': 774, 
   'T': 679, 
   'C': 635, 
   'F': 593, 
   'R': 558, 
   'V': 548, 
   'Y': 505, 
   'N': 472, 
   'I': 436, 
   'K': 413, 
   'G': 394, 
   'A': 375, 
   'H': 332, 
   'Q': 325, 
   'P': 292, 
   'D': 290, 
   'E': 270, 
   'W': 263, 
   'M': 117
})
"""

custom1 = sns.color_palette(palette="husl", n_colors=len(aa_frequency.keys()))
sns.barplot(data=aa_frequency, palette=custom1, width=0.6)
plt.xlabel("Amino Acid")
plt.ylabel("Frequency")
plt.title("Amino Acid Frequency of COVID 19")
plt.show()

poly_peps = covid_aa.split("*")  # split the sequence at any stop codon
print(poly_peps)
print(len(poly_peps))
"""
[
  Seq('IKGLYLPR'), 
  Seq('QTNQLSISCRSVL'), 
  Seq('TNFKICVAVTRLHA'), 
  Seq('CTHAV'), 
  Seq('LITNYCR'), 
  ...

775

"""
poly_peps1 = []

for poly_pep in poly_peps:
    if len(poly_pep) > 20:
        poly_peps1.append([len(poly_pep), poly_pep])

print(poly_peps1[0:5])

"""
[[35, Seq('QDTSNSSIFCRLLTVSSVLQPIISTSRFRPGVTER')], 
[46, Seq('DGEPCPWFQRENTRPTQFACFTGSRRARTWLWRLRGGGLIRGTSTS')], 
[21, Seq('TALCVHQTFGCSNCTSWSCYG')], 
[22, Seq('DTWCPCPSCGRNTSGLPQGSSS')], 
[24, Seq('HLQWGMSKFCISLKFHNQDYSTKG')]] 
"""

# The longest seq
longest = False

for pep in poly_peps1:
    if not longest:
        longest = pep
    if longest[0] < pep[0]:
        longest = pep

print(longest)

# save the longest for further analysis
with open("longest.fasta", "w") as file:
    file.write(f">covid protein\n{longest[1]}")




path = "/home/sam/Documents/projects/biopy/files/longest.fasta"

qov_pro = SeqIO.read(handle=path, format="fasta")

print(qov_pro.seq)

"""
CTIVFKRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVKRHTFSNYQHEETIYNLLKDCPAVAKHDFF
KFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEILVTYNCCDDDYFNKKDWYDFVENPDILRVYANLGERVRQALLKTVQFCDAMRNAGI
VGVLTLDNQDLNGNWYDFGDFIQTTPGSGVPVVDSYYSLLMPILTLTRALTAESHVDTDLTKPYIKWDLLKYDFTEERLKLFDRYFKYWDQTYHPNCVN
CLDDRCILHCANFNVLFSTVFPPTSFGPLVRKIFVDGVPFVVSTGYHFRELGVVHNQDVNLHSSRLSFKELLVYAADPAMHAASGNLLLDKRTTCFSVA
...
"""

blast_res = NCBIWWW.qblast("blastp", "pdb", qov_pro.seq)

qov_pro = SearchIO.read(blast_res, "blast-xml")

print(qov_pro[0:10])
"""
Program: blastp (2.14.1+)
  Query: unnamed (2701)
         protein product
 Target: pdb
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  pdb|7D4F|A  Chain A, RNA-directed RNA polymerase [Sever...
            1      1  pdb|6YYT|A  Chain A, nsp12 [Severe acute respiratory sy...
            2      1  pdb|6XEZ|A  Chain A, RNA-directed RNA polymerase [Sever...
            3      1  pdb|7BW4|A  Chain A, RNA-directed RNA polymerase [Sever...
            4      1  pdb|6XQB|A  Chain A, RNA-directed RNA polymerase [Sever...
            5      1  pdb|7BV1|A  Chain A, RNA-directed RNA polymerase [Sever...
            6      1  pdb|7C2K|A  Chain A, RNA-directed RNA polymerase [Sever...
            7      1  pdb|6M71|A  Chain A, RNA-directed RNA polymerase [Sever...
            8      1  pdb|7ED5|A  Chain A, RNA-directed RNA polymerase [Sever...
            9      1  pdb|7AAP|A  Chain A, Non-structural protein 12 [Severe ...

"""

for pro in qov_pro:
    print(pro.id)
    print(pro.description)
    print(pro[0].evalue)
    print(pro[0].bitscore)
    print(pro[0].aln)

"""
pdb|7D4F|A
Chain A, RNA-directed RNA polymerase [Severe acute respiratory syndrome coronavirus 2]
0.0
1938.7
Alignment with 2 rows and 926 columns
FKRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKT...LQA unnamed
LNRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKT...LQG pdb|7D4F|A
pdb|6YYT|A
Chain A, nsp12 [Severe acute respiratory syndrome coronavirus 2]
0.0
1938.31
Alignment with 2 rows and 925 columns
FKRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKT...VLQ unnamed
LNRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKT...VLQ pdb|6YYT|A
...

"""


seq_id = "pdb|6YYT|A"
id = seq_id.split("|")[1]

print(id)

"""
download the pdb file of the protein
!wget https://files.rcsb.org/download/6YYT.pdb

sam@sam-IdeaPad-L340-15IRH-Gaming:~/Documents/projects$ !wget https://files.rcsb.org/download/6YYT.pdb
wget -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz https://files.rcsb.org/download/6YYT.pdb
--2023-12-23 10:37:12--  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz
           => ‘SRR003265.filt.fastq.gz’
Resolving ftp.1000genomes.ebi.ac.uk (ftp.1000genomes.ebi.ac.uk)... 193.62.193.167
Connecting to ftp.1000genomes.ebi.ac.uk (ftp.1000genomes.ebi.ac.uk)|193.62.193.167|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /vol1/ftp/phase3/data/NA18489/sequence_read ... done.
==> SIZE SRR003265.filt.fastq.gz ... 28919712
==> PASV ... done.    ==> RETR SRR003265.filt.fastq.gz ... done.
Length: 28919712 (28M) (unauthoritative)

SRR003265.filt.fastq.gz                                100%[===========================================================================================================================>]  27.58M  2.06MB/s    in 8.9s    

2023-12-23 10:37:24 (3.09 MB/s) - ‘SRR003265.filt.fastq.gz’ saved [28919712]

--2023-12-23 10:37:24--  https://files.rcsb.org/download/6YYT.pdb
Resolving files.rcsb.org (files.rcsb.org)... 128.6.159.245
Connecting to files.rcsb.org (files.rcsb.org)|128.6.159.245|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: unspecified [application/octet-stream]
Saving to: ‘6YYT.pdb’

6YYT.pdb                                                   [    <=>                                                                                                                     ] 954.91K  1.49MB/s    in 0.6s    

2023-12-23 10:37:25 (1.49 MB/s) - ‘6YYT.pdb’ saved [977832]

FINISHED --2023-12-23 10:37:25--
Total wall clock time: 12s
Downloaded: 2 files, 29M in 9.5s (2.99 MB/s)

"""





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


