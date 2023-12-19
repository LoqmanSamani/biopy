from Bio import Entrez, SeqIO
from Bio.SeqUtils import molecular_weight as mw
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter


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




