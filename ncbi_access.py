from Bio import Entrez, SeqIO



# Set your email address. NCBI requires users to provide an email address when accessing their services.
Entrez.email = "samaniloqman91@gmail.com"



# List of the available Entrez databases
handle = Entrez.einfo()

conn = Entrez.read(handle)

handle.close()

print(conn.keys())

"""
dict_keys(['DbList'])
"""

print(conn['DbList'])  # Print the list of available databases

"""
['pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide', 'structure', 'genome', 'annotinfo',
 'assembly', 'bioproject', 'biosample', 'blastdbinfo', 'books', 'cdd', 'clinvar', 'gap',
 'gapplus', 'grasp', 'dbvar', 'gene', 'gds', 'geoprofiles', 'homologene', 'medgen',
 'mesh', 'nlmcatalog', 'omim', 'orgtrack', 'pmc', 'popset', 'proteinclusters', 'pcassay',
 'protfam', 'pccompound', 'pcsubstance', 'seqannot', 'snp', 'sra', 'taxonomy', 'biocollections', 'gtr']

"""



handle1 = Entrez.esearch(db="nucleotide", term="CRT|[Gene Name] AND Plasmodium falciparum[Organism]", retmax="40")

conn1 = Entrez.read(handle1)

print(conn1)
"""
{
 'Count': '8404', 'RetMax': '40', 'RetStart': '0',
 'IdList': ['2587918588', '2576334895', '2576334893', '2576334891', '2576334889', '2576334887', '2576334885',
            '2576334883', '2576334881', '2576334879', '2576334877', '2576334875', '2576334873', '2576334871',
            '2576334869', '2576334867', '2576334865', '2576334863', '2576334861', '2576334859', '2576334857', 
            '2576334855', '2576334853', '2576334851', '2576334849', '2576334847', '2576334845', '2576334843', 
            '2576334841', '2576334839', '2576334837', '2576334835', '2576334833', '2576334831', '2576334829', 
            '2576334827', '2576334825', '2576334823', '2576334821', '2576334819'],
 'TranslationSet': [{'From': 'Plasmodium falciparum[Organism]', 
 'To': '"Plasmodium falciparum"[Organism]'}], 
 'TranslationStack': [{'Term': 'CRT[All Fields]', 'Field': 'All Fields', 'Count': '287347', 'Explode': 'N'},
                      {'Term': 'Gene[All Fields]', 'Field': 'All Fields', 'Count': '220606680', 'Explode': 'N'}, 
                      {'Term': 'Name[All Fields]', 'Field': 'All Fields', 'Count': '171729592', 'Explode': 'N'},
                       'AND', 'GROUP', 'OR', {'Term': '"Plasmodium falciparum"[Organism]', 'Field': 'Organism',
                       'Count': '268657', 'Explode': 'Y'}, 'AND'],
 'QueryTranslation': 'CRT[All Fields] OR (Gene[All Fields] AND Name[All Fields]) AND "Plasmodium falciparum"[Organism]'
 }

"""




handle1.close()

print(conn1["Count"])

""" 8404 """

print(len(conn1["IdList"]))

""" 40 """

# List of the 40 first ids

ids = conn1["IdList"]
print(ids)


"""
['2587918588', '2576334895', '2576334893', '2576334891', '2576334889', '2576334887',
 '2576334885', '2576334883', '2576334881', '2576334879', '2576334877', '2576334875',
 '2576334873', '2576334871', '2576334869', '2576334867', '2576334865', '2576334863', 
 '2576334861', '2576334859', '2576334857', '2576334855', '2576334853', '2576334851', 
 '2576334849', '2576334847', '2576334845', '2576334843', '2576334841', '2576334839', 
 '2576334837', '2576334835', '2576334833', '2576334831', '2576334829', '2576334827', 
 '2576334825', '2576334823', '2576334821', '2576334819']

"""




handle2 = Entrez.efetch(db="nucleotide", id=ids, rettype="gb")  # Genbank format

records = list(SeqIO.parse(handle=handle2, format="gb"))

handle2.close()

print(records)


"""
[SeqRecord(seq=Seq('GGTGGAGGTTCTTGTCTTGGTAAATGTGCTCATGTGTTTAAACTTATTTTTAAA...AAA'), id='OR483864.1', name='OR483864', description='Plasmodium falciparum isolate PE-26 chloroquine resistance transporter (crt) gene, partial cds', dbxrefs=[]), 
 SeqRecord(seq=Seq('TTATTTCTGTAATTTGATACAAAAAGCTATTGATTATAAAAATAAAGGACAAAA...TTA'), id='OR349195.1', name='OR349195', description='Plasmodium falciparum isolate MSC144 multidrug resistance protein 1 gene, partial cds', dbxrefs=[]), 
 SeqRecord(seq=Seq('TTATTTCTGTAATTTGATACAAAAAGCTATTGATTATAAAAATAAAGGACAAAA...TTA'), id='OR349194.1', name='OR349194', description='Plasmodium falciparum isolate MSC143 multidrug resistance protein 1 gene, partial cds', dbxrefs=[]), 
 ...
"""
print(len(records))
""" 40 """

target = []

for rec in records:
    if rec.name == 'OR349171':
        target = rec  # try to find CRT gene
        break





print(target.id)
print(target.description)
print(target.name)
print(target.seq)

"""
OR349171.1
Plasmodium falciparum isolate MSC119 multidrug resistance protein 1 gene, partial cds
OR349171
TTTCTGTAATTTGATACAAAAAGCTATTGATTATAAAAATAAAGGACAAAAAAGAAGAATTATTGTAAATGCAGCTTTATGGGGATTCAGTCAAAGCGCTC
AATTATTTATTAATAGTTTTGCCTATTGGTTTGGATCCTTCTTAATTAAAAGAGGTACTATATTAGTTGATGACTTTATGAAATCCTTATTTACTTTTATA
TTTACTGGTAGTTATGCTGGAAAATTAATGTCCTTAAAAGGAGATTCAGAAAATGCAAAATTATCATTTGAGAAATATTATCCATTAATGATTAGAAAATC
AAATATTGATGTAAGAGATGATGGTGGAATAAGAATAAATAAAAATTTAATAAAAGGTAAAGTTGATATTAAAGATGTAAATTTCCGTTATATTTCAAGAC
CAAATGTACCTATTTATAAAAATTTATCTTTTACATGTGATAGTAAAAAAACTACAGCAATCGTTGGAGAAACAGGTAGTGGAAAATCAACTTTTATGAAT
CTCTTATTAAGATTTTATGACTTGAAAAATGATCACATTATATTAAAAAATGATATGACAAATTTTCAAGATTATCAAAATAATAATAATAATTCATTGGT
TTTAAAAAATGTAAATGAATTTTCAAACCAATCTGGATCTGCAGAAGATTATACTGTATTTAATAATAATGGAGAAATATTATTAGATGATATTAATATAT
GTGATTATAACTTAAGAGATCTTAGAAACTTATTTTCAATAGTTAGTCAAGAACCCATGTTATTTAATATGTCCATATATGAAAATATCAAATTTGGAAGA
GAAGATGCAACATTGGAAGATGTTAAACGTGTTAGTAAGTTTGCTGCTATAGATGAATTTA

"""




