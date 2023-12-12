from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML  # the output of blast is in XML format



help(NCBIWWW.qblast)


"""
Help on function qblast in module Bio.Blast.NCBIWWW:

qblast(program, database, sequence, url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi', auto_format=None, 
       composition_based_statistics=None, db_genetic_code=None, endpoints=None, entrez_query='(none)', expect=10.0, 
       filter=None, gapcosts=None, genetic_code=None, hitlist_size=50, i_thresh=None, layout=None, lcase_mask=None, 
       matrix_name=None, nucl_penalty=None, nucl_reward=None, other_advanced=None, perc_ident=None, phi_pattern=None, 
       query_file=None, query_believe_defline=None, query_from=None, query_to=None, searchsp_eff=None, service=None, 
       threshold=None, ungapped_alignment=None, word_size=None, short_query=None, alignments=500, alignment_view=None, 
       descriptions=500, entrez_links_new_window=None, expect_low=None, expect_high=None, format_entrez_query=None, 
       format_object=None, format_type='XML', ncbi_gi=None, results_file=None, show_overview=None, megablast=None, 
       template_type=None, template_length=None, username='blast', password=None)
       
    BLAST search using NCBI's QBLAST server or a cloud service provider.
    
    Supports all parameters of the old qblast API for Put and Get.
    
    Please note that NCBI uses the new Common URL API for BLAST searches
    on the internet (http://ncbi.github.io/blast-cloud/dev/api.html). Thus,
    some of the parameters used by this function are not (or are no longer)
    officially supported by NCBI. Although they are still functioning, this
    may change in the future.
    
    The Common URL API (http://ncbi.github.io/blast-cloud/dev/api.html) allows
    doing BLAST searches on cloud servers. To use this feature, please set
    ``url_base='http://host.my.cloud.service.provider.com/cgi-bin/blast.cgi'``
    and ``format_object='Alignment'``. For more details, please see
    https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=CloudBlast
    
    Some useful parameters:
    
     - program        blastn, blastp, blastx, tblastn, or tblastx (lower case)
     - database       Which database to search against (e.g. "nr").
     - sequence       The sequence to search.
     - ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.
     - descriptions   Number of descriptions to show.  Def 500.
     - alignments     Number of alignments to show.  Def 500.
     - expect         An expect value cutoff.  Def 10.0.
     - matrix_name    Specify an alt. matrix (PAM30, PAM70, BLOSUM80, BLOSUM45).
     - filter         "none" turns off filtering.  Default no filtering
     - format_type    "HTML", "Text", "ASN.1", or "XML".  Def. "XML".
     - entrez_query   Entrez query to limit Blast search
     - hitlist_size   Number of hits to return. Default 50
     - megablast      TRUE/FALSE whether to use MEga BLAST algorithm (blastn only)
     - short_query    TRUE/FALSE whether to adjust the search parameters for a
                      short query sequence. Note that this will override
                      manually set parameters like word size and e value. Turns
                      off when sequence length is > 30 residues. Default: None.
     - service        plain, psi, phi, rpsblast, megablast (lower case)
    
    This function does no checking of the validity of the parameters
    and passes the values to the server as is.  More help is available at:
    https://ncbi.github.io/blast-cloud/dev/api.html
    
"""


seq = "ggtaagtcctctagtacaaacacccccaatattgtgatataattaaaattatattcatattctgttgccagaaaaaacacttttaggctatattagagccatcttctttgaagcgttgtc"



blast = NCBIWWW.qblast(
    program="blastn",
    database="nt",
    sequence=seq
)



result1 = NCBIXML.parse(handle=blast)  # parse the results using BLASTXML function

records1 = list(result1)

print(records1)

""" [<Bio.Blast.Record.Blast object at 0x7fe138732350>] """




threshold1 = 10e-11  # e-value threshold to limit the results

count = 0  # number of hits

results = []  # store the results


for record in records1:
    for align in record.alignments:
        for hsp in align.hsps:
            if hsp.expect < threshold1:
                count += 1
                results.append(record)
                print("****Alignment****")
                print("sequence:", align.title)
                print("length:", align.length)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")
                print()




print(count)


"""
****Alignment****
sequence: gi|1853088208|gb|CP054431.1| Chlamydia trachomatis strain CH2_mutant_L2/434/Bu(i) plasmid unnamed
length: 7676
GGTAAGTCCTCTAGTACAAACACCCCCAATATTGTGATATAATTAAAATTATATTCATATTCTGTTGCCAGAAAA...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGTAAGTCCTCTAGTACAAACACCCCCAATATTGTGATATAATTAAAATTATATTCATATTCTGTTGCCAGAAAA...

****Alignment****
sequence: gi|1853087206|gb|CP054433.1| Chlamydia trachomatis strain CH1_mutant_L2/434/Bu(i) plasmid unnamed
length: 7676
GGTAAGTCCTCTAGTACAAACACCCCCAATATTGTGATATAATTAAAATTATATTCATATTCTGTTGCCAGAAAA...
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
GGTAAGTCCTCTAGTACAAACACCCCCAATATTGTGATATAATTAAAATTATATTCATATTCTGTTGCCAGAAAA...
...

number of sequences: 53
"""




# A function which calculates the same

def blast_result(program="blastn", database="nt", sequence=None, threshold=10e-10):

    alignments = []

    if sequence:

        blast = NCBIWWW.qblast(program=program, database=database, sequence=sequence)
        results = NCBIXML.parse(handle=blast)
        records = list(results)


        for record in records:

            alignment = {}
            for align in record.alignments:
                for hsp in align.hsps:
                    if hsp.expect < threshold:
                        alignment["Title"] = align.title
                        alignment["Length"] = align.length
                        alignment["Query"] = hsp.query
                        alignment["Match"] = hsp.match
                        alignment["Subject"] = hsp.sbjct

            alignments.append(alignment)

    else:
        print("There is no input sequence!!!")

    return alignments


results = blast_result(program="blastn", database="nt", sequence=seq, threshold=10e-11)

result1 = results[0]


print(result1["Title"])

print(result1["Length"])

print(result1["Query"])

print(result1["Match"])

print(result1["Subject"])


"""
gi|2577273931|gb|OR270824.1| Chlamydia trachomatis isolate 1CP_CRYP_F06 hypothetical protein gene, complete cds
315
GTGATATAATTAAAATTATATTCATATTCTGTTGCCAGAAAAAACACTTTTAGGCTATATTAGAGCCATCTTCTTTGAAGCGTTGTC
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
GTGATATAATTAAAATTATATTCATATTCTGTTGCCAGAAAAAACACTTTTAGGCTATATTAGAGCCATCTTCTTTGAAGCGTTGTC

"""





