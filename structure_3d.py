from Bio import PDB as pdb
import os
import nglview as nv
import numpy as np


# Download a 3D structure from PDB

pdb1 = pdb.PDBList()

print(pdb1.retrieve_pdb_file("4xp1"))
"""  
Structure exists: '/home/sam/Documents/projects/biopy/xp/4xp1.cif' 
/home/sam/Documents/projects/biopy/files/4xp1.cif
"""

path = "/home/sam/Documents/projects/biopy/files/4xp1.cif"  # path to the downloaded file

print(os.listdir("files"))
""" 
['samp2.fastq', 'samp4.fastq', 'longest.fasta', 'samp1.fasta', 
'samp3.fastq', '4xp1.cif', 'ls_orchid.fasta', 'sample.txt', 
'6YYT.pdb', 'samp5.fastq', 'samp6.fasta']
"""

parser = pdb.MMCIFParser()

structure = parser.get_structure("4xp1", path)
print(structure)


# this function print the entire functions in an object
def clean_dir(obj):
    print([",".join([o for o in dir(obj) if not o.startswith("_")])])


clean_dir(structure)

"""
['add,atom_to_internal_coordinates,center_of_mass,child_dict,child_list,copy,detach_child,
detach_parent,full_id,get_atoms,get_chains,get_full_id,get_id,get_iterator,get_level,
get_list,get_models,get_parent,get_residues,has_id,header,id,insert,internal_to_atom_coordinates,
level,parent,set_parent,transform,xtra']
"""

view = nv.show_biopython(structure)
view.clear_representations()
view.add_ball_and_stick()

# print(view)  # 3D ball and stick structure


view.clear_representations()

view.add_cartoon("protein")

view.add_ball_and_stick("not protein")

# print(view)


for model in structure:  # how many models are there in this structure
    print(f"model {model}")

""" model <Model id=0> """

model = structure[0]  # since we have just one model
for chain in model:
    print(f"chain: {chain}, chain_id: {chain.id}")

"""
chain: <Chain id=A>, chain_id: A
chain: <Chain id=L>, chain_id: L
chain: <Chain id=H>, chain_id: H
chain: <Chain id=B>, chain_id: B
chain: <Chain id=C>, chain_id: C
"""

chain_a = model["A"]

for res in chain_a:
    print(f"Residue Name: {res.resname}, Number: {res.id[1]}")

"""
Residue Name: ASP, Number: 25
Residue Name: GLU, Number: 26
Residue Name: ARG, Number: 27
Residue Name: GLU, Number: 28
Residue Name: THR, Number: 29
...
Residue Name: HOH, Number: 814
Residue Name: HOH, Number: 815
"""


res66 = chain_a[66]

for atom in res66:
    print(atom.name)
"""
N
CA
C
O
CB
CG
CD1
CD2
"""

atom_count = 0

# print every atom in every chain
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print(atom)
                atom_count += 1

"""
<Atom N>
<Atom CA>
<Atom C>
<Atom O>
<Atom CB>
<Atom CG>
<Atom OD1>
<Atom OD2>
<Atom N>
<Atom CA>
...
<Atom O4>
<Atom O5>
<Atom O6>
"""


print(atom_count)
""" 7675 """


struc_dict = pdb.MMCIF2Dict.MMCIF2Dict(path)

print(struc_dict.keys())
"""
dict_keys(['data_', '_entry.id', '_audit_conform.dict_name', '_audit_conform.dict_version', 
'_audit_conform.dict_location', '_database_2.database_id', '_database_2.database_code', 
'_pdbx_database_related.db_name', '_pdbx_database_related.details', '_pdbx_database_related.db_id', 
'_pdbx_database_related.content_type', '_pdbx_database_status.status_code', '_pdbx_database_status.status_code_sf', 
... 

"""

print(struc_dict["_audit_conform.dict_name"])


site_ID = struc_dict['_struct_site_gen.site_id']
site_chain = struc_dict['_struct_site_gen.auth_asym_id']
site_resnum = struc_dict['_struct_site_gen.auth_seq_id']
site_resname = struc_dict['_struct_site_gen.label_comp_id']

cif_binding_residues = []
for bind_id, ch, num, name in zip(site_ID, site_chain, site_resnum, site_resname):
    if bind_id == "AC7":
        print(bind_id, ch, num, name)
        try:
            cif_binding_residues.append(structure[0][ch][int(num)])
        except:
            continue
    else:
        continue




# LDP residues

LDP = None
for res in structure[0].get_residues():
    if res.resname == "LDP":
        LDP = res
        break

print(LDP)
"""
<Residue LDP het=H_LDP resseq=708 icode= >

"""


res_1_ca = structure[0]["A"][56]["CA"]
print(res_1_ca)

print(res_1_ca.coord)  # coordinates of an atom (carbon alpha of the 56's residue)
""" [  3.765 -10.389 -31.972] """


res_2_ca = structure[0]["A"][327]["CA"]
print(res_2_ca.coord)
""" [-12.762   9.25  -31.894] """


diff = res_2_ca.coord - res_1_ca.coord  # distance between the two residues

print(diff)
""" [-16.527       19.639        0.07800102] """


dist = np.sqrt(diff * diff)  # absolute values of the distances
print(dist)
""" [16.527      19.639       0.07800102] """



def near_res(protein, cutoff):

    bind_res = []

    for res in protein[0].get_residues():

        if res == LDP:
            continue

        elif res.id[0].startswith("H"):
            continue

        else:
            alpha_carbon = res['CA']
            distances = []
            for atom in LDP:

                diff_vector = alpha_carbon.coord - atom.coord

                distances.append(np.sqrt(np.sum(diff_vector * diff_vector)))

            if min(distances) < cutoff:
                bind_res.append(res)

    return cutoff, bind_res


cutoff, bind_res = near_res(structure, 10)

print(cutoff)
print(bind_res)

"""
[<Residue GLY het=  resseq=42 icode= >, <Residue PHE het=  resseq=43 icode= >, <Residue ALA het=  resseq=44 icode= >, 
<Residue VAL het=  resseq=45 icode= >, <Residue ASP het=  resseq=46 icode= >, <Residue LEU het=  resseq=47 icode= >, 
<Residue ALA het=  resseq=48 icode= >, <Residue ASN het=  resseq=49 icode= >, <Residue VAL het=  resseq=113 icode= >, 
<Residue VAL het=  resseq=114 icode= >, <Residue LEU het=  resseq=115 icode= >, <Residue ILE het=  resseq=116 icode= >, 
<Residue ALA het=  resseq=117 icode= >, <Residue PHE het=  resseq=118 icode= >, <Residue TYR het=  resseq=119 icode= >, 
<Residue VAL het=  resseq=120 icode= >, <Residue ASP het=  resseq=121 icode= >, <Residue PHE het=  resseq=122 icode= >, 
<Residue TYR het=  resseq=123 icode= >, <Residue TYR het=  resseq=124 icode= >, <Residue ASN het=  resseq=125 icode= >, 
<Residue VAL het=  resseq=126 icode= >, <Residue ILE het=  resseq=127 icode= >, <Residue ILE het=  resseq=128 icode= >,
<Residue CYS het=  resseq=250 icode= >, <Residue PHE het=  resseq=319 icode= >, <Residue SER het=  resseq=320 icode= >, 
...
"""


view = nv.show_biopython(structure)

residues = structure[0].get_residues()
colors = ['0x0000FF' if r not in bind_res else '0xFF0000' for r in residues]

view._set_color_by_residue(colors, component_index=0, repr_index=0)

print(view)


view = nv.show_biopython(structure)


residues = structure[0].get_residues()

colors = ['0x0000FF' if r not in cif_binding_residues else '0xFF0000' for r in residues]
view._set_color_by_residue(colors, component_index=0, repr_index=0)

print(view)