## Find binding sites of ATP from experimental PDB data
import os
from Bio import PDB

base_path = os.path.abspath(os.path.dirname(__file__))

# Parse the PDB file
pdb_parser = PDB.PDBParser(QUIET=True)
structure = pdb_parser.get_structure("protein", os.path.join(base_path, "../Files/pdb/2j9c_chain_A.pdb"))

# Get the model and ligand (ATP)
model = structure[0]  # Select the first model
ligand_name = "ATP"   # ATP residue name

# Identify ATP molecules
atp_residues = [res for res in model.get_residues() if res.get_resname() == ligand_name]

# Compute interacting amino acids with ATP
ns = PDB.NeighborSearch(list(model.get_atoms()))  # Create a neighbor search object
interaction_cutoff = 4.0  # Set the binding site distance threshold (unit: Å)

interacting_residues = set()

for atp in atp_residues:
    print(f"ATP binding site ({atp}) interacting amino acids:")
    for atom in atp.get_atoms():
        nearby_atoms = ns.search(atom.coord, interaction_cutoff)
        nearby_residues = {atom.get_parent() for atom in nearby_atoms if atom.get_parent().get_id()[0] == " "}  # Filter protein residues
        for res in nearby_residues:
            interacting_residues.add((res.get_resname(), res.get_id()[1]))  # Store unique amino acid name and sequence number

# Sort the residues by sequence number
sorted_residues = sorted(interacting_residues, key=lambda x: x[1])

# Mapping of three-letter amino acid codes to one-letter codes
aa_three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# Output the sorted residues
for resname, resid in sorted_residues:
    one_letter_code = aa_three_to_one.get(resname, '?')  # Get the one-letter code or '?' if not found
    print(f"- {resname} ({one_letter_code}) {resid}")
