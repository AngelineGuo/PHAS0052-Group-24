from Bio import PDB

# Parse the PDB file
pdb_parser = PDB.PDBParser(QUIET=True)
structure = pdb_parser.get_structure("protein",
                                     r"C:")

# Get the model and ligand (ATP)
model = structure[0]  # Select the first model
ligand_name = "ATP"  # ATP residue name

# Identify ATP molecules
atp_residues = [res for res in model.get_residues() if res.get_resname().strip() == ligand_name]  # Remove spaces

# Compute interacting amino acids with ATP
ns = PDB.NeighborSearch(list(model.get_atoms()))  # Create a neighbor search object
interaction_cutoff = 5.0  # Set the binding site distance threshold (unit: Å)

# Custom three-letter to one-letter amino acid mapping
one_letter_aa = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Handle common non-standard amino acids
    "HSD": "H", "HSE": "H", "HSP": "H",  # HIS variants
    "MSE": "M",  # Selenomethionine
    "SEP": "S",  # Phosphoserine
    "TPO": "T",  # Phosphothreonine
    "PTR": "Y"  # Phosphotyrosine
}

for atp in atp_residues:
    interacting_residues = set()  # Use a set to store unique interacting residues

    for atom in atp.get_atoms():
        nearby_atoms = ns.search(atom.coord, interaction_cutoff)
        for nearby_atom in nearby_atoms:
            parent_res = nearby_atom.get_parent()
            if parent_res.get_id()[0] == " ":  # Ensure it's a protein residue
                res_name = parent_res.get_resname().strip()  # Remove spaces
                res_id = parent_res.get_id()[1]

                # Convert three-letter code to one-letter code
                res_one_letter = one_letter_aa.get(res_name, "?")  # "?" for unknown residues

                interacting_residues.add((res_name, res_one_letter, res_id))  # Store as tuple

    print(f"ATP binding site ({atp}) interacting amino acids:")
    for res_name, res_one_letter, res_id in sorted(interacting_residues, key=lambda x: x[2]):  # Sort by residue ID
        print(f"- {res_name} ({res_one_letter}) {res_id}")
