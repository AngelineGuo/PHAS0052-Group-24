from Bio import PDB

# Parse the PDB file
pdb_parser = PDB.PDBParser(QUIET=True)
structure = pdb_parser.get_structure("protein", r"C:\Users\YourName\Documents\example.pdb")

# Get the model and ligand (ATP)
model = structure[0]  # Select the first model
ligand_name = "ATP"   # ATP residue name

# Identify ATP molecules
atp_residues = [res for res in model.get_residues() if res.get_resname() == ligand_name]

# Compute interacting amino acids with ATP
ns = PDB.NeighborSearch(list(model.get_atoms()))  # Create a neighbor search object
interaction_cutoff = 5.0  # Set the binding site distance threshold (unit: Å)

for atp in atp_residues:
    print(f"ATP binding site ({atp}) interacting amino acids:")
    for atom in atp.get_atoms():
        nearby_atoms = ns.search(atom.coord, interaction_cutoff)
        nearby_residues = {atom.get_parent() for atom in nearby_atoms if atom.get_parent().get_id()[0] == " "}  # Filter protein residues
        for res in nearby_residues:
            print(f"- {res.get_resname()} {res.get_id()[1]}")  # Output amino acid name and sequence number
