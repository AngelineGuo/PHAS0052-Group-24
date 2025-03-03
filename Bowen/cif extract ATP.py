import os
from Bio import PDB

# Set paths
input_folder = r"C:\Users\ybw14\Desktop\大三项目\AlphaFold model 0 result cifs"  # Folder containing CIF files
output_folder = r"C:\Users\ybw14\Desktop\大三项目\Alphafold model 0 binding sites 4.0 A"  # Folder to save TXT files

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Set distance threshold (in Å)
DISTANCE_THRESHOLD = 4.0

# Manual Three-Letter to One-Letter Conversion Dictionary
three_to_one_dict = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

# Initialize CIF parser
cif_parser = PDB.MMCIFParser(QUIET=True)

# Loop through all CIF files in the input folder
for filename in os.listdir(input_folder):
    if filename.endswith(".cif"):
        cif_file_path = os.path.join(input_folder, filename)
        output_file_path = os.path.join(output_folder, f"{filename.replace('.cif', '.txt')}")

        try:
            # Parse structure
            structure = cif_parser.get_structure("protein", cif_file_path)

            # Identify ATP molecule (common residue name: ATP)
            ligand_resname = "ATP"

            # Store all protein atoms and ATP atoms
            all_atoms = []
            atp_atoms = []

            for model in structure:
                for chain in model:
                    for residue in chain:
                        # Collect all protein atoms
                        if PDB.is_aa(residue, standard=True):
                            all_atoms.extend(residue.get_atoms())
                        # Collect ATP ligand atoms
                        elif residue.resname == ligand_resname:
                            atp_atoms.extend(residue.get_atoms())

            # Calculate distances between ATP and protein atoms
            binding_residues = set()
            for atp_atom in atp_atoms:
                for protein_atom in all_atoms:
                    distance = atp_atom - protein_atom
                    if distance <= DISTANCE_THRESHOLD:
                        res_name = protein_atom.get_parent().resname
                        res_id = protein_atom.get_parent().id[1]  # Residue sequence number
                        one_letter_code = three_to_one_dict.get(res_name, "X")  # 'X' for unknown residues
                        binding_residues.add((res_id, one_letter_code))

            # Sort binding sites by residue sequence number
            binding_residues = sorted(binding_residues, key=lambda x: x[0])

            # Save results to a TXT file
            with open(output_file_path, "w") as f:
                f.write("Binding Sites\n")
                f.write("=" * 20 + "\n")
                for res_id, aa in binding_residues:
                    f.write(f"{res_id:<6} {aa}\n")

            print(f"Processed: {filename} -> {output_file_path}")

        except Exception as e:
            print(f"Error processing {filename}: {e}")

print("Processing completed.")
