import os
import gzip
import csv
from Bio.PDB import PDBList
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.Data.IUPACData import protein_letters_3to1

# Use the HTTPS server to avoid old FTP mirror issues
pdb_list = PDBList(server="https://files.rcsb.org/")

# Define current Dir as base path
base_path = os.path.abspath(os.path.dirname(__file__))

# Define output directory for files
files_output_dir = os.path.join(base_path, "..", "Files")
os.makedirs(files_output_dir, exist_ok=True)

# Define output directories for fasta and cif format files

fasta_output_dir = os.path.join(files_output_dir, "fasta")
os.makedirs(fasta_output_dir, exist_ok=True)

cif_output_dir = os.path.join(files_output_dir, "cif")
os.makedirs(cif_output_dir, exist_ok=True)

# Read representatives from a CSV file
csv_file = os.path.join(base_path, "..", "Files", "representatives.csv") # Change to complete file path later
representatives = []

try:
    with open(csv_file, "r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            pdb_chain_combined = row["pdb_id"]
            pdb_id = pdb_chain_combined[:4].lower()  # First 4 characters for the PDB ID
            chain_id = pdb_chain_combined[4:].upper()  # Remaining characters for the chain ID
            representatives.append((pdb_id, chain_id))
except FileNotFoundError:
    print(f"Error: CSV file '{csv_file}' not found.")
    exit(1)
except KeyError as e:
    print(f"Error: Missing expected column in CSV file - {e}")
    exit(1)

# Process each PDB ID and chain
for pdb_id, chain in representatives:
    try:
        # Download mmCIF file
        pdb_file = pdb_list.retrieve_pdb_file(
            pdb_id,
            file_format="mmCif",
            pdir=cif_output_dir,
            overwrite=True
        )
        print(f"Biopython returned path: {pdb_file}")

        # Download PDB file
        # pdb_file_pdb_format = pdb_list.retrieve_pdb_file(
        #     pdb_id,
        #     file_format="pdb",
        #     pdir=pdb_output_dir,
        #     overwrite=True
        # )
        # print(f"Biopython returned path for PDB format: {pdb_file_pdb_format}")

    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        continue

    # Check if the file or gz version exists
    if not os.path.exists(pdb_file):
        if os.path.exists(pdb_file + ".gz"):
            pdb_file += ".gz"
        else:
            print(f"File not found after download: {pdb_file}")
            continue

    # Parse mmCIF
    parser = MMCIFParser(QUIET=True)
    try:
        if pdb_file.endswith(".gz"):
            with gzip.open(pdb_file, "rt") as handle:
                structure = parser.get_structure(pdb_id, handle)
        else:
            structure = parser.get_structure(pdb_id, pdb_file)

        # Extract chain sequence
        found_chain = False
        sequence = []
        for model in structure:
            for chain_obj in model:
                if chain_obj.id == chain:
                    found_chain = True
                    print(f"Processing chain {chain}")
                    for residue in chain_obj.get_residues():
                        if residue.id[0] == " ":  # Exclude non-amino acids
                            residue_name = residue.resname.strip().capitalize()
                            if residue_name in protein_letters_3to1:
                                sequence.append(protein_letters_3to1[residue_name])
                            else:
                                print(f"Skipping unknown residue: {residue_name}")
                    break

        if not found_chain:
            print(f"Chain {chain} not found in PDB {pdb_id}")
            continue

        # Save sequence as FASTA
        if sequence:
            fasta_sequence = "".join(sequence)
            fasta_file = os.path.join(fasta_output_dir, f"{pdb_id}_{chain}.fasta")
            record = SeqRecord(Seq(fasta_sequence), id=f"{pdb_id}_{chain}", description="")
            with open(fasta_file, "w") as out_fasta:
                FastaIO.FastaWriter(out_fasta).write_record(record)
            print(f"FASTA file saved: {fasta_file}")
        else:
            print(f"No sequence found for {pdb_id}_{chain}. Skipping.")

    except Exception as e:
        print(f"Error processing {pdb_id}_{chain}: {e}")