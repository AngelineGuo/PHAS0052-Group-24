import os
import requests
import csv
from Bio.PDB import PDBList
from Bio.Data.IUPACData import protein_letters_3to1

# Define current Dir as base path
base_path = os.path.abspath(os.path.dirname(__file__))

# Define output directory for files
pdb_output_dir = os.path.join(base_path, "..", "Files", "pdb")
os.makedirs(pdb_output_dir, exist_ok=True)

# Read representatives from a CSV file
csv_file = os.path.join(base_path, "..", "Files", "all_proteins.csv")
representatives = []

# Extract the pdb ids and chain ids from the csv file
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
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)

    if response.status_code == 200:
        file_path = os.path.join(pdb_output_dir, f"{pdb_id}_{chain}.pdb")
        with open(file_path, "wb") as file:
            file.write(response.content)
        print(f"File downloaded and saved as {pdb_id}_{chain}.pdb")
    else:
        print(f"Failed to download file. Status code: {response.status_code}")