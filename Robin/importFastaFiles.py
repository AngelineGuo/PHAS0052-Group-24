import os
import csv
import requests
import re

# Define path of Files directory
files_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "Files")

# Define output directory for fasta files
fasta_output_dir = os.path.join(files_dir, "fasta")
os.makedirs(fasta_output_dir, exist_ok=True)

# Read representatives from a CSV file
csv_file = os.path.join(files_dir, "all_proteins.csv")
representatives = []

try:
    with open(csv_file, "r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            pdb_chain_combined = row["pdb_id"]     # Example: 1a2bA
            pdb_id = pdb_chain_combined[:4].lower()  # First 4 characters for the PDB ID
            chain_id = pdb_chain_combined[4:].upper()  # Remaining characters for the chain ID
            representatives.append((pdb_id, chain_id))
except FileNotFoundError:
    print(f"Error: CSV file '{csv_file}' not found.")
    exit(1)
except KeyError as e:
    print(f"Error: Missing expected column in CSV file - {e}")
    exit(1)

# Download and save FASTA files
for pdb_id, chain in representatives:
    fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    try:
        response = requests.get(fasta_url)
        response.raise_for_status()    # Raise an exception for error codes
        fasta_data = response.text

        # Split the response string into paragraphs
        paragraphs = fasta_data.split(">")

        # Extract the required chain sequence
        chain_sequence = ""
        chain_pattern = re.compile(rf"{pdb_id.upper()}_\d+\|[^|]*\bChain(s)?\b[^|]*\b{chain}\b")
        for paragraph in paragraphs:
            if chain_pattern.search(paragraph):
                chain_sequence = f">{paragraph}"
                break

        if chain_sequence:
            # Save the FASTA file
            fasta_file = os.path.join(fasta_output_dir, f"{pdb_id}_{chain}.fasta")
            with open(fasta_file, "w") as file:
                file.write(chain_sequence)
            print(f"FASTA file saved: {fasta_file}")
        else:
            print(f"Chain {chain} not found in FASTA for {pdb_id}")

    except requests.RequestException as e:
        print(f"Error downloading FASTA for {pdb_id}: {e}")
