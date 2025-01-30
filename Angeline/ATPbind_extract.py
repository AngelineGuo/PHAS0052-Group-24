import requests
import os
import csv
from Bio import PDB

### website requests ###

atpbind = 'https://seq2fun.dcmb.med.umich.edu//ATPbind/'
pdb_ext = []
with open('PHAS0052 PDB File - Sheet1.csv') as csvDataFile:
    read = csv.reader(csvDataFile)
    for i in read:
        pdb_ext.append(i[0])

pdb_info = [(entry[:-1].lower(), entry[-1]) for entry in pdb_ext]

print(pdb_info)

# ATPbind submission URL
URL = "https://zhanggroup.org/ATPbind/atpbind.cgi"
EMAIL = "your_email@example.com"  # Replace with your actual email

# Folder to store PDB files
PDB_DIR = "PDB_files"
os.makedirs(PDB_DIR, exist_ok=True)  # Create folder if it doesn't exist

# PDB downloader
pdbl = PDB.PDBList()

# Function to extract chain from downloaded PDB file
def extract_chain(pdb_id, chain_id):
    pdb_path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    
    # Download PDB file
    pdbl.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=PDB_DIR)
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, os.path.join(PDB_DIR, f"pdb{pdb_id}.ent"))

    output_filename = os.path.join(PDB_DIR, f"{pdb_id}_chain_{chain_id}.pdb")
    
    with open(output_filename, "w") as out_f:
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    io = PDB.PDBIO()
                    io.set_structure(chain)
                    io.save(out_f)
                    return output_filename
    return None

# Function to check if ATPbind allows new submissions
def can_submit_job():
    response = requests.post(URL, data={"REPLY-E-MAIL": EMAIL, "TARGET-NAME": "test_job"})
    if "10 jobs already submitted" in response.text:
        print("🚫 ATPbind job limit reached (10 active jobs). Stopping submissions.")
        return False  # Stop submitting
    return True  # Submission allowed

# Submit each extracted chain to ATPbind
for pdb_id, chain in pdb_info:
    if not can_submit_job():
        break  # Stop submitting when the limit is reached

    extracted_chain_file = extract_chain(pdb_id, chain)

    if not extracted_chain_file:
        print(f"❌ Failed to extract Chain {chain} from {pdb_id}, skipping...")
        continue

    with open(extracted_chain_file, "rb") as f:
        response = requests.post(
            URL,
            files={"str_file": (extracted_chain_file, f, "chemical/x-pdb")},
            data={"REPLY-E-MAIL": EMAIL, "TARGET-NAME": pdb_id},
        )

    if response.status_code == 200:
        print(f"✅ {extracted_chain_file} submitted successfully!")
    else:
        print(f"❌ Failed to submit {extracted_chain_file} (HTTP {response.status_code})")