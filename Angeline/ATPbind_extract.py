import requests
import os
from Bio import PDB

# ATPbind URL
URL = "https://zhanggroup.org/ATPbind/atpbind.cgi"
EMAIL = "angelineguoatp@gmail.com"  # Replace with your actual email

# PDB downloader
pdbl = PDB.PDBList()

# Function to extract chain from downloaded PDB file
def extract_chain(pdb_id, chain_id):
    pdb_filename = f"{pdb_id}.pdb"
    pdbl.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=".")
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, f"pdb{pdb_id}.ent")
    
    output_filename = f"{pdb_id}_chain_{chain_id}.pdb"
    with open(output_filename, "w") as out_f:
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    io = PDB.PDBIO()
                    io.set_structure(chain)
                    io.save(out_f)
                    return output_filename
    return None

# Submit each extracted chain to ATPbind
for pdb_id, chain in pdb_info:
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