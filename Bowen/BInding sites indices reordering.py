import pandas as pd
from Bio import PDB, pairwise2, SeqIO


# === Step 1: Parse the PDB file and extract ATP binding sites ===
def extract_atp_binding_sites(pdb_file, ligand_name="ATP", cutoff=4.0):  # Set the binding site distance threshold (unit: Å)
    """
    Parse the PDB structure and find residues interacting with ATP.
    Returns:
        - binding_residues: List of (three-letter name, one-letter name, PDB residue ID)
    """
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("protein", pdb_file)

    model = structure[0]  # Select the first model
    ns = PDB.NeighborSearch(list(model.get_atoms()))  # Create a neighbor search object

    # Three-letter to one-letter amino acid mapping
    one_letter_aa = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "HSD": "H", "HSE": "H", "HSP": "H",  # HIS variants
        "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y"  # Non-standard amino acids
    }

    binding_residues = set()

    # Identify ATP molecules
    atp_residues = [res for res in model.get_residues() if res.get_resname().strip() == ligand_name]

    for atp in atp_residues:
        for atom in atp.get_atoms():
            nearby_atoms = ns.search(atom.coord, cutoff)  # Find residues within the cutoff distance
            for nearby_atom in nearby_atoms:
                parent_res = nearby_atom.get_parent()
                if parent_res.get_id()[0] == " ":  # Filter non-protein molecules
                    res_name = parent_res.get_resname().strip()
                    res_id = parent_res.get_id()[1]
                    res_one_letter = one_letter_aa.get(res_name, "?")
                    binding_residues.add((res_name, res_one_letter, res_id))

    return sorted(binding_residues, key=lambda x: x[2])  # Sort by residue ID


# === Step 2: Extract the amino acid sequence and residue numbers from the PDB structure ===
def extract_pdb_sequence(pdb_file):
    """
    Extract the amino acid sequence from the PDB structure.
    Returns:
        - pdb_sequence: Extracted sequence (one-letter)
        - pdb_residue_numbers: List of residue IDs from the PDB structure
    """
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("protein", pdb_file)

    pdb_seq = []
    pdb_residue_numbers = []

    one_letter_aa = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "HSD": "H", "HSE": "H", "HSP": "H", "MSE": "M", "SEP": "S",
        "TPO": "T", "PTR": "Y"
    }

    for model in structure:
        for chain in model:
            for res in chain.get_residues():
                if res.get_id()[0] == " ":  # Filter non-protein residues
                    res_name = res.get_resname().strip()
                    res_id = res.get_id()[1]
                    res_one_letter = one_letter_aa.get(res_name, "?")
                    pdb_seq.append(res_one_letter)
                    pdb_residue_numbers.append(res_id)

    return "".join(pdb_seq), pdb_residue_numbers


# === Step 3: Read the FASTA sequence ===
def read_fasta_sequence(fasta_file):
    """
    Read the amino acid sequence from a FASTA file.
    Returns:
        - fasta_sequence: Extracted sequence from the FASTA file
    """
    fasta_seq = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_seq = str(record.seq)
    return fasta_seq


# === Step 4: Map PDB residue IDs to FASTA sequence positions ===
def map_pdb_to_fasta(pdb_sequence, pdb_residue_numbers, fasta_sequence):
    """
    Perform sequence alignment to map PDB residue numbers to FASTA sequence positions.
    Returns:
        - pdb_to_fasta_mapping: {PDB residue ID: FASTA sequence position}
    """
    alignments = pairwise2.align.globalxx(pdb_sequence, fasta_sequence, one_alignment_only=True)

    aligned_pdb_seq, aligned_fasta_seq = alignments[0][:2]

    pdb_to_fasta_mapping = {}
    fasta_index = 0  # Index for FASTA sequence

    for pdb_index, (pdb_res, fasta_res) in enumerate(zip(aligned_pdb_seq, aligned_fasta_seq)):
        if pdb_res == "-" or fasta_res == "-":
            continue
        pdb_to_fasta_mapping[pdb_residue_numbers[pdb_index]] = fasta_index + 1  # FASTA sequence is 1-based
        fasta_index += 1

    return pdb_to_fasta_mapping


# === Step 5: Match ATP binding site residues to their FASTA sequence positions ===
def get_fasta_positions(binding_sites, pdb_to_fasta_mapping):
    """
    Get the FASTA sequence positions for PDB ATP binding site residues.
    Returns:
        - fasta_positions: {PDB residue ID: FASTA sequence position}
    """
    fasta_positions = {}
    for res_id in binding_sites:
        fasta_pos = pdb_to_fasta_mapping.get(res_id, "N/A")
        fasta_positions[res_id] = fasta_pos
    return fasta_positions


# === Main Execution ===

pdb_file = r"C:\Users\ybw14\Desktop\大三项目\extracted_chains\extracted_chains\3f5m_chainA.pdb"
fasta_file = r"C:\Users\ybw14\Desktop\大三项目\BioPython 测试\3f5m_A.fasta"

# Extract ATP binding sites from PDB
binding_residues = extract_atp_binding_sites(pdb_file)

# Extract PDB sequence & residue numbers
pdb_sequence, pdb_residue_numbers = extract_pdb_sequence(pdb_file)

# Read the FASTA sequence
fasta_sequence = read_fasta_sequence(fasta_file)

# Align PDB & FASTA sequences to establish mapping
pdb_to_fasta_mapping = map_pdb_to_fasta(pdb_sequence, pdb_residue_numbers, fasta_sequence)

# Get FASTA sequence positions for ATP binding residues
binding_fasta_positions = get_fasta_positions([res[2] for res in binding_residues], pdb_to_fasta_mapping)

# Print final mapping results
print("PDB Residue -> FASTA Sequence Position Mapping:")
for res_name, res_one_letter, res_id in binding_residues:
    fasta_res = binding_fasta_positions[res_id]
    print(f"{res_name} ({res_one_letter}) PDB: {res_id} → FASTA: {fasta_res}")
