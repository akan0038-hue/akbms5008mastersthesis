import subprocess
import os

# Define your mapping: { "filename": residue_index }
# The '0' remains constant as the first chain.
mod_mapping = {
    "pris1.fa": 73,
    "pris2.fa": 72,
    "vir1.fa": 73,
    "vir2.fa": 74,
    "vir3.fa": 73,
    "vir4.fa": 73,
}

fasta_dir = "FASTA"
output_dir = "json_dir/dap"
smiles = "O=C(N[C@H](C(N[C@H](C(N1CCC[C@H]1C(N[C@H](C(N2CCCC[C@H]2C(N[C@H](C(O)=O)c3ccccc3)=O)=O)Cc4ccccc4)=O)=O)CC)=O)CN)c5ccccc5"

os.makedirs(output_dir, exist_ok=True)

for fasta_name, residue_index in mod_mapping.items():
    input_path = os.path.join(fasta_dir, fasta_name)
    base_name = os.path.splitext(fasta_name)[0]
    output_path = os.path.join(output_dir, f"{base_name}.json")
    
    print(f"Processing {fasta_name} at residue {residue_index}...")
    
    cmd = [
        "python", "af3_te_substrate.py",
        "--protein_files", input_path,
        # '0' is the chain, residue_index is the position
        "--te_modification", "0", str(residue_index), smiles, "dt",
        "--max_iterations", "5000",
        "--num_conformers", "200",
        "--save_ccd",
        "--output", output_path
    ]
    
    subprocess.run(cmd, check=True)

print("Batch complete!")
