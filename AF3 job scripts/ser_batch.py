import os
import subprocess

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

fasta_dir = "FASTA/pris"
output_dir = "json_dir/ser"
pipeline_out_dir = "out_dir/ser"
smiles = "O=C(N[C@H](C(N[C@H](C(N1CCC[C@H]1C(N([C@H](C(N2CCCC[C@H]2C(N[C@H](C(O)=O)c3ccccc3)=O)=O)Cc4ccccc4)C)=O)=O)CC)=O)CO)c5ccccc5"

# 1. Create output directories using Python
os.makedirs(output_dir, exist_ok=True)
os.makedirs(pipeline_out_dir, exist_ok=True)

# 2. Run substrate modification batch loop
for fasta_name, residue_index in mod_mapping.items():
    input_path = os.path.join(fasta_dir, fasta_name)
    base_name = os.path.splitext(fasta_name)[0]
    output_path = os.path.join(output_dir, f"{base_name}.json")

    print(f"Processing {fasta_name} at residue {residue_index}...")

    cmd = [
        "python",
        "af3_te_substrate.py",
        "--protein_files",
        input_path,
        "--te_modification",
        "0",
        str(residue_index),
        smiles,
        "ser",
        "--max_iterations",
        "5000",
        "--num_conformers",
        "200",
        "--save_ccd",
        "--output",
        output_path,
    ]

    subprocess.run(cmd, check=True)

print("Batch JSON generation complete!")

# 3. Execute the pipeline command via subprocess
print("Running AF3 PPANT pipeline...")

pipeline_cmd = [
    "python",
    "run_af3_ppant_pipeline.py",
    "--json_dir",
    output_dir,
    "--output_dir",
    pipeline_out_dir,
    "--num_seeds",
    "25",
    "--random_seeds",
    "--partition",
    "gpu",
    "--time",
    "0-05:00:00",
    "--memory",
    "154G",
]

subprocess.run(pipeline_cmd, check=True)

print("Pipeline execution complete!")
