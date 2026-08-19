import os
import subprocess

# --- 1. DEFINE PROTEIN FASTA FILES ---
RAMO_FASTA = "FASTA/ramo_te199.fasta"
ENDURA_FASTA = "FASTA/endura_te_199a.fasta"
CHERS_FASTA = "FASTA/che_nohisortag.fa"

# --- 2. DEFINE SUBSTRATE MODIFICATIONS (SMILES & AF3 Res-Code) ---
RAMO_SMILES = "O[C@H](C(N)=O)[C@H](NC([C@@H](NC(/C=C\\C=C\\CC(C)C)=O)CC(N)=O)=O)C(N[C@H](C1=CC=C(O)C=C1)C(N[C@H](CCCN)C(N[C@H]([C@H](O)C)C(N[C@@H](C2=CC=C(O)C=C2)C(N[C@H](C3=CC=C(O)C=C3)C(N[C@@H]([C@H](C)O)C(N[C@@H](CC4=CC=CC=C4)C(N[C@H](CCCN)C(N[C@@H](C5=CC=C(O)C=C5)C(N[C@H]([C@@H](C)O)C(N[C@@H](C6=CC=C(O)C=C6)C(NCC(N[C@@H](CC(C)C)C(N[C@H](C)C(N[C@@H](C7=CC=C(O)C(Cl)=C7)C(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O"
RAMO_CODE = "ramo"

ENDURA_SMILES = "O[C@H](C)[C@H](NC([C@@H](NC(/C=C\\C=C\\CCCCC(C)C)=O)CC(O)=O)=O)C(N[C@H](C1=CC=C(O)C=C1)C(N[C@H](CCCN)C(N[C@H]([C@@H](O)C)C(N[C@@H](C2=CC=C(O)C=C2)C(N[C@H](C3=CC=C(O)C=C3)C(N[C@@H]([C@H](C)O)C(N[C@@H](CCCNC(N)=O)C(N[C@H](CC4NC(NC4)=N)C(N[C@@H](C5=CC=C(O)C=C5)C(N[C@H](CO)C(N[C@@H](C6=CC(Cl)=C(O)C(Cl)=C6)C(NCC(N[C@@H](CC7NC(NC7)=N)C(N[C@H](C)C(N[C@@H](C8=CC=C(O)C=C8)C(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O"
ENDURA_CODE = "endu"

CHERS_SMILES = "O[C@H](C)[C@H](NC([C@@H](NC(CCCCCC(C)C)=O)CC(N)=O)=O)C(N[C@H](C1=CC=C(O)C=C1)C(N[C@H](CCCN)C(N[C@H]([C@H](O)C)C(N[C@@H](C2=CC=C(O)C=C2)C(N[C@H](C3=CC=C(O)C=C3)C(N[C@@H]([C@H](C)O)C(N[C@@H](CC4=CC=CC=C4)C(N[C@H](CCCN)C(N[C@@H](C5=CC=C(O)C=C5)C(N[C@H]([C@@H](C)O)C(N[C@@H](C6=CC(O)=CC(O)=C6)C(NCC(N[C@@H](C(C)C)C(N[C@H](C)C(N[C@@H](C7=CC=C(O)C(Cl)=C7)C(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O"
CHERS_CODE = "chrs"

SYNTH_SMILES = "OC([C@H](C1=CC=C(O)C(Cl)=C1)NC([C@@H](C)NC([C@H](CC(C)C)NC(CNC([C@H](C2=CC=C(O)C=C2)NC([C@@H]([C@H](C)O)NC([C@H](C3=CC=C(O)C=C3)NC([C@@H](CCCN)NC([C@H](CC4=CC=CC=C4)NC([C@@H](NC([C@@H](C5=CC=C(O)C=C5)NC([C@@H](NC([C@@H]([C@@H](O)C)NC([C@@H](CCCN)NC([C@@H](C6=CC=C(O)C=C6)NC([C@@H](NC([C@@H](NC(CCCCCCC)=O)CC(N)=O)=O)CN)=O)=O)=O)=O)C7=CC=C(O)C=C7)=O)=O)[C@@H](C)O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O"
SYNTH_CODE = "synth"


# --- 3. RUN PIPELINE FUNCTION ---
def run_screen(sub_name, smiles, code, dom_name, fasta):
    target_json_dir = f"json_dir/{sub_name}_on_{dom_name}"
    target_out_dir = f"out_dir/{sub_name}_on_{dom_name}"

    print("=" * 73)
    print(f"[LAUNCHING] Substrate: {sub_name} ---> TE Domain: {dom_name}")
    print("=" * 73)

    os.makedirs(target_json_dir, exist_ok=True)

    # Step 1: Generate the JSON
    cmd_step1 = [
        "python",
        "af3_te_substrate.py",
        "--protein_files",
        fasta,
        "--te_modification",
        "0",
        "91",
        smiles,
        code,
        "--num_conformers",
        "5000",
        "--max_iterations",
        "200",
        "--save_ccd",
        "--output",
        f"{target_json_dir}/{sub_name}_on_{dom_name}.json",
    ]
    subprocess.run(cmd_step1, check=True)

    # Step 2: Submit to Pipeline
    cmd_step2 = [
        "python",
        "run_af3_ppant_pipeline.py",
        "--json_dir",
        target_json_dir,
        "--output_dir",
        target_out_dir,
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
    subprocess.run(cmd_step2, check=True)

    print(f"[QUEUED] {sub_name}_on_{dom_name} successfully setup.\n")


# --- 4. EXECUTE ISOLATED CROSS MATRIX ---

if __name__ == "__main__":
    # Ramoplanin Crosses
    run_screen("ramoplanin", RAMO_SMILES, RAMO_CODE, "enduracin", ENDURA_FASTA)
    run_screen(
        "ramoplanin", RAMO_SMILES, RAMO_CODE, "chersinamycin", CHERS_FASTA
    )

    # Enduracin Crosses
    run_screen("enduracin", ENDURA_SMILES, ENDURA_CODE, "ramoplanin", RAMO_FASTA)
    run_screen(
        "enduracin", ENDURA_SMILES, ENDURA_CODE, "chersinamycin", CHERS_FASTA
    )

    # Chersinamycin Crosses
    run_screen(
        "chersinamycin", CHERS_SMILES, CHERS_CODE, "ramoplanin", RAMO_FASTA
    )
    run_screen(
        "chersinamycin", CHERS_SMILES, CHERS_CODE, "enduracin", ENDURA_FASTA
    )

    # Synthetic Crosses
    run_screen("synthetic", SYNTH_SMILES, SYNTH_CODE, "ramoplanin", RAMO_FASTA)
    run_screen("synthetic", SYNTH_SMILES, SYNTH_CODE, "enduracin", ENDURA_FASTA)
    run_screen(
        "synthetic", SYNTH_SMILES, SYNTH_CODE, "chersinamycin", CHERS_FASTA
    )
