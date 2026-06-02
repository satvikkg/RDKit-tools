#!/usr/bin/env python3
import os, sys, json, time
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools as pt
from roshambo2 import Roshambo2
from confGenV2 import confGEN
from confGenDirect import confGEN

'''
WARNING: Works only with Roshambo 2
'''

query_sdf = "/home/cadd/PROJECTS/GENAI/rein4/roshambo2-based/development-scripts/deucravacitinib.sdf"
confGEN_path = "/home/cadd/PROJECTS/GENAI/rein4/roshambo2-based/development-scripts"

# Get the working directory from the reinvent4 command line argument
if len(sys.argv) < 2:
    sys.stderr.write("Error: Working Directory not supplied by REINVENT4. Please check and set the working directory in the .toml file\n")
    sys.exit(1)

wd = sys.argv[1]
step_count = 1
step_dir = f"{wd}/step{step_count}"
while os.path.exists(step_dir):
    step_count += 1
    step_dir = f"{wd}/step{step_count}"
os.mkdir(step_dir)
os.chdir(step_dir)

# Read SMILES strings supplied from REINVENT4
smiles_list = [line.strip() for line in sys.stdin]
smiles_file = f"step{step_count}.sdf"
smiles_dict = {}

def smiles_to_2d_mol(smiles, mol_name=None):
    """Convert SMILES to RDKit molecule with 2D coordinates"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.stderr.write(f"Warning: Could not parse SMILES: {smiles}\n")
        return None
    Chem.AllChem.Compute2DCoords(mol)
    if mol_name:
        mol.SetProp("_Name", mol_name)
    return mol

writer = Chem.SDWriter(smiles_file)
for smiles_count, smiles in enumerate(smiles_list, start=1):
    smiles_id = f"ReinventGen{step_count}_{smiles_count}"
    smiles_dict[smiles_id] = smiles
    mol = smiles_to_2d_mol(smiles, smiles_id)
    if mol is not None:
        writer.write(mol)
writer.close()

def suppress_fd_stdout(step_dir):
    """Redirect OS-level stdout (fd 1) to a log file to capture all library output."""
    saved_fd = os.dup(1)
    log_path = f"{step_dir}/suppressed_stdout.log"
    log_fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
    os.dup2(log_fd, 1)
    os.close(log_fd)
    sys.stdout = open(log_path, 'a')
    return saved_fd

# Suppress ALL stdout (Python + C-level) for the entire library processing block
saved_fd = suppress_fd_stdout(step_dir)
try:
    confGenOutFile = f"{os.path.splitext(os.path.basename(smiles_file))[0]}-RDconfGen.sdf"
    confGEN(
        input_file=smiles_file,
        numConfs=10,
        maxAttempts=1000,
        pruneRmsThresh=0.1,
        clusterMethod="RMSD",
        clusterThreshold=2.0,
        minimizeIterations=0,
        n_jobs=-2,
        output_file=confGenOutFile
    )

    roshambo2_calculator = Roshambo2(query_sdf, confGenOutFile, color=True, verbosity=1,
                                     remove_Hs_before_color_assignment=False)
    scores = roshambo2_calculator.compute(
        backend='cuda', optim_mode='combination',
        combination_param=0.5, reduce_over_conformers=True,
        tanimoto_threshold=0.1, max_results=1000000, start_mode=1
    )
    for key in scores:
        df = scores[key]
    roshambo2_calculator.write_best_fit_structures()

    # Keep ALL RDKit post-processing inside suppression block
    load_query_sdf = Chem.SDMolSupplier(query_sdf)[0]
    query_molname = load_query_sdf.GetProp("_Name")
    roshambo_out_sdf = f"hits_for_query_{query_molname}.sdf"
    roshambo_out = pt.LoadSDF(roshambo_out_sdf, removeHs=False)
    roshambo_out['score'] = roshambo_out['score'].astype('float')
    filtered = roshambo_out[roshambo_out['name'].isin(smiles_dict.keys())]
    score_dict = dict(zip(filtered['name'], filtered['score']))
    score_dict = {id: score_dict.get(id, 0) for id in smiles_dict.keys()}
    roshambo_scores = list(score_dict.values())

finally:
    # Flush and close the log-redirected stdout
    sys.stdout.flush()
    sys.stdout.close()

    # Build and write JSON directly to saved_fd (the original real stdout)
    # fd 1 stays pointed at the log the entire time — no restore, no race condition
    output = {
        "version": 1,
        "payload": {
            "predictions": roshambo_scores
        }
    }
    json_bytes = (json.dumps(output) + "\n").encode()
    # os.write(saved_fd, json_bytes)
    # os.close(saved_fd)
