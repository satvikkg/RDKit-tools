import sys, os
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import multiprocessing

'''
Usage: python confGenV2.py test.sdf 100 1000 0.1 RMSD 2.0 0 -2
'''

def gen_conformers(
    mol,
    mol_index,
    numConfs=100,
    maxAttempts=1000,
    pruneRmsThresh=0.1,
    useExpTorsionAnglePrefs=True,
    useBasicKnowledge=True,
    enforceChirality=True
):
    ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=numConfs,
        maxAttempts=maxAttempts,
        pruneRmsThresh=pruneRmsThresh,
        useExpTorsionAnglePrefs=useExpTorsionAnglePrefs,
        useBasicKnowledge=useBasicKnowledge,
        enforceChirality=enforceChirality,
        numThreads=0
    )
    return list(ids)

def calc_energy(mol, conformerId, minimizeIts):
    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conformerId)
    ff.Initialize()
    ff.CalcEnergy()
    results = {}
    if minimizeIts > 0:
        results["converged"] = ff.Minimize(maxIts=minimizeIts)
    results["energy_abs"] = ff.CalcEnergy()
    return results

def cluster_conformers(mol, mode="RMSD", threshold=2.0):
    if mode == "TFD":
        dmat = TorsionFingerprints.GetTFDMatrix(mol)
    else:
        dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
    rms_clusters = Butina.ClusterData(dmat, mol.GetNumConformers(), threshold, isDistData=True, reordering=True)
    return rms_clusters

def align_conformers(mol, clust_ids):
    rmslist = []
    AllChem.AlignMolConformers(mol, confIds=clust_ids, RMSlist=rmslist)
    return rmslist

def process_molecule(args):
    mol, mol_name, numConfs, maxAttempts, pruneRmsThresh, clusterMethod, clusterThreshold, minimizeIterations = args

    if mol is None:
        return [], mol_name

    m = Chem.AddHs(mol)
    conformerIds = gen_conformers(m, 0, numConfs, maxAttempts, pruneRmsThresh)

    conformerPropsDict = {}
    for confId in conformerIds:
        props = calc_energy(m, confId, minimizeIterations)
        conformerPropsDict[confId] = props

    rmsClusters = cluster_conformers(m, clusterMethod, clusterThreshold)

    minEnergy = min([conformerPropsDict[confId]["energy_abs"] for confId in conformerIds]) if conformerIds else 0.0
    results = []

    clusterNumber = 0
    for cluster in rmsClusters:
        clusterNumber += 1
        rmsWithinCluster = align_conformers(m, cluster)
        for idx, confId in enumerate(cluster):
            props = conformerPropsDict[confId]
            props["cluster_no"] = clusterNumber
            props["cluster_centroid"] = cluster[0] + 1
            props["rms_to_centroid"] = rmsWithinCluster[idx-1] if idx > 0 else 0.0
            props["energy_delta"] = props["energy_abs"] - minEnergy
            props["molecule_name"] = mol_name
            props["conformer_id"] = confId + 1
            results.append((m, confId, props))

    return results, mol_name

def confGEN():
    # Parse command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python confGenV2.py input.sdf [numConfs] [maxAttempts] [pruneRmsThresh] [clusterMethod] [clusterThreshold] [minimizeIterations] [n_jobs]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    numConfs = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    maxAttempts = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    pruneRmsThresh = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1
    clusterMethod = sys.argv[5] if len(sys.argv) > 5 else "RMSD"
    clusterThreshold = float(sys.argv[6]) if len(sys.argv) > 6 else 2.0
    minimizeIterations = int(sys.argv[7]) if len(sys.argv) > 7 else 0
    n_jobs = int(sys.argv[8]) if len(sys.argv) > 8 else -2
    
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f"{base_name}-RDconfGen.sdf"
    
    if n_jobs == -2:
        n_jobs_used = max(1, multiprocessing.cpu_count() - 1)
    elif n_jobs == -1:
        n_jobs_used = multiprocessing.cpu_count()
    else:
        n_jobs_used = max(1, n_jobs)

    suppl = Chem.SDMolSupplier(input_file)

    mol_data_list = [
        (mol, mol.GetProp("_Name") if mol and mol.HasProp("_Name") else f"Molecule_{i+1}", 
         numConfs, maxAttempts, pruneRmsThresh, clusterMethod, clusterThreshold, minimizeIterations)
        for i, mol in enumerate(suppl)
    ]

    writer = Chem.SDWriter(output_file)
    lock = multiprocessing.Lock()

    with ProcessPoolExecutor(max_workers=n_jobs_used) as executor:
        futures = {executor.submit(process_molecule, args): i for i, args in enumerate(mol_data_list)}
        for f in tqdm(as_completed(futures), total=len(futures), desc="Processing Molecules"):
            res, mol_name = f.result()
            if res:
                with lock:
                    for mol, confId, props in res:
                        mol.ClearProp("*")
                        mol.SetProp("_Name", props["molecule_name"])
                        for key, val in props.items():
                            if key != "molecule_name":
                                mol.SetProp(str(key), str(val))
                        writer.write(mol, confId=confId)

    writer.close()
    print(f"All molecules written to {output_file}")                        

if __name__ == "__main__":
    confGEN()
