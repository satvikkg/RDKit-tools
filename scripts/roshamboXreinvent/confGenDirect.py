import os
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import multiprocessing

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

def confGEN(
    input_file,
    numConfs=100,
    maxAttempts=1000,
    pruneRmsThresh=0.1,
    clusterMethod="RMSD",
    clusterThreshold=2.0,
    minimizeIterations=0,
    n_jobs=-2,
    output_file=None
):
    """
    Generate conformers for molecules in an SDF file.
    
    Parameters:
    -----------
    input_file : str
        Path to input SDF file
    numConfs : int, default=100
        Number of conformers to generate
    maxAttempts : int, default=1000
        Maximum attempts for conformer generation
    pruneRmsThresh : float, default=0.1
        RMSD threshold for pruning conformers
    clusterMethod : str, default="RMSD"
        Clustering method ("RMSD" or "TFD")
    clusterThreshold : float, default=2.0
        Threshold for clustering
    minimizeIterations : int, default=0
        Number of iterations for energy minimization
    n_jobs : int, default=-2
        Number of parallel jobs (-2: all CPUs - 1, -1: all CPUs, >0: specific number)
    output_file : str, optional
        Output file path. If None, uses input filename with "-RDconfGen.sdf" suffix
    
    Returns:
    --------
    str
        Path to output file
    """
    
    if output_file is None:
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
        for f in as_completed(futures):
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
    return output_file
