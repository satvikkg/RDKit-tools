################################################################################################
### This script can be used to calculate the synthetic accessibility of compounds.           ###
### A positive score indicates the compound is ES(Easy to synthesize) whereas a negative     ###
### score indicates that a molecule is HS(Hard to synthesize).                               ###
###                                                                                          ###
### DOI:10.1186/s13321-017-0206-2                                                            ###
### Code reference: https://github.com/lich-uct/syba                                         ###
### Installation: conda create -n syba_ra (Create new environment)                           ###
###               conda install -c rdkit -c lich                                             ###
###                                                                                          ###
################################################################################################

import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools as pt
from syba.syba import SybaClassifier

syba=SybaClassifier()
syba.fitDefaultScore()

#Loading csv file 
df=pt.LoadSDF("dcm_random_200.sdf")
df.shape

#Calculating SYBA synthetic accessibiity score
df['sa_score']=[syba.predict(mol=mol) for mol in df['ROMol']]

pt.WriteSDF(df, "dcm_random_200_SA.sdf", properties=list(df.columns))
print("Synthetic score calculations complete")
