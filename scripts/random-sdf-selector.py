##############################################################################
# Requires rdkit and pandas                                                  #
# This scripts selects a given number of molecules randomly from an sdf file #
##############################################################################

import os
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools as pt

truepath=input("Enter the path for the sdf file: ")
norm=input("Enter the number of random molecules to select: ")
num=int(norm)

base=os.path.basename(truepath)
y=os.path.splitext(base)[0]
x=y+"_random_"+norm+".sdf"

df = pt.LoadSDF(truepath)
randsel=df.sample(n=num)
pt.WriteSDF(randsel,x, properties=list(df.columns))
print(norm+"random selection done from"+truepath)
