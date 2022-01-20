"""
##############################################################################
# This sript can be used to add the file name as a column named as "source". #
# The path for the directory containing the sdf files should be entered.     #  
# Requires RDKIT installation.                                               #
##############################################################################
"""

import os
import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools as pt

truepath=input("Enter the path for the folder containing sdf files: ")
files=os.listdir(truepath)

for f in files:
   name=truepath+f
   x=pt.LoadSDF(str(name),molColName="Frag ID",includeFingerprints=False)
   x["source"]=str(f)
   pt.WriteSDF(x,name,molColName="Frag ID", properties=list(x.columns))
