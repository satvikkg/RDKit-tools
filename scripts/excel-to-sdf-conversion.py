############################################################################################################
###                                                                                                      ###
###        This script can be used to convert excel sheet containing smile into a 2D sdf file.           ###
###        Requires rdkit and pandas.                                                                    ###
###                                                                                                      ###
############################################################################################################

import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools as pt

#Reading the excel file using pandas.
df=pd.read_excel("path for excel file", sheet_name="enter the sheet name")

#Converting exccel file to sdf file using rdkit.
pt.AddMoleculeColumnToFrame(df, smilesCol='smiles')

#Saving dataframe as sdf file.
pt.WriteSDF(df,"path and file name of output",molColName="ID for the molecule", properties=list(df.columns))
print("Conversion of excel to sdf file complete.")
