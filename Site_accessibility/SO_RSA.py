import pandas as pd
import os
import OxiAnalysis as OA
import warnings; warnings.simplefilter('ignore')
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
quantfilepath = sys.argv[1]
print(quantfilepath)
# quant = pd.read_csv("/home/robbe/ionbot/ionbot_0.9.5/PXD022545/mzml_files/QuantifiedPeptides.tsv", sep="\t")
# quant = quant.drop(quant.filter(regex=r"Detection|Gene|Organism|H2O2"), axis=1)
# uniprot = pd.read_csv("/home/robbe/ionbot/full_projects_/PXD022545/PXD022545_first.csv")
# uniprot = uniprot[["matched_peptide", "proteins"]]