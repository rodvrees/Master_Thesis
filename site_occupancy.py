import pandas as pd
import os
import OxiAnalysis as OA
import warnings; warnings.simplefilter('ignore')
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import urllib
import sys
from Bio.PDB import DSSP
from Bio.PDB import PDBParser
from pypdb import *
from Bio.Align import PairwiseAligner
import numpy as np
import os  
from tqdm import tqdm, tqdm_pandas
import time
tqdm.pandas()

def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    if os.path.isfile(outfnm):
        return outfnm
    else:
        try:
            urllib.request.urlretrieve(url, outfnm)
            return outfnm
        except Exception as err:
            print(str(err), file=sys.stderr)
            return None 

def RSA_calc(row):
    peptide = row["Base Sequence"]
    pos = int(row["position"]) -1
    acc = row["UP_acc"]
    AA = peptide[pos]
    p = PDBParser()
    try:
        PDB_ID = Query(acc).search()[:20]
        for i in PDB_ID:
            try:
                if os.path.isfile('/public/compomics2/Robbe/PDBfiles/{}.pdb'.format(i)):
                    print("PDB file found, getting structure...\n")
                    structure = p.get_structure("{}".format(i), "/public/compomics2/Robbe/PDBfiles/{}.pdb".format(i))
                    model = structure[0]
                    try:
                        print("Getting DSSP\n")
                        dssp = DSSP(model, "/public/compomics2/Robbe/PDBfiles/{}.pdb".format(i))
                        DSSPseq = ""
                        for i in range(len(dssp)):
                            a_key = list(dssp.keys())[i]
                            DSSPseq += dssp[a_key][1]
                        aligner = PairwiseAligner()
                        aligner.mode = "local"
                        aligner.gap_score = -1
                        alignments = aligner.align(DSSPseq, peptide)
                        alignment = alignments[0]
                        if alignment.score == float(len(peptide)):
                            index = alignment.aligned[0][0][0] + pos
                            a_key = list(dssp.keys())[index]
                            if dssp[a_key][1] == AA:
                                RSA = dssp[a_key][3]
                                print("RSA found\n")
                                return RSA
                    except Exception as e:
                        print(e)
                        continue

                else:
                    print("downloading")
                    pdb_file = download_pdb(i, '/public/compomics2/Robbe/PDBfiles/')
                    print("download done, getting structure...")
                    structure = p.get_structure("{}".format(i), "/public/compomics2/Robbe/PDBfiles/{}.pdb".format(i))
                    model = structure[0]
                    try:
                        print("getting DSSP")
                        dssp = DSSP(model, "/public/compomics2/Robbe/PDBfiles/{}.pdb".format(i))
                        DSSPseq = ""
                        for i in range(len(dssp)):
                            a_key = list(dssp.keys())[i]
                            DSSPseq += dssp[a_key][1]
                        aligner = PairwiseAligner()
                        aligner.mode = "local"
                        aligner.gap_score = -1
                        alignments = aligner.align(DSSPseq, peptide)
                        alignment = alignments[0]
                        if alignment.score == float(len(peptide)):
                            index = alignment.aligned[0][0][0] + pos
                            a_key = list(dssp.keys())[index]
                            if dssp[a_key][1] == AA:
                                RSA = dssp[a_key][3]
                                print(RSA)
                                return RSA
                    except Exception as e:
                        print(e)
                        continue
            except FileNotFoundError:
                print("File not found\n")
                continue
        return "no RSA returned\n"
    except TypeError:
        print("TypeError\n")
        return np.nan




quant = pd.read_csv("/home/robbe/ionbot/PXD014381/raw_files/QuantifiedPeptides.tsv", sep="\t")
quant = quant.drop(["Gene Names","Organism", "Intensity_QX01860","Intensity_QX01863","Intensity_QX01866","Intensity_QX01869","Intensity_QX01872","Intensity_QX01862","Intensity_QX01865","Intensity_QX01868","Intensity_QX01984","Intensity_QX01874"], axis=1)
quant = quant.drop(list(quant.filter(regex = "Detection")), axis = 1)
uniprot = pd.read_csv("/home/robbe/ionbot/full_projects_/PXD014381/PXD014381_first.csv")
uniprot = uniprot[["matched_peptide", "proteins"]]

cols = ["Intensity_QX01983","Intensity_QX01981_160316090220","Intensity_QX01867","Intensity_QX01870","Intensity_QX01873"]

OA.quantile_transform(quant,cols)
quantprot = quant[~quant["Protein Groups"].isna()]
quantprot['Counts'] = quantprot.groupby(['Base Sequence'])['Sequence'].transform('count')
quantprot = quantprot.dropna(thresh=7, axis=0)
quantprot['median']=quantprot.filter(like='Intensity').apply(lambda x: x.median(), axis=1)
quantprot["Total_peptide_intensity"] = quantprot.groupby(["Base Sequence"])["median"].transform('sum')
quantprot["Site_occupancy"] = quantprot["median"] / quantprot["Total_peptide_intensity"]
quantprot["modifications"] = quantprot["Sequence"].apply(OA.flashLFQmods)
oxpept = quantprot[quantprot["modifications"].isin(OA.modslist)]

pattern = r"(\d+)\|"
oxpept["position"] = oxpept["Sequence"].apply(lambda x: re.findall(pattern=pattern, string=x)[0])
oxpept["Protein Groups"] = oxpept["Protein Groups"].apply(lambda x: re.sub("sp\|", "", string=x))
data = oxpept.merge(uniprot,left_on="Base Sequence", right_on="matched_peptide", how="left").drop_duplicates()
data = data[data["proteins"].notna()]
data["UP_acc"] = data["proteins"].apply(lambda x: x.split("(")[4].replace(")",""))
data["UP_acc"] = data["UP_acc"].apply(lambda x: re.sub("sp\|","", x))
data["ProtPos"] = data["UP_acc"].astype(str) + "|" +  data["position"].astype(str)


data1 = data.iloc[:1000,:]
data2 = data.iloc[1000:2000,:]
data3 = data.iloc[2000:,:]
datalist = [data1, data2, data3]

for i, df in enumerate(datalist):
    df["RSA"] = df.progress_apply(RSA_calc, axis=1)
    df.to_csv("/home/robbe/ionbot/site_occupancy_data/data{}.csv".format(i+5))
    print("data {} finished".format(i+5))