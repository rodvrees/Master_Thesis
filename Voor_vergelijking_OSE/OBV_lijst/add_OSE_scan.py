import pandas as pd
import os
from numpy import nan
from tqdm import tqdm
tqdm.pandas()


scanstocheck = pd.read_csv("scanstocheck.csv")
ionbotPXD002516 = pd.read_csv("ionbotPXD002516.csv")
ionbotPXD004010 = pd.read_csv("ionbotPXD0040106.csv")
fragPXD002516 = pd.read_csv("FragpipePXD002516.tsv", sep="\t")
fragPXD004010 = pd.read_csv("FragpipePXD004010.tsv", sep="\t")
pfindPXD002516 = pd.read_csv("pFindPXD002516.csv", sep="\t")
fragPXD002516["Spectrum_file"] = fragPXD002516["Spectrum"].apply(lambda x: x.split(".")[0])
fragPXD004010["Spectrum_file"] = fragPXD004010["Spectrum"].apply(lambda x: x.split(".")[0])
pfindPXD002516["Spectrum_file"] = pfindPXD002516["File_Name"].apply(lambda x: x.split(".")[0])
fragPXD002516["Scan"] = fragPXD002516["Spectrum"].apply(lambda x: x.split(".")[1])
fragPXD004010["Scan"] = fragPXD004010["Spectrum"].apply(lambda x: x.split(".")[1])

def ionbot_psm(row):
    spectrum = row["Spectrum_file"]
    spectrummgf = spectrum + ".mgf"
    scan = row["scannumber"]

    if row["Project"] == "PXD002516":
        ionbotscanindex = ionbotPXD002516[(ionbotPXD002516["spectrum_file"] == spectrummgf) & (ionbotPXD002516["scan"] == scan)].index
        try:
            peptide = ionbotPXD002516._get_value(ionbotscanindex[0],"matched_peptide")
            modification = ionbotPXD002516._get_value(ionbotscanindex[0], "modifications")
            unexpected = ionbotPXD002516._get_value(ionbotscanindex[0], "unexpected_modification")
        except IndexError:
            return nan

    elif row["Project"] == "PXD004010":
        ionbotscanindex = ionbotPXD004010[(ionbotPXD004010["spectrum_file"] == spectrummgf) & (ionbotPXD004010["scan"] == scan)].index
        try:
            peptide = ionbotPXD004010._get_value(ionbotscanindex[0], "matched_peptide")
            modification = ionbotPXD004010._get_value(ionbotscanindex[0], "modifications")
            unexpected = ionbotPXD004010._get_value(ionbotscanindex[0], "unexpected_modification")
        except IndexError:
            return nan

    psm = peptide + "/" + str(modification) + "/" + str(unexpected)
    return psm

def fragger_psm(row):
    spectrum = row["Spectrum_file"]
    scan = str(row["scannumber"]).zfill(5)
    if row["Project"] == "PXD002516":
        fragscanindex = fragPXD002516[(fragPXD002516["Spectrum_file"] == spectrum) & (fragPXD002516["Scan"] == scan)].index
        try:
            peptide = fragPXD002516._get_value(fragscanindex[0], "Peptide")
            assmod = fragPXD002516._get_value(fragscanindex[0], "Assigned Modifications")
            obsmod = fragPXD002516._get_value(fragscanindex[0], "Observed Modifications")
        except IndexError:
            return nan

    elif row["Project"] == "PXD004010":
        fragscanindex = fragPXD004010[(fragPXD004010["Spectrum_file"] == spectrum) & (fragPXD004010["Scan"] == scan)].index
        try:
            peptide = fragPXD004010._get_value(fragscanindex[0], "Peptide")
            assmod = fragPXD004010._get_value(fragscanindex[0], "Assigned Modifications")
            obsmod = fragPXD004010._get_value(fragscanindex[0], "Observed Modifications")
        except IndexError:
            return nan
    
    psm = peptide + '/' + str(assmod) + "/" + str(obsmod)
    return psm

def pfind_psm(row):
    spectrum = row['Spectrum_file']
    scan = row["scannumber"]

    if row["Project"] == "PXD002516":
        pfindscanindex = pfindPXD002516[(pfindPXD002516["Spectrum_file"] == spectrum) & (pfindPXD002516["Scan_No"] == scan)].index
        try:
            peptide = pfindPXD002516._get_value(pfindscanindex[0], "Sequence")
            mod = pfindPXD002516._get_value(pfindscanindex[0], "Modification")
        except IndexError:
            return nan
        psm = peptide + "/" + str(mod)
    else:
        return nan
    return psm

scanstocheck["ionbot_PSM"] = scanstocheck.progress_apply(ionbot_psm, axis=1)
scanstocheck["fragger_PSM"] = scanstocheck.progress_apply(fragger_psm, axis=1)
scanstocheck["pfind_PSM"] = scanstocheck.progress_apply(pfind_psm, axis=1)

scanstocheck.to_csv("psmstoscanwithpfind.csv", sep="\t")