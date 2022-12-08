import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pyteomics import pepxml
import os


def parsemods(modifications, peptide):
    moddict = {14.015650: "Methyl", 79.966331: "Phospho", 0.984016: "Deamidated", 203.079373: "HexNAc", 42.010565: "Acetyl", 15.994915: "Oxidation"}
    modslist = []
    PTMlist = []
    for mod in modifications:
        mass = mod["variable"]
        if mass in moddict:
            PTMvariablename = moddict[mass]
            position = mod["position"]
            AA = peptide[position-1]
            formatted = AA+str(position)
            PTMlist.append(PTMvariablename)
            modslist.append(formatted)
        else:
            continue
    if len(modslist) == 1:
        return PTMlist[0], modslist[0]
    else:
        return PTMlist, modslist

PXD004010files = ["/public/compomics2/Robbe/PXD004010/Comet_out/{}".format(f) for f in os.listdir("/public/compomics2/Robbe/PXD004010/Comet_out") if "pep.xml" in f]
PXD004010 = pepxml.chain.from_iterable(PXD004010files)
PXD005216files = ["/public/compomics2/Robbe/PXD005216/Comet_out/{}".format(f) for f in os.listdir("/public/compomics2/Robbe/PXD005216/Comet_out") if "pep.xml" in f]
PXD005216 = pepxml.chain.from_iterable(PXD005216files)

with open("Scannumbers.tsv", "a") as f:
    # print("Starting PXD004010")
    # f.write("Peptide\tPTM\tSite\tspectrum\tscannumber\tProject\n")
    # for i, scan in enumerate(PXD004010):
    #     print("Progress: scan "+str(i))
    #     try:
    #         PSM = scan["search_hit"][0]
    #         peptide = PSM["peptide"]
    #         PTM, Site = parsemods(PSM["modifications"], peptide)
    #         scannumber = scan["spectrumNativeID"].split("=")[3]
    #         spectrum = scan["spectrum"]
    #         f.write("{}\t{}\t{}\t{}\t{}\tPXD004010\n".format(peptide, PTM, Site, spectrum, scannumber))
    #     except KeyError:
    #         pass
    print("Starting PXD002516")
    for i, scan in enumerate(PXD005216):
        print("Progress: scan "+str(i))
        try:
            PSM = scan["search_hit"][0]
            peptide = PSM["peptide"]
            PTM, Site = parsemods(PSM["modifications"], peptide)
            scannumber = scan["start_scan"]
            spectrum = scan["spectrum"]
            f.write("{}\t{}\t{}\t{}\t{}\tPXD002516\n".format(peptide, PTM, Site, spectrum, scannumber))
        except KeyError:
            pass 
    f.close()
        
