import sys
import pandas as pd

def rename_files(string):
    if string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\Control_1\\interact.pep.xml":
        return "QX01983"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\Control_2\\interact.pep.xml":
        return 'QX01981_160316090220'
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\Control_3\\interact.pep.xml":
        return "QX01867"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\Control_4\\interact.pep.xml":
        return "QX01870"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\Control_5\\interact.pep.xml":
        return "QX01873"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\Diamide_1\\interact.pep.xml":
        return "QX01862"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\Diamide_2\\interact.pep.xml":
        return "QX01865"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\Diamide_3\\interact.pep.xml":
        return "QX01868"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\Diamide_4\\interact.pep.xml":
        return "QX01984"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\Diamide_5\\interact.pep.xml":
        return "QX01874"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\RA_1\\interact.pep.xml":
        return "QX01860"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\RA_2\\interact.pep.xml":
        return "QX01863"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\RA_3\\interact.pep.xml":
        return "QX01866"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\RA_4\\interact.pep.xml":
        return "QX01869"
    elif string == "C:\\Users\\robbe\\Fragpipe_results\\PXD014381\\RA_2\\interact.pep.xml":
        return "QX01872"

d = []
for fn in sys.argv[1:]:
    data = pd.read_csv(fn, sep="\t")
    #data["is_unlocalized"] = data["modifications"].apply(is_unlocalized)
    data["Full Sequence"] = data["Peptide"] + "|" + data["Assigned Modifications"].astype(str) + "|" + data["Observed Modifications"].astype(str)
    #TODO: #10 Make it so that Full Sequence is either with modifications like above but if no modifications: make it just the Sequence
    # data["Full Sequence"].fillna(data["Peptide"])
    #data['Full Sequence'] = data.apply(lambda row:  row["Peptide"] + row["Assigned Modifications"].astype(str) + row["Observed Modifications"] if np.isnan(row['c']) else row['c'],axis=1)
    data["Spectrum File"] = data["Spectrum File"].apply(rename_files)
    data = data[data["Is Unique"]]
    proteins = pd.read_csv(fn.replace("psm.tsv","peptide.tsv"), sep="\t")
    # proteins = proteins[proteins["is_shared_peptide"]==False]
    proteins = proteins[["Peptide","Protein ID"]].drop_duplicates()
    data = data.merge(proteins,on="Peptide",how="left")
    
    tmp = pd.DataFrame()
    tmp["Scan Retention Time"] = data["Retention"] / 60
    tmp["Precursor Charge"] = data["Charge"]
    tmp["Base Sequence"] = data["Peptide"]
    tmp["Full Sequence"] = data["Full Sequence"]
    tmp["Peptide Monoisotopic Mass"] = data["Calibrated Observed Mass"]
    tmp["Protein Accession"] = data["Protein ID_y"]
    tmp["File Name"] = data["Spectrum File"]

    d.append(tmp)

data = pd.concat(d)
data["Peptide Monoisotopic Mass"] = data.groupby("Full Sequence")["Peptide Monoisotopic Mass"].transform('median')

data.to_csv("flashlfqfrag.tsv",sep="\t",index=False)