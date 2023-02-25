import sys
import pandas as pd

d = []
for fn in sys.argv[1:]:
    data = pd.read_csv(fn, sep="\t")
    #data["is_unlocalized"] = data["modifications"].apply(is_unlocalized)
    data["Full Sequence"] = data["Peptide"] + "|" + data["Assigned Modifications"].astype(str) + "|" + data["Observed Modifications"].astype(str)
    #TODO: #10 Make it so that Full Sequence is either with modifications like above but if no modifications: make it just the Sequence
    # data["Full Sequence"].fillna(data["Peptide"])
    #data['Full Sequence'] = data.apply(lambda row:  row["Peptide"] + row["Assigned Modifications"].astype(str) + row["Observed Modifications"] if np.isnan(row['c']) else row['c'],axis=1)
    data["Spectrum File"] = data["Spectrum"].str.split(".").str[0]
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