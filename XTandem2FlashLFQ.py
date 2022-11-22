import sys
import pandas as pd

d = []
for i, fn in enumerate(sys.argv[1:]):
    print(i, ": first file")
    data = pd.read_csv(fn, sep="\t")
    data = data[data["Confidence"] >= 95]
    #data["is_unlocalized"] = data["modifications"].apply(is_unlocalized)
    data["Full Sequence"] = data["Sequence"]
    #TODO: #10 Make it so that Full Sequence is either with modifications like above but if no modifications: make it just the Sequence
    # data["Full Sequence"].fillna(data["Peptide"])
    #data['Full Sequence'] = data.apply(lambda row:  row["Peptide"] + row["Assigned Modifications"].astype(str) + row["Observed Modifications"] if np.isnan(row['c']) else row['c'],axis=1)
    if i == 0:
        data["Spectrum File"] = "2020-07-14_MS_PG_Proteome_Davide_Casserini_CTRL_01_20200717061234"
    elif i == 1:
        data["Spectrum File"] = "2020-07-14_MS_PG_Proteome_Davide_Casserini_CTRL_02_20200717083719"
    elif i==2:
        data["Spectrum File"] = "2020-07-14_MS_PG_Proteome_Davide_Casserini_CTRL_03_20200717110239"
    elif i==3:
        data["Spectrum File"] = "2020-07-14_MS_PG_Proteome_Davide_Casserini_H2O2_01"
    elif i==4:
        data["Spectrum File"] = "2020-07-14_MS_PG_Proteome_Davide_Casserini_H2O2_02"
    elif i==5:
        data["Spectrum File"] = "2020-07-14_MS_PG_Proteome_Davide_Casserini_H2O2_03"

    
    
    tmp = pd.DataFrame()
    tmp["Scan Retention Time"] = data["RT (min)"] #Check if this is okay, might be /60
    tmp["Precursor Charge"] = data["Charge"].astype(str).str[0]
    tmp["Base Sequence"] = data["Sequence"].str.split("-").str[1]
    tmp["Full Sequence"] = data["Sequence"]
    tmp["Peptide Monoisotopic Mass"] = data["m/z"] * data["Charge"].astype(str).str[0].astype(int)
    tmp["Protein Accession"] = data["Protein(s)"]
    tmp["File Name"] = data["Spectrum File"]

    d.append(tmp)

data = pd.concat(d)
data["Peptide Monoisotopic Mass"] = data.groupby("Full Sequence")["Peptide Monoisotopic Mass"].transform('median')

data.to_csv("flashlfqXtandem.tsv",sep="\t",index=False)