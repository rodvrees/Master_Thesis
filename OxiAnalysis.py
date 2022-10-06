#list of Oxidative modifications
from pyteomics import mass as pymass
db = pymass.Unimod()
modslist = []
for p in range(len(db.mods)):
    for pp in db.mods[p]['specificity']:
        if db.mods[p]['record_id'] in [35, 53, 129, 130, 205, 206, 275, 288, 318, 319, 335, 340, 342, 344, 345, 348, 349, 350, 351, 352, 354, 
        359, 360, 368, 378, 392, 401, 421, 425, 534, 540, 548, 569, 720, 721, 743, 743, 860, 936, 936, 937, 949, 1384, 1914, 1915, 1916, 1917, 1918, 
        1922, 1923, 1924, 1925, 1926, 1927, 1928, 1929]:
            t = db.mods[p]['title']
            t = t.replace("[",":").replace("]",":")
            mod = "[" + str(db.mods[p]['record_id']) + "]" + t + "[" + pp['site'] + "]"
            modslist.append(mod)
modslist.append('[35]oxidation[M]')


def peptidoform_name(row):
    """
    Writes peptidoform names in column.
    Parameter: df row 
    Use with df.apply()
    """
    import re
    matched_peptide = row["matched_peptide"]
    
    def splitatn(strng, sep, pos):
        strng = strng.split(sep)
        return sep.join(strng[:pos]), sep.join(strng[pos:])

    #If no modifications, return None
    if row["modifications"] == "None":
        return matched_peptide
    else:
        #Separate modifications into a list
        modifications = splitatn(row["modifications"],"|",2)
        modifications = [i for i in modifications if i]

    
    poslist = []
    modlist = []
    #For every modification, get position-1 (because ionbot gives nth amino acid, starting from 1)
    for mod in modifications:
        pos = mod.split("|")[0]
        if pos == "0":
            poslist.append("N-TERM")
        elif pos != "x": 
            poslist.append(int(pos)-1)
    
        #Get name of modification
        modi = mod.split("|")[1]
        #Remove modified amino acid from modification name
        pattern = re.compile(pattern= r"\[\D+\]")
        modif = re.sub(pattern, "", modi)
        modlist.append(modif)

    #Moddict {position: modification}
    moddict ={poslist[i]: modlist[i] for i in range(len(poslist))}
    peptidoform_list = []
    #Reconstruct peptidoform
    if "N-TERM" in moddict:
        peptidoform_list.append(moddict["N-TERM"])
    for i, aa in enumerate(matched_peptide):
        if i in moddict:
            peptidoform_list.append(aa)
            peptidoform_list.append(moddict[i])
        else:
            peptidoform_list.append(aa)
    
    peptidoform = "".join(peptidoform_list)
    return peptidoform

#Get position of modification out of ionbot row (use with apply)
def get_positions(row):
    """
    Get position of modification out of ionbot modifications row
    Parameter: df row
    Use with df.apply()
    """
    if str == "None":
            return None
    if str != "":
        lijst = str.split("|")
        result = lijst[0::2]
        if len(result) == 1:
            return result[0]
        else:
            return result
        

#Get modifications out of ionbot row (use with apply)
def get_modification(str):
    """
    Get modifications out of ionbot row
    Parameter: df row
    Use with df.apply()
    """
    if str == "None":
        return None
    lijst = str.split("|")
    if len(lijst) <= 1:
        return None
    else:
        result = lijst[1::2]
        if len(result) == 1:
            return result[0]
        else:
            return result

#Returns boolean whether PSM is oxidatively modified or not
def oxidatively_modified(str):
    """
    Returns boolean value based on whether PSM is oxidatively modified or not
    Parameter: df row
    Use with df.apply()
    """
    for mod in modslist:
        if mod in str:
            return True
        else:
            continue
    return False

#Filters PSMs found in all replicates 
def replicate_filter(df, n_of_replicates):
    """
    Filters dataframe for only PSMs found in n replicates
    Parameters: dataframe, number of replicates
    """
    return df[df.groupby("Peptidoform_name")["spectrum_file"].transform('nunique') == n_of_replicates]

#Makes modratios file of occurence of PTM/total PSMs
def modratios(df):
    """
    Makes modratios dataframe of occurrence of oxPTM/total PSMs
    Parameter: df
    """
    import pandas as pd
    peptidoforms = df["Modification"]
    Modification = []
    Ratios = []
    for mod in modslist:
        Modification.append(mod)
        val = 0
        total = 0
        for i, peptidoform in peptidoforms.items():
            if peptidoform == None:
                total += 1
                continue
            elif mod in peptidoform:
                total += 1
                val += 1
        Ratios.append(val/total)
    dt = {"Modification" : Modification, "Ratios" : Ratios}
    result = pd.DataFrame(dt)
    return result

#Makes Venn diagram of shared and non-shared peptidoforms (figure out how to put this in OxiAnalysis)
def condition_venn(listofdf, listoflabels):
    """
    Make Venn diagram of shared and non-shared peptidoforms between dataframes
    Parameters: list of dataframes, list of dataframe labels
    """
    import matplotlib.pyplot as plt
    import matplotlib_venn as venn
    if len(listofdf) == 2:
        setlist = []
        for df in listofdf:
            s = set(df["Peptidoform_name"])
            setlist.append(s)
        total = len(set().union(*setlist))
        plt.figure(figsize=(10,10))
        diagram = venn.venn2(setlist, set_labels=listoflabels, subset_label_formatter=lambda x: f"{(x/total):1.0%}")
        return diagram
    elif len(listofdf) == 3:
        setlist = []
        for df in listofdf:
            s = set(df["Peptidoform_name"])
            setlist.append(s)
        total = len(set().union(*setlist))
        plt.figure(figsize=(10,10))
        diagram = venn.venn3(setlist, set_labels=listoflabels, subset_label_formatter=lambda x: f"{(x/total):1.0%}")
        return diagram
        

#Returns list and dataframe of PSMs that occur in treatment df but not in control df
def comparelist(treatment_df, control_df):
    """
    Returns list and dataframe of PSMs that occur in treatment df but not in control df
    Parameters: treatment dataframe, control dataframe
    """
    treatmentset = set(treatment_df["Peptidoform_name"])
    controlset = set(control_df["Peptidoform_name"])
    diffset = treatmentset - controlset
    difflist = list(diffset)
    diffdf = treatment_df[treatment_df['Peptidoform_name'].isin(difflist)]
    return difflist, diffdf

#Makes modcounts file of how many PTM occurs in PSM df
def modcounts(df):
    """
    Makes modcounts file of how many times an oxPTM occurs in df
    Parameter: dataframe
    """
    import pandas as pd
    peptidoforms = df["Modification"]
    Modification = []
    countlist = []
    for mod in modslist:
        Modification.append(mod)
        val = 0
        for i, peptidoform in peptidoforms.items():
            if peptidoform == None:
                continue
            elif mod in peptidoform:
                val += 1
        countlist.append(val)
    dt = {"Modification" : Modification, "Counts" : countlist}
    result = pd.DataFrame(dt)
    return result

#Gives relative level of PSMs containing unmodified residues 
def relative_PSM_modification(df):
    """
    Gives relative level of PSMs containing unmodified residues
    Parameter: dataframe
    """
    import pandas as pd
    import re
    amino_acids = ["A","R","N","D","C","Q","E","G","H","I","K","M","F","P","S","T","W","Y","V"] #No L because L not found ==> I = I/L
    ratiolist = []
    for aa in amino_acids:
        filtered = df[df['matched_peptide'].str.contains(aa)]
        n_of_psms = filtered.shape[0]
        n_of_modified_psms = 0
        
        for index, row in filtered.iterrows():
            peptidoform = row['Peptidoform_name']
            pattern = re.compile(r"[A-Z](?=\[)")
            modified_aa_list = re.findall(pattern, peptidoform)
            for i in modified_aa_list:
                if i == aa:
                    n_of_modified_psms += 1
                    break
        modified_psm_ratio = n_of_modified_psms/n_of_psms
        unmodified_psm_ratio = 1 - modified_psm_ratio
        if unmodified_psm_ratio < 0:
            ratiolist.append(0) 
        else:
            ratiolist.append(unmodified_psm_ratio)           
    dt = {"Amino acid" : amino_acids, "Relative level of PSMs containing unmodified residue" : ratiolist}
    df = pd.DataFrame(dt)
    return df

#Returns pie chart with distribution of non-modified methionines, singly oxidized methionines (Met Sulfoxide) and doubly oxidized methionines (Met Sulfones)
def methionine_overview(df, ax = None):
    """
    Returns pie chart with distribution of 
    non-modified methionines, singly oxidized methionines (Met Sulfoxide) and 
    doubly oxidized methionines (Met Sulfones)
    Parameters: dataframe, ax name (default: ax, you can probably leave out this parameter)
    """
    import matplotlib.pyplot as plt
    ax = ax or plt.gca()
    filtered = df[df["matched_peptide"].str.contains("M")]
    non_modified = 0
    met_sulfoxide = 0
    met_sulfone = 0
    for index, row in filtered.iterrows():
                peptidoform = row['Peptidoform_name']
                modifications = row['Modification']
                if modifications == None:
                    non_modified += 1
                elif type(modifications) == str:
                    if modifications == "[35]oxidation[M]":
                        met_sulfoxide += 1
                    elif modifications == "[425]Dioxidation[M]":
                        met_sulfone += 1
                elif type(modifications) == list:
                    for mod in modifications:
                        if modifications == "[35]oxidation[M]":
                            met_sulfoxide += 1
                        elif modifications == "[425]Dioxidation[M]":
                            met_sulfone += 1
    labels = []
    sizes = []
    if non_modified != 0:
        sizes.append(non_modified)
        labels.append("Met")
    if met_sulfoxide != 0:
        sizes.append(met_sulfoxide)
        labels.append("Met Sulfoxide")
    if met_sulfone != 0:
        sizes.append(met_sulfone)
        labels.append("Met Sulfone")

    return ax.pie(sizes, labels=labels, autopct='%1.1f%%', shadow= True, startangle = 90)