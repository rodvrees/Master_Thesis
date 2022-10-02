#Write peptidoform names in column
def peptidoform_name(row):
    import re
    matched_peptide = row["matched_peptide"]
    def splitatn(strng, sep, pos):
        strng = strng.split(sep)
        return sep.join(strng[:pos]), sep.join(strng[pos:])

    if row["modifications"] == "None":
        return matched_peptide
    else:
        modifications = splitatn(row["modifications"],"|",2)
        modifications = [i for i in modifications if i]

    
    poslist = []
    modlist = []
    for mod in modifications:
        pos = mod.split("|")[0]
        if pos != "x":
            poslist.append(int(pos)-1)
        modi = mod.split("|")[1]
        pattern = re.compile(pattern= r"\[\D\]")
        modif = re.sub(pattern, "", modi)
        modlist.append(modif)

    moddict ={poslist[i]: modlist[i] for i in range(len(poslist))}
    peptidoform_list = []
    for i, aa in enumerate(matched_peptide):
        if i in moddict:
            peptidoform_list.append(aa)
            peptidoform_list.append(moddict[i])
        else:
            peptidoform_list.append(aa)
    
    peptidoform = "".join(peptidoform_list)
    return peptidoform

