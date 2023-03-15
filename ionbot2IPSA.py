import argparse
import os
import numpy as np
import re


def format_modification(string):
    formatted_modslist = []
    if string == np.nan:
        return np.nan
    elif string == "semi_tryptic":
        return np.nan
    else:
        lijst = string.split("|")
        positions = lijst[0::2]
        mods = lijst[1::2]
        for pos, mod in zip(positions, mods):
            pos = str(pos).replace("0","1")
            modpipe = mod.replace(":", "|")
            # modded_AA = re.match(pattern=r"\[\D|N-TERM\]", string=mod)[0]
            # modification = re.match(pattern=r"\[\d\](.*)\[\D|N-TERM\]")[1]
            formattedmod = modpipe + ":" + pos
            # formattedmod = modification + " " + modded_AA + ":" + pos
            formatted_modslist.append(formattedmod)
    return ";".join(formatted_modslist)




#Define CL arguments
argParser = argparse.ArgumentParser()
argParser.add_argument("input", help= "Input ionbot result file name (.csv)")
argParser.add_argument("-o", "--output", metavar="DIR", help= "Full path to output directory")

#Parse CL arguments
args = argParser.parse_args()
input = vars(args)["input"]
output = vars(args)["output"]

if not output:
    output = os.path.dirname(__file__)

if not os.path.isfile(input):
    raise FileNotFoundError("Not a valid ionbot result file!")

os.chdir(output)
with open(input) as ifl:
    with open('outputfile.csv', "w") as ofl:
        ofl.write("Scan,Sequence,Charge,Modification\r\n")
        for line in ifl.readlines()[1:]:
            linelist = line.split(",")
            scan = linelist[2]
            seq = linelist[9]
            charge = linelist[7]
            modificationraw = linelist[10]
            if not "x|" in modificationraw:  #Filter out PSMs with unlocalized modifications
                formatted_mod = format_modification(modificationraw)
                ofl.write("{},{},{},{}\r\n".format(scan, seq, charge, formatted_mod))

ifl.close()
ofl.close()







