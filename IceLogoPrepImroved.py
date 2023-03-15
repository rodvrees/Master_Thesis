import warnings; warnings.simplefilter('ignore')
import requests as r
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from time import sleep
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import os
import logging
import argparse
from tqdm import tqdm

session = r.Session()
retry = Retry(connect=3, backoff_factor=0.5)
adapter = HTTPAdapter(max_retries=retry)
session.mount('http://', adapter)
session.mount('https://', adapter)
logging.basicConfig(level=logging.INFO)


#Define CL arguments
argParser = argparse.ArgumentParser()
argParser.add_argument("input", help= "Input file name")
argParser.add_argument("-o", "--output", metavar="DIR", help= "Full path to output directory")
argParser.add_argument("-n", "--name", help="Output file name (Default: IceLogoSeqs.txt)", default="IceLogoSeqs.txt")
argParser.add_argument("-a", "--AA", help="center AA; if variable, leave empty")

#Parse arguments
args = argParser.parse_args()
input = vars(args)["input"]
output = vars(args)["output"]
name = vars(args)["name"]
AA = vars(args)["AA"]
if AA not in ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","X","Y", None]:
    raise ValueError("Not a valid amino acid")
    

if not output:
    output = os.path.dirname(__file__)

logging.info("IceLogo Seq preparation!")
logging.info("Input file: {}".format(input))
logging.info("Output directory: {}".format(output))
logging.info("Writing IceLogo sequences to {}".format(output))

os.chdir(output)
with open(input) as fn:
    with open(name, "w") as of:
        for line in tqdm(fn.readlines()[1:]):
            linelist = line.split("\t")
            protein = linelist[0]
            position = int(linelist[1])

            #Get Uniprot sequence for protein
            baseUrl="http://www.uniprot.org/uniprot/"
            currentUrl=baseUrl+protein+".fasta"
            sleep(1)
            response = session.post(currentUrl)
            cData=''.join(response.text)

            Seq=StringIO(cData)
            pSeq=SeqIO.parse(Seq,'fasta')
            
            for record in pSeq:
                Uniprotseq = record.seq
                if AA:
                    if Uniprotseq[position-1] == AA:
                        cheat = 10*"X" + str(Uniprotseq) + 10*"X"
                        newpos = position + 9
                
                        prefix = cheat[newpos-8:newpos]
                        suffix = cheat[newpos+1:newpos+9]
                        mod_res = cheat[newpos]
                        alignment = prefix + mod_res + suffix
                        of.write(alignment + "\n")

                else:
                    cheat = 10*"X" + str(Uniprotseq) + 10*"X"
                    newpos = position + 9
            
                    prefix = cheat[newpos-8:newpos]
                    suffix = cheat[newpos+1:newpos+9]
                    mod_res = cheat[newpos]
                    alignment = prefix + mod_res + suffix
                    of.write(alignment + "\n")
        of.close()
    fn.close()


                










