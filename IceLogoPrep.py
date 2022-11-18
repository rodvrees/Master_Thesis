from pypdb import *
from Bio.PDB import PDBParser
import requests as r
from Bio import SeqIO
from io import StringIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from time import sleep
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import pandas as pd
import os
import OxiAnalysis as OA
import urllib
import sys
import logging
import importlib
importlib.reload(OA)
import warnings; warnings.simplefilter('ignore')
session = requests.Session()
retry = Retry(connect=3, backoff_factor=0.5)
adapter = HTTPAdapter(max_retries=retry)
session.mount('http://', adapter)
session.mount('https://', adapter)
logging.basicConfig(level=logging.INFO)

mastermod = pd.read_csv("/home/robbe/ionbot/mastersets/PXD012477_modifications.csv")
g = mastermod.groupby(["protein", "uniprot_id", "unexpected_modification", "position"])["#PSMs"].sum().to_frame().reset_index()
g10 = g[g["#PSMs"] > 10]
g10 = g10[g10["unexpected_modification"].isin(OA.modslist)]
logging.info('Making alignments...')

os.chdir("/home/robbe/ionbot/SE_for_IceLogo")
logging.info('Writing sequence files to {}'.format(os.getcwd()))
for mod in OA.modslist:
    alignmentlist = []
    with open("{}.txt".format(mod), "w") as f:
        #f.write("------------RSA values for {}------------\n".format(mod))
        logging.info('Writing alignments for proteins modified with {}'.format(mod))
        moddf = g10[g10["unexpected_modification"] == mod]
        #get modified AA for current modification
        pattern = re.compile(r"\[(.*?)\]")
        AA = re.findall(pattern= pattern, string= mod)[1]
        #Remove peptidoforms from contaminant proteins
        moddf = moddf[~moddf["uniprot_id"].str.startswith("sp|")]
        #get list of uniprot_ids that carry current mod
        protein_list = moddf["uniprot_id"].tolist()
        logging.info('Trying to write sequence alignments for {} proteins...'.format(len(protein_list)))
        #get list of modification position
        poslist = moddf["position"].tolist()
        #combine together, e.g. P105486|65 means that protein P105486 carried the modification at site 65
        pro_poslist = []
        for i in range(len(protein_list)):
            pro_pos = protein_list[i] + "|" + str(poslist[i])
            pro_poslist.append(pro_pos)

        if len(pro_poslist) == 0:
            logging.warning("No proteins found with this modification, skipping...")

        else:  
            check_list = []
            # PDBnotfound = 0
            # PDBfilenotfound = 0
            # noDSSP = 0
            # mismatch = 0
            for i, pr in enumerate(pro_poslist):
                logging.info("{:.2f}% finished".format((i/len(pro_poslist))*100))
                #Get protein and position
                protein = pr.split("|")[0]
                pos = int(pr.split("|")[1])
                #Search pdb accession id's for the protein
                # found_pdbs = Query(protein).search()
                
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
                try:
                    if Uniprotseq[pos-1] == AA:
                        cheat = 10*"X" + str(Uniprotseq) + 10*"X"
                        newpos = pos + 9
                
                    try:
                
                        prefix = cheat[newpos-8:newpos]
                        suffix = cheat[newpos+1:newpos+9]
                        mod_res = cheat[newpos]
                        alignment = prefix + mod_res + suffix
                        if alignment not in alignmentlist:
                            f.write(alignment+"\n")
                    except ValueError:
                        continue
                except IndexError:
                    continue
                        
                    # try:
                
                    #     prefix = cheat[pos-8:pos]
                    #     suffix = cheat[pos+1:pos+9]
                    #     mod_res = cheat[pos]
                    #     alignment = prefix + mod_res + suffix
                    #     f.write(str(alignment) + "\n")
                    # except ValueError:
                    #     continue
        
        #TODO #9 don't do RSA file for each mod, but one for all with headers
        if len(pro_poslist) != 0:
            logging.info("Done writing alignments for residues modified with {}".format(mod))

            # #Print all RSA's to file                    
            # with open("{}_RSA.txt".format(mod),'w') as f:
            #     for i in check_list:
            #         f.write(str(i))
            #         f.write("\n")
    logging.info("Writing output file finished!")
                

