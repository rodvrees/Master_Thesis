from pypdb import *
from Bio.PDB import PDBParser
import requests as r
from Bio import SeqIO
from io import StringIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import pandas as pd
import os
import OxiAnalysis as OA
import urllib
import sys
import logging
import importlib
importlib.reload(OA)
import warnings; warnings.simplefilter('ignore')

logging.basicConfig(level=logging.INFO)
def seqmatch(alignment, pos): #TODO do sanity check met if nieuw gemapte AA == doel AA of check nog even hoe beter op te lossen
    pos = pos - 1

    for i, (a, b) in enumerate(zip(alignment[0], alignment[1])):
        try:
            if alignment.seqA[pos] == alignment.seqB[pos]:
                PDBseq2B = alignment.seqB[0:pos+1]
                countB = PDBseq2B.count("-")
                PDBseq2A = alignment.seqA[0:pos+1]
                countA = PDBseq2A.count("-")
                PDBseq2 = alignment.seqB.replace('-', "")
                if countA > countB:
                    need_added = countA - countB
                    PDBpos = pos + need_added
                elif countA < countB:
                    need_substracted = countB - countA
                    PDBpos = pos - need_substracted
                else:
                    PDBpos = pos
                #PDBpos = pos - count #get index of the AA out of this
                return PDBpos
            else:
                return None
        except IndexError:
                return None

def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None

mastermod = pd.read_csv("/home/robbe/ionbot/mastersets/PXD012477_modifications.csv")
mastermod = mastermod[mastermod["Accession"] != 6657]
g = mastermod.groupby(["protein", "uniprot_id", "unexpected_modification", "position"])["#PSMs"].sum().to_frame().reset_index()
g10 = g[g["#PSMs"] > 10]
g10[g10["unexpected_modification"].isin(OA.modslist)]
logging.info('RSA Calculation started!')

os.chdir("/home/robbe/ionbot/RSA_files")
logging.info('Writing result files to {}'.format(os.getcwd()))

# for each mod, get the peptidoforms that had at least 10 PSMs over all datasets
for mod in OA.modslist:
    logging.info('Calculating RSA for residues modified with {}'.format(mod))
    moddf = g10[g10["unexpected_modification"] == mod]
    #get modified AA for current modification
    pattern = re.compile(r"\[(.*?)\]")
    AA = re.findall(pattern= pattern, string= mod)[1]
    #Remove peptidoforms from contaminant proteins
    moddf = moddf[~moddf["uniprot_id"].str.startswith("sp|")]
    #get list of uniprot_ids that carry current mod
    protein_list = moddf["uniprot_id"].tolist()
    logging.info('Trying to calculate RSA values for {} proteins...'.format(len(protein_list)))
    #get list of modification position
    poslist = moddf["position"].tolist()
    #combine together, e.g. P105486|65 means that protein P105486 carried the modification at site 65
    pro_poslist = []
    for i in range(len(protein_list)):
        pro_pos = protein_list[i] + "|" + str(poslist[i])
        pro_poslist.append(pro_pos)

    if len(pro_poslist) == 0:
        logging.warning("No residues found with this modification, skipping...")

    else:  
        check_list = []
        PDBnotfound = 0
        PDBfilenotfound = 0
        noDSSP = 0
        mismatch = 0
        for i, pr in enumerate(pro_poslist):
            logging.info("{:.2f}% finished".format((i/len(pro_poslist))*100))
            #Get protein and position
            protein = pr.split("|")[0]
            pos = int(pr.split("|")[1])
            #Search pdb accession id's for the protein
            found_pdbs = Query(protein).search()
            
            #Get Uniprot sequence for protein
            baseUrl="http://www.uniprot.org/uniprot/"
            currentUrl=baseUrl+protein+".fasta"
            response = r.post(currentUrl)
            cData=''.join(response.text)

            Seq=StringIO(cData)
            pSeq=SeqIO.parse(Seq,'fasta')
            
            for record in pSeq:
                Uniprotseq = record.seq
            #first search result in PDB for Uniprot Accession gives PDB ID
            if found_pdbs != None:
                PDB_ID = found_pdbs[0]
                #dictionary to convert PDB file sequence to one-letter code
                d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
                #Download PDB file
                pdb_file = download_pdb(PDB_ID, '/home/robbe/ionbot/PDB_files')

                #For large structures, PDB file might not be available, and for those, DSSP can't be used
                if pdb_file != None:
                    # run parser
                    parser = PDBParser(QUIET=True)
                    structure = parser.get_structure('struct', pdb_file)    

                    # iterate each model, chain, and residue
                    # printing out the sequence for each chain
                    #Together forms PDBseq
                    for model in structure:
                        PDBseqlist = []
                        for chain in model:
                            
                            for residue in chain:
                                if residue.resname in d3to1:
                                    PDBseqlist.append(d3to1[residue.resname])
                                    PDBseq = "".join(PDBseqlist)
                    
                    #PDB sequence and Uniprot sequence often don't fully correspond. Because ionbot uses Uniprot seq to localize modification site, but DSSP algorithm uses PDB seq, 
                    #The modification site for Uniprot seq first needs to mapped onto the site for PDB seq, to do this, we align the sequences
                    alignments = pairwise2.align.globalxx(Uniprotseq, PDBseq)

                    #seqmatch gives the index of the modification site on the PDB sequence
                    m = seqmatch(alignments[0], pos)
                    model = structure[0]
                    
                    try:
                        dssp = DSSP(model, pdb_file)
                    
                        a_key = list(dssp)
                        #DSSP sometimes has unknown residues include in the sequence (denoted as X), if these are removed, you are left with the same as the PDBseq
                        for i in a_key:
                            if i[1] == "X":
                                a_key.remove(i)
                        #if matched amino acid matches the modified amino acid (i.e. alignment was succesful), append RSA to list
                        if m != None:
        
                            if a_key[m][1] == AA:
                                RSA = a_key[m][3]
                                check_list.append(RSA)
                        else:
                            mismatch += 1
                    except Exception:
                        noDSSP +=1
                        continue
                else:
                    PDBfilenotfound += 1
            else:
                PDBnotfound += 1
    
    #TODO #9 don't do RSA file for each mod, but one for all with headers
    if len(pro_poslist) != 0:
        logging.info("Done calculating RSA values for residues modified with {}".format(mod))
        logging.info("{} RSA values succesfully calculated ({:.2f}%)".format(len(check_list), ((len(check_list)/len(pro_poslist))*100)))
        logging.warning("{} PDB Query's were unsuccesful".format(PDBnotfound))
        logging.warning("{} PDB files could not be downloaded".format(PDBfilenotfound))
        logging.warning("{} RSA values could not be calculated by DSSP".format(noDSSP))
        logging.warning("{} modification sites could not be succesfully mapped from the Uniprot sequence on the PDB sequence".format(mismatch))
        logging.info("Writing {}_RSA output file...".format(mod))
        #Print all RSA's to file                    
        with open("{}_RSA.txt".format(mod),'w') as f:
            for i in check_list:
                f.write(str(i))
                f.write("\n")
        logging.info("Writing output file finished!")
            

