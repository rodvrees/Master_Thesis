# Master thesis Robbe Devreese
2022-2023
#### *ELUCIDATING THE PROTEIN OXIDATION LANDSCAPE WITH COMPREHENSIVE MODIFICATION SEARCHING*
This github repository contains the code and analysis work done during the Research Internship and Master Dissertation to compile the master thesis.
A short description of the folders/files/scripts in root directory (see README.md files in folders for further information on contents):
- Ageing_study: this folder contains analyses of the four datasets on ageing research.
- Alzheimers: this folder contains analyses on the AD CSF dataset
- Site_accessibility: this folder contains the python script used to calculate RSA values for modified protein sites, and the analysis notebook thereof
- Voor_vergelijking_OSE/From_scratch: this folder contains the analyses notebooks of the search engine comparison dataset

- [IPSAmodfilewriter.py](https://github.com/rodvrees/Master_Thesis/blob/main/IPSAmodfilewriter.py) is a python script which converts the downloaded Unimod database into a IPSA compatible modification file for bulk data visualisation.

- [IceLogoPrepImproved.py](https://github.com/rodvrees/Master_Thesis/blob/main/IceLogoPrepImroved.py) is a python script which takes an input file (.tsv file with two columns: Uniprot protein accession and modified AA position within this protein sequence) and outputs a file in a format which is compatible with IceLogo. This script was in the end not used in the final master dissertation. 

-  [MSFragger2FlashLFQ.py](https://github.com/rodvrees/Master_Thesis/blob/main/MSFragger2FlashLFQ.py) is a python script that converts MSFragger identifications to a FlashLFQ compatible identification file. Adapted from [ionbot2FlashLFQ.py](https://github.com/sdgroeve/ionbot.quant/blob/main/ionbot2FlashLFQ.py), written by prof. Degroeve.

-  [OxiAnalysis.py](https://github.com/rodvrees/Master_Thesis/blob/main/OxiAnalysis.py) is a small python module with functions written that were used in various analyses.

-  [ionbot2IPSA.py](https://github.com/rodvrees/Master_Thesis/blob/main/ionbot2IPSA.py) is a python script which takes an ionbot identification file (.csv) as input and outputs an IPSA compatible data file which can be used for bulk data visualisation. 

- [ionbotsearch.nf](https://github.com/rodvrees/Research_Internship/blob/main/ionbotsearch.nf) is a Nextflow script which can be used to
    1. Download RAW files from the given PRIDE accession code(s)
    2. Convert these RAW files to MGF files using ThermoRawFileParser
    3. Search these MGF files with ionbot using the appropriate configuration file

- [noftpnoproblem.nf](https://github.com/rodvrees/Master_Thesis/blob/main/noftpnoproblem.nf) is identical to [ionbotsearch.nf](https://github.com/rodvrees/Research_Internship/blob/main/ionbotsearch.nf), but does not include the first process. Instead of downloading RAW files from PRIDE, the user can provide a path to raw files.
