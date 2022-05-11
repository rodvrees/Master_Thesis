# Research Internship Robbe Devreese
2021-2022
#### *Search engine-based elucidation of the protein modification landscape under oxidative stress*
This github repository contains preliminary work done during the Research Internship in preparation for next year's master thesis.
A short description of the files/scripts:
- Uniprot: this folder contains csv files with SwissProt entries annotated with oxidation modifications.
    - [Human_Uniprot_DB.csv](https://github.com/rodvrees/Research_Internship/blob/main/Uniprot/Human_Uniprot_DB.csv) and [Mouse_Uniprot_DB.csv](https://github.com/rodvrees/Research_Internship/blob/main/Uniprot/Mouse_Uniprot_DB.csv) contain SwissProt entries with a 'Modified Residue' corresponding to an oxidative modification, for human and mural proteins, respectively.
    - [Human_ptm_oxidation.csv](https://github.com/rodvrees/Research_Internship/blob/main/Uniprot/Human_ptm_oxidation.csv) and [Mouse_ptm_oxidation.csv](https://github.com/rodvrees/Research_Internship/blob/main/Uniprot/Mouse_ptm_oxidation.csv) contain SwissProt entries for which the 'Post-translational modification' section contains the term 'oxidation', for human and mural proteins, respectively. These were collected because the “Modified Residue” filter does not include important sequence-independent annotations on some post-translational modifications. The same goes for [Human_ptm_oxidative_stress.csv](https://github.com/rodvrees/Research_Internship/blob/main/Uniprot/Human_ptm_oxidative_stress.csv) and [Mouse_ptm_oxidative_stress.csv](https://github.com/rodvrees/Research_Internship/blob/main/Uniprot/Mouse_ptm_oxidative_stress.csv), but with the term 'oxidative stress'.
- Configfiles: this folder contains four ionbot configuration files needed for ionbot searching, dependent on what search needs to be performed.
- [Pride_download.nf](https://github.com/rodvrees/Research_Internship/blob/main/PRIDE_download.nf) is a Nextflow script which can be used to download RAW files from the given PRIDE accession code(s). This script has since been incorporated as one of the processes in [ionbotsearch.nf](https://github.com/rodvrees/Research_Internship/blob/main/ionbotsearch.nf).
- [ionbotsearch.nf](https://github.com/rodvrees/Research_Internship/blob/main/ionbotsearch.nf) is a Nextflow script which can be used to
    1. Download RAW files from the given PRIDE accession code(s)
    2. Convert these RAW files to MGF files using ThermoRawFileParser
    3. Search these MGF files with ionbot using the appropriate configuration file
- [Oxidative_modifications.xlsx](https://github.com/rodvrees/Research_Internship/blob/main/Oxidative_modifications.xlsx) is an excel file which contains oxidative modifications collected from the literature, with their Unimod and Uniprot accession codes, if available, plus some additional information.
- [Verified_oxmods.xlsx](https://github.com/rodvrees/Research_Internship/blob/main/Verified_oxmods.xlsx) is an excel file containing verified oxidative modifications from the literature.
