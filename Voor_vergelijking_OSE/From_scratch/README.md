# Search engine comparison
This folder contains the analysis notebooks in which the 4 search engines Comet, MSFraggger, ionbot and pFind were compared on the same MS dataset.
- [Analysis.ipynb](https://github.com/rodvrees/Master_Thesis/blob/main/Voor_vergelijking_OSE/From_scratch/Analysis.ipynb): main analysis, contains Upsetplots, confidence metrics comparison and main PSM difference plot
- [CheckSimils.ipynb](https://github.com/rodvrees/Master_Thesis/blob/main/Voor_vergelijking_OSE/From_scratch/CheckSimils.ipynb): contains the Venn diagram showing overlap between deamidation ambiguity spectra and overlap between spectra assigned a propionamide mod by ionbot and pFind
- [Formatter.ipynb](https://github.com/rodvrees/Master_Thesis/blob/main/Voor_vergelijking_OSE/From_scratch/Formatter.ipynb): is used to format the PSMs made by the 4 SEs in a comparable way
- [IorL.ipynb](https://github.com/rodvrees/Master_Thesis/blob/main/Voor_vergelijking_OSE/From_scratch/IorL.ipynb) shows examples of PSMs that differ in sequence but not in modification are usually due to I->L changes in peptide sequence, coming from different proteins in the same protein family
- [OpenSEcomp.ipynb](https://github.com/rodvrees/Master_Thesis/blob/main/Voor_vergelijking_OSE/From_scratch/OpenSEcomp.ipynb) shows the comparison between the three OSEs
- [Oxmods.ipynb](https://github.com/rodvrees/Master_Thesis/blob/main/Voor_vergelijking_OSE/From_scratch/Oxmods.ipynb) shows the plot of the oxPTMs identified by ionbot which were missed by Comet 
- [Preprocess.ipynb](https://github.com/rodvrees/Master_Thesis/blob/main/Voor_vergelijking_OSE/From_scratch/Preprocess.ipynb) is the notebook where the data from the 4 SEs is preprocessed before formatting and analysis
- [Replacements.ipynb](https://github.com/rodvrees/Master_Thesis/blob/main/Voor_vergelijking_OSE/From_scratch/Replacements.ipynb) shows common modification replacements between SEs
- [deamidation_check.ipynb](https://github.com/rodvrees/Master_Thesis/blob/main/Voor_vergelijking_OSE/From_scratch/deamidation_check.ipynb) shows the nature of PSMs with deamidation ambiguity
