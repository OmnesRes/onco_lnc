# onco_lnc
This code will reproduce the Cox regression analyses in OncoLnc, www.oncolnc.org, along with creating convenient data structures for the expression data and clinical data, which are used to create the database backend of OncoLnc.

This repository already contains the clinical data necessary, but the expression data will need to be downloaded from https://tcga-data.nci.nih.gov/tcga/ and http://mitranscriptome.org/.

For the mRNA Cox regressions, the "rsem.genes.normalized_results" files from TCGA will need to be placed in the correct folders in tcga_data

For the miRNA Cox regressions, counting.py, located in the mirna folder, will need to be run on the "isoform.quantification" files from TCGA.  Make sure the processed files get placed in the correct folders for each cancer.

For the MiTranscriptome beta analysis, processing.py, located in the lncrna folder, will need to be run on mitranscriptome.expr.counts.tsv from http://mitranscriptome.org/.  Warning, this script may require over 6GB of RAM.  Place the resulting files in the correct locations for each cancer.

Once the expression files are in the correct locations you should be able to successfully run the desired Cox regressions.  Go to the cancer of interest and run the cox_regression.py file from the command line.  This code was run with Python 2.7, rpy2, and NumPy.  I don't know if it will run on other Python versions.

The resulting file is already included in the repository.  You can rename this file and check if you get the same results with your analysis.

The clinical data in this repository was downloaded Jan. 5th and 6th, 2016.  If you update the clinical data the code will likely still run.  However, the clinical data parsing is complex, and most of the cancers only required a partial merge of the clinical files to extract 100% of the data.  If you download new clinical data I cannot guarantee that a partial merge will be sufficient for perfect parsing of the files.  Full merges were performed for BRCA and UCEC.  If you want to guarantee maximum extraction of clinical data from newly downloaded data, I would recommend taking a look at those scripts, but they may be overly complex.  Feel free to write your own clinical data parser!  Just keep in mind that clinical data from TCGA can be spread across up to four files, and the data will be nonredundant with no clear order as to which file is the most recent.

Please report any problems to omnesresnetwork@gmail.com


