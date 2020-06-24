# A model for miRNA-mediated repression, implementing translation-related features
This repository contains the code and parameters for the miRNA-mediated repression model, used in our paper [TBD].

# Prerequisites
The model was implemented in MATLAB R2019a, and works only on linux, as it augments several linux-based tools.

# Running the model r


The main function, used to predict repression, is predict.m. It accepts a target site and outputs its predicted repression value.

# Input
"predict" requires the following input:

1. RNA [char] 
The complete nucleotide sequence of the RNA, including 5'UTR, ORF and 3'UTR.

2. miRNA [char]
