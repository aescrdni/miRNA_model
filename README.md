# A model for miRNA-mediated repression, implementing translation-related features
This repository contains the code and parameters for the miRNA-mediated repression model, used in our [Bioinformatics paper](https://doi.org/10.1093/bioinformatics/btaa1021). The code is freely available for academic usage.

# Requirements
The model was implemented in MATLAB R2019a, and works only on linux, as it incorporates several linux-based tools.

# Running the model
The main function, used to predict repression, is ```predict.m```. It accepts an ENST number and a miRNA sequence, and outputs their canonical binding sites and predicted repression value (optionally including the calculated features as well). Before running the code, the gene data file should be downloaded from [here](https://www.cs.tau.ac.il/~tamirtul/miRNA_Data.zip) and unzipped to the "Data" folder.

# Input

```[total_repression,feature_table] = predict(mRNA_ID,miRNA_seq,with_features)```

```predict``` requires the following input:

1. ```mRNA_ID``` [char]<br />
The ENST of the mRNA, e.g. ```'ENST00000376838'```.

2. ```miRNA_seq``` [char]<br />
The complete nucleotide sequence of the miRNA.

3. ```with_features``` [double]<br />
0/1 index for including calculated features in the output or only the predicted repression

# Output
```total_repression``` is the predicted total log2 fold change of mRNA level, as a result of miRNA-mediated repression.

```feature_table``` depends on ```with_features```:

If ```with_features == 0```, ```feature_table``` includes canonical binding sites in the 3'UTR and the last third of the ORF, along with the per-site predicted repression.

If ```with_features == 1```, ```feature_table``` includes all per-site features as well.

# Custom prediction
Users can run custom mRNAs and data, if their mRNA doesn't exist in the model database, or if they wish to include custom conservation/RiboSeq data. See predict_site manual.md for more details.

# Credits
Authors: Shaked Bergman, Alon Diament and Tamir Tuller.

The model incorporates programs and data from multiple sources; please see the paper for full credits.
