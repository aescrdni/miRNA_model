# A model for miRNA-mediated repression, implementing translation-related features
This repository contains the code and parameters for the miRNA-mediated repression model, used in our paper [TBD].

# Requirements
The model was implemented in MATLAB R2019a, and works only on linux, as it incorporates several linux-based tools.

# Running the model
The main function, used to predict repression, is predict.m. It accepts a target site (along with some additional data) and outputs its predicted repression value - possibly with the calculated features.

predict_example.m is an example script of running the model, using example data (10 target sites from the validation dataset) stored in input_example.mat

# Input
"predict" requires the following input:

1. RNA [char]<br />
The complete nucleotide sequence of the mRNA, including 5'UTR, ORF and 3'UTR.

2. ORF/UTR3_start [double]<br />
The first coordinate of the ORF/UTR3 region in RNA (e.g. if the length of the 5'UTR is 50nt, ORF_start = 51)

3. miRNA [char]<br />
The complete nucleotide sequence of the miRNA.

4. RNA_start [double]<br />
The first coordinate of the target site in RNA, i.e. the one paired with position 7 in the miRNA for 6mer/7mer-A1 sites, or the one paired with position 8 in the miRNA for 7mer-m8/8mer sites.

5. seed_type [char]<br />
The canonical seed type, can be one of 6mer/7mer-A1/7mer-m8/8mer

6. phastcons20/100 [double vector of size 1x(RNA length)]<br />
A vector of the PhastCons scores for each coordinate of the RNA, based on multiple alignment of hg38 with 19 mammals/99 vertebrates, respectively (as appears in the UCSC Genome Browser's database)

7. phylops20/100 [double vector of size 1x(RNA length)]<br />
Similar to phastcons20/100, based on PhyloP scores

8. riboseq [1x3 double vector]<br />
A 1x3 vector, containing average RiboSeq score for the first, second and last third of the ORF - as appears in the paper. Estimated Ribo-Seq values for the training and validation gene sets appear in Table S2.

9. with_features [double]<br />
0/1 index for including calculated features in the output or only the predicted repression

# Output
The output type depends on "with_features":

If with_features == 0, output is the log2 fold change of mRNA levels (as a result of the miRNA-mediated repression via the input target site).

If with_features == 1, output is a 1-row table, containing all calculated features as well as the log2 fold change value.

# Finding target sites

# Credits
Authors: Shaked Bergman, Alon Diament and Tamir Tuller.

The model incorporates programs and data from multiple sources; please see the paper for full credits.
