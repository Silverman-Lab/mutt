paper is behind paywall so included here in repo. Assessed the distribution of ( 10.1136/gutjnl-2015-311004 ) 106 previously amplicon-sequenced primary sclerosing cholangitis and/or inflammatory bowel disease (PSC/IBD) faecal samples - (PSC n=18, Crohn’s disease (CD) n=29, ulcerative colitis (UC) n=13, PSC–CD n=20 and PSC–UC n=26) 10 against
a background of non-disease-associated microbiome variation as
observed in the Flemish Gut Flora Project cohort (FGFP; n=1,120). On the basis of genus-level inter-individual differences in microbiome composition, the combined patient group was observed to separate from a set of 66 FGFP
healthy controls (mHCs) 10 that were matched by age, gender and
body mass index with PSC or PSC–IBD diagnosed individuals. 224 individuals with UC or CD (Chron's Disease) were then independently recruited as a validation cohort. 

code for QMP: https://github.com/raeslab/QMP/
counts obtained from: http://raeslab.org/software/QMP2/ this includes counts of 161 samples. 

Measured covariates:
Serum alkaline phosphatase (U/L), Faecal calprotectin (µg/g), Serum CRP (mg/L)	Average faecal cell count (cells/g), Moisture content (%), Predicted DMM Enterotype,	 RMP Observed Richness (N genera)

This study did not have counts, relative abundances, batch information, or the ASV-Sequences. However, the authors supplied the QMP matrix. It was from the QMP matrix in which we transformed to obtain the relative abundance matrix. 

parse script does some weird stuff: 
-requires an NCBI entrez api key because it uses the taxize package to get the corresponding taxa name from the relevant database we need to convert to for MLSCALE study. 
-it then manually corrects some that Maxwell looked up by hand as it cant classify some of the corresponding taxa

Sequencing data is not publicly available -- Jeroen Raes study.

all data collected by Maxwell Konnaris

parse.R finished 04/09/25


