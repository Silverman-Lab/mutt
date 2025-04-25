This and '2019_contijoch_eLife_5dat16S' are from the same paper doi; https://doi.org/10.7554/eLife.40553

parse_2019_contijoch_eLife_metagenomic return one output for shotgun metagenomics data. This data is not used to generate main results in the paper, but was used only for detection/estimation of non-bacterial DNA in the samples (they found only a minority of samples had any non-bacterial sequences identified, and for the sample that did, it was low abundance).

parse_2019_contijoch_eLife_16S will return list with five elements where each element is a list that includes scale, count, metadata, and taxonomy for each source data based on 16S rRNA amplicon sequencing. To find the data used to generate the main results in the paper, go to '2019_contijoch_eLife_5dat16S'.