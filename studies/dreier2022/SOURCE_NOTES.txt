paper: https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-022-02451-y
code: https://github.com/biologger/htqpcr_ngs_data

Cheese from different locations.

Data collection, parser, and reprocessed by Maxwell.

Author correspondence:

Dear Maxwell,

I am happy to answer your questions. I am pleased to support the objective of enhancing the reproducibility of microbiome analyses. I have gained valuable insights over the years and would make several improvements today to ensure reproducibility and accessibility to data.

1) 
All samples were subjected to sequencing on a single chip, precluding the likelihood of batch effects. The cheese samples were of a uniform type, namely Raclette du Valais PDO, produced from raw cow's milk, and collected from 21 distinct cheesemakers. As the samples were PDO, their production was undertaken in accordance with the specified guidelines, but there is known variation between the cheese dairies. Furthermore, the starter cultures, the milk, and the production environment were different, which indicates that these are sources of variation between samples. However, detailed metadata from the cheese makers was not collected, as they typically maintain the confidentiality of their precise recipes. The name of the manufacturer or the exact location could not be provided, but the region in Switzerland where the cheeses were produced is documented in the biosample records on NCBI: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA786903
Example from the Full (text) summary from NCBI:
https://www.ncbi.nlm.nih.gov/biosample?Db=biosample&DbFrom=bioproject&Cmd=Link&LinkName=bioproject_biosample&LinkReadableName=BioSample&ordinalpos=1&IdsFromResult=786903

2: Metagenome or environmental sample from food fermentation metagenome
Identifiers: BioSample: SAMN23733339; Sample name: RdVEH20; SRA: SRS11245430
Organism: food fermentation metagenome
Attributes:
    /isolation source="cheese"
    /collection date="2018-03-29"
    /geographic location="Switzerland: Simplon-Dorf"
    /latitude and longitude="missing"
Description:
Raclette du Valais sample S20
Accession: SAMN23733339         ID: 23733339

2)

If I understand correctly, you are referring to the "copydict.json" file.

Unfortunately, I no longer have access to the pipeline output, but I can assist you in obtaining the connection from the IDs to the labels. Within the csv file, we find analogous information: ASV_ID as ID, the ASV sequence as Sequence, and the taxonomy hierarchy from the classification using the DAIRYdb-IDTAXA pipeline (https://github.com/marcomeola/DAIRYdb) with the following columns: "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species".

As demonstrated in the Jupyter Notebook used for the analysis (https://github.com/biologger/htqpcr_ngs_data/blob/main/htqpcr_ngs_comparison_R.ipynb), we never worked on the ASV level (because we did not assume strain level resolution using the V1-V2 rRNA partial sequence), we worked on the species level. We aggregated the ASVs to the species level (after filtering) for the comparison with the qPCR data by using the function:

ngs_count = hf().asv_to_speciestable(rac_data)

This functions similarly to the tax_glom function in phyloseq if you are familiar with this package:
https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/tax_glom

We implemented this by taking only the Species column in the csv file into consideration and sums all the counts for each ASV with the same Species value in a new dataframe with a species name. We then additionally added a conversion of the old lactobacilli taxonomy to the new one introduced in 2020 (https://doi.org/10.1128/AEM.02116-15).

These species names were then mapped to the copydict.json file.

I have created a JSON file that maps the species names used in the paper and in copydict.json back to the ASVs and the 16S rRNA gene copy number (if it was available). I will attach this to the email (species_to_asv_dict.json).The code to reproduce it in Python is here: (run the script in the hqpcr_ngs_data directory of the original repo).


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import pandas as pd

def newtax_dict():
    # dictionary for new taxonomic classification of Lactobacillus group
    with open('rawdata/new_species_names_dict.json') as f:
        for line in f:
            ntax_dict = json.loads(line)
    return ntax_dict

def new_species_labels(specieslist):
    ntax_dict = newtax_dict()
    newlabel = []
    for label in specieslist:
        newl = " ".join(label.split("_"))
        if newl in ntax_dict:
            newl = ntax_dict[newl]
            if 'subsp.' in newl:
                newl = newl.split(' subsp.')[0]
        newlabel.append(newl)
    return newlabel
def asv_to_speciestable(df):
    # Combine ASVs to species
    species = list(set(df["Species"]))
    spec_dict = {}
    for spec in species:
        q = df[df["Species"] == spec]
        spec_dict.update({spec: q.iloc[:, 9:].sum()})

    ngs_count = pd.DataFrame.from_dict(spec_dict).T
    ngs_count.index = new_species_labels(ngs_count.index)
    return ngs_count

# all paths are relative to the main directory of the original repo 
# load copy number dictionary
with open("rawdata/copydict.json", "r") as f:
    for line in f:
        copy_dict = json.loads(line)

# load raw reads data
raw_reads_df = pd.read_csv("rawdata/V18-22-21_ASV_counts_table.csv", header=0, index_col=0)

# convert ASV to species table (aggregate ASVs to species)
df = asv_to_speciestable(raw_reads_df)

# replace underscores with spaces
asv_dict = raw_reads_df["Species"].str.split("_").str.join(" ").to_dict()

# load new taxonomic classification dictionary
ntax_dict = newtax_dict()

# take unique species names with spaces instead of underscores
species_list = raw_reads_df["Species"].str.split("_").str.join(" ").unique()

# create a dictionary of species to ASVs and copy numbers
species_to_asv_dict = {}
for species in species_list:
    if species in ntax_dict.keys():
        # convert the names to new taxonomic classification
        species = ntax_dict[species]

    for k,v in asv_dict.items():
        if species in v:
            if species in species_to_asv_dict.keys():
                species_to_asv_dict[species]["asvs"].append(k)
            else:
                if species in copy_dict.keys():
                    species_to_asv_dict.update({species: {"copy_number": copy_dict[species]}})

                else:
                    species_to_asv_dict.update({species: {"copy_number": "unknown"}})
                species_to_asv_dict[species].update({"asvs": [k]})

# save the dictionary
with open("species_to_asv_dict.json", "w") as f:
    json.dump(species_to_asv_dict, f)
# print the species to ASV mapping
for k, v in species_to_asv_dict.items():
    print(f"{k}:")
    print(v["asvs"])


 3)

The reads in the SRA archive are the raw data that were solely primer trimmed using cutadapt. There was no Host removal step.

The Method described in “Reference sequence alignments” only applies to the sequences used in Supplementary Figure S1: the alignment of the 16S rRNA gene V1-V2 partial sequences.

The representative genomes are:

https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1590,60520,1589&reference_only=true

L. pentosus: GCF_003641185.1
L. paraplantarum: GCF_029025825.1
L. plantarum: GCF_009913655.1 



I hope this helps and feel free to ask if you need more information. 



Best regards,
Matthias


