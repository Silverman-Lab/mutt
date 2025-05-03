Dear Maxwell,

We intended to put these types of files in an organized way in a public Dryad directory but sadly that's not operational. For expediency I am giving you the files directly from my laptop, as I used them for the analysis (therefore, the file names are all in the original shorthand instead). I give you a key of what's what below. I also attached the R script I used to integrate these data for the publication. Again, this R script is not unfortunately well-annotated as a public document, but it may be helpful. I can easily imagine some further confusion or clarification will arise - feel welcome to reach out to me and I'll try to help you get going.

Greetings,
Derek


>16S OTU table, taxonomy, and metadata respectively:
515_blackblue_all_Zotutab_20181206.txt
515_blackblue_all_zotus.tax
515_blackblue_metadata.txt

>ITS OTU table, taxonomy, and metadata respectively:


AgITS1blackblue_READ1_271_all_Zotutab.txt
AgITS1blackblue_READ1_271_all_zotus.tax
ITS_blackblue_metadata_v2.txt

>Metagenome count tables at various taxonomic limits:
BacteriaGenusRaw.txt
FungiGenusRaw.txt
BacteriaFamilyRaw.txt
FungiFamilyRaw.txt
BacteriaPhylumRaw.txt

>Metagenome total reads for each sample (needed for any normalization of metagenome data)
total_seq.txt

>Metagenome reads assigned to the plant chromosome
chromCountsNew.txt

>Unassignable metagenome reads per sample
numUnaligned.txt

>The rather messy R script I used to integrate these data and prepare the figures, not written to be read easily or followed, but still perhaps useful:
correct_16S_load.R

>A metadata table of the metagenomes that connects the IDs from ENA to the IDs used throughout these other files
Metadata_PRJEB31530.xlsx



