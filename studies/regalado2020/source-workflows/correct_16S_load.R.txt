#correlation

source("/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/scripts/R/functions/microbiome_custom_functions.R")

library(vegan)
library(reshape)
library(ggplot2)
library(RColorBrewer)

normPLANT=TRUE
addPLANT=FALSE

normTOTAL=FALSE
addOTHER=TRUE

add_plant_to_LOAD=FALSE

removePLASTIDSmetagenome=TRUE

LOAD_CORRECTION="COMMON"  #TOTAL OR COMMON

#FOR NORMALIZATION DIAGRAMS
SUBTRACT=FALSE
SQRT=FALSE

MAKE_TOTAL_METABREAKDOWN_PLOT=TRUE
if(MAKE_TOTAL_METABREAKDOWN_PLOT==TRUE){
	normTOTAL=addOTHER=add_plant_to_LOAD=TRUE
	normPLANT=addPLANT=removePLASTIDSmetagenome=addEMPTY=FALSE
}
	
MAKE_SUPPLEMENTARY_NORM=FALSE
if(MAKE_SUPPLEMENTARY_NORM==TRUE){
	#NORM BY TOTAL READS
		norm_fig_by_total=TRUE	
		normTOTAL=addOTHER=add_plant_to_LOAD=TRUE
		addEMPTY=normPLANT=addPLANT=removePLASTIDSmetagenome=FALSE
		norm_fig_by_microbes=FALSE
	#NORM BY PLANT READS
		#norm_fig_by_plant=FALSE
		#norm_fig_by_microbes=FALSE
		#addEMPTY=normPLANT=addPLANT=TRUE
		#addOTHER=normTOTAL=add_plant_to_LOAD=removePLASTIDSmetagenome=FALSE
		#norm_fig_by_plant=TRUE	
		#normTOTAL=addOTHER=add_plant_to_LOAD=TRUE
		#addEMPTY=normPLANT=addPLANT=removePLASTIDSmetagenome=FALSE
}
	

#16S
#otutab="/Users/dlundberg/Documents/abt6/Regalado2017/16S_blackblue_platesV4/515_blackblue_all_zotutab_raw.txt"
otutab="/Users/dlundberg/Documents/abt6/Regalado2017/16S_blackblue_platesV4/515_blackblue_all_Zotutab_20181206.txt"
#tax="/Users/dlundberg/Documents/abt6/Regalado2017/16S_blackblue_platesV4/515_blackblue_all_zotus.tax"
tax="/Users/dlundberg/Documents/abt6/Regalado2017/16S_blackblue_platesV4/515_blackblue_all_zotus_20181207.tax"
metadata="/Users/dlundberg/Documents/abt6/Regalado2017/16S_blackblue_platesV4/515_blackblue_metadata.txt"
#load OTU table
	read.countTable(file=otutab)->otureads
#load metadata
	as.matrix(read.table(metadata, sep="\t"))->metadata
#convert otu_sample names to S numbers
	as.matrix(metadata[match(colnames(otureads), metadata[,1]),11])->colnames(otureads)
#load taxonomy
	as.matrix(read.table(file=tax, sep="\t"))->taxonomy
#remove contaminants from OTU table
	remove_contaminant_taxa(otureads, taxonomy, keywords=c("Chloroplast", "Mitochondria", "Archaea"))->otureads
	
#rename f:Cytophagaceae to f:Hymenobacteraceae
	gsub("f:Cytophagaceae,g:Hymenobacter", "f:Hymenobacteraceae,g:Hymenobacter", taxonomy)->taxonomy
	gsub("f:Pseudomonadaceae,g:Cellvibrio", "f:Cellvibrionaceae,g:Cellvibrio", taxonomy)->taxonomy

#filter by taxonomy to ___ level
	otureads<-filter_by_taxonomyic_level(countTable=otureads, taxonomy=taxonomy, keywords=c("f"), last_filter=TRUE)
	#otureads<-filter_by_taxonomyic_level(countTable=otureads, taxonomy=taxonomy, keywords=c("p"), last_filter=TRUE)
	
#threshold samples with <1000 reads
	otureads=otureads[,which(colSums(otureads)>=1000)]
#normalize table by dividing by column sums
	normalize100(otureads)->otureads

	
if(MAKE_TOTAL_METABREAKDOWN_PLOT==TRUE){
	mystery<-rep(0.1, ncol(otureads))
	rbind(otureads, mystery)->otureads
	allfungi<-rep(0.1, ncol(otureads))
	rbind(otureads, allfungi)->otureads
}
if(MAKE_SUPPLEMENTARY_NORM==TRUE){
	mystery<-rep(0.1, ncol(otureads))
	rbind(otureads, mystery)->otureads
	allfungi<-rep(0.1, ncol(otureads))
	rbind(otureads, allfungi)->otureads
}		


#ITS
ITStab="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/ITS1_blackblue_plates/countTables/AgITS1blackblue_READ1_271_all_Zotutab.txt"
ITStax="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/ITS1_blackblue_plates/countTables/AgITS1blackblue_READ1_271_all_zotus.tax"
ITSmetadata="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/ITS1_blackblue_plates/countTables/ITS_blackblue_metadata_v2.txt"
#load ITS table
	read.countTable(file=ITStab)->ITSreads
#load ITSmetadata
	as.matrix(read.table(ITSmetadata, sep="\t"))->ITSmetadata
#convert ITS_sample names to S numbers
	as.matrix(ITSmetadata[match(colnames(ITSreads), ITSmetadata[,1]),11])->colnames(ITSreads)
#load taxonomy
	as.matrix(read.table(file=ITStax, sep="\t"))->ITStaxonomy
#remove contaminants from ITS table
	remove_contaminant_taxa(ITSreads, ITStaxonomy, keywords=c("Chloroplast", "Mitochondria", "Archaea"))->ITSreads
#filter by taxonomy to family level
	ITSreads<-filter_by_taxonomyic_level(countTable=ITSreads, taxonomy=ITStaxonomy, keywords=c("f"), last_filter=TRUE)
#threshold samples with <1000 reads
	ITSreads=ITSreads[,which(colSums(ITSreads)>=1000)]
#normalize table by dividing by column sums
	normalize100(ITSreads)->ITSreads
	
	#do all the tech reps form part of all the common samples?


#load total reads
totalseq="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/count_tables/total_seq.txt"
read.table(totalseq, comment.char = "",  colClasses = "character")->totalseq
totalseq[2:nrow(totalseq),]->totalseq
as.numeric(totalseq[,2])->totalseq[,2]
c("sample", "total")->names(totalseq)

#plant load reads
#RepSamples="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/count_tables/RepSamples.txt"
#read.table(RepSamples, comment.char = "",  colClasses = "character")->RepSamples


#plant load reads
#chromCounts1="/Users/dlundberg/Documents/abt6/Regalado2017/count_tables/chromCounts.txt"
#read.countTable(chromCounts1)->chromCounts1
#chromCounts1[grep("188smplDtSt_S", rownames(chromCounts1)),]->chromCounts1
#gsub("188smplDtSt_", "", names(chromCounts1))->names(chromCounts1)
chromCounts2="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/count_tables/chromCountsNew.txt"
read.countTable(chromCounts2)->chromCounts2
chromCounts2[grep("S", rownames(chromCounts2)),]->chromCounts2
chromCounts2->chromCounts


#mystery load reads
mystery="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/count_tables/numUnaligned.txt"
read.countTable(mystery)->mystery
mystery[grep("188smplDtSt_S", rownames(mystery)),]->mystery
gsub("188smplDtSt_", "", names(mystery))->names(mystery)



mystery->mystery2
chromCounts->chromCounts2
for(t in 1:length(mystery2)){
	mystery[t]/totalseq[t,2]->mystery2[t]
	chromCounts[t]/totalseq[t,2]->chromCounts2[t]
}
#plot(mystery2, chromCounts2)

	
	
	
#next load the fungi family table
	metatab_F="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/count_tables/FungiFamilyRaw.txt"
	read.MetacountTable(metatab_F)->metareads_F
	metareads_F[setdiff(1:nrow(metareads_F),match(c("Bacteria", "TotalSeq"), rownames(metareads_F))),]->metareads_F
	#normalize table by dividing by total reads
	totalseq[match(totalseq[,1], colnames(metareads_F)),2]->total
	colSums(metareads_F)->allfungi
	for(c in 1:ncol(metareads_F)){
		(100*(metareads_F[,c])/total[c])->metareads_F[,c]
	}
	
############################################################################################################
#Fungi table centrifuge
	#read.countTable(file="/Users/dlundberg/Documents/abt6/Regalado2017/count_tables/FungiCentrifugeABS.txt")->metareads_F
	#for(c in 1:ncol(metareads_F)){
	#	(100*metareads_F[,c])/sum(metareads_F[,c])->metareads_F[,c]
	#}




#bacterial metagenome FAMILY
#first load the metagenome bacteria family table
	metatab_B="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/count_tables/BacteriaFamilyRaw.txt"
	#metatab_B="/Users/dlundberg/Documents/abt6/Regalado2017/count_tables/BacteriaGenusRaw.txt"  ###change back
	read.MetacountTable(metatab_B)->metareads_B
		as.matrix(metareads_B[c("Bacteria", "TotalSeq"),])->Total_bacteria
		Total_bacteria["Bacteria",]/Total_bacteria["TotalSeq",]->Total_bacteria_true_ratio		
	metareads_B[setdiff(1:nrow(metareads_B),match(c("Bacteria", "TotalSeq"), rownames(metareads_B))),]->metareads_B
		chromCounts[match(names(chromCounts), colnames(Total_bacteria))]->Total_plant
		Total_plant/Total_bacteria["TotalSeq",]->Total_plant_true_ratio		
	







if(MAKE_TOTAL_METABREAKDOWN_PLOT==TRUE){
	rbind(metareads_B, mystery, allfungi)->metareads_B	
}
if(MAKE_SUPPLEMENTARY_NORM==TRUE){
	rbind(metareads_B, mystery, allfungi)->metareads_B	
}	

	
	
	
	if("subsample"=="Subsample"){
		metatab_B="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/count_tables/BacteriaFamilyABS200kSubSample188Dtdt.txt"
		t(read.MetacountTable(metatab_B))->metareads_B
		#normalize table by taking sum of all taxa and dividing by the total
	}
	
	if(normTOTAL==TRUE){    #normalize by total reads
		totalseq[match(totalseq[,1], colnames(metareads_B)),2]->total
		chromCounts[match(names(chromCounts), colnames(metareads_B))]->plant
		for(c in 1:ncol(metareads_B)){
			(100*(metareads_B[,c])/total[c])->metareads_B[,c]
		}
	}
	if(normPLANT==TRUE){    #normalize by plant
		chromCounts[match(names(chromCounts), colnames(metareads_B))]->plant
		for(c in 1:ncol(metareads_B)){
			(100*(metareads_B[,c])/plant[c])->metareads_B[,c]
		}
		
		
	}
	


	
	
	
#bacterial metagenome PHYLUM
if("phylum"=="NO"){
	metatab_B="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/count_tables/BacteriaPhylumRaw.txt"
	read.MetacountTable(metatab_B)->metareads_B	
	#normalize table by taking sum of all taxa and dividing by the total
	totalseq[match(totalseq[,1], colnames(metareads_B)),2]->total
	for(c in 1:ncol(metareads_B)){
		(100*(metareads_B[,c])/total[c])->metareads_B[,c]
	}
	gsub("_<phylum>", "", rownames(metareads_B))->rownames(metareads_B)
}
if("subsample"=="sFFubsample"){
	otureads="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/count_tables/BacteriaFamilyABS200kSubSample188Dtdt.txt"
	t(read.MetacountTable(otureads))->otureads
	#normalize table by taking sum of all taxa and dividing by the total
	for(c in 1:ncol(otureads)){
		(100*(otureads[,c])/sum(otureads[,c]))->otureads[,c]
	}
}

	

#combine the tech rep samples for all tables
combine_techreps<-function(countTable){
	techreps1=c("S23","S24","S41","S50","S65","S66","S69","S70","S82","S83","S84","S85","S86","S87","S90","S96","S117","S118","S119","S131","S136","S144", "S141","S142")
	techreps2=c("115newAB","115newAB","150newAB","150newAB","113AB","113AB","124AB","124AB","206AB","206AB","156AB","156AB","177AB","177AB","174AB","174AB","39AB","166AB","166AB","39AB","21AB","21AB","33AB","33AB")
	cbind(techreps1, techreps2)->techreps
	for(pair in unique(techreps[,2])){
		which(match(colnames(countTable),techreps[which(techreps[,2]==pair),1], nomatch=0)>0)->matching_columns
		print(pair)
		print(length(matching_columns))
		if(length(matching_columns)>1){
			rowSums(countTable[, matching_columns])->newvalue
			countTable[, setdiff(1:ncol(countTable), matching_columns)]->countTable
			cbind(newvalue, countTable)->countTable
			pair->colnames(countTable)[1]
		}
	}
	return(countTable)
}
combine_techreps(otureads)->otureads
combine_techreps(ITSreads)->ITSreads
combine_techreps(metareads_B)->metareads_B
combine_techreps(metareads_F)->metareads_F
rbind(Total_plant_true_ratio, Total_plant_true_ratio)->Total_plant_true_ratio
combine_techreps(Total_plant_true_ratio)->Total_plant_true_ratio
Total_plant_true_ratio[1,]->Total_plant_true_ratio
Total_plant_true_ratio[grep("AB", names(Total_plant_true_ratio))]/2->Total_plant_true_ratio[grep("AB", names(Total_plant_true_ratio))]

#first look only at shared samples between all 3 tables
intersect(colnames(metareads_B), colnames(otureads))->common_samples
intersect(common_samples, colnames(ITSreads))->common_samples
metareads_B[,match(common_samples, colnames(metareads_B))]->metareads_B
metareads_F[,match(common_samples, colnames(metareads_F))]->metareads_F
otureads[,match(common_samples, colnames(otureads))]->otureads
ITSreads[,match(common_samples, colnames(ITSreads))]->ITSreads
totalseq[match(common_samples, totalseq[,1]),]->totalseq


#BEGIN LOAD CORRECTION BY SHARED TAXA
#calculate load adjustment for ITS and OTU tables based on shared taxa.
#find shared fungal taxa between raw metagenome and ITS
#make shared and unique tables for 16S, ITS, and metagenome.
	normalize100(otureads)->otureads 
	intersect(rownames(metareads_F), rownames(ITSreads))->common_ITS_taxa	
	setdiff(rownames(ITSreads),rownames(metareads_F))->unique_ITS_taxa	
	ITSreads[match(common_ITS_taxa, rownames(ITSreads)),]->ITSreads_common	
	ITSreads[match(unique_ITS_taxa, rownames(ITSreads)),]->ITSreads_unique	
		setdiff(rownames(metareads_F),common_ITS_taxa)->unique_metareads_F_taxa	
		metareads_F[match(common_ITS_taxa, rownames(metareads_F)),]->metareads_F_common	
		rownames(metareads_F_common)<-common_ITS_taxa	
		metareads_F[match(unique_metareads_F_taxa, rownames(metareads_F)),]->metareads_F_unique	
	intersect(rownames(metareads_B), rownames(otureads))->common_16S_taxa	
	setdiff(rownames(otureads),rownames(metareads_B))->unique_16S_taxa	
	otureads[match(common_16S_taxa, rownames(otureads)),]->otureads_common	
	otureads[match(unique_16S_taxa, rownames(otureads)),]->otureads_unique	
		setdiff(rownames(metareads_B),common_16S_taxa)->unique_metareads_B_taxa	
		metareads_B[match(common_16S_taxa, rownames(metareads_B)),]->metareads_B_common	
		rownames(metareads_B_common)<-common_16S_taxa	
		metareads_B[match(unique_metareads_B_taxa, rownames(metareads_B)),]->metareads_B_unique
#First we calculate the ratio of unique taxa to common taxa 
	colSums(ITSreads_unique)/colSums(ITSreads_common)->unique_common_ratio_ITS
	colSums(otureads_unique)/colSums(otureads_common)->unique_common_ratio_16S
#Then we scale the amplicon common tables so that the sum of shared reads matches the sum of shared reads
#for the metagenome table
#load adjustment for common fungi is % common taxa should represent of the TOTAL reads
	#LOAD CORRECTION BY COMMON FUNGI
	if(LOAD_CORRECTION=="COMMON"){
		load_adjustment_common_fungi=colSums(metareads_F_common)/100    #percent common families are of metagenome total
		for(c in 1:length(load_adjustment_common_fungi)){
			load_adjustment_common_fungi[c]*ITSreads_common[,c]->ITSreads_common[,c]
			load_adjustment_common_fungi[c]*ITSreads_unique[,c]->ITSreads_unique[,c]
		}
	}
	#LOAD CORRECTION BY TOTAL FUNGI
	if(LOAD_CORRECTION=="TOTAL"){
		load_adjustment_total_fungi=colSums(metareads_F)/100    #percent common families are of metagenome total
		for(c in 1:length(load_adjustment_total_fungi)){
			load_adjustment_total_fungi[c]*ITSreads_common[,c]->ITSreads_common[,c]
			load_adjustment_total_fungi[c]*ITSreads_unique[,c]->ITSreads_unique[,c]
		}
	}
	#LOAD CORRECTION BY COMMON BACTERIA
	#if(LOAD_CORRECTION=="COMMON"){
	#	load_adjustment_common_bacteria=colSums(metareads_B_common)/100  #percent common families are of metagenome total
	#	for(c in 1:length(load_adjustment_common_bacteria)){
	#		load_adjustment_common_bacteria[c]*otureads_common[,c]->otureads_common[,c]
	#		load_adjustment_common_bacteria[c]*otureads_unique[,c]->otureads_unique[,c]
	#	}
	#}
	#LOAD CORRECTION BY COMMON BACTERIA
	if(LOAD_CORRECTION=="COMMON"){
		  #percent common families are of metagenome total
		load_adjustment_common_bacteria=colSums(metareads_B_common)/100 
		for(c in 1:length(load_adjustment_common_bacteria)){
				colSums(metareads_B_common)[c]/sum(otureads_common[,c])->adjustment_factor
				adjustment_factor*otureads_common[,c]->otureads_common[,c]
				adjustment_factor*otureads_unique[,c]->otureads_unique[,c]
		}
	}
		
	#LOAD CORRECTION BY TOTAL BACTERIA
	if(LOAD_CORRECTION=="TOTAL"){
		load_adjustment_total_bacteria=colSums(metareads_B)/100 
		#calculate load adjustment for ITS and OTU tables based on total bacteria.
		for(c in 1:length(load_adjustment_total_bacteria)){
			load_adjustment_total_bacteria[c]*otureads_common[,c]->otureads_common[,c]
			load_adjustment_total_bacteria[c]*otureads_unique[,c]->otureads_unique[,c]
		}
	}	

	colSums(otureads_common)->sum_common_otus
	colSums(otureads_unique)->sum_unique_otus
	colSums(ITSreads_common)->sum_common_ITS
	colSums(ITSreads_unique)->sum_unique_ITS

#sort(rowMeans((metareads_B_unique)))

#metareads_B_
	colSums(metareads_B_unique)->sum_unique_metareads_B
	match(rownames(metareads_B_unique), c("Chroococcidiopsidaceae","Gloeobacteraceae","Gloeomargaritaceae","Aphanizomenonaceae","Chlorogloeopsidaceae","Fortieaceae","Hapalosiphonaceae","Nostocaceae","Rivulariaceae","Scytonemataceae","Stigonemataceae","Symphyonemataceae","Tolypothrichaceae","Aphanothecaceae","Chroococcaceae","Cyanobacteriaceae","Microcystaceae","Coleofasciculaceae","Cyanothecaceae","Desertifilaceae","Gomontiellaceae","Microcoleaceae","Oscillatoriaceae","Dermocarpellaceae","Hydrococcaceae","Hyellaceae","Xenococcaceae","Spirulinaceae","Acaryochloridaceae","Chamaesiphonaceae","Leptolyngbyaceae","Merismopediaceae","Prochloraceae","Pseudanabaenaceae","Schizotrichaceae","Synechococcaceae"),
			 nomatch=0)->cyanos
	####REMOVES CYANOBACTERIA UNIQUE TO METAGENOME
if(removePLASTIDSmetagenome==TRUE){
	colSums(metareads_B_unique[which(cyanos==0),])->sum_unique_metareads_B
	colSums(metareads_B_unique[which(cyanos>0),])->sum_unique_metareads_B_cyanos
	sum_unique_metareads_B_cyanos[1:length(sum_unique_metareads_B_cyanos)]<-0
}

#metareads_F_common
	colSums(metareads_F_unique)->sum_unique_metareads_F

#combine unique sums and common all taxa for total table, scaled for load by the common ones
	rbind(otureads_common, sum_unique_otus)->otureads_LOAD
		colSums(otureads_LOAD)->sum_all_otu
	rbind(ITSreads_common, sum_unique_ITS)->ITSreads_LOAD
		colSums(ITSreads_LOAD)->sum_all_ITS
	rbind(metareads_B_common, sum_unique_metareads_B)->metareads_B_LOAD
		colSums(metareads_B_LOAD)->sum_all_metareads_B
	rbind(metareads_F_common, sum_unique_metareads_F)->metareads_F_LOAD
		colSums(metareads_F_LOAD)->sum_all_metareads_F


	
#add to the 16S table a taxa that is the sum of all fungi and vice versa
	#rbind(otureads_LOAD, sum_all_ITS)->otureads_LOAD
	#rbind(ITSreads_LOAD, sum_all_otu)->ITSreads_LOAD
	#rbind(metareads_B_LOAD, sum_all_metareads_F)->metareads_B_LOAD
	#rbind(metareads_F_LOAD, sum_all_metareads_B)->metareads_F_LOAD
	
	#THIS WAS USED AT ISME TO LOOK AT ONLY FUNGI PLUS BACTERIA
		#rbind(sum_all_metareads_B, sum_all_metareads_F)->metareads_B_LOAD
	
#load metadata
metadata=read.table(file="/Users/dklu0001/Documents/20220220_Germany/DOCUMENTS/abt6/Regalado2017/count_tables/sample_infoFinal.txt", sep="\t")


#Make relative abundance version of metareads tables
otureads_LOAD->otureads_ReAb
ITSreads_LOAD->ITSreads_ReAb
metareads_B_LOAD->metareads_B_ReAb
metareads_F_LOAD->metareads_F_ReAb
normalize100(metareads_B_LOAD)->metareads_B_ReAb
normalize100(metareads_F_LOAD)->metareads_F_ReAb
normalize100(otureads_LOAD)->otureads_ReAb
normalize100(ITSreads_LOAD)->ITSreads_ReAb


			
if(addOTHER==TRUE){
	add.other.taxa(otureads_LOAD, subtract=SUBTRACT)->otureads_LOAD
	add.other.taxa(ITSreads_LOAD, subtract=SUBTRACT)->ITSreads_LOAD
	add.other.taxa(metareads_B_LOAD, subtract=SUBTRACT)->metareads_B_LOAD
	add.other.taxa(metareads_F_LOAD, subtract=SUBTRACT)->metareads_F_LOAD
	if(add_plant_to_LOAD==TRUE){
		100*Total_plant_true_ratio[match(colnames(metareads_B_LOAD), names(Total_plant_true_ratio))]->plant
		metareads_B_LOAD["other",]-plant->metareads_B_LOAD["other",]
		rbind(plant, metareads_B_LOAD)->metareads_B_LOAD
	}
}


if(addPLANT==TRUE){
	add.plant(otureads_LOAD)->otureads_LOAD
	add.plant(ITSreads_LOAD)->ITSreads_LOAD
	add.plant(metareads_B_LOAD)->metareads_B_LOAD
	add.plant(metareads_F_LOAD)->metareads_F_LOAD
}




#order taxa in amplicon tables by that in metagenome tables (most abundant phyla first)
	order(rowSums(metareads_B_ReAb), decreasing=FALSE)->ordering_vector
		metareads_B_ReAb[ordering_vector,]->metareads_B_ReAb
		match(rownames(metareads_B_ReAb), rownames(otureads_ReAb), nomatch=0)->shared	
		setdiff(1:nrow(otureads_ReAb), shared)->notshared
		otureads_ReAb[c(shared, notshared),]->otureads_ReAb
	order(rowSums(metareads_F_ReAb), decreasing=FALSE)->ordering_vector
		metareads_F_ReAb[ordering_vector,]->metareads_F_ReAb
		match(rownames(metareads_F_ReAb), rownames(ITSreads_ReAb), nomatch=0)->shared	
		setdiff(1:nrow(ITSreads_ReAb), shared)->notshared
		ITSreads_ReAb[c(shared, notshared),]->ITSreads_ReAb
	###maybe the below is the same as ReAB?   check anyway if necessary
	order(rowSums(metareads_B_LOAD), decreasing=FALSE)->ordering_vector
		metareads_B_LOAD[ordering_vector,]->metareads_B_LOAD
		match(rownames(metareads_B_LOAD), rownames(otureads_LOAD), nomatch=0)->shared	
		setdiff(1:nrow(otureads_LOAD), shared)->notshared
		otureads_LOAD[c(shared, notshared),]->otureads_LOAD
	order(rowSums(metareads_F_LOAD), decreasing=FALSE)->ordering_vector
		metareads_F_LOAD[ordering_vector,]->metareads_F_LOAD
		match(rownames(metareads_F_LOAD), rownames(ITSreads_LOAD), nomatch=0)->shared	
		setdiff(1:nrow(ITSreads_LOAD), shared)->notshared
		ITSreads_LOAD[c(shared, notshared),]->ITSreads_LOAD



#order taxa for graphing
c("plant","sum_unique_metareads_B_cyanos","other","empty","low_abundance", "positive_blank_space", "negative_blank_space", "sum_unique_metareads_F",
	"sum_unique_metareads_B", "sum_unique_ITS","sum_unique_otus", 
	"sum_all_otu", "sum_all_ITS", "sum_all_metareads_B", "sum_all_metareads_F")->toporder
	if(addPLANT==TRUE){toporder[toporder!="plant"]->toporder}
c("Ceratobasidiaceae", "Sphingomonas", "Pseudomonas", "Sphingomonadaceae", "Pseudomonadaceae", "Sphingobacteriaceae", "Oxalobacteraceae","Hymenobacteraceae",
	"Flavobacteriaceae","Cytophagaceae","Comamonadaceae","Microbacteriaceae","Methylobacteriaceae")->bottomorder
	if(addPLANT==TRUE){c("plant", bottomorder)->bottomorder}
	bottomorder[length(bottomorder):1]->bottomorder

#c("sum_all_metareads_B")->bottomorder
topOrder(ITSreads_ReAb, toporder, bottomorder)->ITSreads_ReAb
topOrder(otureads_ReAb, toporder, bottomorder)->otureads_ReAb
topOrder(metareads_B_ReAb, toporder, bottomorder)->metareads_B_ReAb
topOrder(metareads_F_ReAb, toporder, bottomorder)->metareads_F_ReAb
topOrder(ITSreads_LOAD, toporder, bottomorder)->ITSreads_LOAD
topOrder(otureads_LOAD, toporder, bottomorder)->otureads_LOAD
topOrder(metareads_B_LOAD, toporder, bottomorder)->metareads_B_LOAD
topOrder(metareads_F_LOAD, toporder, bottomorder)->metareads_F_LOAD


order(metareads_B_LOAD["Pseudomonadaceae",], decreasing=TRUE)->ordervector
#order(metareads_B_LOAD["Proteobacteria",], decreasing=TRUE)->ordervector
#order(metareads_B_ReAb["Proteobacteria",], decreasing=TRUE)->ordervector
#order(metareads_F_ReAb["Ceratobasidiaceae",], decreasing=TRUE)->ordervector


#order(metareads_B_LOAD["Pseudomonadaceae",], decreasing=TRUE)->ordervector
#order(metareads_B_LOAD["sum_all_metareads_B",], decreasing=TRUE)->ordervector
#order(metareads_B_LOAD["other",], decreasing=TRUE)->ordervector


ITSreads_ReAb[,ordervector]->ITSreads_ReAb
otureads_ReAb[,ordervector]->otureads_ReAb
metareads_B_ReAb[,ordervector]->metareads_B_ReAb
metareads_F_ReAb[,ordervector]->metareads_F_ReAb
ITSreads_LOAD[,ordervector]->ITSreads_LOAD
otureads_LOAD[,ordervector]->otureads_LOAD
metareads_B_LOAD[,ordervector]->metareads_B_LOAD
metareads_F_LOAD[,ordervector]->metareads_F_LOAD



#CALCULATE HISTOGRAM REPRESENTING THE TAXA DIFFERENT BETWEEN LOAD AND RELATIVE ABUNDANCE METHODS
#HIGHER IN LOAD THAN IN REAB

find_difference_taxa_tables(metareads_B_ReAb, otureads_ReAb)->meta_B_ReAb_minus_otureads_ReAb
topOrder(meta_B_ReAb_minus_otureads_ReAb, toporder, bottomorder)->meta_B_ReAb_minus_otureads_ReAb
find_difference_taxa_tables(metareads_F_ReAb, ITSreads_ReAb)->meta_F_ReAb_minus_ITSreads_ReAb



library(reshape2)


currentreads=meta_B_ReAb_minus_otureads_ReAb
currentreads=metareads_B_ReAb
currentreads=otureads_ReAb
currentreads=meta_F_ReAb_minus_ITSreads_ReAb
currentreads=ITSreads_ReAb
#currentreads=otureads_ReAb
#currentreads=metareads_B_ReAb

#currentreads=meta_F_ReAb_minus_ITSreads_ReAb

currentreads=metareads_F_ReAb
currentreads=metareads_B_ReAb
#currentreads=metareads_B_LOAD
#currentreads=ITSreads_ReAb
#currentreads=meta_B_ReAb_minus_otureads_ReAb
#currentreads=metareads_B_ReAb
#currentreads=meta_F_ReAb_minus_ITSreads_ReAb

#metareads_B_LOAD->metareads_B_LOAD_PLANT
#colnames(currentreads)->columnorder

#metareads_B_LOAD_PLANT=currentreads
#metareads_B_LOAD_TOTAL[,match(columnorder,colnames(metareads_B_LOAD_TOTAL))]->metareads_B_LOAD_TOTAL
#metareads_B_LOAD_PLANT[,match(columnorder,colnames(metareads_B_LOAD_PLANT))]->metareads_B_LOAD_PLANT

#
	#
		#
			#
			currentreads=otureads_LOAD
			#
		#
	#
#


	#if(1==1){currentreads=currentreads[grep("sum_unique_metareads_B", rownames(currentreads), invert=TRUE),]}

#order taxa for graphing
c("empty", "mystery", "other", "plant", "allfungi")->toporder
c("HI")->bottomorder

#c("sum_all_metareads_B")->bottomorder
topOrder(currentreads, toporder, bottomorder)->currentreads


#normalize100(otureads_LOAD[rownames(otureads_LOAD)!="plant",])->otureads_LOAD

#currentreads=otureads_LOAD
#currentreads=metareads_B_ReAb

 
#colSums(currentreads[c(grep("eae", rownames(currentreads)), grep("unique", rownames(currentreads)),grep("allfungi", rownames(currentreads))),])->allmicrobes
#colSums(currentreads[c(grep("other", rownames(currentreads)), grep("plant", rownames(currentreads))),])->all_non_bacteria

currentreads[,order(currentreads["Pseudomonadaceae",], decreasing=TRUE)]->currentreads
currentreads[,order(currentreads["Pseudomonas",], decreasing=TRUE)]->currentreads
#currentreads[,order(allmicrobes, decreasing=TRUE)]->currentreads

COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
ALPHA_COLUMN_IN_COLORSLIST=7

if(MAKE_TOTAL_METABREAKDOWN_PLOT==TRUE){
	currentreads=metareads_B_LOAD
	c("empty", "mystery", "other", "plant", "allfungi")->toporder
	c("sum_unique_metareads_B")->bottomorder
	COLUMN_IN_COLORSLIST=4    # 3 is colors, 4 is black / white / green / purple
	colSums(currentreads[c(grep("eae", rownames(currentreads)), grep("unique", rownames(currentreads)),grep("allfungi", rownames(currentreads))),])->allmicrobes
	colSums(currentreads[c(grep("other", rownames(currentreads)), grep("plant", rownames(currentreads))),])->all_non_bacteria
	currentreads[,order(allmicrobes, decreasing=TRUE)]->currentreads
	colSums(currentreads[c(grep("eae", rownames(currentreads)), grep("unique", rownames(currentreads)),grep("allfungi", rownames(currentreads))),])->allmicrobes
	ALPHA_COLUMN_IN_COLORSLIST=8
}




if(MAKE_SUPPLEMENTARY_NORM==TRUE){
	currentreads=metareads_B_LOAD
	c("other","empty","low_abundance", "positive_blank_space", "negative_blank_space", "mystery", "allfungi")->toporder
	
	
	colSums(currentreads[c(grep("eae", rownames(currentreads)), grep("unique", rownames(currentreads)),grep("allfungi", rownames(currentreads))),])->allmicrobes
	colSums(currentreads[c(grep("other", rownames(currentreads)), grep("plant", rownames(currentreads))),])->all_non_bacteria
	currentreads[,order(allmicrobes, decreasing=TRUE)]->currentreads
	
	
	if(norm_fig_by_plant==TRUE){
		c("empty","other","low_abundance", "positive_blank_space", "negative_blank_space", "mystery", "allfungi")->toporder
		c("sum_unique_metareads_B_cyanos","sum_unique_metareads_B", "plant")->bottomorder	
		renorm_plant<-function(matrix){
			for(c in 1:ncol(matrix)){ matrix[,c]/matrix["plant",c]->matrix[,c] }
				return(matrix)
		}
		renorm_plant(currentreads)->currentreads
		add.empty.taxa(currentreads)->currentreads
		normalize100(currentreads)->currentreads
	}
	if(norm_fig_by_microbes==TRUE){
		currentreads=metareads_B_LOAD
		c("empty","other","low_abundance", "positive_blank_space", "negative_blank_space", "mystery", "plant", "allfungi")->toporder
		c("sum_unique_metareads_B_cyanos","sum_unique_metareads_B")->bottomorder	
		colSums(currentreads[c(grep("eae", rownames(currentreads)), grep("unique", rownames(currentreads)),grep("allfungi", rownames(currentreads))),])->allmicrobes
		
		rbind(currentreads[c("other","mystery", "plant"),], allmicrobes)->currentreads	
		renorm_microbes<-function(matrix){
			for(c in 1:ncol(matrix)){ matrix[,c]/allmicrobes[c]->matrix[,c] }
				return(matrix)
		}
		renorm_microbes(currentreads)->currentreads
		add.empty.taxa(currentreads)->currentreads
		normalize100(currentreads)->currentreads
	}
	
	COLUMN_IN_COLORSLIST=4    # 3 is colors, 4 is black / white / green / purple

}



date=format(Sys.Date(), format="%Y%m%d")
	

#currentreads[,order(currentreads["Pseudomonadaceae",], decreasing=TRUE)]->currentreads


topOrder(currentreads, toporder, bottomorder)->currentreads


counter=1
rep(0, ncol(currentreads))->currentreads["other",]
#APPLY CLR AND ADD OTHER


#counter+1->counter


melt(currentreads)->histDL
histDL[,c(2, 1, 3)]->histDL
cbind(rep(1, nrow(histDL)))->histDL$alphas
c("sample", "taxa", "abundance", "alphas")->names(histDL)


factor(histDL$taxa, levels=rownames(currentreads))->histDL$taxa
factor(histDL$sample, levels=colnames(currentreads))->histDL$sample


#SET COLOR SCHEME
rownames(currentreads)->taxa_color_pairs

source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"blue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Pseudomonas"),2]
	"red"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Sphingomonas"),2]

taxa_color_pairs[match(histDL$taxa, taxa_color_pairs[,1]),2]->histDL$colors
taxa_color_pairs[match(histDL$taxa, taxa_color_pairs[,1]),3]->histDL$alphas
as.numeric(histDL$alphas)->histDL$alphas



library(ggplot2)
library(reshape2)
library(RColorBrewer)
options(stringsAsFactors = FALSE)

#PLAY
#ggplot(histDL,aes(x=histDL$sample,y=histDL$abundance, fill=histDL$taxa)) + 

#ggplot(histDL,aes(x=histDL$sample,y=histDL$abundance, fill=histDL$taxa,alpha=histDL$alphas)) + 
ggplot(histDL,aes(x=histDL$sample,y=histDL$abundance, fill=histDL$taxa)) + 
  geom_bar(position="stack", stat="identity", width=1.0) + 
  scale_fill_manual(values=histDL$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("") + ylab("") + guides(fill=FALSE) +
  theme_classic() + theme(legend.position="none") + 
	scale_alpha_continuous(range = c(0.3, 1)) 
  

#FIGURE
#png("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/supp_normBACTERIA_LOAD.png", width = 3200, height = 2400, units = "px", res = 400)
pdf(paste("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/LOAD_ANIMATION_", date, counter, ".pdf", sep="", collapse=""), width = 5, height = 2, useDingbats=FALSE)
ggplot(histDL,aes(x=histDL$sample,y=histDL$abundance, fill=histDL$taxa,alpha=histDL$alphas)) + 
  geom_bar(position="fill", stat="identity", width=1.0) + 
 scale_fill_manual(values=histDL$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("") + ylab("") + guides(fill=FALSE) + scale_alpha_continuous(range = c(0.3, 1)) + 
	 theme_classic() + theme(legend.position="none")
  # theme_minimal() + 
dev.off()
png(paste("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/TOTALREADS_BREAKDOWN_", date, ".png", sep="", collapse=""), width = 3200, height = 2400, units = "px", res = 400)
ggplot(histDL,aes(x=histDL$sample,y=histDL$abundance, fill=histDL$taxa)) + 
  geom_bar(position="fill", stat="identity", width=1.0) + 
 scale_fill_manual(values=histDL$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("") + ylab("") + guides(fill=FALSE) + scale_alpha_continuous(range = c(0.3, 1)) + 
	 theme_classic() + theme(legend.position="none")
  # theme_minimal() + 
dev.off()

alpha=histDL$alphas



#LEGEND
par(mar = c(.2, .2, .2, .2)) 
png("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/Fig3_legend.png", width = 1000, height = 800, units = "px", res = 400)
plot(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), rep(10, 10), type="h",
	col=c("#3182BD","#de2d26","#a5d1c8","#9f5f2e","#d8853a","#ad84bc","#5f4298","#dd9e99","#6e9d35","#b2cde2"), 
	xlim=c(-1, 11), lwd=5)
dev.off()



#scatterplot for correlation 
#currentreadsX=metareads_B_ReAb
#currentreadsY=otureads_ReAb

currentreadsX=metareads_F_ReAb
currentreadsY=ITSreads_ReAb


#currentreadsX=metareads_B_LOAD[which(rownames(metareads_B_LOAD)!="other"),]
#currentreadsY=otureads_LOAD[which(rownames(otureads_LOAD)!="other"),]


#currentreadsX=matrix_format(metareads_F_ReAb[nrow(metareads_F_ReAb):(nrow(metareads_F_ReAb)-1),])
#currentreadsY=matrix_format(ITSreads_ReAb[nrow(ITSreads_ReAb):(nrow(ITSreads_ReAb)-1),])

#currentreadsX=metareads_B_ReAb[nrow(metareads_B_ReAb):(nrow(metareads_B_ReAb)-4),]
#currentreadsY=otureads_ReAb[nrow(otureads_ReAb):(nrow(otureads_ReAb)-4),]


intersect(rownames(currentreadsX), rownames(currentreadsY))->common
setdiff(rownames(currentreadsX), rownames(currentreadsY))->addtoY
setdiff(rownames(currentreadsY), rownames(currentreadsX))->addtoX
rbind(rep(0, ncol(currentreadsX)), currentreadsX)->currentreadsX
rbind(rep(0, ncol(currentreadsY)), currentreadsY)->currentreadsY
rownames(currentreadsX)[1]<-addtoX
rownames(currentreadsY)[1]<-addtoY

c("plant","sum_unique_metareads_B_cyanos","other","empty","low_abundance", "positive_blank_space", "negative_blank_space", "sum_unique_metareads_F",
	"sum_unique_metareads_B", "sum_unique_ITS","sum_unique_otus", 
	"sum_all_otu", "sum_all_ITS", "sum_all_metareads_B", "sum_all_metareads_F")->toporder
topOrder(currentreadsX, toporder, bottomorder)->currentreadsX
topOrder(currentreadsY, toporder, bottomorder)->currentreadsY

rownames(currentreadsX)->taxa_color_pairs

melt(currentreadsX)->histDL1
histDL1[,c(2, 1, 3)]->histDL1
c("sample1", "taxa1", "abundance1")->names(histDL1)
factor(histDL1$taxa1, levels=rownames(currentreadsX))->histDL1$taxa1

melt(currentreadsY)->histDL2
histDL2[,c(2, 1, 3)]->histDL2
c("sample2", "taxa2", "abundance2")->names(histDL2)
factor(histDL2$taxa2, levels=rownames(currentreadsY))->histDL2$taxa2

cbind(histDL1, histDL2)->histDL3


#ONLY EXTRACT THE "Ceratobasidiaceae"
histDL3[grep("Ceratobasidiaceae", histDL3[,2]),]->histDL3


COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
ALPHA_COLUMN_IN_COLORSLIST=7
source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
taxa_color_pairs[match(histDL3$taxa1, taxa_color_pairs[,1]),2]->histDL3$colors
taxa_color_pairs[match(histDL3$taxa1, taxa_color_pairs[,1]),3]->histDL3$alphas
as.numeric(histDL3$alphas)->histDL3$alphas

correlation0=round(cor(histDL3$abundance1, histDL3$abundance2, method="pearson"), 4)^2
correlation=round(cor(sqrt(histDL3$abundance1), sqrt(histDL3$abundance2), method="pearson"), 4)^2
correlation4=round(cor(sqrt(sqrt(histDL3$abundance1)), sqrt(sqrt(histDL3$abundance2)), method="pearson"), 4)^2
round(cor(sqrt(histDL3$abundance1), sqrt(histDL3$abundance2), method="pearson"), 4)^2
round(cor(histDL3$abundance1, histDL3$abundance2, method="pearson"), 4)^2



#PLAY VERSION
ggplot(histDL3, aes(x=abs(sqrt(abundance1)), y=abs(sqrt(abundance2)), color=colors, alpha=alphas)) + geom_point(size=2) + 
	theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
	scale_color_manual(values=DerekPalette) + theme(legend.position="none") + 
	coord_fixed(ratio=1/1) +  
	scale_x_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)^2, limits = c(0,10)) + 
	scale_y_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)^2, limits = c(0,10)) + 
	annotate(geom="text", size=5, x = 4, y = 9, label = deparse(bquote(R^2~"="~.(correlation))), parse=T)

ggplot(histDL3, aes(x=abs(sqrt(abundance1)), y=abs(sqrt(abundance2)), color=taxa1)) + geom_point(alpha = 0.5, size=2) + 
	theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
	#scale_color_manual(values=DerekPalette) + theme(legend.position="none") + 
	scale_color_manual(values=histDL3$colors) + theme(legend.position="none") + 
	coord_fixed(ratio=1/1) +  
	scale_x_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)^2, limits = c(0,10)) + 
	scale_y_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)^2, limits = c(0,10)) + 
	annotate(geom="text", size=5, x = 4, y = 5, label = deparse(bquote(R^2~"="~.(correlation))), parse=T)


#FIGURE VERSION
#png("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/Correlation_scatter_OTUvsMetaB_ReAb_TOP5.png", width = 3200, height = 3200, units = "px", res = 400)
pdf(paste("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/Correlation_scatter_ITSvsMetaF_ReAb_2root_", date, ".pdf", sep="", collapse=""), width = 5, height = 2, useDingbats=FALSE)
ggplot(histDL3, aes(x=abs(sqrt(histDL3$abundance1)), y=abs(sqrt(histDL3$abundance2)), color=histDL3$taxa1, alpha=histDL3$alphas)) + geom_point(shape=16, size=.5) + 
	theme_minimal() + theme(axis.text.x = element_text(size=4),axis.text.y = element_text(size=4)) +  
	scale_color_manual(values=histDL3$colors) + theme(legend.position="none") + 
	coord_fixed(ratio=1/1) +  
	scale_x_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)^2, limits = c(0,10)) + 
	scale_y_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)^2, limits = c(0,10)) + 
	annotate(geom="text", size=3, x = 4, y = 9, label = deparse(bquote(R^2~"="~.(correlation))), parse=T) + 
	theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
	theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
	theme(panel.grid.minor = element_blank())
dev.off()
#FIG VERSION
pdf(paste("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/Correlation_scatter_ITSvsMetaF_ReAb_4root_", date, ".pdf", sep="", collapse=""), width = 5, height = 2, useDingbats=FALSE)
ax=c(0, 0.5, 5, 20, 50, 100)
ggplot(histDL3, aes(x=abs(sqrt(sqrt(histDL3$abundance1))), y=abs(sqrt(sqrt(histDL3$abundance2))), color=histDL3$taxa1, alpha=histDL3$alphas)) + geom_point(shape=16, size=.5) + 
	theme_minimal() + theme(axis.text.x = element_text(size=4),axis.text.y = element_text(size=4)) +  
	scale_color_manual(values=histDL3$colors) + theme(legend.position="none") + 
	coord_fixed(ratio=1/1) +  
	scale_x_continuous(name = "", breaks = sqrt(sqrt(ax)), labels = ax, limits = c(0,sqrt(sqrt(100))))+ 
	scale_y_continuous(name = "", breaks = sqrt(sqrt(ax)), labels = ax, limits = c(0,sqrt(sqrt(100))))+ 
	annotate(geom="text", size=3, x = 4, y = 9, label = deparse(bquote(R^2~"="~.(correlation4))), parse=T) + 
	theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
	theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
	theme(panel.grid.minor = element_blank())
dev.off()


#png("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/Correlation_scatter_OTUvsMetaB_ReAb_TOP5.png", width = 3200, height = 3200, units = "px", res = 400)
pdf(paste("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/Correlation_scatter_ITSvsMetaF_ReAb_notSQRT_", date, ".pdf", sep="", collapse=""), width = 5, height = 2, useDingbats=FALSE)
ggplot(histDL3, aes(x=abs(histDL3$abundance1), y=abs(histDL3$abundance2), color=histDL3$taxa1, alpha=histDL3$alphas)) + geom_point(shape=16, size=.5) + 
	theme_minimal() + theme(axis.text.x = element_text(size=4),axis.text.y = element_text(size=4)) +  
	scale_color_manual(values=histDL3$colors) + theme(legend.position="none") + 
	coord_fixed(ratio=1/1) +  
	scale_x_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10)*10, labels = c(0, 2, 4, 6, 8, 10)*10, limits = c(0,100)) + 
	scale_y_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10)*10, labels = c(0, 2, 4, 6, 8, 10)*10, limits = c(0,100)) + 
	annotate(geom="text", size=3, x = 4, y = 9, label = deparse(bquote(R^2~"="~.(correlation))), parse=T) + 
	theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
	theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
	theme(panel.grid.minor = element_blank())
dev.off()
       


###########################################
########################################### Pseudo vs Sphingo
currentreadsX=otureads_LOAD[which(rownames(otureads_LOAD)!="other"),]
currentreadsY=otureads_LOAD[which(rownames(otureads_LOAD)!="other"),]

#currentreadsX=otureads_ReAb[which(rownames(otureads_ReAb)!="other"),]
#currentreadsY=otureads_ReAb[which(rownames(otureads_ReAb)!="other"),]


currentreadsX=metareads_B_LOAD[which(rownames(metareads_B_LOAD)!="other"),]
currentreadsY=metareads_B_LOAD[which(rownames(metareads_B_LOAD)!="other"),]

#currentreadsX=metareads_B_ReAb[which(rownames(metareads_B_ReAb)!="other"),]
#currentreadsY=metareads_B_ReAb[which(rownames(metareads_B_ReAb)!="other"),]

library(compositions)

CLR=t(clr(t(metareads_B_ReAb[which(rownames(metareads_B_ReAb)!="other"),])))
currentreadsX=CLR[1:nrow(CLR),]
currentreadsY=CLR[1:nrow(CLR),]


melt(currentreadsX)->histDL1
histDL1[,c(2, 1, 3)]->histDL1
c("sample1", "taxa1", "abundance1")->names(histDL1)
factor(histDL1$taxa1, levels=rownames(currentreadsX))->histDL1$taxa1

melt(currentreadsY)->histDL2
histDL2[,c(2, 1, 3)]->histDL2
c("sample2", "taxa2", "abundance2")->names(histDL2)
factor(histDL2$taxa2, levels=rownames(currentreadsY))->histDL2$taxa2

histDL1[histDL1[,"taxa1"]=="Pseudomonadaceae",]->histDL1
histDL2[histDL2[,"taxa2"]=="Sphingomonadaceae",]->histDL2

cbind(histDL1, histDL2)->histDL3



rownames(currentreadsX)->taxa_color_pairs
COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
ALPHA_COLUMN_IN_COLORSLIST=7
source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
taxa_color_pairs[match(histDL3$taxa1, taxa_color_pairs[,1]),2]->histDL3$colors
taxa_color_pairs[match(histDL3$taxa1, taxa_color_pairs[,1]),3]->histDL3$alphas
as.numeric(histDL3$alphas)->histDL3$alphas

correlation=round(cor(sqrt(histDL3$abundance1), sqrt(histDL3$abundance2), method="pearson"), 2)^2
correlation4=round(cor(sqrt(sqrt(histDL3$abundance1)), sqrt(sqrt(histDL3$abundance2)), method="pearson"), 2)^2

round(cor(histDL3$abundance1, histDL3$abundance2, method="pearson"), 2)^2

plot(histDL3$abundance1, histDL3$abundance2)

library(compositions)
round(cor(histDL3$abundance1, histDL3$abundance2, method="pearson"), 2)^2


#PLAY VERSION
ggplot(histDL3, aes(x=abs(sqrt(abundance1)), y=abs(sqrt(abundance2)))) + geom_point(alpha = 0.5, size=2) + 
	theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
	scale_color_manual(values=DerekPalette) + theme(legend.position="none") + 
	coord_fixed(ratio=1/1) +  
	scale_x_continuous(name = "", breaks = c(0, 2, 4), labels = c(0, 2, 4)^2, limits = c(0,10)) + 
	scale_y_continuous(name = "", breaks = c(0, 2, 4), labels = c(0, 2, 4)^2, limits = c(0,10)) + 
	annotate(geom="text", size=5, x = 3, y = 3, label = deparse(bquote(R^2~"="~.(correlation))), parse=T)

#fig version scatter
pdf("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/OTUreads_PseudoVsSphingo_LOAD.pdf", width = 5, height = 2, useDingbats=FALSE)
ggplot(histDL3, aes(x=abs(sqrt(histDL3$abundance1)), y=abs(sqrt(histDL3$abundance2)), color="black", alpha=histDL3$alphas)) + geom_point(shape=16, size=.5, color="black") + 
	theme_classic() + theme(axis.text.x = element_text(size=4),axis.text.y = element_text(size=4)) +  
	scale_color_manual(values=histDL3$colors) + theme(legend.position="none") + 
	coord_fixed(ratio=1/1) +  
	scale_x_continuous(name = "", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)^2, limits = c(0,8)) + 
	scale_y_continuous(name = "", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)^2, limits = c(0,8)) + 
	annotate(geom="text", size=3, x = 4, y = 3, label = deparse(bquote(R^2~"="~.(correlation))), parse=T) + 
	theme(axis.line = element_line(colour = "black", size = .15), axis.ticks = element_line(colour = "black", size = .15))
dev.off()
ax=c(0, 0.5, 5, 20, 50, 100)
pdf("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/OTUreads_PseudoVsSphingo_LOAD_4root.pdf", width = 5, height = 2, useDingbats=FALSE)
ggplot(histDL3, aes(x=abs(sqrt(sqrt(histDL3$abundance1))), y=abs(sqrt(sqrt(histDL3$abundance2))), color="black", alpha=histDL3$alphas)) + geom_point(shape=16, size=.5, color="black") + 
	theme_minimal() + theme(axis.text.x = element_text(size=4),axis.text.y = element_text(size=4)) +  
	scale_color_manual(values=histDL3$colors) + theme(legend.position="none") + 
	coord_fixed(ratio=1/1) +  
	scale_x_continuous(name = "", breaks = sqrt(sqrt(ax)), labels = ax, limits = c(0,sqrt(sqrt(100))))+ 
	scale_y_continuous(name = "", breaks = sqrt(sqrt(ax)), labels = ax, limits = c(0,sqrt(sqrt(100))))+ 
	annotate(geom="text", size=3, x = 4, y = 9, label = deparse(bquote(R^2~"="~.(correlation4))), parse=T) + 
	theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
	theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
	theme(panel.grid.minor = element_blank())
dev.off()
	
myplot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))



#PLAY VERSION BOXPLOT
names(histDL1)->names(histDL2)
rbind(histDL1, histDL2)->histDL3
rownames(currentreadsX)->taxa_color_pairs
source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
taxa_color_pairs[match(histDL3$taxa1, taxa_color_pairs[,1]),2]->DerekPalette
as.character(histDL3[,1])->histDL3[,1]
ggplot(histDL3, aes(taxa1, sqrt(abundance1), colour = taxa1)) + theme_minimal() + 	geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + 
	scale_color_manual(values=c("#de2e26", "#3182BD")) + 
	scale_y_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)^2, limits = c(0,10)) 
	
	geom_boxplot() + 
	
	boxplot(sqrt(ReAb["Pseudomonadaceae",greater]), sqrt(ReAb["Pseudomonadaceae",less]),
		sqrt(ReAb["Sphingomonadaceae",greater]), sqrt(ReAb["Sphingomonadaceae",less]), outline=FALSE, ylim=c(0, 10))
	points(jitter(rep(1, length(ReAb["Pseudomonadaceae",greater])), amount=.2),sqrt(ReAb["Pseudomonadaceae",greater]),bg="#000000", pch=21)
	points(jitter(rep(2, length(ReAb["Pseudomonadaceae",greater])), amount=.2),sqrt(ReAb["Pseudomonadaceae",less]), bg="#000000", pch=21)
	points(jitter(rep(3, length(ReAb["Sphingomonadaceae",greater])), amount=.2),sqrt(ReAb["Sphingomonadaceae",greater]),bg="#ffff00", pch=21)
	points(jitter(rep(4, length(ReAb["Sphingomonadaceae",greater])), amount=.2),sqrt(ReAb["Sphingomonadaceae",less]), bg="#ffff00", pch=21)
	axis(side=1, at=c(1, 2, 3, 4), labels=c("Ps > median", "Ps < median", "Sp >", "Sp <"))


	
png("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/Correlation_scatter_pseudoVSsphingoReAb.png", width = 3200, height = 3200, units = "px", res = 400)
ggplot(histDL3, aes(x=abs(sqrt(abundance1)), y=abs(sqrt(abundance2)))) + geom_point(alpha = 0.5, size=4) + 
	theme_minimal() + theme(axis.text.x = element_text(size=30),axis.text.y = element_text(size=30)) +  
	scale_color_manual(values=DerekPalette) + theme(legend.position="none") + 
	coord_fixed(ratio=1/1) +  
	scale_x_continuous(name = "", breaks = c(0, 2, 4), labels = c(0, 2, 4)^2, limits = c(0,4)) + 
	scale_y_continuous(name = "", breaks = c(0, 2, 4), labels = c(0, 2, 4)^2, limits = c(0,4)) + 
	annotate(geom="text", size=15, x = 2, y = 3, label = deparse(bquote(R^2~"="~.(correlation))), parse=T)
dev.off()

png("/Users/dlundberg/Documents/abt6/Regalado2017/limits_to_adjust_16S_by_metagenome/Correlation_scatter_pseudoVSsphingoReAb.png", width = 3200, height = 3200, units = "px", res = 400)
ggplot(histDL3, aes(x=abs(sqrt(abundance1)), y=abs(sqrt(abundance2)))) + geom_point(alpha = 0.5, size=4) + 
	theme_minimal() + theme(axis.text.x = element_text(size=30),axis.text.y = element_text(size=30)) +  
	scale_color_manual(values=DerekPalette) + theme(legend.position="none") + 
	coord_fixed(ratio=1/1) +  
	scale_x_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4,6,8,10)^2, limits = c(0,10)) + 
	scale_y_continuous(name = "", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4,6,8,10)^2, limits = c(0,10)) + 
	annotate(geom="text", size=15, x = 5, y = 8, label = deparse(bquote(R^2~"="~.(correlation))), parse=T)
dev.off()
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################


#compare Sphingomonas vs. Pseudomonas for relative abundance vs load

metareads_B_ReAb[c("Pseudomonadaceae","Sphingomonadaceae"),]->ReAb
metareads_B_LOAD[c("Pseudomonadaceae","Sphingomonadaceae"),]->LOAD

otureads_ReAb[c("Pseudomonadaceae","Sphingomonadaceae"),]->ReAb
otureads_LOAD[c("Pseudomonadaceae","Sphingomonadaceae"),]->LOAD











par(mfrow=c(1,2))
margins=par(mar=c(5.1, 4.1, 4.1, 2.1))
	
	plot(sqrt(metareads_B_ReAb["Pseudomonadaceae",]), sqrt(metareads_B_ReAb["Sphingomonadaceae",]), xlab="Pseudo", ylab="Sphingo", xlim=c(0,10), ylim=c(0,10), xaxt="n", yaxt="n")
		axis(side=1, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		
	plot(sqrt(metareads_B_LOAD["Pseudomonadaceae",]), sqrt(metareads_B_LOAD["Sphingomonadaceae",]), xlab="Pseudo", ylab="Sphingo", xlim=c(0,10), ylim=c(0,10), xaxt="n", yaxt="n")
		axis(side=1, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		#cor(sqrt(metareads_B_LOAD["Pseudomonadaceae",]), sqrt(metareads_B_LOAD["Sphingomonadaceae",]))
		#cor.test(sqrt(metareads_B_LOAD["Pseudomonadaceae",]), sqrt(metareads_B_LOAD["Sphingomonadaceae",]),alternative="greater")


par(mfrow=c(1,2))
margins=par(mar=c(5.1, 4.1, 4.1, 2.1))

	plot(sqrt(otureads_ReAb["Pseudomonadaceae",]), sqrt(otureads_ReAb["Sphingomonadaceae",]), main="Relative Abundance", xlab="Pseudo", ylab="Sphingo", xlim=c(0,10), ylim=c(0,10), xaxt="n", yaxt="n")
		axis(side=1, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		ReAbcor<-cor.test(sqrt(otureads_ReAb["Pseudomonadaceae",]), sqrt(otureads_ReAb["Sphingomonadaceae",]),alternative="greater")
		text(7, 7, round(ReAbcor[[4]],4))

	plot(sqrt(otureads_LOAD["Pseudomonadaceae",]), sqrt(otureads_LOAD["Sphingomonadaceae",]), main="LOAD", xlab="Pseudo", ylab="Sphingo", xlim=c(0,10), ylim=c(0,10), xaxt="n", yaxt="n")
		axis(side=1, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		#cor(sqrt(metareads_B_LOAD["Pseudomonadaceae",]), sqrt(metareads_B_LOAD["Sphingomonadaceae",]))
		LOADcor<-cor.test(sqrt(otureads_LOAD["Pseudomonadaceae",]), sqrt(otureads_LOAD["Sphingomonadaceae",]),alternative="greater")
		text(7, 7, round(LOADcor[[4]],4))

	plot(sqrt(metareads_B2["Pseudomonadaceae",]), sqrt(metareads_B2["Sphingomonadaceae",]), main="LOAD", xlab="Pseudo", ylab="Sphingo", xaxt="n", yaxt="n")
		axis(side=1, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		#cor(sqrt(metareads_B_LOAD["Pseudomonadaceae",]), sqrt(metareads_B_LOAD["Sphingomonadaceae",]))
		LOADcor<-cor.test(sqrt(metareads_B2["Pseudomonadaceae",]), sqrt(metareads_B2["Sphingomonadaceae",]),alternative="greater")
		text(7, 7, round(LOADcor[[4]],4))


metareads_B2

## 23 Nov 2018
par(mfrow=c(2,2))
margins=par(mar=c(2, 2, 2, 2))
	plot(1:ncol(otureads_ReAb), otureads_ReAb["Pseudomonadaceae",], type="h", xaxt="n", yaxt="n", xlab="", ylab="", main="Pseudomonadaceae")
	plot(1:ncol(otureads_ReAb), otureads_ReAb["Sphingomonadaceae",], type="h", xaxt="n", yaxt="n", xlab="", ylab="", main="Sphingomonadaceae")
	
	1:ncol(otureads_ReAb)->xaxis
	otureads_ReAb["Pseudomonadaceae",]->yaxis
		plot(xaxis, yaxis, type="p", xaxt="n", yaxt="n", xlab="", ylab="")
		abline(lm(yaxis~xaxis), col="red") # regression line (y~x) 
		
	1:ncol(otureads_ReAb)->xaxis
	otureads_ReAb["Sphingomonadaceae",]->yaxis
		plot(xaxis, yaxis, type="p", xaxt="n", yaxt="n", xlab="",ylab="")
				abline(lm(yaxis~xaxis), col="red") # regression line (y~x) 
		
par(mfrow=c(2,2))
margins=par(mar=c(2, 2, 2, 2))
	plot(1:ncol(otureads_LOAD), otureads_LOAD["Pseudomonadaceae",], type="h", xaxt="n", yaxt="n", xlab="", ylab="", main="Pseudomonadaceae")
	plot(1:ncol(otureads_LOAD), otureads_LOAD["Sphingomonadaceae",], type="h", xaxt="n", yaxt="n", xlab="", ylab="", main="Sphingomonadaceae")
	
	1:ncol(otureads_LOAD)->xaxis
	otureads_LOAD["Pseudomonadaceae",]->yaxis
		plot(xaxis, yaxis, type="p", xaxt="n", yaxt="n", xlab="", ylab="")
		abline(lm(yaxis~xaxis), col="red") # regression line (y~x) 
		
	1:ncol(otureads_LOAD)->xaxis
	otureads_LOAD["Sphingomonadaceae",]->yaxis
		plot(xaxis, yaxis, type="p", xaxt="n", yaxt="n", xlab="",ylab="")
				abline(lm(yaxis~xaxis), col="red") # regression line (y~x) 
	



which(LOAD["Pseudomonadaceae",]>median(LOAD["Pseudomonadaceae",]))->greaterLOAD
which(LOAD["Pseudomonadaceae",]<median(LOAD["Pseudomonadaceae",]))->lessLOAD
which(ReAb["Pseudomonadaceae",]>median(ReAb["Pseudomonadaceae",]))->greaterReAb
which(ReAb["Pseudomonadaceae",]<median(ReAb["Pseudomonadaceae",]))->lessReAb

wilcox.test(sqrt(LOAD["Pseudomonadaceae",greaterLOAD]), sqrt(LOAD["Pseudomonadaceae",lessLOAD]), alternative = c("greater"))
wilcox.test(sqrt(LOAD["Sphingomonadaceae",greaterLOAD]), sqrt(LOAD["Sphingomonadaceae",lessLOAD]), alternative = c("less"))
wilcox.test(sqrt(ReAb["Pseudomonadaceae",greaterLOAD]), sqrt(ReAb["Pseudomonadaceae",lessLOAD]), alternative = c("greater"))
wilcox.test(sqrt(ReAb["Sphingomonadaceae",greaterLOAD]), sqrt(ReAb["Sphingomonadaceae",lessLOAD]), alternative = c("less"))

t.test(sqrt(LOAD["Pseudomonadaceae",greaterLOAD]), sqrt(LOAD["Pseudomonadaceae",lessLOAD]), alternative = c("greater"))
t.test(sqrt(LOAD["Sphingomonadaceae",greaterLOAD]), sqrt(LOAD["Sphingomonadaceae",lessLOAD]), alternative = c("less"))
t.test(sqrt(ReAb["Pseudomonadaceae",greaterReAb]), sqrt(ReAb["Pseudomonadaceae",lessReAb]), alternative = c("greater"))
t.test(sqrt(ReAb["Sphingomonadaceae",greaterReAb]), sqrt(ReAb["Sphingomonadaceae",lessReAb]), alternative = c("less"))





boxplot(sqrt(ReAb["Pseudomonadaceae",greater]), sqrt(ReAb["Pseudomonadaceae",less]),
		sqrt(ReAb["Sphingomonadaceae",greater]), sqrt(ReAb["Sphingomonadaceae",less]), outline=FALSE, ylim=c(0, 10))
	points(jitter(rep(1, length(ReAb["Pseudomonadaceae",greater])), amount=.2),sqrt(ReAb["Pseudomonadaceae",greater]),bg="#000000", pch=21)
	points(jitter(rep(2, length(ReAb["Pseudomonadaceae",greater])), amount=.2),sqrt(ReAb["Pseudomonadaceae",less]), bg="#000000", pch=21)
	points(jitter(rep(3, length(ReAb["Sphingomonadaceae",greater])), amount=.2),sqrt(ReAb["Sphingomonadaceae",greater]),bg="#ffff00", pch=21)
	points(jitter(rep(4, length(ReAb["Sphingomonadaceae",greater])), amount=.2),sqrt(ReAb["Sphingomonadaceae",less]), bg="#ffff00", pch=21)
	axis(side=1, at=c(1, 2, 3, 4), labels=c("Ps > median", "Ps < median", "Sp >", "Sp <"))


	
par(mfrow=c(1,2))
margins=par(mar=c(5,5,1,1))
plot(sqrt(metareads_B_LOAD["Pseudomonadaceae",]), sqrt(metareads_B_LOAD["Sphingomonadaceae",]), 
		 xlab="Pseudo", ylab="Sphingo", xlim=c(0,10), ylim=c(0,10), xaxt="n", yaxt="n")
		axis(side=1, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		abline(lm(sqrt(metareads_B_LOAD["Sphingomonadaceae",])~sqrt(metareads_B_LOAD["Pseudomonadaceae",])), col="red")

#par(mfrow=c(1,1))	
margins=par(mar=c(5,5,1,1))
plot(sqrt(metareads_B_ReAb["Pseudomonadaceae",]), sqrt(metareads_B_ReAb["Sphingomonadaceae",]), 
		 xlab="Pseudo", ylab="Sphingo", xlim=c(0,10), ylim=c(0,10), xaxt="n", yaxt="n")
		axis(side=1, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		abline(lm(sqrt(metareads_B_ReAb["Sphingomonadaceae",])~sqrt(metareads_B_ReAb["Pseudomonadaceae",])), col="red")
		
		

#16S ADJUSTED VERSION	
par(mfrow=c(1,2))
margins=par(mar=c(5,5,1,1))
plot(sqrt(otureads_LOAD["Pseudomonadaceae",]), sqrt(otureads_LOAD["Sphingomonadaceae",]), 
		 xlab="Pseudo", ylab="Sphingo", xlim=c(0,10), ylim=c(0,10), xaxt="n", yaxt="n")
		axis(side=1, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		abline(lm(sqrt(metareads_B_LOAD["Sphingomonadaceae",])~sqrt(metareads_B_LOAD["Pseudomonadaceae",])), col="red")

#par(mfrow=c(1,1))	
margins=par(mar=c(5,5,1,1))
plot(sqrt(otureads_ReAb["Pseudomonadaceae",]), sqrt(otureads_ReAb["Sphingomonadaceae",]), 
		 xlab="Pseudo", ylab="Sphingo", xlim=c(0,10), ylim=c(0,10), xaxt="n", yaxt="n")
		axis(side=1, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
		abline(lm(sqrt(metareads_B_ReAb["Sphingomonadaceae",])~sqrt(metareads_B_ReAb["Pseudomonadaceae",])), col="red")
		
				
		

		
		
margins=par(mar=c(0,3,1,1))
boxplot(sqrt(LOAD["Sphingomonadaceae",greaterLOAD]), sqrt(LOAD["Sphingomonadaceae",lessLOAD]), outline=FALSE, ylim=c(0, 4), xaxt="n", yaxt="n")
	points(jitter(rep(1, length(LOAD["Sphingomonadaceae",greaterLOAD])), amount=.2),sqrt(LOAD["Sphingomonadaceae",greaterLOAD]),bg=alpha("#ffff00", 0.5), pch=21)
	points(jitter(rep(2, length(LOAD["Sphingomonadaceae",greaterLOAD])), amount=.2),sqrt(LOAD["Sphingomonadaceae",lessLOAD]), bg=alpha("#ffff00",0.5), pch=21)
	#axis(side=1, at=c(1, 2), labels=c("> median Pseudomonas", "< median Psuedomonas"))
	axis(side=2, at=c(0, 2, 4), labels=c(0, 2, 4)^2, las=2)

margins=par(mar=c(0,3,1,1))
boxplot(sqrt(ReAb["Sphingomonadaceae",greaterReAb]), sqrt(ReAb["Sphingomonadaceae",lessReAb]), outline=FALSE, ylim=c(0, 8),xaxt="n", yaxt="n")	
	points(jitter(rep(1, length(ReAb["Sphingomonadaceae",greaterReAb])), amount=.2),sqrt(ReAb["Sphingomonadaceae",greaterReAb]),bg=alpha("#ffff00", 0.5), pch=21)
	points(jitter(rep(2, length(ReAb["Sphingomonadaceae",greaterReAb])), amount=.2),sqrt(ReAb["Sphingomonadaceae",lessReAb]), bg=alpha("#ffff00",0.5), pch=21)
	#axis(side=1, at=c(1, 2), labels=c("> median Pseudomonas", "< median Psuedomonas"))
	axis(side=2, at=c(0, 4, 8), labels=c(0, 4, 8)^2, las=2)

margins=par(mar=c(1,3,0,1))
boxplot(sqrt(LOAD["Pseudomonadaceae",greaterLOAD]), sqrt(LOAD["Pseudomonadaceae",lessLOAD]), outline=FALSE, ylim=c(0, 10), xaxt="n",yaxt="n")
	points(jitter(rep(1, length(LOAD["Pseudomonadaceae",greaterLOAD])), amount=.2),sqrt(LOAD["Pseudomonadaceae",greaterLOAD]),bg=alpha("#000000", 0.5), pch=21)
	points(jitter(rep(2, length(LOAD["Pseudomonadaceae",greaterLOAD])), amount=.2),sqrt(LOAD["Pseudomonadaceae",lessLOAD]), bg=alpha("#000000", 0.5), pch=21)
	#axis(side=1, at=c(1, 2), labels=c("> median", "< median"))
	axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
	abline(h=sqrt(median(LOAD["Pseudomonadaceae",])), col="red")
	
margins=par(mar=c(1,3,0,1))
boxplot(sqrt(ReAb["Pseudomonadaceae",greaterReAb]), sqrt(ReAb["Pseudomonadaceae",lessReAb]), outline=FALSE, ylim=c(0, 10),xaxt="n", yaxt="n")
	points(jitter(rep(1, length(ReAb["Pseudomonadaceae",greaterReAb])), amount=.2),sqrt(ReAb["Pseudomonadaceae",greaterReAb]),bg=alpha("#000000", 0.5), pch=21)
	points(jitter(rep(2, length(ReAb["Pseudomonadaceae",greaterReAb])), amount=.2),sqrt(ReAb["Pseudomonadaceae",lessReAb]), bg=alpha("#000000", 0.5), pch=21)
	#axis(side=1, at=c(1, 2), labels=c("> median", "< median"))
	axis(side=2, at=c(0, 5, 10), labels=c(0, 5, 10)^2, las=2)
	abline(h=sqrt(median(ReAb["Pseudomonadaceae",])), col="red")

	






	#first look only at shared taxa
intersect(meta_names, otu_names)->common_taxa
metareads[match(common_taxa, meta_names),]->metareads_common
otureads[match(common_taxa, otu_names),]->otureads_common
	#sqrt transform the tables
#log10(1000*metareads_common+1)->metareads_common
#log10(1000*otureads_common+1)->otureads_common
sqrt(sqrt(metareads_common))->metareads_common
sqrt(sqrt(otureads_common))->otureads_common
	#plot samples against each other
plot(0, 0, xlim=c(0, max(otureads_common)), ylim=c(0, max(metareads_common)))
for(c in 1:ncol(otureads_common)){
	points(otureads_common[,c], metareads_common[,c])
}


	#next look at taxa unique to the metagenome, except all bacteria and total seq)
setdiff(meta_names, otu_names)->meta_unique_taxa
setdiff(meta_unique_taxa, c("Bacteria", "TotalSeq"))->meta_unique_taxa
metareads[match(meta_unique_taxa, meta_names),]->metareads_unique
sqrt(sqrt(metareads_unique))->metareads_unique
	#plot samples in new color
for(c in 1:ncol(metareads_unique)){
	points(rep(0, length=length(metareads_unique[,c])),metareads_unique[,c], col="red")
}



	#next look at taxa unique to the 16S
setdiff(otu_names, meta_names)->otu_unique_taxa
otureads[match(otu_unique_taxa, otu_names),]->otureads_unique
sqrt(sqrt(otureads_unique))->otureads_unique
	#plot samples in new color
for(c in 1:ncol(otureads_unique)){
	points(otureads_unique[,c],rep(0, length=length(otureads_unique[,c])), col="blue")
}


	#make correlations by sample. test metagenomes vs otu, same sample
matrix(nrow=ncol(otureads_common),ncol=4)->correlations
rep("metaVotu", ncol(otureads_common))->correlations[,4]
for(c in 1:ncol(otureads_common)){
	sample(1:ncol(otureads_common))[c]->d
	#print(c(c, d))
	common_samples[c]->correlations[c,1]
	common_samples[c]->correlations[c,2]
	cor(metareads_common[,c], otureads_common[,c], method="pearson")->correlations[c,3]
}
correlations->meta_vs_otu_sample_correlations
	#make correlations by sample. Test metagenome vs metagenome, random sample
matrix(nrow=ncol(otureads_common),ncol=4)->correlations
rep("metaVmetaR", ncol(otureads_common))->correlations[,4]
for(c in 1:ncol(otureads_common)){
	sample(1:ncol(otureads_common))[c]->d
	print(c(c, d))
	common_samples[c]->correlations[c,1]
	common_samples[d]->correlations[c,2]
	cor(metareads_common[,c], metareads_common[,d], method="pearson")->correlations[c,3]
}
correlations->meta_vs_meta_sample_correlations
	#make correlations by sample. Test otu vs otu, random sample
matrix(nrow=ncol(otureads_common),ncol=4)->correlations
rep("otuVotuR", ncol(otureads_common))->correlations[,4]
for(c in 1:ncol(otureads_common)){
	sample(1:ncol(otureads_common))[c]->d
	print(c(c, d))
	common_samples[c]->correlations[c,1]
	common_samples[d]->correlations[c,2]
	cor(otureads_common[,c], otureads_common[,d], method="pearson")->correlations[c,3]
}
correlations->otu_vs_otu_sample_correlations
	#make correlations by sample. Test meta vs otu, random sample
matrix(nrow=ncol(otureads_common),ncol=4)->correlations
rep("metaVotuR", ncol(otureads_common))->correlations[,4]
for(c in 1:ncol(otureads_common)){
	sample(1:ncol(otureads_common))[c]->d
	print(c(c, d))
	common_samples[c]->correlations[c,1]
	common_samples[d]->correlations[c,2]
	cor(metareads_common[,c], otureads_common[,d], method="pearson")->correlations[c,3]
}
correlations->meta_vs_otu_random_correlations

pearson_sample_correlations=rbind(meta_vs_otu_sample_correlations, meta_vs_meta_sample_correlations, otu_vs_otu_sample_correlations,meta_vs_otu_random_correlations)

as.data.frame(pearson_sample_correlations)->pearson_sample_correlations
as.numeric(as.character(pearson_sample_correlations[,3]))->pearson_sample_correlations[,3]

c("metasample", "otusample","correlation", "comparison")->colnames(pearson_sample_correlations)
factor(pearson_sample_correlations$comparison, levels=c("metaVotu", "metaVmetaR", "metaVotuR", "otuVotuR"))->pearson_sample_correlations$comparison

#make a ggplot version of the correlations

library(ggplot2)

#Graph 
p <- ggplot(pearson_sample_correlations, aes(comparison, correlation))
p + geom_boxplot()
p + geom_boxplot() + geom_jitter(width = 0.2)



.76 otureads vs metareads
.78 random otureads
.88 random metareads


#make correlations by taxa. Dooes the rank of taxa abundances across samples correlate with taxa abundances in sample by another method?
vector(length=nrow(otureads_common))->correlations
for(c in 1:nrow(otureads_common)){
	sample(1:nrow(otureads_common))[c]->d
	print(c(c, d))
	cor(metareads_common[c,], metareads_common[d,], method="pearson")->correlation
	correlation->correlations[c]
}
#plot(sort(correlations), ylim=c(0, 1), main="taxa")
points(sort(correlations), col="blue")



	(otureads_unique[,c],rep(0, length=length(otureads_unique[,c])), col="blue")
}









#augment each table with those bacteria that are only found in one but not the other so that each table contains the union

setdiff(meta_names, otu_names)->meta_unique_taxa
meta_reads



otu_names
