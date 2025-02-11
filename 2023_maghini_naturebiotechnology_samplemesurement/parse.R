library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)

parse_2023_maghini_naturebio_metagenomic <- function(YourFolderPaths = NULL){

        # save all data in one folder and change the path to the directory to that path
        base_path <-YourFolderPaths

        # merge count data with taxonomic level into a list, if you want to call the data for genus, use tax_list$genus for example.
        tax_level <- c("class","family","genus","kingdom","order","phylum")
        file_path <- paste0(path = base_path, pattern ="/bracken_",tax_level,"_reads.txt")
        tax_list <- lapply(file_path, function(file) {
        read.table(file,header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
        })
        names(tax_list) <- tax_level

        # merge qPCR data into one table
        plate_list <- list.files(path = base_path, pattern ="qPCR_plate", full.names = TRUE)
        qPCR <- plate_list %>% 
        map_dfr(read_csv)%>% # create NA if there is no data in the column in one of the files
        mutate(ID = gsub("_","-",SampleName))%>% # match the ID in the count table
        dplyr::select(Plate,ID,Donor,Condition,PCR_Replicate,Replicate,logCopyNumber,CopyNumber) %>%
        filter(!Condition %in% c("B1","B2","B3","B4")) #have no idea what they are but they are not donor samples so I remove them

        # there are few samples (water and buffer) that are not in the taxonomic data, remove them
        filtersample <- setdiff(qPCR$ID,colnames(tax_list$family)) 
        qPCR <- qPCR %>% filter(!ID %in% filtersample)

        # remove non-donor samples from the taxonomic data
        filtersample <- setdiff(colnames(tax_list$family),qPCR$ID)
        count <- lapply(tax_list, function(x) x[,-which(colnames(x) %in% filtersample)])

        # Note: There are two PCR_replicate, so that in the qPCR sample, the total obs is 410, but the unique obs is 205 in the taxonomic data.
        # I am not sure which replicate to remove, so I will keep the first one but you can change to another based on the PCR_replicate column.
        scale <- qPCR %>% filter(PCR_Replicate == "Rep1") %>% dplyr::select(ID,logCopyNumber,CopyNumber)

        # exrtact metadata
        metadata <- qPCR %>% dplyr::select(ID,Donor,Condition,Replicate) %>% distinct()

        all.equal(nrow(scale),nrow(metadata))
        all.equal(nrow(scale),ncol(count$family)) 

        warning("NOTE1:The data is ready to be used, since the author upload count data with different taxonomic level, 
                I put them in a list, so you can call the count data based on the taxonomic level such as count$genus; 
                if you need the taxa name, just call the rownames(count$genus) to get the taxa name.")

        warning("NOTE2: Below are the information for metadata$condition:
                N = no preservative
                O = Omni preservative
                Z = Zymo preservative
                F = Frozen (-80C) right away
                H = hot temperature (40C)
                R = room temperature (23C))")

    return(list(scale=scale,count=count,metadata=metadata))

}

