# Function: SampleTrackeR
# Version: 1.0

SampleTrackeR <- function(sample_plate_layout,
                          read_count_table, 
                          stm_compositions,
                          read.threshold = 1,
                          fraction.threshold = 1) {

####################################################################################################
  
  # verifying required arguments -------------------------------------------------------------------
  
  if (missing(sample_plate_layout))
    stop(paste0("Argument \"sample_plate_layout\" missing ... stopping SampleTrackeR script."),
         call. = FALSE) 
  
  if (missing(read_count_table))
    stop(paste0("Argument \"read_count_table\" missing ... stopping SampleTrackeR script."),
         call. = FALSE)
 
  if (missing(stm_compositions))
    stop(paste0("Argument \"stm_compositions\" missing ... stopping SampleTrackeR script."),
         call. = FALSE)
  
  # importing and formatting tab-delimited input files ---------------------------------------------
  
  input.search.path = getwd()
  
  ## sample_plate_layout
  
  if (file.exists(paste0(input.search.path, "/", sample_plate_layout))) {
    sample_plate_layout.df <- read.delim(paste0(input.search.path, "/", sample_plate_layout),
                                         header = TRUE, 
                                         stringsAsFactors = FALSE, 
                                         sep = "\t", 
                                         check.names = FALSE) %>% as.data.frame()
    } else {
      stop(paste0("File \"", paste0(input.search.path, "/", sample_plate_layout), "\" not found ... stopping SampleTrackeR script."),
           call. = FALSE) 
    }
  
  if (isTRUE(any(names(sample_plate_layout.df) == 'libID')) &
      isTRUE(any(names(sample_plate_layout.df) == 'stmID')) &
      isTRUE(any(names(sample_plate_layout.df) == 'row')) &
      isTRUE(any(names(sample_plate_layout.df) == 'column'))) {
    sample_plate_layout.df %<>% dplyr::select(libID, stmID, row, column)
    } else {
    stop(paste0("Mandatory column(s) missing in \"sample_plate_layout\" file ... stopping SampleTrackeR script."),
         call. = FALSE)
    }
  
  if(!is.numeric(sample_plate_layout.df$row) | !is.numeric(sample_plate_layout.df$column)) {
    stop(paste0("Row and column identifiers should be numeric in \"sample_plate_layout\" file ... stopping SampleTrackeR script."), 
         call. = FALSE)
  }
  
  ## read_count_table
  
  if (file.exists(paste0(input.search.path, "/", read_count_table))) {
    read_count_table.df <- read.delim(paste0(input.search.path, "/", read_count_table),
                                      header = TRUE, 
                                      stringsAsFactors = FALSE, 
                                      sep = "\t", 
                                      check.names = FALSE) %>% as.data.frame()
  } else {
    stop(paste0("File \"", paste0(input.search.path, "/", read_count_table), "\" not found ... stopping SampleTrackeR script."), call. = FALSE)
  }
  
  if (isTRUE(any(names(read_count_table.df) == 'otuID'))) {
    read_count_table.df %<>% tidyr::gather(variable, value, -otuID) %>% `colnames<-`(c("controlID", "libID", "reads"))
  } else {
    stop(paste0("Mandatory \"otuID\" column missing in \"read_count_table\" file ... stopping SampleTrackeR script."), call. = FALSE)
  }
  
  ## stm_compositions
  
  if (file.exists(paste0(input.search.path, "/", stm_compositions))) {
    stm_compositions.df <- read.delim(paste0(input.search.path, "/", stm_compositions),
                                      header = TRUE, 
                                      stringsAsFactors = FALSE, 
                                      sep = "\t", 
                                      check.names = FALSE) %>% as.data.frame()
    } else {
      stop(paste0("File \"", paste0(input.search.path, "/", stm_compositions), "\" not found ... stopping SampleTrackeR script."), call. = FALSE)
    }
  
  if (isTRUE(any(names(stm_compositions.df) == 'stmID')) & 
      isTRUE(any(names(stm_compositions.df) == 'controlID')) &
      isTRUE(any(names(stm_compositions.df) == 'value'))) {
    stm_compositions.df %<>% dplyr::select(stmID, controlID, value)
  } else {
    stop(paste0("Mandatory column(s) missing in \"stm_compositions\" file ... stopping SampleTrackeR script."), call. = FALSE)
  }
  
####################################################################################################
  
  # identifying STMs for all samples ---------------------------------------------------------------
  
  ## inititating empty stm.map data frame
  
  stm.map <- data.frame("libID" = as.character(), 
                        "stmID" = as.character(), 
                        "reads" = as.numeric(),
                        "fraction" = as.numeric())
  
  ## filling out stm.map by looping through all samples and STMs
  
  for ( i in as.vector(unique(read_count_table.df$libID))) { # looping samples
    
    for (j in as.vector(unique(sample_plate_layout.df$stmID))) { # looping STMs
      
      # extracting read counts of sample i
      observed <- read_count_table.df[read_count_table.df$libID == i, ][,c(1, 3)]

      # extracting composition of STM j
      expected <- stm_compositions.df[stm_compositions.df$stmID == j, ][,c(2, 3)] 
      
      # combining observed and expected data frames
      merged <- merge(expected, observed, by = "controlID", all.x = TRUE)
      merged$reads[is.na(merged$reads)] <- 0
      
      # obtaining summed read count of controls present in STM j in sample i
      reads <- sum(merged[merged$value == 1 & merged$reads >= read.threshold, ]$reads) 
      
      # obtaining fraction controls of STMs j present in sample i
      fraction <- nrow(merged[merged$value == 1 & merged$reads >= read.threshold, ]) / nrow(merged[merged$value == 1, ])
      
      # updating stm.map
      stm.map <- rbind(stm.map, data.frame("libID" = i, "stmID" = j, "reads" = reads, "fraction" = fraction))
    }
  }
  
  # filtering identified STMs based on fraction.threshold
  stm.map.filtered <- stm.map %>% dplyr::filter(., fraction >= fraction.threshold)
  
  # extracting majority STMs
  majority.stms <- stm.map.filtered %>%
    dplyr::group_by(libID) %>%
    dplyr::mutate(max = max(reads)) %>%
    dplyr::filter(reads == max) %>%
    dplyr::select(-max) %>%
    `colnames<-`(c("libID", "majority.stm_ID", "majority.stm_reads", "majority.stm_prevalence")) %>%
    as.data.frame()
  
  # extracting minority STMs
  minority.stms <- stm.map.filtered %>%
    dplyr::group_by(libID) %>%
    dplyr::mutate(max = max(reads)) %>%
    dplyr::filter(reads < max) %>%
    dplyr::select(-max) %>%
    `colnames<-`(c("libID", "minority.stm_ID", "minority.stm_reads", "minority.stm_prevalence")) %>%
    as.data.frame()
    
####################################################################################################
  
  # Part 1 ::: quality control for sample swaps ----------------------------------------------------
  
  # Part 1.1 ::: tabular output --------------------------------------------------------------------
  
  df1 <- base::merge(read_count_table.df %>%
                       dplyr::mutate_if(is.factor, as.character) %>%
                       dplyr::filter(controlID %in% as.vector(unique(stm_compositions.df$controlID))) %>%
                       dplyr::group_by(libID) %>%
                       dplyr::summarize(total.control.reads = sum(reads)),
                     read_count_table.df %>% 
                       dplyr::group_by(libID) %>%
                       dplyr::summarize(total.reads = sum(reads)),
                     by = "libID")
  
  df1 <- base::merge(majority.stms, df1, by = "libID", all.y = TRUE) %>% 
    base::merge(sample_plate_layout.df, ., by = "libID") %>% 
    dplyr::select(libID, 
                  stmID, 
                  total.reads, 
                  total.control.reads, 
                  majority.stm_ID, 
                  majority.stm_reads, 
                  majority.stm_prevalence) %>%
    `colnames<-`(c("libID", 
                   "expected.stm_ID", 
                   "total.reads", 
                   "total.control.reads", 
                   "majority.stm_ID", 
                   "majority.stm_reads", 
                   "majority.stm_prevalence")) %>% 
    dplyr::mutate_if(is.factor, as.character)
  
  notadded <- stm_compositions.df %>% 
    dplyr::group_by(stmID) %>% 
    dplyr::summarize(max = max(value)) %>%
    dplyr::filter(max == 0 ) %>%
    dplyr::select(stmID) %>% unlist(use.names = FALSE)
  
  df1$qc.passed <- "-"
  df1$qc.passed[df1$expected.stm_ID == df1$majority.stm_ID] <- "yes"
  df1$qc.passed[df1$expected.stm_ID != df1$majority.stm_ID] <- "NO"
  if(length(notadded) > 0) {
    df1$expected.stm_ID[df1$expected.stm_ID == notadded] <- "-"
    }
  df1$majority.stm_ID[is.na(df1$majority.stm_ID)] <- "-"
  df1$majority.stm_reads[is.na(df1$majority.stm_reads)] <- "-"
  df1$majority.stm_prevalence[is.na(df1$majority.stm_prevalence)] <- "-"
  
  # Part 1.2 ::: grahical output -------------------------------------------------------------------
  
  ## setting up dataframe
  df2 <- base::merge(read_count_table.df %>% 
                       dplyr::mutate_if(is.factor, as.character) %>%
                       dplyr::filter(controlID %in% as.vector(unique(stm_compositions.df$controlID))) %>%
                       dplyr::group_by(libID) %>%
                       dplyr::summarise(total.read.count = sum(reads)),
                     read_count_table.df %>% 
                       dplyr::mutate_if(is.factor, as.character) %>%
                       dplyr::filter(controlID %in% as.vector(unique(stm_compositions.df$controlID))),
                     by = "libID") %>% 
    dplyr::mutate(proportion = reads / total.read.count * 100) %>%
    dplyr::select(libID, controlID, proportion) %>%
    dplyr::right_join(., sample_plate_layout.df, by = "libID") %>%
    dplyr::left_join(., majority.stms %>% dplyr::mutate_if(is.factor, as.character), by = "libID") %>%
    dplyr::mutate(libID = factor(libID,levels = unique(libID))) %>%
    dplyr::mutate(majority.stm_ID = factor(majority.stm_ID,levels = unique(stmID))) %>%
    as.data.frame()
  
  if(length(notadded) > 0) {
    df2 %<>% dplyr::filter(., !grepl(notadded, stmID))
  }
  
  ## generating plot
  theme_set(theme_gray(base_size = 14))
  plot1 <- ggplot()+
    geom_bar(data = df2,
             mapping=aes(x = factor(1), y = proportion, fill = controlID),
             stat = "identity", width = 1, 
             colour = NA, size = 0.25) +
    coord_polar(theta="y") +
    facet_grid(majority.stm_ID ~ libID) +
    theme(axis.text=element_blank())+
    theme(panel.grid.major=element_blank())+
    theme(panel.grid.minor=element_blank())+
    theme(axis.ticks=element_blank())+
    xlab("") + ylab("")+
    theme(legend.title = element_blank())+
    theme(panel.spacing = unit(0.1, "lines"))+
    theme(strip.text.x = element_text(angle = 90))+
    theme(strip.text.y = element_text(angle = 0, hjust = 1))
  
  ## generating message
  if (nrow(df1[df1$qc.passed == "NO", ]) > 0) {
    message(paste0("\nBased on the majority STMs, it appears that ", nrow(df1[df1$qc.passed == "NO", ]), " samples may have been mislabeled/swapped.\n ... The expected STMs will be updated based on the identified majority_STM for subsequent assessment of between-sample carry-over ... \n"))
  } else {
    message("Based on the identified majority STMs, no samples appear to have been swapped/mislabeled.")
  }

####################################################################################################
  
  # Part 2 ::: quality control for sample cross-contamination --------------------------------------
  
  # Part 2.1 ::: tabular output --------------------------------------------------------------------
  
  contamination.map <- merge(majority.stms[,c(1,2)], 
                             minority.stms[,c(1,2)], 
                             all = TRUE, 
                             by = "libID") %>% na.omit()
  
  if (nrow(contamination.map) >= 1) { # if no minority STMs are found this portion will be skipped #
    
    carryover.map <- data.frame("libID" = as.character(),
                                "majority.stm_ID" = as.character(),
                                "minority.stm_ID" = as.character(),
                                "majority.stm_reads" = as.numeric(),
                                "minority.stm_reads" = as.numeric(),
                                "distinguishing.controls" = as.numeric(),
                                "carryover_estimated.percent" = as.numeric(),
                                "carryover_number.reads" = as.numeric())
    
    for (i in 1:nrow(contamination.map)) {
      
      # obtaining composition majority STM
      major.mix <- stm_compositions.df[grep(contamination.map[i,2], stm_compositions.df$stmID), ][,-c(1)] %>%
        `colnames<-`(c("controlID", "major.expected"))
      
      # obtaining composition minority STM
      minor.mix <- stm_compositions.df[grep(contamination.map[i,3], stm_compositions.df$stmID), ][,-c(1)] %>%
        `colnames<-`(c("controlID", "minor.expected"))
      
      # obtaining reads for sample
      reads <- read_count_table.df[grep(contamination.map[i,1], read_count_table.df$libID), ] %>% 
        dplyr::select(controlID, reads)
      
      # generating merged data frame
      dat <- base::merge(reads, major.mix, by = "controlID") %>% merge(., minor.mix, by = "controlID")
     
      # summing reads for majority STMs
      majority.stm_reads <- sum(dat[dat$major.expected == 1, ]$reads)
      
      # summing reads for minority STMs
      minority.stm_reads <- sum(dat[dat$minor.expected == 1, ]$reads)
      
      # extracting spike-in controls only in majority STM
      dat1 <- dat[dat$major.expected == 1 & dat$minor.expected == 0, ]
     
      # extracting spike-in controls only in minority STM
      dat2 <- dat[dat$major.expected == 0 & dat$minor.expected == 1, ]
      
      # combining both sets #
      dat <- rbind(dat1,dat2)
      
      # counting number of distinguishing controls, not present in majority STM
      distinguishing.controls <- nrow(dat2) 
      
      # normalizing reads counts
      dat$reads.normalized <- dat$reads / sum(dat$reads) * 100
      
      # subsetting to reads in minority STM
      dat <- dat[dat$minor.expected == 1, ]
      
      # summing normalized reads for minority STM as a measure of carry-over
      carryover_estimated.percent <- sum(dat$reads.normalized)
      
      # obtaining reads in carry-over
      carryover_number.reads <- sum(dat2$reads) # carryover_number.reads
      
      # updating carryover.map
      carryover.map <- rbind(carryover.map,
                           data.frame("libID" = contamination.map[i,1],
                                      "majority.stm_ID" = stm_compositions.df[grep(contamination.map[i,2], stm_compositions.df$stmID), ][1,1],
                                      "minority.stm_ID" = stm_compositions.df[grep(contamination.map[i,3], stm_compositions.df$stmID), ][1,1],
                                      "majority.stm_reads" = majority.stm_reads,
                                      "minority.stm_reads" = minority.stm_reads,
                                      "distinguishing.controls" = distinguishing.controls,
                                      "carryover_estimated.percent" = carryover_estimated.percent,
                                      "carryover_number.reads" = carryover_number.reads))
    }
    
    df3 <- merge(carryover.map, sample_plate_layout.df %>% dplyr::select(libID), all.y = TRUE)
    
    message(paste0("Based on the identified minority STMs, ", nrow(carryover.map), " potential cases of between-sample carry-over have been identified.\n"))
    
    # Part 2.2 ::: graphical output ----------------------------------------------------------------
    
    temp <- sample_plate_layout.df
    temp$majority.stm_ID <- temp$stmID
    temp$minority.stm_ID <- temp$stmID
    
    crosstalkmap <- merge(df3[,c(2, 3, 6, 7)], temp[,c(5, 3, 4)], by = "majority.stm_ID")
    colnames(crosstalkmap)[5] <- "row.end"
    colnames(crosstalkmap)[6] <- "column.end"
    
    crosstalkmap <- merge(crosstalkmap, temp[,c(6, 3, 4)], by = "minority.stm_ID")
    colnames(crosstalkmap)[7] <- "row.start"
    colnames(crosstalkmap)[8] <- "column.start"
                            
    plot2 <- ggplot() + 
      geom_curve(data = crosstalkmap, 
                 mapping = aes(x = column.start, xend = column.end,
                               y = row.start, yend = row.end,
                               size = carryover_estimated.percent,
                               color = as.factor(distinguishing.controls)),
                 alpha = 0.75, curvature = 0.09, ncp = 100, arrow = grid::arrow(length = grid::unit(0.03, "npc"))) +
      scale_size(range = c(0.1, 2), limits = c(min(crosstalkmap$carryover_estimated.percent), ceiling(max(crosstalkmap$carryover_estimated.percent)))) +
      xlab("") + ylab("") +
      scale_y_reverse(limits = c(max(crosstalkmap$row.end, crosstalkmap$row.start) + 0.9, min(crosstalkmap$row.end, crosstalkmap$row.start) - 0.9), expand = c(0, 0)) +
      scale_x_continuous(limits = c(min(crosstalkmap$column.end, crosstalkmap$column.start) - 0.9, max(crosstalkmap$column.end, crosstalkmap$column.start) + 0.9), expand = c(0, 0)) +
      guides(colour = guide_legend(override.aes = list(size=2)))
    
    } else {
    message("No minority STMs were identified in any of the samples.\n ...Evaluation of between-sample carry-over was not performed ...\n")
    }

####################################################################################################

  if (nrow(contamination.map) >= 1) {
    return(list("tab1" = df1,
                "plot1" = plot1,
                "tab2" = df3,
                "plot2" = plot2))
    } else {
    return(list("tab1" = df1,
                "plot1" = plot1))
    }
}