Data from paper: https://www.nature.com/articles/s41587-023-01754-3


### **NOTE:**
1. Download and unzip the file.
2. Copy the folder path and use it in the function:
   `parse_2023_maghini_naturebio_metagenomic(YourFolderPaths = PATH)`
3. The author uploaded count data at different taxonomic levels, which I merged into a list. As a result, this function does not output taxa directly. If you wish to use the genus-level counts, access them with `count$genus`. To extract taxa names, use `rownames(count$genus)`.
4. Metadata/modeling information from author:
    + D## = the donor number.
    + N = no preservative
    + O = Omni preservative
    + Z = Zymo preservative
    + F = Frozen (-80C) right away
    + H = hot temperature (40C)
    + R = room temperature (23C)
    + R1 = replicate 1
    + R2 = replicate 2
    + R3 = replicate 3
    + Note: In the GEE model, no other covariates was included in the analysis. (Quto "No covariates, it was really the simplest thing. As we say in the stats methods section, we used an ""unadjusted regression model with fixed effects for the seven conditions and an exchangeable correlation structure between the participant-level clusters to account for repeat measurements."  So, for instance, if we were modeling the log of the total abundance, the model would have just been log_total_abundance ~ C(Sample_Type).")

