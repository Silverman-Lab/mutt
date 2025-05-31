paper:https://www.nature.com/articles/ismej2016117
Flow cytometry data (.fcs format) are available on the FlowRepository archive under repository ID FR-FCM-ZZNA and the Dryad Digital Repository (http://dx.doi.org/10.5061/dryad.m1c04). Sequences are available on the NCBI Sequence Read Archive (SRA) under accession number SRP066190.
Case study data used in MEASURING THE BIODIVERSITY OF MICROBIAL COMMUNITIES BY FLOW CYTOMETRY by Ruben Props, Pieter Monsieurs, Mohamed Mysara, Lieven Clement and Nico Boon

Data provided by Nico Boon, email: <Nico.Boon@UGent.be>.
See Materials and Methods of the paper for detailed information about the data.
The file metadata folder contains additional information that was collected concurrently with the flow cytometry analysis:
1. Denoising strategy that was used on the flow cytometry data. Explanation for this can be found in SI figure 1 of the paper.

2. Extra metadata
- Cycle: the number of the survey (1 or 2)
- Reactor phase: phase of the reactor operation
- Conductivity, pH and temperature of the ecosystem
- Timepoint of sampling relative to the start-up of the reactors for each respective cycle
- D2.tax/D1.tax/D0.tax: Hill species diversities based on the 16S rRNA gene amplicon sequencing. Further description of these parameters is available in the materials and methods section of the paper.
- D2.fcm/D1.fcm/D0.fcm: Hill diversities calculated from the flow cytometry data. Calculation is described in materials and methods of the paper.
- D2.fcm.sd/D1.fcm.sd/D0.fcm.sd: standard deviations on the above diversity estimates based on triplicate measurements.
- Cell.density: total cell density of the community as calculated by flow cytometry
- Cell.density.sd: standard deviation on the total cell density measurements based on triplicate measurements.

<<<<<<< Updated upstream
^ The above is great, but there is no way to connect this information to the count data or the sequencing files because the sample ID is not synonymous. So that extra metadata is unfortunately not useful at the moment. Hopefully we can retain this information.

Maxwell: I recollected all information, used the SRA file and supplemental table 2 to pair the sequences to the flow data through the sample title column. Nico Boon gave me Dr. Props email on 4/25/25 although I never reached out because there was no need. Was going to ask about the A in the Sample_name column, but I just dropped it. I wrote the parse.R and cleaned the data.  
=======
>>>>>>> Stashed changes

