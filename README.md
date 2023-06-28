# Microbiome
Analysis of fungal ITS1 amplicon sequence data in R


### What fungal species are in the air and what are their weather determinants?
The air contains a multitude of microorganisms - bacteria, fungi, viruses - that we breathe in all day. Some of these organisms cause economically important diseases on agricultural crops. This data set includes paired end 250b Illumina amplicon sequencing data of the ITS1 region - the universal barcode for fungi. Air samples were collected over a range of agricultural fields over two different years. Below is a rough explanation on each of the .R files included.

For more background information and to access the publication associated with this data set, see [here](https://apsjournals.apsnet.org/doi/abs/10.1094/PHYTOFR-10-22-0108-R).

### R files:
1. create_phyloseq.R - takes the amplicon sequence data, metadata, and taxonomy files to create a phyloseq object
2. rarefaction_of_phyloseq.R - preliminary data exploration and rarefaction of data (N.B. some issues with rarefying are addressed in the code)
3. alpha_diversity_analysis.R - plotting alpha diversity metrics across data sets; statistical analysis (linear models and ANOVAs) to determine whether there are differences in Location and how alpha diversity changes over time
4. beta_diversity_analysis.R - ordinations and PermANOVAs to investigate community level differences
5. funguild_analysis.R - using the FunGUILD database to investigate fungal plant pathogens of interest
