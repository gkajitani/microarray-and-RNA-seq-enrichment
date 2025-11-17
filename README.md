# GEO Microarray and RNA-SEQ-Alingment, followed by differential expression and enrichment analysis
Pipeline used for Kajitani et al., 2025, "Altered pathways in Cockayne syndrome: Involvement of MAPK, PI3K-Akt, extracellular matrix, inflammation, and neuronal signaling"

The article analyzed transcriptome datasets generated using either microarrays or RNA-Seq. Both types of data were obtained from Gene Expression Omnibus (GEO). Microarray data was
processed using R, while RNA-Seq data was processed using a pipeline based on Hisat2 alignment followed by read mapping and counting through featureCounts.

The RNA-Seq pipeline is used to obtain data from NCBI SRA (Sequence Read Archive), convert it to a .fastq format, remove adapter sequence and sequences with low quality scores.
Reads are then aligned to a reference genome, generating a .sam file. This file is then converted into a .bam file, used to get the number of raw read counts aligned in each gene.
Read counts table is then read in R and read count numbers are normalized via DESeq2 package, also used to identify differentially expressed genes, along with their respective log2 fold change.
R packages clusterprofiler and enrichplot are then used for enrichment analysis (Over-representation analysis (ORA) and Gene set enrichment analysis (GSEA) and graph plotting.


##############
# Microarray #
##############

# Main packages, available at bioconductor, scripts based on GEO2R:

- Biobase;
- GEOquery;
- limma;

##############
# RNA-seq #
##############

# Main softwares:


- SRA Tools for obtaining .sra data and converting to fastq:
  - Avaiable at https://github.com/ncbi/sra-tools;
  - Check github page to learn how to install and use;
    
- Trimmomatic for adapter and low-quality sequence trimming:
  - Available at https://github.com/usadellab/Trimmomatic;
  - Manual available at http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf;
  - Used for both single-end and paired-end sequences;
 
- Hisat2 for read alignment:
  - Available at http://daehwankimlab.github.io/hisat2/
  - Check github page to learn how to install and use, as well as index download links;
 
- Samtools to convert from .sam to .bam:
  - Software and manual available at https://www.htslib.org/download/ or https://github.com/samtools/samtools;

- Subread to get read counts per gene:
  - Software and manual available at https://github.com/ShiLab-Bioinformatics/subread;
  - Will output a .txt read count per sample. Combine those into a .csv table containing all samples and gene symbols.
 
###########
# Using R #
###########

# Main packages, available at bioconductor:

- DESeq2;
- clusterprofiler;
- enrichplot;
