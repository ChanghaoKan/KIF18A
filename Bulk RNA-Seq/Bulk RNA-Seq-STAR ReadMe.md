RNA-seq Analysis Pipeline Script
This bash script processes RNA-seq raw sequencing data into gene count data using the STAR alignment method. 
Users are required to provide their own sequencing data and reference files. Please cite if using this script:
Kan, C. et al. KIF18A in HCC Metastasis: Dual Role in Anoikis Resistance and Chromosome Instability. bioRxiv, 2024.2012.2020.629627, doi. 10.1101/ 2024.12.20.629627 (2024) https://www.biorxiv.org/content/10.1101/2024.12.20.629627v1.full.
Below is a list of all the packages and tools utilized by the script, along with brief descriptions of their roles.
Tools and Packages Used
The script depends on the following tools and packages, most of which are installed within a Conda environment:
Conda Environment Packages:
fastqc: Performs quality control checks on raw sequencing reads to assess data quality.

multiqc: Aggregates multiple FastQC reports into a single, easy-to-read summary.

trim-galore: Trims low-quality bases and adapters from sequencing reads to improve alignment accuracy.

star: Aligns sequencing reads to a reference genome and generates gene counts for downstream analysis.

samtools: Provides utilities for manipulating and analyzing BAM files produced during alignment.

subread: Offers tools for read counting and other sequence processing tasks.

r-base (version 4.4.0): A programming language and environment for statistical computing and visualization.

bioconductor-deseq2: An R package for differential gene expression analysis using count data.

Download Tool:
wget: A command-line utility used to download sequencing data files from specified URLs.

Text Editing Tool:
nano: A simple text editor used to create and edit files, such as urls.txt for storing data download links.

Linux Command-Line Tools:
find: Searches for files and directories within the filesystem.

sed: A stream editor for filtering and transforming text in files or pipelines.

awk: A text-processing tool for extracting and manipulating data from structured files.

These tools and packages collectively enable the script to handle the entire RNA-seq workflow, from data retrieval and quality control to alignment, counting, and preparation for differential expression analysis.

