# Ok-Seq_Processing

This repository contains scripts for processing OkSeq data.  

 (1) process_okseq_SE_PE.py and run_process_okseq.sh can be used to process FASTQ files to obtain text files contain OkSeq read density on W and C strand for every 1kb.

 (2) pre_process_gene_table_file.py is the script that filters a list of genes (obtained from UCSC genome browser) for only protein coding genes and adds fpkm values from a separate file of RNA-Seq data (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89413)

 (3) add_okazaki_data_to_site_list.py is the script that will then take the gene list and add columns for the OkSeq data -100 to 100 kb around TSS and TTS for each gene.

 (4) make_okseq_plots.py will then plot this data as strand bias plots, as in Figs1 and 2 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6320713/, as well as gene normalized strand bias plots and calculation of p-values for comparison of strand bias plots.


Input data and results are also included.  (2 large files must be unzipped before running.)  
The input data is in the folder okseq_data.  Subfolders and their contents are as follows:

* gene_name_aliases: protein_coding_gene.txt - this file was obtained from HGNC (https://www.genenames.org) and lists gene symbols and their aliases.  This file is used when mapping gene names from the RNA-Seq data (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89413) for RPE-1 cells to gene names from the table of all human genes (GRCh37/hg19) obtained from UCSC table browser (https://genome.ucsc.edu/cgi-bin/hgTables).    

* raw_files/rpe: rpe_edu_2.txt - this file contains the Watson and Crick stand reads for every 1kb bin.  It is the result of the initial processing as described in (1) above.

* rna_seq_RPE: GSE89413_2016-10-30-NBL-cell-line-STAR-fpkm.txt - the RNA-seq FPKM results per gene (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89413)

* strand_bias_plots: the results from running the script in step (4) above.  

