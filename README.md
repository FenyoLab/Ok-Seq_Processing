# Ok-Seq_Processing

This repository contains scripts for processing OkSeq data.  
(1) process_okseq_SE_PE.py and run_process_okseq.sh can be used to process FASTQ files to obtain text files contain OkSeq read density on W and C strand for every 1kb.

(2) pre_process_gene_table_file.py is the script that filters a list of genes (obtained from UCSC genome browser) for only protein coding genes and adds fpkm values from a separate file of RNA-Seq data (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89413)

(3) add_okazaki_data_to_site_list.py is the script that will then take the gene list and add columns for the OkSeq data -100 to 100 kb around TSS and TTS for each gene.

(4) make_strand_bias_plots.py will then plot this data as in Figs1 and 2 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6320713/

Input data and results are also included.  2 large files must be unzipped before running.
