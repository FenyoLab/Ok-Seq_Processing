# Input is a set of genomic positions in table format (csv)
# This program will look up the raw W and C reads +/- 100 for each genomic position
# This W/C data can come from the following sources:
# It will add this data to the table and output a new table.
# for example, if 'TSS' is the name of the column containing position data,
# new columns will be named like: distr_c_100_from_TSS_1, distr_c_100_from_TSS_2, etc,
# and there will be 201 columns:
# dist_c_100_from_TSS_1 is -100 position,
# dist_c_100_from_TSS_101 is TSS position,
# dist_c_201_from_TSS is +100 position

import pandas as pd
import numpy as np
import sys
import os

if(len(sys.argv) >= 3):
    data_dir = sys.argv[1]
    input_txt_file = sys.argv[2]
    sites_file = sys.argv[3]
    fpkm_prefix = sys.argv[4]
    results_dir = sys.argv[5]
else:
    print("Missing command line input.  Attempting to run with default settings.")
    data_dir = "/Users/sarahkeegan/okseq_data"
    input_txt_file = "/raw_files/rpe/rpe_edu_2.txt"
    sites_file = "hgTables_filtered_with_fpkm_first.txt"
    fpkm_prefix = "RPE1"
    results_dir = "strand_bias_plots"

length_around_site = 200

allowed_chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14',
                'chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']

# READ IN TXT FILES OF OKAZAKI W AND C READ COUNTS
print("Reading in W/C Okazaki data...")

all_data = pd.read_table(data_dir + '/' + input_txt_file)

#correct small inconsistency in RPE raw files
all_data.loc[all_data.chr > 1, 'pos'] = all_data.loc[all_data.chr > 1, 'pos']-1
index_data = all_data.copy()
index_data['row_index'] = range(len(index_data))
index_data.drop(['w','c'], axis=1, inplace=True)
index_data.set_index(['chr','pos'], inplace=True)
all_data.index = range(len(all_data))

# READ IN SITE TABLE FILE
print("Processing sites file...")

# col file names for w/c data in gene table
col_str_list_w = {'TSS': [], 'TTS': []}
col_str_list_c = {'TSS': [], 'TTS': []}
for j in range(1, length_around_site + 1 + 1):
    for site in ['TSS', 'TTS']:
        col_str_list_w[site].append('dist_w_100_from_' + site + '_' + str(j))
        col_str_list_c[site].append('dist_c_100_from_' + site + '_' + str(j))

gene_sites = pd.read_table(data_dir + '/' + sites_file)

#rename fpkm col for the data source we are processing
gene_sites['fpkm']=gene_sites[fpkm_prefix + '_fpkm']

#for each gene, add W/C data as columns in the file, for all raw data files
gene_sites.index = range(len(gene_sites))

w_rows_arr={}
c_rows_arr={}
for site in ['TSS','TTS']:
    w_rows_arr[site] = np.zeros(shape=(len(gene_sites),length_around_site+1))
    c_rows_arr[site] = np.zeros(shape=(len(gene_sites), length_around_site + 1))

prev_chrom = 0
for i, row in enumerate(gene_sites.iterrows()):
    if(i % 1000 == 0):
        print(i)
        #if(i>0):break

    row = row[1]
    strand = row['hg19.knownGene.strand']
    chr_ = row['hg19.knownGene.chrom']
    chr_ = int(chr_.lstrip('chr'))
    if(chr_ != prev_chrom):
        # work with only this chrom drop all else
        cur_wc_data = all_data[all_data['chr'] == chr_]
        first_index = cur_wc_data.first_valid_index()
        last_index = cur_wc_data.last_valid_index()
        first_pos = cur_wc_data.loc[first_index, 'pos']
        last_pos = cur_wc_data.loc[last_index, 'pos']
        prev_chrom = chr_

    site_pos={}
    if(strand == '+'):
        site_pos['TSS'] = round(row['hg19.knownGene.txStart'],-3)
        site_pos['TTS'] = round(row['hg19.knownGene.txEnd'], -3)
    else:
        site_pos['TSS'] = round(row['hg19.knownGene.txEnd'],-3)
        site_pos['TTS'] = round(row['hg19.knownGene.txStart'],-3)

    for site in ['TSS','TTS']:
        #around TSS/TTS
        if (site_pos[site] >= first_pos and site_pos[site] <= last_pos):
            base_pos = index_data.loc[(chr_,site_pos[site])]['row_index']
            st = int(base_pos-length_around_site/2)
            end = int(base_pos+length_around_site/2)
            if(st >= first_index and end < last_index):
                w_rows_arr[site][i]=list(cur_wc_data.loc[st:end]['w'])
                c_rows_arr[site][i]=list(cur_wc_data.loc[st:end]['c'])

    #sum TSS to TTS
    sum_st = round(row['hg19.knownGene.txStart'],-3)
    sum_end = round(row['hg19.knownGene.txEnd'], -3)
    if (sum_st >= first_pos and sum_end <= last_pos):
        st = index_data.loc[(chr_, sum_st)]['row_index']
        end = index_data.loc[(chr_, sum_end)]['row_index']
        if (st >= first_index and end < last_index):
            gene_sites.loc[i,'w_sum_TSS_to_TTS'] = np.sum(cur_wc_data.loc[st:end]['w'])
            gene_sites.loc[i,'c_sum_TSS_to_TTS'] = np.sum(cur_wc_data.loc[st:end]['c'])

print("Finished. Concating and saving data files...")
w_rows_df={}
c_rows_df={}
for site in ['TSS','TTS']:
    w_rows_df[site] = pd.DataFrame(index=range(len(gene_sites)), columns=col_str_list_w[site], data=w_rows_arr[site])
    c_rows_df[site] = pd.DataFrame(index=range(len(gene_sites)), columns=col_str_list_c[site], data=c_rows_arr[site])
#add w_rows and c_rows as additional columns
gene_sites = pd.concat([gene_sites, w_rows_df['TSS'], c_rows_df['TSS'], w_rows_df['TTS'], c_rows_df['TTS']], axis=1)

file_name=os.path.split(input_txt_file)[1]
file_prefix = 'sites_table_with_distributions-' + file_name[:-4]
gene_sites.to_csv(data_dir+'/'+results_dir+'/'+file_prefix + '.csv')
