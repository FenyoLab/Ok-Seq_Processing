# Performs all steps necessary to process OkSeq sequencing files (fastq PE)
# Result is density of reads mapped to W and C strands for every 1KB

import pandas as pd
import numpy as np
import sys

if(len(sys.argv) >= 4):
    bedGraph_fwd=sys.argv[1]
    bedGraph_rev=sys.argv[2]
    output_file=sys.argv[3]
    binsize=sys.argv[4]
else:
    bedGraph_fwd="/Users/sarahkeegan/Downloads/aligned.fwd.200.bedGraph"
    bedGraph_rev = "/Users/sarahkeegan/Downloads/aligned.rev.200.bedGraph"
    output_file="/Users/sarahkeegan/Downloads/FT194.200.txt"
    binsize=200
    #print("Error: bad input.")
    #exit(0)

chr_lengths = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067,
  'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
  'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248,
  'chr19': 59128983, 'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566,}
chrom_list = range(1, 22+1, 1)
chrom_list = list( map(str, chrom_list) )
chrom_list.extend(['X','Y'])

def convert_bedGraph_to_txt(infile_fwd, infile_rev, outfile, bs, fw):
    df_w = pd.read_table(infile_fwd, names=['chr', 'start', 'end', 'cov'], header=None, low_memory=False)
    df_c = pd.read_table(infile_rev, names=['chr', 'start', 'end', 'cov'], header=None, low_memory=False)

    output_dicts = [{}, {}]
    for df_i, cur_df in enumerate([df_w, df_c]):
        print(df_i)
        prev_chrom = ''
        for i, row in enumerate(cur_df.iterrows()):
            row = row[1]
            chr_ = row['chr']
            if (chr_ != prev_chrom):
                print(chr_)
                prev_chrom = chr_
            st_pos = int(row['start'])
            end_pos = int(row['end'])
            cov = row['cov']

            for pos in range(st_pos + bs, end_pos + bs, bs):
                output_dicts[df_i][chr_ + '_' + str(pos)] = cov

    f = open(outfile, mode='w')
    f.write("chr\tpos\tw\tc\n")
    for chrom in chrom_list:
        end = int(round(chr_lengths['chr' + chrom], -3))
        while (chr_lengths['chr' + chrom] < end):
            end = end - bs
        for pos in range(bs, end + bs, bs):
            f.write(chrom + "\t" + str(pos) + '\t' + str(output_dicts[0]['chr' + chrom + '_' + str(pos)]) +
                '\t' + str(output_dicts[1]['chr' + chrom + '_' + str(pos)]) + '\n')

import os
(fp,fn)=os.path.split(output_file)

convert_bedGraph_to_txt(bedGraph_fwd, bedGraph_rev, output_file, int(binsize))















