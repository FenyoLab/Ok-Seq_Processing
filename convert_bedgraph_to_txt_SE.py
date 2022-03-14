# Performs all steps necessary to process OkSeq sequencing files (fastq PE)
# Result is density of reads mapped to W and C strands for every 1KB

import pandas as pd
import numpy as np
import sys

if(len(sys.argv) >= 4):
    bedGraph=sys.argv[1]
    output_file=sys.argv[2]
    binsize=sys.argv[3]
else:
    bedGraph="/Users/sarahkeegan/Downloads/perc_100.fwd.bedGraph"
    output_file="/Users/sarahkeegan/Downloads/perc_100.fwd.txt"
    binsize=500
    #print("Error: bad input.")
    #exit(0)

chr_lengths = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067,
  'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
  'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248,
  'chr19': 59128983, 'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566,}
chrom_list = range(1, 22+1, 1)
chrom_list = list( map(str, chrom_list) )
chrom_list.extend(['X','Y'])

def convert_bedGraph_to_txt_SE(infile, bs, fw):

    fw.write(infile+"\n")
    fw.flush()
    df = pd.read_csv(infile, sep='\t', header=None, names=['chr', 'start', 'end', 'signal'])
    df_to_save=pd.DataFrame()

    for chr in chr_lengths.keys():
        fw.write(chr+"\n")
        fw.flush()

        df_chr = df[df['chr'] == chr].copy()
        df_chr = df_chr[df_chr['signal'] > 0]
        df_chr.index = range(len(df_chr))

        chr_i = chr.lstrip('chr')
        if (chr_i == 'X'):
            chr_i = 23.0
        elif (chr_i == 'Y'):
            chr_i = 24.0
        else:
            chr_i = float(chr_i)

        # expand into dictionary
        data = {}
        #k_len = len(df_chr)
        for k, row in df_chr.iterrows():
            for l in range(row['start'] + 1, row['end'] + 1, 1):
                data[l] = row['signal']

        # sum read count over every 1000 base pairs
        cur_len = int(chr_lengths[chr] / bs) * bs
        data_kb = np.empty([int(cur_len / bs), 3], dtype=float)
        arr_i = 0
        arr_i_len = int(cur_len / bs)
        for kb_pos in range(0, cur_len, bs):
            if (arr_i % 10000 == 0):

                fw.write(str(arr_i) + ' ' + str(arr_i_len) + "\n")
                fw.flush()

            sum = 0.0
            for k in range(kb_pos + 1, kb_pos + bs + 1, 1):
                if (k in data):
                    sum += data[k]
                # else read count is 0 at this position

            data_kb[arr_i] = [chr_i, float(kb_pos + bs), sum]
            arr_i += 1

        # save to file
        to_save_temp = pd.DataFrame(data_kb, columns=['chr', 'pos', 'signal'])
        to_save_temp['pos'] = to_save_temp['pos'].astype('int')

        if (chr_i == 23.0):
            to_save_temp['chr'] = 'X'
        elif (chr_i == 24.0):
            to_save_temp['chr'] = 'Y'
        else:
            to_save_temp['chr'] = str(int(chr_i))

        df_to_save = df_to_save.append(to_save_temp, sort=True)
    return df_to_save

import os
(fp,fn)=os.path.split(output_file)
if(fp==""):
    logfile="convert_"+fn+".log.txt"
else:
    logfile=fp+'/'+"convert_"+fn+".log.txt"
with open(logfile,'w') as fw:
    df_to_save_fwd=convert_bedGraph_to_txt_SE(bedGraph, int(binsize), fw)

    #
    df_to_save_fwd.to_csv(output_file, index=False, sep='\t', columns=['chr', 'pos', 'signal'])















