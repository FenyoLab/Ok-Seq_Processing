# Performs all steps necessary to process OkSeq sequencing files (fastq PE)
# Result is density of reads mapped to W and C strands for every 1KB

import pandas as pd
import numpy as np
import sys

if(len(sys.argv) >= 4):
    w_file=sys.argv[1]
    c_file=sys.argv[2]
    output_file=sys.argv[3]
else:
    print("Error: bad input.")
    exit(0)

df_w=pd.read_csv(w_file, sep='\t', low_memory=False)
df_c=pd.read_csv(c_file, sep='\t', low_memory=False)
df_w['w']=df_w['signal']
df_w['c']=df_c['signal']

df_w.to_csv(output_file, index=False, sep='\t', columns=['chr', 'pos', 'w', 'c'])















