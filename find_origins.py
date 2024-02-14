# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 18:47:39 2017

@author: sarahkeegan, davidfenyo

this program will use the derivative of smoothed watson and crick strand signal to find origins based on okazaki
sequencing data: DER-SCORE = -1 * DER(W) * DER(C), where DER-SCORE > 0, an origin region is identified

input is tab delimited text files with position and watson and crick signals per 1kb

output is txt file with origin locations, origin midpoint

"""

import pandas as pd
import numpy as np
import sys
import os

if len(sys.argv) >= 3:
    filename = sys.argv[1]          # input file
    output_dir = sys.argv[2]        # output directory
else:
    print("Missing command line input.  Attempting to run with default settings.")
    filename = "/raw_files/rpe_edu_2.txt"
    output_dir = "/raw_files/results/origins"

############################

window = 'hanning'
width = 60  # 80

if window in ['hanning', 'hamming', 'bartlett', 'blackman']:
    w = eval('np.'+window+'(2*width+1)')
else:
    window = 'flat'
    w = np.ones(2*width+1)

peaks = {}
peaks_width = {}

print(f" {filename}: Finding origins...")

# read file contents into pandas DataFrame object, input columns: chr, pos, w, c
# each row is the signal values averaged over 1Kb
df = pd.read_table(filename, low_memory=False)

# add rfd
df['rfd']=df['c']-df['w']/(df['c']+df['w'])
df['rfd']=np.nan_to_num(df['rfd'])

# add perc_c
df['perc_c']=df['c']/(df['w']+df['c'])
df['perc_c']=np.nan_to_num(df['perc_c'])

# convolves values from c and w columns with the window function (for smoothing)
df['c_conv']=np.convolve(w/w.sum(),df['c'],'same')
df['w_conv']=np.convolve(w/w.sum(),df['w'],'same')

# calculate the 'derivative' of the crick/watson smoothed signal at each point
# the diff() function calculates the first order difference, i.e. for each point, diff[i] = x[i+1]-x[i]
# (this is the change in x, and the change in y is always 1)
df['c_diff']=df['c_conv'].diff()
df['w_diff']=df['w_conv'].diff()

# at initiation sites, crick should be increasing so if it is decreasing set to 0
df.loc[df['c_diff']<0,'c_diff']=0

# at initiation sites, watson should be decreasing so if it is increasing set to 0
df.loc[df['w_diff']>0, 'w_diff']=0

# calculate our 'Score' by multipying c and w derivatives (add minus sign to make it positive)
df['cw_diff_prod']=-1*df['c_diff']*df['w_diff']
df.loc[df['cw_diff_prod'] == -0, 'cw_diff_prod'] = 0

# loop through each 1kb position and look for peaks by finding max cw_diff_prod between zeros values
df.index = range(len(df))
df.fillna(0, inplace=True)
height = 0        # stores running max score found for current peak
height_idx = 0    # index of the current max score (cw_diff_prod)
start_idx = 0     # start index of the current 'peak'
previous_idx = 0

for idx in df.index:
    if idx > 0 and idx % 100000 == 0:
        print(idx)

    # idx is the index position in the DataFrame
    if df['cw_diff_prod'].iloc[idx] > 0:
        # if Score greater than 0, we are in a peak, continue to find max
        if height == 0: # at the beginning of a peak
            start_idx = previous_idx
        if height < df['cw_diff_prod'].iloc[idx]:
            # looking for max height in the peak
            height = df['cw_diff_prod'].iloc[idx]
            height_idx = idx

    else:  # if Score is 0, we are between peaks or at the end of a peak
        if height > 0:  # height has a value, we are at the end of a peak
            # store the index of the max height; the width is the start/end of the 0-values on either side of the peak
            peaks[(df['chr'].iloc[height_idx],
                   df['pos'].iloc[height_idx])] = height
            peaks_width[(df['chr'].iloc[height_idx],
                         df['pos'].iloc[height_idx])] = (df['pos'].iloc[start_idx],
                                                         df['pos'].iloc[idx])

        height = 0
        height_idx = 0
        start_idx = 0

    previous_idx = idx

# save origin locations and scores to output file
print("Saving origin locations...")
df_to_save = pd.DataFrame()
for chr_, peak in peaks:
    (start_pos1, end_pos1) = peaks_width[(chr_, peak)]
    df_to_save = df_to_save.append(dict(chr=chr_,
                                        pos=peak,
                                        height=peaks[(chr_, peak)],
                                        start=start_pos1,
                                        end=end_pos1,
                                        ), ignore_index=True)
fileroot = os.path.splitext(os.path.split(filename)[1])[0]
df_to_save.to_csv(output_dir + '/' + fileroot + '-' + window + str(width) + '-origins.txt', sep='\t')

