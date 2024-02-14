# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 18:47:39 2017

@author: sarahkeegan

this program will calculate the following scores for origin finding for each 1kb position of input file
input file is w/c okazaki read density for every 1kb

these "scores" are output for each origin:
 RFD
 Percent Crick
 OEM 20-40 and 20-60
 Efficiency 20-40 and 20-60
"""

import pandas as pd
import numpy as np
import sys
import os


if len(sys.argv) >= 3:
    filename = sys.argv[1]    # input file
    output_dir = sys.argv[2]  # output directory
    width = int(sys.argv[3])  # smoothing width
else:
    print("Missing command line input.  Attempting to run with default settings.")
    filename = '/Users/snk218/Dropbox (NYU Langone Health)/mac_files/fenyolab/data_and_results/huang/raw_files/rpe_edu_1.txt'
    output_dir = '/Users/snk218/Dropbox (NYU Langone Health)/mac_files/fenyolab/data_and_results/okazaki_origins/origins/test/'
    width = 60

print(filename)

window='hanning'
if window in ['hanning', 'hamming', 'bartlett', 'blackman']:
    w=eval('np.'+window+'(2*width+1)')
else:
    window='flat'
    w=np.ones(2*width+1)

print("Calculating scores...")

# read file contents into pandas DataFrame object, input columns: chr, pos, w, c
# each row is the signal values averaged over 1Kb
df = pd.read_table(filename, low_memory=False)
df.index = range(len(df))

# add rfd
df['rfd'] = df['c'] - df['w'] / (df['c'] + df['w'])
df['rfd'] = np.nan_to_num(df['rfd'])

# smoothed rfd, and delta rfd
width_ = 15
df[f'w_sm{width_}'] = np.convolve(df['w'], np.ones(width_), 'same')
df[f'c_sm{width_}'] = np.convolve(df['c'], np.ones(width_), 'same')
df[f'rfd_sm{width_}'] = (df[f'c_sm{width_}'] - df[f'w_sm{width_}']) / (df[f'c_sm{width_}'] + df[f'w_sm{width_}'])
df[f'rfd_sm{width_}'] = np.nan_to_num(df[f'rfd_sm{width_}'])
df[f'delta_rfd_sm{width_}'] = df[f'rfd_sm{width_}'].diff()

# add perc_c
#df['perc_c'] = df['c'] / (df['w'] + df['c'])
#df['perc_c'] = np.nan_to_num(df['perc_c'])

#add DER-SCORE
# convolves values from c and w columns with the window function (for smoothing)
df['c_conv'] = np.convolve(w / w.sum(), df['c'], 'same')
df['w_conv'] = np.convolve(w / w.sum(), df['w'], 'same')

# calculate the 'derivative' of the crick/watson smoothed signal at each point
# the diff() function calculates the first order difference, i.e. for each point, diff[i] = x[i+1]-x[i]
# (this is the change in x, and the change in y is always 1)
df['c_diff_ini'] = df['c_conv'].diff()
df['w_diff_ini'] = df['w_conv'].diff()

df['c_diff_term'] = df['c_conv'].diff()
df['w_diff_term'] = df['w_conv'].diff()

# at initiation sites, crick should be increasing so if it is decreasing set to 0
df.loc[df['c_diff_ini'] < 0, 'c_diff_ini'] = 0

# at initiation sites, watson should be decreasing so if it is increasing set to 0
df.loc[df['w_diff_ini'] > 0, 'w_diff_ini'] = 0

# calcuate 'Score' by multipying c and w derivatives (add minus sign to make it positive)
df['der_score_ini'] = -1 * df['c_diff_ini'] * df['w_diff_ini']
df.loc[df['der_score_ini'] == -0, 'der_score_ini'] = 0

# at termination sites, crick should be decreasing so if it is increasing set to 0
df.loc[df['c_diff_term'] > 0, 'c_diff_term'] = 0

# at termination sites, watson should be increasing so if it is decreasing set to 0
df.loc[df['w_diff_term'] < 0, 'w_diff_term'] = 0

# calcuate 'Score' by multipying c and w derivatives (add minus sign to make it positive)
df['der_score_term'] = -1 * df['c_diff_term'] * df['w_diff_term']
df.loc[df['der_score_term'] == -0, 'der_score_term'] = 0

# eff_20_40 and eff_20_60
# percentage crick averaged over window of width 20 or 40 (width_) on rhs minus lhs, starting at 20 from center (offset)
#offset=20
#for width_ in [20,40]:
#    w1=np.concatenate(((1./width_)*np.ones(width_), np.zeros(offset*2), (-1/width_)*np.ones(width_))) # (convolve will flip array)
#    df['eff_'+str(offset)+'_'+str(offset+width_)]=np.convolve(df['perc_c'], w1, 'same')
#    df['eff_'+str(offset)+'_'+str(offset+width_)]=np.nan_to_num(df['eff_'+str(offset)+'_'+str(offset+width_)])

# oem3 - go out 20 from center (offset) to L and R, average next 20 or 40 (width_) for wL/wR and cL/cR
# oem = wL / (wL+cL) - wR / (wR+cR)
offset=0 #20
for width_ in [10,20,40]:
    window_L=np.concatenate((np.zeros(offset+width_), np.zeros(offset), np.ones(width_))) # (convolve will flip array)
    window_R=np.concatenate((np.ones(width_), np.zeros(offset), np.zeros(offset+width_)))
    df['wL']=np.convolve(df['w'], window_L, 'same')
    df['cL']=np.convolve(df['c'], window_L, 'same')
    df['wR']=np.convolve(df['w'], window_R, 'same')
    df['cR']=np.convolve(df['c'], window_R, 'same')
    df['oem_'+str(offset)+'_'+str(offset+width_)]=df['wL']/(df['wL']+df['cL']) - df['wR']/(df['wR']+df['cR'])
    df['oem_'+str(offset)+'_'+str(offset+width_)] = np.nan_to_num(df['oem_'+str(offset)+'_'+str(offset+width_)])

# save origin locations and scores to output file
fileroot = os.path.splitext(os.path.split(filename)[1])[0]

df.fillna(0, inplace=True)
print("Saving scores...")
df.to_csv(output_dir + '/' + fileroot+'-'+window+str(width) + '-scores.txt', sep='\t',
        #columns=('chr', 'pos', 'w', 'c', 'der_score', 'rfd', 'perc_c', 'oem_20_40', 'oem_20_60', 'eff_20_40', 'eff_20_60'))
        columns = ('chr', 'pos',
                   'w', 'c',
                   'w_conv',
                   'c_conv',
                   'der_score_ini',
                   'der_score_term',
                   'rfd',
                   'rfd_sm15',
                   'delta_rfd_sm15',
                   'oem_0_10', 'oem_0_20', 'oem_0_40'))

# break into individual bedGraph files - make crick minus for display
df['st'] = df['pos']-1000
df['w'] = -1 * df['w']
df['der_score_ini'] = np.log10(df['der_score_ini'] * 1000000 + 1)
df['der_score_term'] = -1 * np.log10(df['der_score_term'] * 1000000 + 1)
df.loc[df['der_score_term'] == -0, 'der_score_term'] = 0
df.fillna(0, inplace=True)
df.replace([np.inf, -np.inf], 0, inplace=True)
for col in ['w', 'c',
            'w_conv',
            'c_conv',
            'der_score_ini',
            'der_score_term',
            'rfd',
            'rfd_sm15',
            'delta_rfd_sm15',
            'oem_0_10', 'oem_0_20', 'oem_0_40']:
    print(col)
    print(np.nanpercentile(df[col], q=[10,25,50,75,90]))
    df.to_csv(f"{output_dir}/{fileroot}-{col}.bedGraph",
              sep='\t', header=False, index=False, columns=('chr', 'st', 'pos', col))
