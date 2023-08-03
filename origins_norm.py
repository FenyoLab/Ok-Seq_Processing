import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import seaborn as sns
from scipy.stats import gaussian_kde
import os

if len(sys.argv) >= 3:
    filename = sys.argv[1]  # input file
    output_dir = sys.argv[2]  # output directory
else:
    print("Missing command line input.  Attempting to run with default settings.")
    filename = "/Users/snk218/Dropbox (NYU Langone Health)/mac_files/fenyolab/data_and_results/huang/raw_files/results/origins/rpe_edu_1-hanning60-origins.txt"
    output_dir = "/Users/snk218/Dropbox (NYU Langone Health)/mac_files/fenyolab/data_and_results/huang/raw_files/results/origins"

# read origins file
df = pd.read_csv(filename, sep='\t', low_memory=False, index_col=0)
print(len(df))

fileroot = os.path.splitext(os.path.split(filename)[1])[0]

id = []
chr_coor = []
heights = []
bins = []

# log normalize DER score (height) to view data distribution
#df['height'] = scipy.stats.boxcox(df['height'], lmbda=0.1, alpha=None)
df['height-log'] = np.log10(df['height'])
#df['height-log'] = (df['height'])

# snsplot of normalized origins
#sns.kdeplot(df['height-log'])
#plt.savefig(output_dir + '/' + fileroot + '-kde.png')

# calculate kde gaussian without plot to find fwhm
kde = gaussian_kde(df['height-log'])
x = np.linspace(df['height-log'].min(), df['height-log'].max(), len(df))
y = kde(x)

plt.plot(x, y, color='black')

halfmax = y.max()/2
maxpos = y.argmax()
leftpos = (np.abs(y[:maxpos] - halfmax)).argmin()
rightpos = (np.abs(y[maxpos:] - halfmax)).argmin() + maxpos
fullwidthhalfmax = x[rightpos] - x[leftpos]
print(x[leftpos])

plt.axvline(x[leftpos])

plt.savefig(output_dir + '/' + fileroot + '-kde.png')
plt.clf()

# use position of max height to define origin +/-10kb on either side
df['pos-up'] = df['pos'] + 10000
df['pos-down'] = df['pos'] - 10000

df_ = df[df['height-log'] >= x[leftpos]]
print(len(df_))

# calculate kde gaussian without plot to find fwhm
kde = gaussian_kde(df_['height-log'])
x = np.linspace(df_['height-log'].min(), df_['height-log'].max(), len(df_))
y = kde(x)
plt.plot(x, y, color='black')
plt.savefig(output_dir + '/' + fileroot + '-kde-filtered.png')
plt.clf()

df_.index = range(len(df_))
df_.to_csv(output_dir + '/' + fileroot + '-log10-filtered.txt', sep='\t', index=True)


