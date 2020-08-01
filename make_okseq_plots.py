# make strand bias plots from data tables

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib.pyplot as plt

matplotlib.rc('pdf', fonttype=42)
matplotlib.rc('ps', fonttype=42)
matplotlib.rc('font', size=6)          # controls default text sizes
matplotlib.rc('axes', titlesize=6)     # fontsize of the axes title
matplotlib.rc('axes', labelsize=6)    # fontsize of the x and y labels
matplotlib.rc('xtick', labelsize=6)    # fontsize of the tick labels
matplotlib.rc('ytick', labelsize=6)    # fontsize of the tick labels
matplotlib.rc('legend', fontsize=6)    # legend fontsize
matplotlib.rc('figure', titlesize=6)  # fontsize of the figure title
matplotlib.rc('axes',linewidth=0.6)
matplotlib.rc('xtick.major',width=0.6)
matplotlib.rc('ytick.major',width=0.6)

mult=1
fig_size_width=1.88 * mult
fig_size_height=1.78 * mult
fig_size_panel_height=5.55
fig_size_panel_width=6.25

base_dir = "/Users/sarahkeegan/okseq_data/strand_bias_plots"
length_around_site = 200
ext='pdf'

rep_file_name = "rpe_edu_2"

def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def get_wc_plot_data(data_df, site):

    plus_strand_df = data_df[data_df['hg19.knownGene.strand'] == '+']
    w_data = plus_strand_df.loc[:,'dist_w_100_from_' + site + '_1':'dist_w_100_from_' + site + '_' + str(length_around_site + 1)].mean()
    c_data = plus_strand_df.loc[:,'dist_c_100_from_' + site + '_1':'dist_c_100_from_' + site + '_' + str(length_around_site + 1)].mean()
    w_data = np.asarray(w_data)
    c_data = np.asarray(c_data)
    ave_perc_c_plus = (c_data / (w_data + c_data))*100

    minus_strand_df = data_df[data_df['hg19.knownGene.strand'] == '-']
    w_data = minus_strand_df.loc[:,'dist_w_100_from_' + site + '_1':'dist_w_100_from_' + site + '_' + str(length_around_site + 1)].mean()
    c_data = minus_strand_df.loc[:,'dist_c_100_from_' + site + '_1':'dist_c_100_from_' + site + '_' + str(length_around_site + 1)].mean()
    w_data = np.asarray(w_data)
    c_data = np.asarray(c_data)
    ave_perc_w_minus = (w_data / (w_data + c_data))*100
    ave_perc_w_minus_flip = np.flip(ave_perc_w_minus, axis=0)

    if (len(data_df[data_df['hg19.knownGene.strand'] == '-']) == 0):
        return ave_perc_c_plus
    if (len(data_df[data_df['hg19.knownGene.strand'] == '+']) == 0):
        return ave_perc_w_minus_flip
    else:
        return np.mean([ave_perc_c_plus, ave_perc_w_minus_flip], axis=0)

def remove_zeros(df_to_rem):
    df_to_rem = df_to_rem[(df_to_rem['c_sum_TSS_to_TTS'] != 0) & (df_to_rem['w_sum_TSS_to_TTS'] != 0)]
    return df_to_rem.copy()

def load_data():
    # load file with data
    df = pd.read_table(base_dir + '/' + "sites_table_with_distributions-" + rep_file_name + ".csv", sep=',', low_memory=False)

    #remove genes with n/a for any cols
    col_list_w = []
    col_list_c = []
    for site in ['TSS','TTS']:
        for i in range(1, length_around_site + 1 + 1, 1):
            col_list_w.append('dist_w_100_from_' + site + '_' + str(i))
            col_list_c.append('dist_c_100_from_' + site + '_' + str(i))

    df.dropna(axis=0, subset=col_list_w,inplace=True)
    df.dropna(axis=0, subset=col_list_c, inplace=True)

    print("Loaded 'sites_table_with_distributions-" + rep_file_name + ".csv'")
    print(len(df))

    return df.copy()

def make_plot(data_to_plot, colors_to_plot, plot_title, filename, sub_dir='',xlim=(-50, 50),xticks=[-40,-20,0,20,40],
              ylim=(29, 71),yticks=[30,40,50,60,70],yline=50,plot_yline=True,legend_labels=[],xtick_labels=[],
              xlabel='Distance to TSS (kb), Txn L-->R',ylabel='% of forks moving L-->R', legend_loc='lower right'):
    plt.figure(figsize=(fig_size_width,fig_size_height), dpi=300)

    for i,data in enumerate(data_to_plot):
        color=colors_to_plot[i]
        if(len(legend_labels)>0):
            lab = legend_labels[i]
        else:
            lab = ''
        plt.plot(range(int(-1 * length_around_site / 2), int(length_around_site / 2 + 1)), data, '-', color=color,
                 label=lab, linewidth=1.25,)

    if (len(ylim) == 2):
        plt.ylim(ylim)
        if(len(yticks)>0):
            plt.yticks(yticks)

    if (len(xlim) == 2):
        plt.xlim(xlim)
        if (len(xticks) > 0):
            plt.xticks(xticks)
        if(len(xtick_labels) > 0):
            plt.xticklabels(xtick_labels)

    plt.title(plot_title, x=0.02, y=0.9, ha='left', va='top')
    if(plot_yline==True):
        plt.axhline(y=yline, color='black',linewidth=0.8, dashes=[2, 4])
    plt.axvline(linestyle='--', color='black',linewidth=0.8, dashes=[2, 4])
    if (len(legend_labels) > 0):
        if (legend_loc == 'lower right'):
            legend = plt.legend(loc='lower right', fancybox=False, bbox_to_anchor=(1.03, -0.035))
        elif (legend_loc == 'upper right'):
            legend = plt.legend(loc='upper right', fancybox=False, bbox_to_anchor=(1.03, 1.035))
        else:
            legend = plt.legend(loc=legend_loc)
        frame = legend.get_frame()
        frame.set_linewidth(0.6)
        frame.set_edgecolor('black')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout(pad=2)
    plt.savefig(base_dir + '/' + sub_dir + '/'+filename+'.'+ext, transparent=True)
    plt.clf()
    plt.close()

def make_panels_plot(data_to_plot, plot_title, filename, sub_dir='',yline=50,line_colors=['black',], legend_labels=[],
                     vert_plot=True, xlim=(-50, 50), xticks=[-40, -20, 0, 20, 40], ylim=(29, 71), yticks=[30, 40, 50, 60, 70],
                     xlabel='Distance to TSS (kb), Txn L-->R', ylabel='% of forks moving L-->R',legend_loc='lower right',
                     fig_size=None):
    if(vert_plot):
        if(fig_size):
            fig, axarr = plt.subplots(nrows=len(data_to_plot), ncols=1, sharex=True,figsize=(fig_size_width, fig_size))
        else:
            fig, axarr = plt.subplots(nrows=len(data_to_plot),ncols=1, sharex=True,
                                      figsize=(fig_size_width, fig_size_panel_height))
    else:
        #horizontal plot
        if(fig_size):
            fig, axarr = plt.subplots(ncols=len(data_to_plot),nrows=1, sharey=True, figsize=(fig_size,fig_size_height))
        else:
            fig, axarr = plt.subplots(ncols=len(data_to_plot), nrows=1, sharey=True,
                                      figsize=(fig_size_panel_width, fig_size_height))

    for i,data_list in enumerate(data_to_plot):
        for j,data in enumerate(data_list):
            if(len(legend_labels) > 0 and len(legend_labels[i]) > 0):
                lab=legend_labels[i][j]
            else:
                lab=''

            axarr[i].plot(range(int(-1 * length_around_site / 2), int(length_around_site / 2 + 1)), data, '-',
                          color=line_colors[j],label=lab,linewidth=0.75,)

            if(len(ylim)==2):
                axarr[i].set_ylim(ylim)
                if (len(yticks) > 0):
                    if(vert_plot):
                        axarr[i].set_yticks(yticks)
                    else:
                        if (i == len(axarr) - 1): axarr[i].set_yticks(yticks)
            if(len(xlim)==2):
                axarr[i].set_xlim(xlim)
                if (len(xticks) > 0):
                    if(vert_plot):
                        if (i == len(axarr) - 1): axarr[i].set_xticks(xticks)
                    else:
                        axarr[i].set_xticks(xticks)

            axarr[i].set_title(plot_title[i],x=0.02, y=0.9, ha='left', va='top')
            axarr[i].axhline(y=yline, color='black',linewidth=0.8, dashes=[2, 4])
            axarr[i].axvline(color='black',linewidth=0.8, dashes=[2, 4])
            if(i == len(axarr)-1): axarr[i].set_xlabel(xlabel)
            axarr[i].set_ylabel(ylabel)
            if (len(legend_labels) > 0 and len(legend_labels[i]) > 0):
                if(legend_loc=='lower right'):
                    legend = axarr[i].legend(loc='lower right', fancybox=False, bbox_to_anchor=(1.03, -0.035))
                elif(legend_loc == 'upper right'):
                    legend = axarr[i].legend(loc='upper right', fancybox=False, bbox_to_anchor=(1.03, 0.93))
                else:
                    legend = axarr[i].legend(loc=legend_loc)
                frame = legend.get_frame()
                frame.set_linewidth(0.6)
                frame.set_edgecolor('black')
    if (vert_plot):
        plt.tight_layout(pad=2,h_pad=1)
    else:
        plt.tight_layout(pad=2, w_pad=1)
    fig.savefig(base_dir + '/' + sub_dir + '/' + filename + '.' + ext)
    fig.clf()

def make_panels_plot_with_data(data_dfs, filter_col, quartile_col, filename, sub_dir='', xlim=(-50, 50),
                               xticks=[-40, -20, 0, 20, 40],ylim=(29,71),yticks=[30, 40, 50, 60, 70],
                               yline=50, line_colors=['black',], plot_grad=False,grad_filename='',
                               grad_ylim=(-0.5,2.0),grad_yticks=[], col_to_title={},
                               xlabel='Distance to TSS (kb), Txn L->R',legend_labels=[], grad_legend_loc='lower right',
                               save_grad_data=False, save_grad_data_names=[],single_plot=False):

    plot_data_list = [[],[],[],[]]
    to_save_list = []
    to_save_names=[]
    plot_data_grad_list = [[],[],[],[]]
    plot_title_list = []
    q_list = ['4', '3', '2', '1']

    # filter first the dfs
    for data_df_i, data_df in enumerate(data_dfs):
        if (filter_col != ''):
            filter_data = data_df[filter_col].dropna()
            med = np.median(filter_data)
            data_df_filt = data_df[data_df[filter_col] > med]
        else:
            data_df_filt = data_df

        # second break up by quartiles
        qt_data = data_df_filt[quartile_col].dropna()
        cutoffs = np.percentile(qt_data, [25, 50, 75])
        print(filename, ' ', cutoffs)
        data_limits = [[cutoffs[2], max(qt_data)+1],
                       [cutoffs[1], cutoffs[2]],
                       [cutoffs[0], cutoffs[1]],
                       [-1, cutoffs[0]]]

        q1_zero = False
        cutoff_df = data_df_filt[(data_df_filt[quartile_col]>=data_limits[3][0]) & (data_df_filt[quartile_col]<data_limits[3][1])]
        if(len(cutoff_df)==0):
            q1_zero=True

        # then add to data list
        for i, cur_limits in enumerate(data_limits):
            if(q1_zero):
                cutoff_df = data_df_filt[data_df_filt[quartile_col] > cur_limits[0]]
                cutoff_df = cutoff_df[cutoff_df[quartile_col] <= cur_limits[1]]
            else:
                cutoff_df = data_df_filt[data_df_filt[quartile_col] >= cur_limits[0]]
                cutoff_df = cutoff_df[cutoff_df[quartile_col] < cur_limits[1]]

            print(len(cutoff_df))

            if(save_grad_data):
                to_save_list.append(cutoff_df)
                if(len(save_grad_data_names) > 0):
                    to_save_names.append(save_grad_data_names[data_df_i] + " Q" + q_list[i] + " " + quartile_col)
                else:
                    to_save_names.append("Q" + q_list[i] + " " + quartile_col)

            to_plot = get_wc_plot_data(cutoff_df, site)

            plot_data_list[i].append(to_plot)
            if(plot_grad):
                plot_data_grad_list[i].append(smooth(np.gradient(to_plot),3))

            if(data_df_i == 0):
                if (col_to_title != {}):
                    col_str = col_to_title[quartile_col]
                else:
                    col_str = quartile_col

                if(filter_col != ''):
                    plot_title_list.append("Q" + q_list[i] + " " + col_str + "\n"+filter_col+" > median\nn=" + str(len(cutoff_df)))
                else:
                    plot_title_list.append("Q" + q_list[i] + " " + col_str + "\nn=" + str(len(cutoff_df)))
    if (filter_col != ''):
        title_ = col_str + "\n" + filter_col + " > median\n"
    else:
        title_ = col_str

    if(single_plot):
        plot_data_list_ = [plot_data_list[0][0],plot_data_list[1][0],plot_data_list[2][0],plot_data_list[3][0]]
        make_plot(plot_data_list_,colors_to_plot=line_colors,plot_title=title_,filename=filename,sub_dir=sub_dir,
                  xlim=xlim,xticks=xticks,ylim=ylim,yticks=yticks,yline=yline,xlabel=xlabel, legend_labels=legend_labels)
    else:
        make_panels_plot(plot_data_list, plot_title=plot_title_list, filename=filename, sub_dir=sub_dir, xlim=xlim, xticks=xticks,
                         ylim=ylim,yticks=yticks, yline=yline, line_colors=line_colors, xlabel=xlabel, legend_labels=legend_labels)
    if(plot_grad):
        if(single_plot):
            plot_data_grad_list_ = [plot_data_grad_list[0][0], plot_data_grad_list[1][0], plot_data_grad_list[2][0], plot_data_grad_list[3][0]]
            make_plot(plot_data_grad_list_,colors_to_plot=line_colors,plot_title=title_,filename=grad_filename,sub_dir=sub_dir,
                      xlim=xlim, ylim=grad_ylim, yline=0, yticks=grad_yticks,xlabel=xlabel,xticks=xticks,
                      legend_labels=legend_labels, legend_loc=grad_legend_loc,ylabel='% increase in L->R forks')
        else:
            make_panels_plot(plot_data_grad_list, plot_title=plot_title_list, filename=grad_filename, sub_dir=sub_dir,
                         xlim=xlim, ylim=grad_ylim, yline=0,line_colors=line_colors, yticks=grad_yticks, xlabel=xlabel,
                         xticks=xticks,legend_labels=legend_labels, legend_loc=grad_legend_loc,ylabel='% increase in L->R forks')


#######################
#plot strand bias plots
df_rpe_edu_2 = load_data()
df_rpe_edu_2['length'] = np.abs(df_rpe_edu_2['hg19.knownGene.txEnd'] - df_rpe_edu_2['hg19.knownGene.txStart'])
df_rpe_edu_2['volume'] = df_rpe_edu_2['fpkm'] * df_rpe_edu_2['length']
site='TSS' # 'TTS'

#Fig 1
#1d
w_data = df_rpe_edu_2.loc[:,'dist_w_100_from_' + site + '_1':'dist_w_100_from_' + site + '_' + str(length_around_site + 1)].mean()
c_data = df_rpe_edu_2.loc[:,'dist_c_100_from_' + site + '_1':'dist_c_100_from_' + site + '_' + str(length_around_site + 1)].mean()
w_data = np.asarray(w_data)
c_data = np.asarray(c_data)
ave_perc_c = (c_data / (w_data + c_data))*100
make_plot([ave_perc_c,],colors_to_plot=['black'],plot_title='',filename='Fig1d')

#1e
plus_strand_df = df_rpe_edu_2[df_rpe_edu_2['hg19.knownGene.strand'] == '+']
w_data = plus_strand_df.loc[:,'dist_w_100_from_' + site + '_1':'dist_w_100_from_' + site + '_' + str(length_around_site + 1)].mean()
c_data = plus_strand_df.loc[:,'dist_c_100_from_' + site + '_1':'dist_c_100_from_' + site + '_' + str(length_around_site + 1)].mean()
w_data = np.asarray(w_data)
c_data = np.asarray(c_data)
ave_perc_c_plus = (c_data / (w_data + c_data))*100

minus_strand_df = df_rpe_edu_2[df_rpe_edu_2['hg19.knownGene.strand'] == '-']
w_data = minus_strand_df.loc[:,'dist_w_100_from_' + site + '_1':'dist_w_100_from_' + site + '_' + str(length_around_site + 1)].mean()
c_data = minus_strand_df.loc[:,'dist_c_100_from_' + site + '_1':'dist_c_100_from_' + site + '_' + str(length_around_site + 1)].mean()
w_data = np.asarray(w_data)
c_data = np.asarray(c_data)
ave_perc_c_minus = (c_data / (w_data + c_data))*100
make_plot([ave_perc_c_plus,ave_perc_c_minus],colors_to_plot=['blue','magenta'],plot_title='',filename='Fig1e')

#1f
ave_perc_w_minus = (w_data / (w_data + c_data))*100
ave_perc_w_minus_flip = np.flip(ave_perc_w_minus, axis=0)
to_plot=np.mean([ave_perc_c_plus, ave_perc_w_minus_flip], axis=0)
make_plot([to_plot,],colors_to_plot=['black'],plot_title='',filename='Fig1f')

#Fig 2
#2a
make_panels_plot_with_data([df_rpe_edu_2,], filter_col='', quartile_col='fpkm', filename='Fig2a',
                           line_colors=['black',], plot_grad=False, )
#2e
make_panels_plot_with_data([df_rpe_edu_2,], filter_col='fpkm', quartile_col='length', filename='Fig2e',
                           line_colors=['black',], plot_grad=False, ylim=(29,81),yticks=[30, 40, 50, 60, 70,80],)
#2g
make_panels_plot_with_data([df_rpe_edu_2,], filter_col='', quartile_col='volume', filename='Fig2g',
                           line_colors=['black',], plot_grad=False, ylim=(29,81),yticks=[30, 40, 50, 60, 70,80],)



