# make strand bias plots from data tables

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib.pyplot as plt
import seaborn as sns
from skimage import exposure
from scipy import stats

def set_rcs_normal():
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

def set_rcs_heatmaps():
    matplotlib.rc('pdf', fonttype=42)
    matplotlib.rc('ps', fonttype=42)
    matplotlib.rc('font', size=12)  # controls default text sizes
    matplotlib.rc('axes', titlesize=12)  # fontsize of the axes title
    matplotlib.rc('axes', labelsize=12)  # fontsize of the x and y labels
    matplotlib.rc('xtick', labelsize=12)  # fontsize of the tick labels
    matplotlib.rc('ytick', labelsize=12)  # fontsize of the tick labels
    matplotlib.rc('legend', fontsize=12)  # legend fontsize
    matplotlib.rc('figure', titlesize=14)  # fontsize of the figure title
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

set_rcs_normal()

#plot_line_width=1
#plot_line_alpha=1
fig_size_width=1.88
fig_size_height=1.78
fig_size_panel_height=5.55
fig_size_panel_width=6.25

base_dir = "."
length_around_site = 200
ext='pdf'
site='TSS'
rep_file_name = ""

def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def get_wc_plot_data(data_df):

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

def get_heatmap_data(data_df, around_mid=50):
    heatmap_data = np.zeros(shape=(len(data_df), around_mid * 2 + 1))
    data_df = data_df.copy()
    data_df.index = range(len(data_df))
    for i, row in enumerate(data_df.iterrows()):
        row = row[1]
        strand = row['hg19.knownGene.strand']

        st_pos = str(int(101-around_mid))
        end_pos = str(int(101+around_mid))

        if (strand == '+'):
            w_data = np.asarray(
                data_df.loc[i, 'dist_w_100_from_' + site + '_' + st_pos:'dist_w_100_from_' + site + '_' + end_pos])
            c_data = np.asarray(
                data_df.loc[i, 'dist_c_100_from_' + site + '_' + st_pos:'dist_c_100_from_' + site + '_' + end_pos])
            perc_c = c_data / (w_data + c_data)
            heatmap_data[i] = perc_c
        else:
            w_data = np.asarray(
                data_df.loc[i, 'dist_w_100_from_' + site + '_' + st_pos:'dist_w_100_from_' + site + '_' + end_pos])
            c_data = np.asarray(
                data_df.loc[i, 'dist_c_100_from_' + site + '_' + st_pos:'dist_c_100_from_' + site + '_' + end_pos])
            perc_w = w_data / (w_data + c_data)
            perc_w = np.flip(perc_w, axis=0)
            heatmap_data[i] = perc_w
    heatmap_data[np.isnan(heatmap_data)] = 0
    heatmap_data[np.isinf(heatmap_data)] = 0
    return 100 * heatmap_data

def load_data(sites_rm_na=('TSS','TTS')):
    # load file with data
    df = pd.read_table(base_dir + '/' + "sites_table_with_distributions-" + rep_file_name + ".csv", sep=',', low_memory=False)

    #remove genes with n/a for any cols
    col_list_w = []
    col_list_c = []
    for cur_site in sites_rm_na:
        for i in range(1, length_around_site + 1 + 1, 1):
            col_list_w.append('dist_w_100_from_' + cur_site + '_' + str(i))
            col_list_c.append('dist_c_100_from_' + cur_site + '_' + str(i))

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
                               yline=50, line_colors=['black',], plot_grad=False, grad_filename='',
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

            to_plot = get_wc_plot_data(cutoff_df)

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

def make_heatmap_plot_with_data(data_dfs, filter_col, sort_col, filename, sub_dir='', around_mid=50,
                                range_limit=(30,70), color='default',
                                range_limit_ticks=(30, 40, 50, 60, 70)):
    set_rcs_heatmaps()

    if(len(data_dfs)==1):
        data_df=data_dfs[0]
        # filter
        if (filter_col != ''):
            filter_data = data_df[filter_col].dropna()
            med = np.median(filter_data)
            data_df_filt = data_df[data_df[filter_col] > med]
            title_pre = 'above median '+filter_col+', '
        else:
            data_df_filt = data_df
            title_pre=''

        data_df_filt_sorted = data_df_filt.sort_values(by=sort_col)
        plot_data = get_heatmap_data(data_df_filt_sorted, around_mid=around_mid)
    else:
        # must subtract the 2 data sets
        plot_data_to_sub=[]
        for data_df in data_dfs:
            # filter
            if (filter_col != ''):
                filter_data = data_df[filter_col].dropna()
                med = np.median(filter_data)
                data_df_filt = data_df[data_df[filter_col] > med]
                title_pre = 'above median ' + filter_col + ', '
            else:
                data_df_filt = data_df
                title_pre = ''
            data_df_filt_sorted = data_df_filt.sort_values(by=sort_col)
            plot_data = get_heatmap_data(data_df_filt_sorted, around_mid=around_mid)

            plot_data_to_sub.append(plot_data)
        plot_data = plot_data_to_sub[0]-plot_data_to_sub[1]

    #flip the plot data
    plot_data = np.flip(plot_data, axis=0)

    fig, ax_ = plt.subplots(1, figsize=(6, 6))

    if(color=='default'):
        sns.heatmap(plot_data, ax=ax_)
    else:
        sns.heatmap(plot_data, cmap=color, ax=ax_)

    ax_.set_title(title_pre + sort_col + ' sorted\nn=' + str(len(data_df_filt)), loc='left')
    ax_.tick_params(labelleft='off', left=False)

    if (around_mid == 50):
        bounds = 50
        ts=[10,30,50,70,90]
        tls=[-40,-20,0,20,40]
        ax_.set_xticks(ts)
        ax_.set_xticklabels(tls)
    elif(around_mid == 100):
        bounds = 100
        ts = [20, 40, 60, 80, 100, 120, 140, 160, 180]
        tls = [-80, -60, -40,-20, 0, 20, 40, 60, 80]
        ax_.set_xticks(ts)
        ax_.set_xticklabels(tls)

    if(len(range_limit) == 0):
        # do not rescale intensity, range is 0-100% for strand bias
        print("Saving heatmap:", filename)
        plt.tight_layout()
        fig.savefig(base_dir + '/' + sub_dir + '/' + filename + '.' + ext)
        fig.clf()
    else:
        # rescale to given limit
        # k,limit in enumerate(range_limits):
        # output rescaled version of same plot...
        rescl_plot_data = exposure.rescale_intensity(plot_data, in_range=(range_limit[0],range_limit[1]), out_range=(range_limit[0],range_limit[1]))

        fig, ax_ = plt.subplots(1, figsize=(6, 6))
        if(color == 'default'):
            sns.heatmap(rescl_plot_data, ax=ax_, cbar_kws={'ticks':range_limit_ticks,})
        else:
            sns.heatmap(rescl_plot_data, cmap=color, ax=ax_, cbar_kws={'ticks': range_limit_ticks, })

        ax_.set_title(title_pre + sort_col + ' sorted\nn=' + str(len(data_df_filt)), loc='left')
        ax_.tick_params(labelleft='off', left=False)

        if (around_mid == 50):
            ts = [10, 30, 50, 70, 90]
            tls = [-40, -20, 0, 20, 40]
            ax_.set_xticks(ts)
            ax_.set_xticklabels(tls)
        elif(around_mid == 100):
            ts = [20, 40, 60, 80, 100, 120, 140, 160, 180]
            tls = [-80, -60, -40, -20, 0, 20, 40, 60, 80]
            ax_.set_xticks(ts)
            ax_.set_xticklabels(tls)

        limit_str = str(range_limit[0])+'_'+str(range_limit[1])
        print("Saving heatmap:", filename + '_' + str(limit_str))
        plt.tight_layout()
        fig.savefig(base_dir + '/' + sub_dir + '/' + filename + '_' + str(limit_str) + '.' + ext)
        fig.clf()
    set_rcs_normal()

def adjust_len(in_arr, new_len):
    new_arr=[]
    old_len = len(in_arr)
    for pos in range(new_len):
        in_arr_pos = round((pos/new_len)*old_len)
        if(in_arr_pos > (old_len-1)):
            in_arr_pos = old_len-1
        new_arr.append(in_arr[int(in_arr_pos)])
    return new_arr

def make_TSS_to_TTS_plots(df_list, colors=['black',], fn_ext='', min_len=5000, max_len=100000,
                          norm_gene_len=50, append_amt=20,xlim=(),
                          ylim=(29,71),xticks=[],yticks=[30, 40, 50, 60, 70],
                          xtick_labels=[], yline=50, der_ylim=(-2.0, 3.0), der_yticks=[], der_yline=0):
    fig_sb = plt.figure(figsize=(fig_size_width, fig_size_height), dpi=300)
    ax_sb = fig_sb.gca()
    fig_der = plt.figure(figsize=(fig_size_width, fig_size_height), dpi=300)
    ax_der = fig_der.gca()

    for df_i,df in enumerate(df_list):

        df.dropna(axis=0, subset=['fpkm', ], inplace=True)
        #df['length'] = np.abs(df['hg19.knownGene.txEnd'] - df['hg19.knownGene.txStart'])

        df['txStart'] = np.round(df['hg19.knownGene.txStart'], -3)
        df['txEnd'] = np.round(df['hg19.knownGene.txEnd'], -3)
        df['length2'] = np.abs(df['txEnd'] - df['txStart']) + 1000 #make sure to include both TSS and TTS

        fpkm_med = np.median(df['fpkm'])
        df = df[df['fpkm']>fpkm_med]

        df = df[df['length2'] >= min_len]
        df = df[df['length2']<=max_len].copy()

        df_plus_len = len(df[df['hg19.knownGene.strand']=='+'])
        df_minus_len = len(df[df['hg19.knownGene.strand'] == '-'])

        plus_strand_c_data_matrix = np.zeros(shape=(df_plus_len, norm_gene_len+(append_amt*2)))
        plus_strand_w_data_matrix = np.zeros(shape=(df_plus_len, norm_gene_len + (append_amt * 2)))
        minus_strand_c_data_matrix = np.zeros(shape=(df_minus_len, norm_gene_len + (append_amt * 2)))
        minus_strand_w_data_matrix = np.zeros(shape=(df_minus_len, norm_gene_len + (append_amt * 2)))

        print(len(df))
        df.index = range(len(df))
        plus_i=0
        minus_i=0
        for i, row in enumerate(df.iterrows()):
            row = row[1]
            gene_len = int(row['length2']/1000)
            strand = row['hg19.knownGene.strand']

            #read in TSS and TTS data:
            w_data_TSS = np.asarray(df.loc[i,'dist_w_100_from_TSS_1':'dist_w_100_from_TSS_201'])
            c_data_TSS = np.asarray(df.loc[i,'dist_c_100_from_TSS_1':'dist_c_100_from_TSS_201'])

            w_data_TTS = np.asarray(df.loc[i,'dist_w_100_from_TTS_1':'dist_w_100_from_TTS_201'])
            c_data_TTS = np.asarray(df.loc[i,'dist_c_100_from_TTS_1':'dist_c_100_from_TTS_201'])

            if (strand == '+'):
                #data over gene
                c_data_arr_gene = c_data_TSS[100:100 + gene_len]  # pull out data from start to end of gene
                c_data_arr_gene_norm = adjust_len(c_data_arr_gene, norm_gene_len)

                w_data_arr_gene = w_data_TSS[100:100 + gene_len]  # pull out data from start to end of gene
                w_data_arr_gene_norm = adjust_len(w_data_arr_gene, norm_gene_len)

                #data upstream of TSS
                c_data_arr_us = c_data_TSS[100 - append_amt:100]
                w_data_arr_us = w_data_TSS[100 - append_amt:100]

                # data downstr of TTS
                c_data_arr_ds = c_data_TTS[101:101 + append_amt]
                w_data_arr_ds = w_data_TTS[101:101 + append_amt]

                c_data_arr_full = np.append(c_data_arr_us, c_data_arr_gene_norm)
                c_data_arr_full = np.append(c_data_arr_full, c_data_arr_ds)

                w_data_arr_full = np.append(w_data_arr_us, w_data_arr_gene_norm)
                w_data_arr_full = np.append(w_data_arr_full, w_data_arr_ds)

                plus_strand_c_data_matrix[plus_i] = c_data_arr_full
                plus_strand_w_data_matrix[plus_i] = w_data_arr_full
                plus_i += 1

            else:
                # data within gene (TSS to TTS) - needs to be normalized
                data_arr_gene = c_data_TTS[100:100 + gene_len]  # pull out data from start to end of gene
                data_arr_gene = np.flip(data_arr_gene, axis=0)
                data_arr_gene_norm = adjust_len(data_arr_gene, norm_gene_len)

                c_data_arr_gene = c_data_TTS[100:100 + gene_len]  # pull out data from start to end of gene
                c_data_arr_gene = np.flip(c_data_arr_gene, axis=0)
                c_data_arr_gene_norm = adjust_len(c_data_arr_gene, norm_gene_len)

                w_data_arr_gene = w_data_TTS[100:100 + gene_len]  # pull out data from start to end of gene
                w_data_arr_gene = np.flip(w_data_arr_gene, axis=0)
                w_data_arr_gene_norm = adjust_len(w_data_arr_gene, norm_gene_len)

                #data upstream of TSS
                c_data_arr_us = c_data_TSS[101:101 + append_amt]
                c_data_arr_us = np.flip(c_data_arr_us, axis=0)

                w_data_arr_us = w_data_TSS[101:101 + append_amt]
                w_data_arr_us = np.flip(w_data_arr_us, axis=0)

                #data downstream of TTS
                c_data_arr_ds = c_data_TTS[100 - append_amt:100]
                c_data_arr_ds = np.flip(c_data_arr_ds, axis=0)

                w_data_arr_ds = w_data_TTS[100 - append_amt:100]
                w_data_arr_ds = np.flip(w_data_arr_ds, axis=0)

                c_data_arr_full = np.append(c_data_arr_us, c_data_arr_gene_norm)
                c_data_arr_full = np.append(c_data_arr_full, c_data_arr_ds)

                w_data_arr_full = np.append(w_data_arr_us, w_data_arr_gene_norm)
                w_data_arr_full = np.append(w_data_arr_full, w_data_arr_ds)

                minus_strand_c_data_matrix[minus_i] = c_data_arr_full
                minus_strand_w_data_matrix[minus_i] = w_data_arr_full
                minus_i += 1

        #plus strand
        c_gene_data = plus_strand_c_data_matrix.mean(axis=0)
        w_gene_data = plus_strand_w_data_matrix.mean(axis=0)
        ave_perc_c_plus = (c_gene_data / (w_gene_data + c_gene_data)) * 100

        # minus strand
        c_gene_data = minus_strand_c_data_matrix.mean(axis=0)
        w_gene_data = minus_strand_w_data_matrix.mean(axis=0)
        ave_perc_w_minus = (w_gene_data / (w_gene_data + c_gene_data)) * 100

        to_plot=np.mean([ave_perc_c_plus, ave_perc_w_minus], axis=0)

        #plot data ----------------------------------------------------
        ax_sb.plot(range(-1*append_amt, norm_gene_len+append_amt), to_plot, '-', color=colors[df_i], linewidth=0.5,)

        #plot gradient of data
        to_plot_grad = smooth(np.gradient(to_plot),2)
        ax_der.plot(range(-1 * append_amt, norm_gene_len + append_amt), to_plot_grad, '-', color=colors[df_i], linewidth=0.5, )

    ### SAVING PLOTS ###
    plot_title = "n=" + str(len(df))

    xlabel = 'Gene (kb), Txn L-->R'
    ylabel = '% of forks moving L-->R'
    filename = 'TSS_to_TTS_Norm_Len'

    ### sb - save plot ###
    if (len(ylim) == 2):
        ax_sb.set_ylim(ylim)
        if (len(yticks) > 0):
            ax_sb.set_yticks(yticks)

    if (len(xlim) == 2):
        ax_sb.set_xlim(xlim)
        if (len(xticks) > 0):
            ax_sb.set_xticks(xticks)
        if (len(xtick_labels) > 0):
            ax_sb.set_xticklabels(xtick_labels)

    ax_sb.set_title(plot_title, x=0.02, y=0.9, ha='left', va='top')
    ax_sb.axhline(y=yline, color='black', linewidth=0.8, dashes=[2, 4])
    ax_sb.axvline(x=0, linestyle='--', color='black', linewidth=0.8, dashes=[2, 4])
    ax_sb.axvline(x=norm_gene_len, linestyle='--', color='black', linewidth=0.8, dashes=[2, 4])

    ax_sb.set_xlabel(xlabel)
    ax_sb.set_ylabel(ylabel)
    fig_sb.tight_layout(pad=2)
    fig_sb.savefig(base_dir + '/' + filename + '_' + fn_ext + '.' + ext, transparent=True)
    fig_sb.clf()

    ### DER - save plot ###
    der_ylabel = '% increase in L-->R forks'
    filename = 'TSS_to_TTS_Norm_Len_grad'

    if (len(der_ylim) == 2):
        ax_der.set_ylim(der_ylim)
        if (len(der_yticks) > 0):
            ax_der.set_yticks(der_yticks)

    ax_der.set_title(plot_title, x=0.02, y=0.9, ha='left', va='top')
    ax_der.axhline(y=der_yline, color='black', linewidth=0.8, dashes=[2, 4])
    ax_der.axvline(x=0,linestyle='--', color='black', linewidth=0.8, dashes=[2, 4])
    ax_der.axvline(x=norm_gene_len, linestyle='--', color='black', linewidth=0.8, dashes=[2, 4])

    ax_der.set_xlabel(xlabel)
    ax_der.set_ylabel(der_ylabel)
    fig_der.tight_layout(pad=2)
    fig_der.savefig(base_dir + '/' + filename + '_'+fn_ext+'.' + ext, transparent=True)
    fig_der.clf()

def calculate_pvalue(data_df1, data_df2, filter_col, quartile_col, stat_range):
    # compares data sets using kruskal wallis test.
    # if data_df2 is empty, quartile_col MUST be given b/c the Q1 will be compared to Q2, Q2 to Q3, and Q3 to Q4

    # if data_df2 is given, then the comparison is between the 2 data sets,
    # and if quartile_col is also given, then the comparison between the 2 data sets will be done for each quartile
    # if quartile_col == '', then data will not be broken up into quartiles

    # if one tuple given for stat_range, ave of this range will be used for comparison
    # if 2 tuples given, difference between ave over range of 2nd and 1st tuple will be used for comparison

    # will filter first if filter_col != ''

    if(data_df2==''):
        df_list = [data_df1,]
        data_to_compare = []
        if(quartile_col == ''):
            print("Only one data set given, quartile_col must be set.")
            return 0
    else:
        if (quartile_col != ''):
            data_to_compare = [[], [], [], []]  # 4 quartiles
        else:
            data_to_compare = []
        df_list = [data_df1, data_df2]

    for cur_df in df_list:
        # first, filter the dfs
        if (filter_col != ''):
            filter_data = cur_df[filter_col].dropna()
            med = np.median(filter_data)
            data_df_filt = cur_df[cur_df[filter_col] > med]
        else:
            data_df_filt = cur_df

        # second break up by quartiles
        if(quartile_col != ''):
            qt_data = data_df_filt[quartile_col].dropna()
            cutoffs = np.percentile(qt_data, [25, 50, 75])
            data_limits = [[cutoffs[2], max(qt_data) + 1],
                           [cutoffs[1], cutoffs[2]],
                           [cutoffs[0], cutoffs[1]],
                           [-1, cutoffs[0]]]

            q1_zero = False
            cutoff_df = data_df_filt[(data_df_filt[quartile_col] >= data_limits[3][0]) & (data_df_filt[quartile_col] < data_limits[3][1])]
            if (len(cutoff_df) == 0):
                q1_zero = True

            # then add to data list
            for i, cur_limits in enumerate(data_limits):
                if (q1_zero):
                    cutoff_df = data_df_filt[data_df_filt[quartile_col] > cur_limits[0]]
                    cutoff_df = cutoff_df[cutoff_df[quartile_col] <= cur_limits[1]]
                else:
                    cutoff_df = data_df_filt[data_df_filt[quartile_col] >= cur_limits[0]]
                    cutoff_df = cutoff_df[cutoff_df[quartile_col] < cur_limits[1]]

                sb_by_gene = get_heatmap_data(cutoff_df, around_mid=100)
                if(data_df2==''):
                    data_to_compare.append(sb_by_gene)
                else:
                    data_to_compare[i].append(sb_by_gene)
        else:
            sb_by_gene = get_heatmap_data(data_df_filt, around_mid=100)
            data_to_compare.append(sb_by_gene)

    if(data_df2==''):
        return_vals = {}
        qlens = {}
        for i, cur_q in enumerate(['Q3-Q4', 'Q2-Q3', 'Q1-Q2']): # order in data_to_compare is ['Q4','Q3','Q2','Q1']
            # calculate stats for each quartile comparison
            if (len(stat_range) == 1):
                ave_over_range1 = data_to_compare[i+1][:, stat_range[0][0] + 100:stat_range[0][1] + 100 + 1].mean(axis=1)
                ave_over_range2 = data_to_compare[i][:, stat_range[0][0] + 100:stat_range[0][1] + 100 + 1].mean(axis=1)

            elif (len(stat_range) == 2):
                ave_over_range1 = data_to_compare[i+1][:, stat_range[1][0] + 100:stat_range[1][1] + 100 + 1].mean(axis=1) - \
                                  data_to_compare[i+1][:, stat_range[0][0] + 100:stat_range[0][1] + 100 + 1].mean(axis=1)
                ave_over_range2 = data_to_compare[i][:, stat_range[1][0] + 100:stat_range[1][1] + 100 + 1].mean(axis=1) - \
                                  data_to_compare[i][:, stat_range[0][0] + 100:stat_range[0][1] + 100 + 1].mean(axis=1)
            qlens = [len(data_to_compare[i+1]), len(data_to_compare[i])]
            s, p = stats.kruskal(ave_over_range1, ave_over_range2)
            return_vals[cur_q] = (p, ave_over_range1.mean() - ave_over_range2.mean(), qlens)
    else:
        # compare the data in the list
        if (quartile_col != ''):
            return_vals = {}
            qlens={}
            for i, cur_q in enumerate(['Q4','Q3','Q2','Q1']):
                #calculate stats for this quartile
                if(len(stat_range) == 1):
                    ave_over_range1 = data_to_compare[i][0][:,stat_range[0][0]+100:stat_range[0][1]+100+1].mean(axis=1)
                    ave_over_range2 = data_to_compare[i][1][:,stat_range[0][0]+100:stat_range[0][1]+100+1].mean(axis=1)

                elif(len(stat_range) == 2):
                    ave_over_range1 = data_to_compare[i][0][:,stat_range[1][0]+100:stat_range[1][1]+100+1].mean(axis=1) - data_to_compare[i][0][:,stat_range[0][0]+100:stat_range[0][1]+100+1].mean(axis=1)
                    ave_over_range2 = data_to_compare[i][1][:,stat_range[1][0]+100:stat_range[1][1]+100+1].mean(axis=1) - data_to_compare[i][1][:,stat_range[0][0]+100:stat_range[0][1]+100+1].mean(axis=1)
                qlens = [len(data_to_compare[i][0]), len(data_to_compare[i][1])]
                s, p = stats.kruskal(ave_over_range1, ave_over_range2)
                return_vals[cur_q] = (p, ave_over_range1.mean() - ave_over_range2.mean(), qlens)
        else:
            if (len(stat_range) == 1):
                ave_over_range1 = data_to_compare[0][:, stat_range[0][0] + 100:stat_range[0][1] + 100 + 1].mean(axis=1)
                ave_over_range2 = data_to_compare[1][:, stat_range[0][0] + 100:stat_range[0][1] + 100 + 1].mean(axis=1)

            elif (len(stat_range) == 2):
                ave_over_range1 = data_to_compare[0][:, stat_range[1][0] + 100:stat_range[1][1] + 100 + 1].mean(axis=1) - data_to_compare[0][:, stat_range[0][0] + 100:stat_range[0][1] + 100 + 1].mean(axis=1)
                ave_over_range2 = data_to_compare[1][:, stat_range[1][0] + 100:stat_range[1][1] + 100 + 1].mean(axis=1) - data_to_compare[1][:, stat_range[0][0] + 100:stat_range[0][1] + 100 + 1].mean(axis=1)
            qlens=[len(data_to_compare[0]), len(data_to_compare[1])]
            s, p = stats.kruskal(ave_over_range1, ave_over_range2)
            return_vals = (p, ave_over_range1.mean() - ave_over_range2.mean(), qlens)

    return return_vals




