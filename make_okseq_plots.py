# make strand bias plots from data tables
import plot_functions as pf
import numpy as np
import sys

if(len(sys.argv) >= 3):
    data_dir = sys.argv[1]
    genomic_loci = sys.argv[2]
    file_name=sys.argv[3]
else:
    print("Missing command line input.  Attempting to run with default settings.")
    data_dir = "/Users/sarahkeegan/okseq_data/strand_bias_plots"
    genomic_loci = "TSS" #change to 'TTS' to make same plots around TTS sites
    file_name = "sites_table_with_distributions-rpe_edu_2.csv"

make_heatmaps=False
make_panels_plots=False

pf.base_dir = data_dir
pf.site=genomic_loci
pf.rep_file_name = file_name

# Load data from table, add columns for gene length and transcriptional volume
df_gene_table = pf.load_data()
df_gene_table['length'] = np.abs(df_gene_table['hg19.knownGene.txEnd'] - df_gene_table['hg19.knownGene.txStart'])
df_gene_table['volume'] = df_gene_table['fpkm'] * df_gene_table['length']

# PMID: 30598550, Fig 1d
# Plotting C/(C+W) read counts (strand bias) around gene TSS/TTS sites, data averaged over all genes
w_data = df_gene_table.loc[:,'dist_w_100_from_' + pf.site + '_1':'dist_w_100_from_' + pf.site + '_' + str(pf.length_around_site + 1)].mean()
c_data = df_gene_table.loc[:,'dist_c_100_from_' + pf.site + '_1':'dist_c_100_from_' + pf.site + '_' + str(pf.length_around_site + 1)].mean()
w_data = np.asarray(w_data)
c_data = np.asarray(c_data)
ave_perc_c = (c_data / (w_data + c_data))*100
pf.make_plot([ave_perc_c,],colors_to_plot=['black'],plot_title='',filename='Fig1d')

# PMID: 30598550, Fig 1e
# Plotting strand bias for Watson and Crick strand genes, data averaged over all W and C genes separately
plus_strand_df = df_gene_table[df_gene_table['hg19.knownGene.strand'] == '+']
w_data = plus_strand_df.loc[:,'dist_w_100_from_' + pf.site + '_1':'dist_w_100_from_' + pf.site + '_' + str(pf.length_around_site + 1)].mean()
c_data = plus_strand_df.loc[:,'dist_c_100_from_' + pf.site + '_1':'dist_c_100_from_' + pf.site + '_' + str(pf.length_around_site + 1)].mean()
w_data = np.asarray(w_data)
c_data = np.asarray(c_data)
ave_perc_c_plus = (c_data / (w_data + c_data))*100

minus_strand_df = df_gene_table[df_gene_table['hg19.knownGene.strand'] == '-']
w_data = minus_strand_df.loc[:,'dist_w_100_from_' + pf.site + '_1':'dist_w_100_from_' + pf.site + '_' + str(pf.length_around_site + 1)].mean()
c_data = minus_strand_df.loc[:,'dist_c_100_from_' + pf.site + '_1':'dist_c_100_from_' + pf.site + '_' + str(pf.length_around_site + 1)].mean()
w_data = np.asarray(w_data)
c_data = np.asarray(c_data)
ave_perc_c_minus = (c_data / (w_data + c_data))*100
pf.make_plot([ave_perc_c_plus,ave_perc_c_minus],colors_to_plot=['blue','magenta'],plot_title='',filename='Fig1e')

# PMID: 30598550, Fig 1f
# Averaging both Watson and Crick strand genes together
# For Crick strand genes, we calculate s.b. in direction of transcription by flipping the order of the data and,
# since we are now reading in the opposite direction, calculating W/(W+C)
ave_perc_w_minus = (w_data / (w_data + c_data))*100
ave_perc_w_minus_flip = np.flip(ave_perc_w_minus, axis=0)
to_plot=np.mean([ave_perc_c_plus, ave_perc_w_minus_flip], axis=0)
pf.make_plot([to_plot,],colors_to_plot=['black'],plot_title='',filename='Fig1f')

if(make_panels_plots):
    # PMID: 30598550, Fig 2a, e, g
    # 2a Plot of strand bias of genes separated into quartiles based on RNA-Seq FPKM values
    # By setting plot_grad=True, the function will also plot the derivative plots of the strand bias
    pf.make_panels_plot_with_data([df_gene_table,], filter_col='', quartile_col='fpkm', filename='Fig2a',
                               line_colors=['black',], plot_grad=True, grad_ylim=(-1,3.25),
                                  grad_filename='Fig2a_derivative')

    #2e Plot of strand bias of genes with RNA-Seq FPKM values > median, separated into quartiles based on gene length
    pf.make_panels_plot_with_data([df_gene_table,], filter_col='fpkm', quartile_col='length', filename='Fig2e',
                                  line_colors=['black',], ylim=(29,81),yticks=[30, 40, 50, 60, 70,80],)

    #2g Plot of strand bias of genes separated into quartiles based on transcriptional volume (FPKM * gene length)
    pf.make_panels_plot_with_data([df_gene_table,], filter_col='', quartile_col='volume', filename='Fig2g',
                                line_colors=['black',], ylim=(29,81),yticks=[30, 40, 50, 60, 70, 80],
                                plot_grad=True, grad_ylim=(-1,3.25), grad_filename='Fig2g_derivative')

if(make_heatmaps):
# PDF files were too large to upload to GitHub, jpg files are present in the repository
# this script will create the PDF files
    # PMID: 30598550, Fig 2b, f, h
    # 2b Heatmap showing same data as in 2a, above - each row is the s.b. for a gene.  Rows sorted by FPKM.
    # Note: for better visualization & smaller size, we converted the pdf file produced by this function to a jpg file
    # using Apple's Previw App, i.e. if using a Mac, right click on the file, Open with "Preview.app" and save as JPG at Best quality
    pf.make_heatmap_plot_with_data([df_gene_table,], filter_col='', sort_col='fpkm', color='YlGnBu_r', filename='Fig2b')

    # 2f Heatmap showing same data as in 2e, above - each row is the s.b. for a gene.  Genes filtered by FPKM>median.
    # Rows sorted by gene length.
    pf.make_heatmap_plot_with_data([df_gene_table,], filter_col='fpkm', sort_col='length', color='YlGnBu_r', filename='Fig2f', around_mid=100)

    # 2h Heatmap showing same data as in 2g, above - each row is the s.b. for a gene.  Rows sorted by volume (FPKM*length)
    pf.make_heatmap_plot_with_data([df_gene_table,], filter_col='', sort_col='volume', color='YlGnBu_r', filename='Fig2g')

# NORM SB plot, and its derivative plot
pf.make_TSS_to_TTS_plots([df_gene_table,],)

# p-values as in PMID: 30598550, Fig 2g:
ret_vals = pf.calculate_pvalue(data_df1=df_gene_table, quartile_col='volume', stat_range=[[-50,-30],[1,10]])

# ret_vals is a dictionary.  in this case, the keys are the quartiles that are compared,
# and the value is a 3-tuple: (p-value, change in efficiency, lengths of the 2 data sets compared):
# {'Q3-Q4': (1.7249844446624105e-91, -9.47402780052795, [4541, 4542]),
#  'Q2-Q3': (2.3932082214793766e-57, -6.462627166090093, [4541, 4541]),
#  'Q1-Q2': (1.9599109843378173e-61, -5.95357630785935, [4541, 4541])}
# (compare to Fig. 2g, delta-efficiency and p-values)
print(ret_vals)

# *NOTE*: to compare 2 data sets, e.g. HU-treated vs. control, as in PMID: 30598550, Fig. 4a, the function can be called as below.
# The additional data set (df_gene_table2) will need to be loaded in.
#ret_vals = pf.calculate_pvalue(data_df1=df_gene_table, quartile_col='volume', stat_range=[[-50,-30],[1,10]], data_df2=df_gene_table2)





