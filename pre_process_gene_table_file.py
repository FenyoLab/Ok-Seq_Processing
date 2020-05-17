# for gene table file, filters for protein coding genes,
# selects gene with the canonical splice site (if avail), removes duplicates
# option to add fpkm data if given

import pandas as pd

add_fpkm_data=True
duplicate_filter = 'first' #'none' # 'canonical'
allowed_chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14',
                'chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
data_dir = "/Users/sarahkeegan/okseq_data"
sites_file = "hgTables.txt"
dupl_subset=['hg19.knownGene.chrom','hg19.kgXref.geneSymbol']
#['hg19.kgXref.geneSymbol'] #['hg19.knownGene.chrom','hg19.kgXref.geneSymbol']

sites_df = pd.read_table(data_dir + '/' + sites_file)
sites_df.fillna({'hg19.kgXref.refseq':"", 'hg19.knownCanonical.chromStart':"none",'hg19.kgXref.geneSymbol':""},inplace=True)

print("Loaded sites file.  Length = ", len(sites_df))

#only looking at chr1-22
sites_df = sites_df[sites_df['hg19.knownGene.chrom'].isin(allowed_chrs)].copy()
print("Filtered for chr1-22.  Length = ", len(sites_df))

#filter only protein coding, remove any that have refseq != NM_...
sites_df = sites_df[sites_df['hg19.kgXref.refseq'].str.startswith('NM_')].copy()
print("Removed all genes that are not protein coding.  Length = ", len(sites_df))

sites_df_unique = sites_df.drop_duplicates(subset=dupl_subset, keep='first')
print("Removed duplicates using 'drop_duplicates'.  Length = ", len(sites_df_unique))

sites_df_canonical = sites_df[sites_df['hg19.knownCanonical.chromStart'] != 'none'].copy()
print("Removed non-canonical genes.  Length = ", len(sites_df_canonical))

if(duplicate_filter=='canonical'):
    dupl_list = sites_df_canonical.duplicated(subset=dupl_subset, keep=False)
    duplicates = sites_df_canonical[dupl_list]
    duplicates.to_csv(data_dir + '/' + "duplicates_after_rem_non-canonical.txt", sep='\t')
    sites_df_canonical_unique = sites_df_canonical.drop_duplicates(subset=dupl_subset, keep='first')
    print("Removed duplicates from canonical using 'drop_duplicates'.  Length = ", len(sites_df_canonical_unique))
    df_to_save = sites_df_canonical_unique
elif(duplicate_filter == 'first'):
    df_to_save = sites_df_unique
else: #duplicate_filter == 'none'
    df_to_save = sites_df

df_to_save = df_to_save.copy()
df_to_save.index = range(len(df_to_save))
suffix=''

if(add_fpkm_data):
    print("Adding fpkm data...")

    # open fpkm table for RPE
    fpkm_data = {}
    fpkm_data['RPE1'] = pd.read_table(data_dir + '/rna_seq_RPE/GSE89413_2016-10-30-NBL-cell-line-STAR-fpkm.txt')
    fpkm_data['RPE1'].fillna({'RPE1': 0,}, inplace=True)
    fpkm_data['RPE1'] = fpkm_data['RPE1'].filter(items=['GeneID', 'RPE1'])

    #####################
    # get gene alias list
    gene_alias_file = pd.read_table(data_dir + '/gene_name_aliases/protein-coding_gene.txt')
    #HUGO_gene_with_protein_product.txt', low_memory=False)  #protein-coding_gene.txt')

    main_row = 'Approved Symbol' #'symbol'  # 'Approved Symbol'
    row1 = 'Synonyms' #'alias_symbol'  # 'Synonyms'
    row2 = 'Previous Symbols' #'prev_symbol'  # 'Previous Symbols'
    split_by = ',' #''|'  # ','
    gene_alias_file.fillna({row1: "", row2: ""}, inplace=True)

    # for each gene in fpkm data table, get gene name and all aliases
    # and set each of them to the fpkm value (using dict)

    cell_line = 'RPE1'

    cur_fpkm_data = fpkm_data[cell_line]
    fpkm_dict = {}
    no_replace = {}
    for i, row in enumerate(cur_fpkm_data.iterrows()):
        row = row[1]
        fpkm_dict[row['GeneID'].strip(' ').upper()] = row[cell_line]
        no_replace[row['GeneID'].strip(' ').upper()]=1

    for i, row in enumerate(gene_alias_file.iterrows()):
        row=row[1]
        gene_name_list=[]
        gene_name_list.append(row[main_row].strip(' ').upper())

        for alias_row_name in [row1,row2]:
            alias_list = row[alias_row_name]
            if(alias_list != ""):
                alias_list = alias_list.split(split_by)
                for alias in alias_list:
                    if(alias != ""):
                        gene_name_list.append(alias.strip(' ').upper())
        alias_with_fpkm=''
        the_fpkm=-1
        for alias in gene_name_list:
            if(alias in fpkm_dict):
                the_fpkm = fpkm_dict[alias]
                alias_with_fpkm = alias
                break

        if(the_fpkm != -1):
            for alias in gene_name_list:
                if(not alias in fpkm_dict):
                    fpkm_dict[alias] = the_fpkm
                else:
                    if(alias != alias_with_fpkm):
                        pass

    #for each gene in sites table, set the fpkm
    genes_missing_fpkm_data={}
    for i, row in enumerate(df_to_save.iterrows()):
        row=row[1]
        sym = row['hg19.kgXref.geneSymbol'].strip(' ').upper()
        if(sym in fpkm_dict):
            if (cell_line == 'FT194'):
                gene_len = row['hg19.knownGene.txEnd']-row['hg19.knownGene.txStart']
                df_to_save.at[i, cell_line + '_fpkm'] = fpkm_dict[sym]/gene_len
            else:
                df_to_save.at[i,cell_line+'_fpkm']=fpkm_dict[sym]
        else:
            genes_missing_fpkm_data[sym]=1

    print("Number of genes from sites table with no fpkm data found (",cell_line,"): ", len(genes_missing_fpkm_data))

    #drop genes with no fpkm data...? NO!
    #df_to_save.dropna(subset=['RPE_fpkm'], inplace=True)
    #print("Removed any genes with no fpkm data.  Length = ", len(df_to_save))
    suffix += 'with_fpkm_'

df_to_save.to_csv(data_dir + '/' + "hgTables_filtered_"+suffix+duplicate_filter+".txt", sep='\t')






