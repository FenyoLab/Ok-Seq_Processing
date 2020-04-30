# Performs all steps necessary to process OkSeq sequencing files (fastq PE)
# Result is density of reads mapped to W and C strands for every 1KB

import subprocess
import pandas as pd
import sys
import os

if(len(sys.argv) >= 7):
    fastq_r1 = sys.argv[1]
    fastq_r2 = sys.argv[2]
    samfile = sys.argv[3]
    logfile = sys.argv[4]
    final_outfile = sys.argv[5]
    fraglengthlimit = int(sys.argv[6]) #if set to 0 save all frag lengths
else:
    print("Error: bad input.")
    exit(0)

picard_cmd = "/gpfs/share/apps/picard/2.18.11/libs/picard.jar"
genome_file = "human.hg19.genome"
chr_lengths = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067,
  'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
  'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248,
  'chr19': 59128983, 'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566,}
chrom_list = range(1, 22+1, 1)
chrom_list = list( map(str, chrom_list) )
chrom_list.extend(['X','Y'])

def convert_bedGraph_to_txt(infile_fwd, infile_rev, outfile):
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

            for pos in range(st_pos + 1000, end_pos + 1000, 1000):
                output_dicts[df_i][chr_ + '_' + str(pos)] = cov

    f = open(outfile, mode='w')
    f.write("chr\tpos\tw\tc\n")
    for chrom in chrom_list:
        end = int(round(chr_lengths['chr' + chrom], -3))
        if (chr_lengths['chr' + chrom] < end):
            end = end - 1000
        for pos in range(1000, end + 1000, 1000):
            f.write(
                chrom + "\t" + str(pos) + '\t' + str(output_dicts[0]['chr' + chrom + '_' + str(pos)] / 1000.0) + '\t' +
                str(output_dicts[1]['chr' + chrom + '_' + str(pos)] / 1000.0) + '\n')

def convert_bedpe_to_bed(input_file, output_file_fwd, output_file_rev):
    lengthlist_fwd = []
    lengthlist_rev = []
    with open(input_file) as f:
        with open(output_file_fwd, 'w') as fw_fwd:
            with open(output_file_rev, 'w') as fw_rev:
                for line in f:
                    sl = line.strip().split('\t')
                    if (len(sl) == 10):

                        if (sl[8] == '+' and sl[9] == '-'):
                            if int(sl[2]) > int(sl[1]) and int(sl[5]) > int(sl[4]) and int(sl[5]) > int(sl[1]):
                                if (fraglengthlimit > 0):
                                    if (int(sl[5]) - int(sl[1]) <= fraglengthlimit):
                                        lengthlist_fwd.append(int(sl[5]) - int(sl[1]))
                                        fw_fwd.write(sl[0] + '\t' + sl[1] + '\t' + sl[5] + '\t' + sl[6] + '\t' + sl[7] + '\t' + sl[8] + '\n')
                                else:
                                    lengthlist_fwd.append(int(sl[5]) - int(sl[1]))
                                    fw_fwd.write(sl[0] + '\t' + sl[1] + '\t' + sl[5] + '\t' + sl[6] + '\t' + sl[7] + '\t' + sl[8] + '\n')

                        if (sl[8] == '-' and sl[9] == '+'):
                            if int(sl[2]) > int(sl[1]) and int(sl[5]) > int(sl[4]) and int(sl[2]) > int(sl[4]):
                                if (fraglengthlimit > 0):
                                    if (int(sl[2]) - int(sl[4]) <= fraglengthlimit):
                                        lengthlist_rev.append(int(sl[2]) - int(sl[4]))
                                        fw_rev.write(sl[0] + '\t' + sl[4] + '\t' + sl[2] + '\t' + sl[6] + '\t' + sl[7] + '\t' + sl[8] + '\n')
                                else:
                                    lengthlist_rev.append(int(sl[2]) - int(sl[4]))
                                    fw_rev.write(sl[0] + '\t' + sl[4] + '\t' + sl[2] + '\t' + sl[6] + '\t' + sl[7] + '\t' + sl[8] + '\n')
                    else:
                        print(line)

    pd.DataFrame(data=lengthlist_fwd, columns=['length']).to_csv(input_file[:-6] + '_fraglengths_fwd.txt', sep='\t')
    pd.DataFrame(data=lengthlist_rev, columns=['length']).to_csv(input_file[:-6] + '_fraglengths_rev.txt', sep='\t')

def run_cmd(cmd_str, f):
    f.write("cmd = " + cmd_str + '\n')
    f.flush()
    output = subprocess.check_output(cmd_str, shell=True)
    f.write("output = " + output.decode("utf-8") + '\n')
    return output

with open(logfile,'w') as fw:
    fw.write("Input parameters:\n")
    fw.write("fastq R1 = " + fastq_r1 + "\n")
    fw.write("fastq R2 = " + fastq_r2 + "\n")
    fw.write("samfile = "+samfile+'\n')
    fw.write("logfile = " + logfile+'\n')
    fw.write("final output file = " + final_outfile +'\n')
    fw.write("frag len lim = " + str(fraglengthlimit)+'\n')
    fw.write('\n\n')
    fw.flush()

    try:
        step_trim=True
        step_align=True
        step1=True
        step2=True
        step3=True
        step4=True
        step5=True
        step6=True

        output_dir = os.path.split(fastq_r1)[0] + '/'
        trimmed_r1 = fastq_r1.rstrip(".fastq.gz") + "_val_1.fq.gz"  # should fix this.
        trimmed_r2 = fastq_r2.rstrip(".fastq.gz") + "_val_2.fq.gz"  # should fix this.
        bamfile = samfile[:-4] + ".bam"
        sortedbamfile = bamfile[:-4] + ".sorted.bam"
        rmdupfile = bamfile[:-4] + ".rmdup.bam"
        metricsfile = bamfile[:-4] + ".rmdup-metrics.txt"
        namesortedfile = rmdupfile[:-4] + ".namesorted.bam"
        bedpefile = bamfile[:-4] + ".bedpe"
        bedfile_fwd = bamfile[:-4] + ".fwd.bed"
        bedfile_rev = bamfile[:-4] + ".rev.bed"
        bedfile_fwd_sorted = bedfile_fwd[:-4] + ".sorted.bed"
        bedfile_rev_sorted = bedfile_rev[:-4] + ".sorted.bed"
        bamfile_fwd = bedfile_fwd_sorted[:-11] + ".bam"
        bamfile_rev = bedfile_rev_sorted[:-11] + ".bam"
        sortedbamfile_fwd = bamfile_fwd[:-4] + ".sorted.bam"
        sortedbamfile_rev = bamfile_rev[:-4] + ".sorted.bam"
        bedGraph_fwd = bamfile_fwd[:-4] + '.bedGraph'
        bedGraph_rev = bamfile_rev[:-4] + '.bedGraph'

        # Step 'trim'
        if(step_trim):
            fw.write("Fastqc/Trim Step: running fastqc and Trimgalore.\n")
            run_cmd("fastqc " + fastq_r1, fw)
            run_cmd("fastqc " + fastq_r2, fw)
            run_cmd("trim_galore --paired --fastqc --output_dir " + output_dir + " " + fastq_r1 + " " + fastq_r2, fw)
            fw.write('\n')

        # Step 'align'
        if(step_align):
            fw.write("Alignment step: running bowtie2.\n")
            run_cmd("bowtie2 -p 8 -x hg19 -1 " + trimmed_r1 + " -2 " + trimmed_r2 + " > " + samfile, fw)
            fw.write('\n')

        # Step 1, extract only paired reads that have mapq score >= 30, and convert to bam
        if(step1):
            fw.write("Step 1: Extract only paired reads that have mapq score >= 30, and convert sam to bam file.\n")
            run_cmd("samtools view -q 30 -f 0x2 -b " + samfile + " > " + bamfile, fw)
            fw.write('\n')

        # Step 2, Picard remove duplicates:
        if(step2):
            fw.write("Step 2: Picard remove duplicates.\n")
            run_cmd("samtools sort " + bamfile + " > " + sortedbamfile, fw)
            run_cmd("java -jar " + picard_cmd + " MarkDuplicates INPUT=" + sortedbamfile + " OUTPUT=" + rmdupfile +
                    " METRICS_FILE=" + metricsfile + " REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate", fw)
            fw.write('\n')

        # Step 3, use Picard to sort by read name, then convert bam to bedpe using "bedtools bamtobed"
        if(step3):
            fw.write("Step 3: Use Picard to sort by read name and then convert bam to bedpe using 'bedtools bamtobed'.\n")
            run_cmd("java -jar " + picard_cmd + " SortSam I=" + rmdupfile + " O=" + namesortedfile + " SORT_ORDER=queryname", fw)
            run_cmd("bedtools bamtobed -i " + namesortedfile + " -bedpe -mate1 > " + bedpefile, fw)

        # Step 4, loop through the bedpe file and combine R1/R2 pairs to a single fragment;
        # split the fwd and rev strand mapped fragments into separate files
        if(step4):
            fw.write("Step 4: Combine R1/R2 pairs in bedpe file to single fragment and split fragments mapped to the fwd/rev strand into separate bed files.\n")
            fw.flush()
            convert_bedpe_to_bed(bedpefile, bedfile_fwd, bedfile_rev)

        # Step 5, sort and then convert back to bam files (bedtools bedtobam) and run deepTools bamCoverage,
        # this will create the histogram of coverage for every 1KB
        # Note: the genome file must be present in the dir as 'human.hg19.genome',
        # and for bamCoverage, an index file (.bai) must be created (samtools index)
        if(step5):
            fw.write("Step 5: Sort bed files by position and then convert back to bam format, create index files and then run deepTools to get bedgraph files.\n")

            #sort by position
            run_cmd("sort -k 1,1 " + bedfile_fwd + " > " + bedfile_fwd_sorted, fw)
            run_cmd("sort -k 1,1 " + bedfile_rev + " > " + bedfile_rev_sorted, fw)

            #bed to bam
            run_cmd("bedToBam -i " + bedfile_fwd_sorted + " -g " + genome_file + " > " + bamfile_fwd, fw)
            run_cmd("bedToBam -i " + bedfile_rev_sorted + " -g " + genome_file + " > " + bamfile_rev, fw)

            #sort bam file
            run_cmd("samtools sort " + bamfile_fwd + " > " + sortedbamfile_fwd, fw)
            run_cmd("samtools sort " + bamfile_rev + " > " + sortedbamfile_rev, fw)

            #make index
            run_cmd("samtools index " + sortedbamfile_fwd, fw)
            run_cmd("samtools index " + sortedbamfile_rev, fw)

            #convert to bedGraph
            run_cmd("bamCoverage -b " + sortedbamfile_fwd + " -o " + bedGraph_fwd + " --binSize 1000 --outFileFormat bedgraph", fw)
            run_cmd("bamCoverage -b " + sortedbamfile_rev + " -o " + bedGraph_rev + " --binSize 1000 --outFileFormat bedgraph", fw)

        # Step 6, read in result from deepTools bamCoverage: one file for W and one for C
        # convert to a single txt file that places W and C in sep columns, cov for every 1kb
        if(step6):
            fw.write("Step 6: Convert bedGraph files to final txt output file with columns for W and C coverage every 1kb.\n")
            fw.flush()
            convert_bedGraph_to_txt(bedGraph_fwd, bedGraph_rev, final_outfile)

    except subprocess.CalledProcessError as e:
        fw.write("Exception occured (CalledProcessError): \n")
        fw.write("return_code="+str(e.returncode)+'\n')
        fw.write("output="+e.output.decode("utf-8")+'\n')

    fw.write("END")














