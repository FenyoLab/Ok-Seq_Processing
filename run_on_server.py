#

import sys
import os
import glob

if(len(sys.argv) >= 2):
    base_dir=sys.argv[1]
else:
    print("Error no input parameter for start directory.")
    exit(1)

#get all dir names in base_dir
for root, dirs, files in os.walk(base_dir, topdown=False):
    for name in dirs:
        cur_dir = os.path.join(root, name)

        #get the fastq.gz R1 and R2 files
        r1_file = glob.glob(cur_dir + '/*R1_001.fastq.gz')[0]
        r1_file = os.path.split(r1_file)[1]
        r2_file = glob.glob(cur_dir + '/*R2_001.fastq.gz')[0]
        r2_file = os.path.split(r2_file)[1]

        print(f"sbatch run_process_okseq.sh '{cur_dir}' '{base_dir}' '{r1_file}' '{r2_file}' 'aligned.sam' '{name}_log.txt' '{name}_out.txt'")
        print("")
        os.system(f"sbatch run_process_okseq.sh '{cur_dir}' '{base_dir}' '{r1_file}' '{r2_file}' 'aligned.sam' '{name}_log.txt' '{name}_out.txt'")

