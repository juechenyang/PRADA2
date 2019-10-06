"""
created by Juechen Yang at 2019-05-17

"""
import os, sys
target_dir = '/home/CBBI/yangj8/TARGET'
srx_sub = [x[0] for x in os.walk(target_dir) if 'log' not in x[0]][1:]
fastq_dirs = list()
for x in srx_sub:
    for file in os.listdir(x):
        if file.endswith("fastq"):
            fastq_dirs.append(x)
            break
import pandas as pd
fastq_dirs = pd.Series(data=fastq_dirs)
fastq_dirs = list(fastq_dirs.sort_values())

for x in fastq_dirs[155:]:
    id = x.split('/CBBI/yangj8/TARGET/')[-1].split('/')[0]
    sys.stdout.write('processing ' + id + '++++++++++++++++++++++++++++++++\n')
    fastq1 = ''
    fastq2 = ''
    for file in os.listdir(x):
        if file.endswith('1.fastq'):
            fastq1 = os.path.join(x, file)
        elif file.endswith('2.fastq'):
            fastq2 = os.path.join(x, file)
    command = 'python prada2.py --read1 ' + fastq1 + \
              ' --read2 '+ fastq2 + ' --runThread 12' + \
              ' --outdir /home/CBBI/yangj8/prada_results/' + id \
              + ' --fusion --rsem'
    os.system(command)
