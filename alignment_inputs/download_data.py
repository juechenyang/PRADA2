"""
created by Juechen Yang at 2/5/20

"""

import os
import subprocess
from subprocess import Popen
#specify the root dir
root_dir = os.path.dirname(os.path.realpath(__file__))
gtf_file = os.path.join(root_dir, 'gencode.v22.annotation.gtf.gz')
fasta_file = os.path.join(root_dir, 'GRCh38.d1.vd1.fa.tar.gz')

#try to download from GDC
try:
    # execute to download the gtf file
    print "dowloading gtf file"
    Process = Popen(['wget', '-O', gtf_file, '-P', root_dir,
                     'https://api.gdc.cancer.gov/data/25aa497c-e615-4cb7-8751-71f744f9691f'],
                    stdout=subprocess.PIPE)
    print Process.communicate()

    # execute to download the FASTA file
    print "dowloading FASTA file"
    Process = Popen(['wget', '-O', fasta_file, '-P', root_dir,
                     'https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834'],
                    stdout=subprocess.PIPE)
    print Process.communicate()
except:
    print "dowloading failed, please check the connection"

#extraction process

#extract gtf file
extract_process = Popen(['gunzip', gtf_file], stdout=subprocess.PIPE)
print extract_process.communicate()

#extract fasta file
extract_process = Popen(['tar', '-zxvf', fasta_file, '--directory', root_dir], stdout=subprocess.PIPE)
print extract_process.communicate()

#remove the fasta tar
os.remove(fasta_file)
