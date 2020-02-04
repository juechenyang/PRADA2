"""
created by Juechen Yang at 2/22/19

"""

import pandas as pd
import os
from subprocess import Popen
from DirTools import check_dir
import argparse
#specify the base dir
base_dir = os.path.dirname(os.path.realpath(__file__))

#bash scripts dir
bash_dir = os.path.join(base_dir, "bash_scripts")


parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--case_result_dir')
requiredNamed.add_argument('--reference_fasta')
requiredNamed.add_argument('--annotation_gtf')
requiredNamed.add_argument('--reference_dir')
args = parser.parse_args()
case_result_dir = args.case_result_dir
reference_fasta = args.reference_fasta
annotation_gtf = args.annotation_gtf
reference_dir = args.reference_dir

rsem_reference = os.path.join(reference_dir, 'rsem_ref', 'ref')
check_dir(rsem_reference)
bam = os.path.join(case_result_dir, 'bam_results', 'Aligned.toTranscriptome.out.bam')

#specify output dir
outdir = os.path.join(case_result_dir, 'rsem_results', 'rsem')
check_dir(outdir)

#create path for scripts
rsem_script = os.path.join(bash_dir, 'rsem.sh')
rsem_build_script = os.path.join(bash_dir, 'rsem_build_ref.sh')

#run rsem
Process = Popen(['bash ' + rsem_script + ' %s %s %s' % (bam, rsem_reference, outdir)], shell=True)
Process.communicate()
if Process.returncode != 0:
    print 'rsem index error, trying to rebulid'
    Process = Popen(['bash ' + rsem_build_script + ' %s %s %s' % (annotation_gtf, reference_fasta, rsem_reference)], shell=True)
    Process.communicate()
    print 'index build completed redo rsem'
    Process = Popen(['bash ' + rsem_script + ' %s %s %s' % (bam, rsem_reference, outdir)], shell=True)
    Process.communicate()