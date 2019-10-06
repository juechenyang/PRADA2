"""
created by Juechen Yang at 2/22/19

"""

import pandas as pd
import os
from subprocess import Popen
from DirTools import check_dir
import argparse
base_dir = os.path.dirname(os.path.abspath(__file__))
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
bash_file = os.path.join(base_dir, 'rsem.sh')
outdir = os.path.join(case_result_dir, 'rsem_results', 'rsem')
rsem_script = os.path.join(base_dir, 'rsem.sh')
rsem_build_script = os.path.join(base_dir, 'rsem_build_ref.sh')
check_dir(outdir)
Process = Popen(['bash ' + rsem_script + ' %s %s %s' % (bam, rsem_reference, outdir)], shell=True)
Process.communicate()
if Process.returncode != 0:
    print 'rsem index error, trying to rebulid'
    Process = Popen(['bash ' + rsem_build_script + ' %s %s %s' % (annotation_gtf, reference_fasta, rsem_reference)], shell=True)
    Process.communicate()
    print 'index build completed redo rsem'
    Process = Popen(['bash ' + rsem_script + ' %s %s %s' % (bam, rsem_reference, outdir)], shell=True)
    Process.communicate()





#calculate Pearson
# GDC_GeneExpDir = os.path.join('..', 'gdc_gene_exp', case)
# for file in os.listdir(GDC_GeneExpDir):
#     if 'FPKM.txt' in file:
#         gdc_gene_exp = pd.read_csv(os.path.join(GDC_GeneExpDir, file), sep='\t',
#                                    names=['gene_id', 'exp'])
# rsem_gene_exp = pd.read_csv(os.path.join(outdir, '..', 'rsem.genes.results'),
#                             sep='\t')
# rsem_gene_exp[['gene_id', 'gene_name']] = rsem_gene_exp['gene_id'].str.split('_', expand=True)[[0,1]]
# rsem_gene_exp = rsem_gene_exp[['gene_id', 'gene_name', 'FPKM']]
# merged = pd.merge(gdc_gene_exp, rsem_gene_exp, how='outer', on=['gene_id'])
# corr_matrix = merged[['exp', 'FPKM']].corr()
# ten_cases_corr.append([case, corr_matrix.iloc[0,1]])
# ten_cases_corr = pd.DataFrame(data=ten_cases_corr, columns=['case', 'pearson_corr'])
# ten_cases_corr.to_csv('gene_exp_corr_FPKM-UQ.csv', index=None)