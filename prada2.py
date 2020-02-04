"""
created by Juechen Yang at 2019-03-29

"""
#import required modules
import argparse
from subprocess import Popen
from DirTools import check_dir
import os
import shutil
import time

#specify the base dir
base_dir = os.path.dirname(os.path.realpath(__file__))

#bash scripts dir
bash_dir = os.path.join(base_dir, "bash_scripts")

#get the path of reference fasta and annotation gtf
reference_fasta = os.path.join(base_dir, 'alignment_inputs', 'GRCh38.d1.vd1.fa')
annotation_gtf = os.path.join(base_dir, 'alignment_inputs', 'gencode.v22.annotation.gtf')



parser = argparse.ArgumentParser(description='''
************************************************************************************************
Basic description:

Integrated Pipeline for RNA Sequencing Data Alignment
**Command**:
python prada-preprocess-v2018.py --read1 XXX.fastq --read2 XXX.fastq

************************************************************************************************
    ''',formatter_class=argparse.RawTextHelpFormatter
    )
parser.add_argument('--mode', default='concise',help='identify what mode user want to use. Default is "concise"')
parser.add_argument('--qc', default=False, action='store_true', help='Whether to enable RNA-SeQC module in the pipeline, default is false')
parser.add_argument('--markdup', default=False, action='store_true', help='Whether to enable RNA-SeQC module in the pipeline, default is false')
parser.add_argument('--kpinterf', default=False, action='store_true', help='Whether to keep intermediate files, default is false')
parser.add_argument('--outdir', default=os.getcwd(), help='specify which directory to store the output files, default is current working directory')
parser.add_argument('--runThread', default=16, help='decide the number of thread used by STAR, default is 16')
parser.add_argument('--fusion', default=False, action='store_true', help='Whether to perform fusion analysis')
parser.add_argument('--rsem', default=False, action='store_true', help='Whether to perform rsem gene expression analysis')
parser.add_argument('--ref_fa', default=reference_fasta, help='specify which reference genome is going to be used, default is GRCh38.d1.vd1.fa in prada home directory')
parser.add_argument('--ano_gtf', default=annotation_gtf, help='specify which annotation gtf is going to be used, default is gencode.v22.annotation.gtf in prada home directory')
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--read1', help='required option, which is used to specify fastq read1. \nPlease use comma to seprate each fastq file if your sample contains multiple files')
requiredNamed.add_argument('--read2', help='required option, which is used to specify fastq read2. \nPlease use comma to seprate each fastq file if your sample contains multiple files')
args = parser.parse_args()

#collect all parameters to variables
mode = args.mode
kp = args.kpinterf
qc = args.qc
outdir = args.outdir
read1 = args.read1
read2 = args.read2
thread = args.runThread
markdup = args.markdup
fusion = args.fusion
rsem = args.rsem
reference_fasta = args.ref_fa
annotation_gtf = args.ano_gtf


#check if the input fastq files are pair-end reads
if not read1 or not read2:
    print 'missing required read '
    import sys
    sys.exit()
read1_list = read1.split(',')
read2_list = read2.split(',')

#finalize output dir
case_dir = os.path.join(outdir,read1_list[0].split('/')[-1][:-6])

#specify alignment reference file dir and genome index dir
reference_file_dir = os.path.join(base_dir, 'reference_files')
genome_index_dir = os.path.join(reference_file_dir, 'genomeIndices')

#specify critical output for alignment
out_bam = os.path.join(case_dir, 'bam_results', 'Aligned.sortedByCoord.out.bam')
out_bai = os.path.join(case_dir, 'bam_results', 'Aligned.sortedByCoord.out.bam.bai')

#if bam and bai files both exists, then we consider the alignment has been done
if os.path.isfile(out_bam) and os.path.isfile(out_bai):
    print 'from python: skipped the alignment process'

    #fusion script to get all fusions
    if fusion:
        print 'from python: Start doing Fusion Detection at %s' % time.ctime()
        fusion_py = os.path.join(base_dir, 'prada-fusion.py')
        os.system('python ' + fusion_py + ' --case_result_dir ' + case_dir +
                  ' --reference_fasta ' + reference_fasta + ' --annotation_gtf ' + annotation_gtf + ' --reference_dir '
                  + reference_file_dir)
    #rsem scripts to get all gene expression
    if rsem:
        print 'from python: Start doing gene expression calculation at %s' % time.ctime()
        rsem_py = os.path.join(base_dir, 'GeneExpAnalysis.py')
        os.system('python ' + rsem_py + ' --case_result_dir ' + case_dir +
                  ' --reference_fasta ' + reference_fasta + ' --annotation_gtf ' + annotation_gtf + ' --reference_dir '
                  + reference_file_dir)
else:
    # check if genome index is ready
    if not os.path.exists(genome_index_dir):
        print 'genome index has not been built, trying to rebuild'
        check_dir(genome_index_dir)
        Process = Popen(['bash ', bash_dir, '/check_index.sh %s %s %s %s' % (
        genome_index_dir, reference_fasta, annotation_gtf, 8)],
                        shell=True)
        print Process.communicate()

    # create directories to store output results
    check_dir(os.path.join(case_dir, 'bam_results'))
    check_dir(os.path.join(case_dir, 'inter_results'))
    if qc:
        check_dir(os.path.join(case_dir, 'QC'))

    for read in read1_list+read2_list:
        #convert all to a lower case format for verification
        read_lower_form = read.lower()

        #if any of the read is not a fastq format, the program will exit
        if not (read_lower_form.endswith('.fastq') or read_lower_form.endswith('.fq')):
            print 'your input is not a fastq file'
            import sys
            sys.exit()
    if len(read1_list) != len(read2_list):
        print 'your input read1 and read2 is not equal in pairs, please check'
        import sys
        sys.exit()


    print 'from python: Start doing Alignment at %s' % time.ctime()
    read1 = ','.join(read1_list)
    read2 = ','.join(read2_list)
    Process = Popen(['bash ' + bash_dir+ '/concise_mode_v1.0.sh %s %s %s %s %s %s %s %s %s %s' % (read1, read2,
                   case_dir, qc, kp, thread, rsem, read1_list[0].split('/')[-1][:-6], markdup, genome_index_dir)], shell=True)
    res = Process.communicate()

    if Process.returncode != 0:
        print 'error in genome index, trying to rebuild index'
        shutil.rmtree(genome_index_dir, ignore_errors=True)
        check_dir(genome_index_dir)
        Process = Popen(['bash ' + bash_dir + '/check_index.sh %s %s %s %s' % (genome_index_dir, reference_fasta, annotation_gtf, 8)],
                        shell=True)
        print Process.communicate()

        print 'redo alignment'
        Process = Popen(['bash ' + bash_dir + '/concise_mode_v1.0.sh %s %s %s %s %s %s %s %s %s %s' % (
        read1, read2, case_dir, qc, kp, thread, rsem, read1_list[0].split('/')[-1][:-6], markdup, genome_index_dir)],
                        shell=True)
        print Process.communicate()


    if fusion:
        print 'from python: Start doing Fusion Detection at %s' % time.ctime()
        fusion_py = os.path.join(base_dir, 'prada-fusion.py')
        os.system('python ' + fusion_py + ' --case_result_dir ' + case_dir +
                  ' --reference_fasta ' + reference_fasta + ' --annotation_gtf ' + annotation_gtf + ' --reference_dir '
                  + reference_file_dir)
    if rsem:
        print 'from python: Start doing gene expression calculation at %s' % time.ctime()
        rsem_py = os.path.join(base_dir, 'GeneExpAnalysis.py')
        os.system('python ' + rsem_py + ' --case_result_dir ' + case_dir +
                  ' --reference_fasta ' + reference_fasta + ' --annotation_gtf ' + annotation_gtf + ' --reference_dir '
                  + reference_file_dir)

