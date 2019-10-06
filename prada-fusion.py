import time
import os
from ioprada import getMatching, getCandidates
import pandas as pd
import pybedtools
from FastaTools import parsingFasta, get_max_tx_seq, blastFastas
import shutil
from DirTools import check_dir
import pysam
from GTFtools import export_to_bed, check_gtf_tsv
from JunctionTools import filter_junction
from Filters import filterOrientation, sameGeneFilter, readThroughFilter, overlap_exon_filter
import argparse

# set the option for reading data to avoid warning messages
pd.options.mode.chained_assignment = None
base_dir = os.path.dirname(os.path.abspath(__file__))
parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--case_result_dir', help='please specify the name of this case')
requiredNamed.add_argument('--reference_fasta')
requiredNamed.add_argument('--annotation_gtf')
requiredNamed.add_argument('--reference_dir')
args = parser.parse_args()
fasta_file = args.reference_fasta
annotation_file = args.annotation_gtf
case_result_dir = args.case_result_dir
reference_dir = args.reference_dir

# get reference dir and annotation file

outdata_dir = os.path.join(case_result_dir, 'fusion_outputs')
fasta_repo = os.path.join(reference_dir, 'fasta_repo')
intermediate_file_dir = os.path.join(reference_dir, 'intermediates')
check_dir(intermediate_file_dir)
check_dir(fasta_repo)
check_dir(outdata_dir)
gtf_tsv = os.path.join(reference_dir, 'gtf.tsv')

# define some thresholds
junctionSpanning_th = 1
discordant_th = 2
position_consistency_th = 0.1
similarity_evalue_th = 0.001
transcriptional_allelic_ratio = 0.01

# note if gtf has not been converted to a tsv file, then this gtf parser is needed
#*************************GTF Parser******************************
check_gtf_tsv(gtf_tsv, annotation_file)
#*****************************************************************

# export exon's annotation into a bed file
print 'start loading gtf at %s' % time.ctime()
gtf = pd.read_csv(gtf_tsv, sep='\t')
export_to_bed(gtf, intermediate_file_dir=intermediate_file_dir, lincRNA=False)


print 'start parsing fasta at %s' % time.ctime()
# get all transcripts
tx = gtf[gtf.feature == 'transcript']
tx['length'] = tx.end - tx.start + 1
# exons
exs = gtf[gtf.feature == 'exon']
exs['length'] = exs.end - exs.start + 1

# get genome reference sequence
seqDB = parsingFasta(fasta_file)





# case_repo = ['TCGA-CF-A47T-01A', 'TCGA-AO-A03U-01B-21R-A10J-07', 'TCGA-27-1835-01A-01R-1850-01'
#             , 'TCGA-76-4925-01A-01R-1850-01', 'TCGA-FG-7643-01A-11R-2090-07', 'TCGA-P5-A72U-01A-31R-A32Q-07',
#              'TCGA-92-8065-01A', 'TCGA-CH-5765-01A', 'TCGA-CH-5768-01A', 'TCGA-V1-A9OY-01A']
# #get case info
# # for i in range(0, len(case_repo)):
# case_name = case_repo[i]
case = case_result_dir.split('/')[-1]
print 'case is ' + case
alignment_result_dir = os.path.join(case_result_dir,'bam_results')
# get Alignment bam file
bam_file = os.path.join(alignment_result_dir, 'Aligned.sortedByCoord.out.bam')
#get junction file
junc_file = os.path.join(alignment_result_dir, 'Chimeric.out.junction')
unmapped_file = os.path.join(alignment_result_dir, 'Unmapped.out.mate1')
samfile = pysam.AlignmentFile(bam_file, 'rb')
print 'start loading exon bed at %s' %time.ctime()
exonBed = pybedtools.example_bedtool(os.path.join(base_dir,intermediate_file_dir, 'exon.bed'))
#get junction dataframe
columns_name = ['donor_chr', 'donor_fb', 'donor_strand', 'acceptor_chr', 'acceptor_fb', \
                'acceptor_strand', 'junc_type', 'lj_repeatedLen', 'rj_repeatedLen', 'rd_name', \
                'donor_start', 'cigar_firstSeg', 'acceptor_start', 'cigar_secondSeg']
junctionDF = pd.read_csv(junc_file, names=columns_name, low_memory=False, sep='\t')
junctionDF = filter_junction(junctionDF)
#construct bed format for donor and acceptor
donorTable = junctionDF[['donor_chr', 'donor_start', 'donor_end', 'juncID', 'donor_strand', 'donor_fb']]
acceptorTable = junctionDF[['acceptor_chr', 'acceptor_start', 'acceptor_end', 'juncID', 'acceptor_strand', 'acceptor_fb']]
print 'bed file matching start at %s' %time.ctime()
#produce junction dataframe that matched to exon
donorMatching = getMatching(donorTable,'donor',exonBed, intermediate_file_dir)
acceptorMatching = getMatching(acceptorTable, 'acceptor',exonBed, intermediate_file_dir)
result = junctionDF[['juncID', 'junc_type', 'rd_name']]
result = pd.merge(result, donorMatching, how='left', on=['juncID'])
result = pd.merge(result, acceptorMatching, how='left', on=['juncID'])
# ******** sameGeneFilter ********
result = sameGeneFilter(result)
#switch to 1-base coordinates
result['donorgene_start'] = result['donorgene_start'] + 1
result['acceptorgene_start'] = result['acceptorgene_start'] + 1
# ******** orientation filer ********
orientedResult = filterOrientation(result)
# ******** read-through filter ********
orientedResult = readThroughFilter(orientedResult)
# ******** exon overlap filter ********
orientedResult['overlapped'] = orientedResult.apply(lambda x: overlap_exon_filter(x['donorgene_start'],
                            x['donorgene_end'], x['acceptorgene_start'], x['acceptorgene_end'], x['donorgene_chromosome'], x['acceptorgene_chromosome']), axis=1)
orientedResult = orientedResult[orientedResult['overlapped'] == 'No']
#get candidate fusion pairs
print 'data reshape start at %s' % time.ctime()
fusionIdentifier = ['donor_gene_id', 'acceptor_gene_id', 'donor_gene', 'acceptor_gene', 'donor_chr', 'acceptor_chr', 'donorgene_strand', 'acceptorgene_strand', 'donorgene_type', 'acceptorgene_type']
cc = getCandidates(orientedResult, fusionIdentifier, samfile, junctionSpanning_th, discordant_th)
readsResult = cc
# ******** position consistency filter ********
numOfDiscordant = '#discordant_read_pairs'
discordantReads = 'discordant_reads'
numOfJunctionSpanning = '#junctionSpanning_reads'
junctionSpanningReads = 'junctionSpanning_reads'
numberConsistent = '#Consistent'
position_consistency = 'position_consistency'
readsResult[numberConsistent] = readsResult[numberConsistent].astype(float)
readsResult[position_consistency] = (readsResult[numberConsistent]/readsResult[numOfDiscordant]).round(4)
readsResult = readsResult[(readsResult[position_consistency] >= position_consistency_th)]
# ******** PFI filter ********
total_reads_at_junction = 'total_reads_at_junction'
junctionRatio = 'Percentage Fusion Index(PFI)'
readsResult[junctionRatio] = readsResult[numOfJunctionSpanning] / readsResult[total_reads_at_junction]
readsResult = readsResult[(readsResult[junctionRatio] >= transcriptional_allelic_ratio)]
print len(readsResult)
if len(readsResult)==0:
    print "no fusion found"
    import sys
    sys.exit()
# ******** sequence similarity filter ************
shutil.rmtree(fasta_repo, ignore_errors=True)
check_dir(fasta_repo)
donor_fasta = 'donor_fasta'
acceptor_fasta = 'acceptor_fasta'
readsResult[donor_fasta] = readsResult.apply(lambda x: get_max_tx_seq(x['donor_gene_id'], seqDB, tx, exs, fasta_repo), axis=1)
readsResult[acceptor_fasta] = readsResult.apply(lambda x: get_max_tx_seq(x['acceptor_gene_id'], seqDB, tx, exs, fasta_repo), axis=1)
seq_similarity = 'sequence similarity'
print 'blast similarity start at %s' % time.ctime()
readsResult[seq_similarity] = readsResult.apply(lambda x: blastFastas(x[donor_fasta], x[acceptor_fasta]), axis=1)
readsResult[seq_similarity] = readsResult[seq_similarity].astype(float)
readsResult = readsResult[(readsResult[seq_similarity] >= similarity_evalue_th)]
readsResult = readsResult.drop(columns=[donor_fasta, acceptor_fasta])
discodantSeries = readsResult[discordantReads]
junctionSpanningSeries = readsResult[junctionSpanningReads]
numdisSeries = readsResult[numOfDiscordant]
numjunspanSeries = readsResult[numOfJunctionSpanning]
readsResult = readsResult.drop(columns=[discordantReads, junctionSpanningReads, numOfDiscordant, numOfJunctionSpanning])
readsResult.insert(10, numOfDiscordant, numdisSeries)
readsResult.insert(10, numOfJunctionSpanning, numjunspanSeries)
readsResult[junctionSpanningReads] = junctionSpanningSeries
readsResult[discordantReads] = discodantSeries
check_dir(outdata_dir)
readsResult.to_csv(os.path.join(outdata_dir, case + '.fullResult.csv'), index=None)
print 'finished successfully at %s. Output file is %s' % (time.ctime(), case + '.fullResult.csv')


frame_detection_py = os.path.join(base_dir, 'FrameInterface.py')
os.system('python ' + frame_detection_py + ' --case_result_dir ' + case_result_dir +
                  ' --reference_fasta ' + fasta_file + ' --annotation_gtf ' + annotation_file + ' --reference_dir '
                  + reference_dir)




# dd = scan_csv_to_df(outdata_dir)
# dd.to_csv('new_benchmark.csv', index=None)
# pp = orientedResult[(orientedResult.donor_gene == 'CA13') & (orientedResult.acceptor_gene == 'CA3')]
#pp = orientedResult[(orientedResult.donor_gene == 'BCL7A') & (orientedResult.acceptor_gene == 'KCNIP3')]
# pp=pd.DataFrame()
# pp.to_csv('naa30_tert_afterorientation.csv', index=None)
# all_ntrk_annotation = gtf[gtf.gene_name == 'NTRK3']
#
# unmappedDB = SeqIO.parse(unmapped_file, 'fastq')
# unmapped_reads = list()
# for item in unmappedDB:
#     unmapped_reads.append(item.id)
# un_se = pd.Series(data=unmapped_reads)
#
#
# ttresult = pd.DataFrame()
# orientedResult = orientedResult.drop_duplicates(['rd_name', 'donor_gene_id', 'acceptor_gene_id'])
# for name, group in orientedResult.groupby(['donor_gene_id', 'acceptor_gene_id']):
#     junction_sp = group[group.junc_type != -1]
#     if len(junction_sp) < 11:
#         continue
#     if junction_sp.donorgene_strand.unique()[0] == '-':
#         donorBreak = junction_sp.drop_duplicates(['rd_name']).donor_fb.mode().min() + 1
#     else:
#         donorBreak = junction_sp.drop_duplicates(['rd_name']).donor_fb.mode().max() - 1
#     if junction_sp.acceptorgene_strand.unique()[0] == '+':
#         acceptorBreak = junction_sp.drop_duplicates(['rd_name']).acceptor_fb.mode().min() + 1
#     else:
#         acceptorBreak = junction_sp.drop_duplicates(['rd_name']).acceptor_fb.mode().max() - 1
#     allreads = list()
#     for read in samfile.fetch(junction_sp.donor_chr.iloc[0], donorBreak - 1, donorBreak):
#         allreads.append(read.qname)
#     for read in samfile.fetch(junction_sp.acceptor_chr.iloc[0], acceptorBreak - 1, acceptorBreak):
#         allreads.append(read.qname)
#     allreads = pd.Series(data=allreads, index=None)
#     junctionSpanningGroup = junction_sp[~(junction_sp.rd_name.isin(allreads))]
#     junctionSpanningGroup['in_unmapped'] = junctionSpanningGroup['rd_name'].isin(un_se)
#     ttresult = ttresult.append(junctionSpanningGroup)
#
# all_etv6 = gtf[gtf.gene_name == 'ETV6']
# all_etv6.to_csv('all_etv6_annotation.csv', index=None)



# tt = readsResult[(readsResult.donor_gene == 'TMPRSS2') & (readsResult.acceptor_gene == 'ERG')]
# all_result = scan_csv_to_df(outdata_dir)
# all_result.to_csv('new_benchmark.csv', index=None)

# donor_freq = readsResult.groupby('donor_gene').size().reset_index(name='num')
# donor_freq = donor_freq[donor_freq.num >= 3]
# acceptor_freq = readsResult.groupby('acceptor_gene').size().reset_index(name='num')
# acceptor_freq = acceptor_freq[acceptor_freq.num >= 3]
#
# haha = readsResult[~((readsResult.donor_gene.isin(donor_freq.donor_gene) | readsResult.acceptor_gene.isin(acceptor_freq.acceptor_gene)))]

# def fetchRead(chr ,start, end, samfile):
#     read_name = 'qname'
#     is_read1 = 'is_read1'
#     is_read2 = 'is_read2'
#     is_unmapped = 'is_unmapped'
#     mate_is_unmapped = 'mate_is_unmapped'
#     reference_start = 'pos'
#     reference_end = 'aend'
#     reference_length = 'alen'
#     next_reference_start = 'mpos'
#     mapping_quality = 'mapping_quality'
#     is_reverse = 'is_reverse'
#     is_secondary = 'is_secondary'
#     reference_id = 'tid'
#     next_reference_id = 'rnext'
#     cigarstring = 'cigarstring'
#     test = list()
#     for read in samfile.fetch(chr, start-1, end):
#         innerList = list()
#         innerList.append(read.qname)
#         innerList.append(read.mapping_quality)
#         innerList.append(read.is_read1)
#         innerList.append(read.is_read2)
#         innerList.append(read.is_unmapped)
#         innerList.append(read.mate_is_unmapped)
#         innerList.append(read.is_reverse)
#         innerList.append(read.is_secondary)
#         innerList.append(read.reference_id)
#         innerList.append(read.next_reference_id)
#         innerList.append(read.reference_start)
#         innerList.append(read.reference_end)
#         innerList.append(read.next_reference_start)
#         innerList.append(read.reference_length)
#         innerList.append(read.cigarstring)
#         test.append(innerList)
#     return pd.DataFrame(data=test, columns=[read_name, mapping_quality, is_read1, is_read2, is_unmapped, mate_is_unmapped,
#                         is_reverse, is_secondary, reference_id, next_reference_id, reference_start, reference_end, next_reference_start, reference_length, cigarstring])
#
# a = fetchRead('chr14', 57399883, 57399883,samfile)
# b = fetchRead('chr5',1282624, 1282624,samfile)
# a.to_csv('naa30_fetched_at_break.csv', index=None)
# b.to_csv('tert_fetched_at_ break.csv', index=None)
# all_tert_annotation = gtf[gtf.gene_name == 'TERT']
# all_tert_annotation.to_csv('all_tert_annotation.csv', index=None)
# # ETV6_NTRK3 = ETV6.append(NTRK3, ignore_index=True)
# # ETV6_NTRK3.to_csv('ETV6_NTRK3.csv', index=None)
# #
# # bcl7a = fetchRead('chr12', 122035427, 122035427,samfile)
# # kcnip3 = fetchRead('chr2',95374296, 95374296,samfile)
# # bcl7a_kcnip3 = bcl7a.append(kcnip3, ignore_index=True)
# # pp = readsResult[readsResult.donor_gene == 'BCL7A']
# # overlapped_reads = (pp[ol_jsp_donor].iloc[0]+ '*' + pp[ol_jsp_acceptor].iloc[0]).split('*')
# # readsFeature = bcl7a_kcnip3[bcl7a_kcnip3.read_name.isin(overlapped_reads)]
# # readsFeature.to_csv('BCL7A_KCNIP3.csv', index=None)
# Group = gtf.groupby('gene_id')
# out = list()
# for name, group in Group:
#     inner = list()
#     if len(group) == 0:
#         continue
#     txs = group[group.feature == 'transcript']
#     txs['tag'] = txs.tag.fillna('unknown')
#     inner.append(name)
#     inner.append(len(txs[txs['tag'].str.contains('principal', regex=False)]))
#     out.append(inner)
# yes = pd.DataFrame(data=out, columns=['gene_id', 'count'])

