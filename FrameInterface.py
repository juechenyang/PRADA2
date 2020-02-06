"""
created by Juechen Yang at 1/31/19

"""
import os, argparse
import pandas as pd
pd.options.mode.chained_assignment = None
from FrameDetection import FrameDetection, is_shifting
from TranscriptAnalysis import annotate_transcript
from DirTools import check_dir


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
fusion_dir = os.path.join(case_result_dir, 'fusion_outputs')

gtf_tsv = os.path.join(reference_dir, 'gtf.tsv')
gtf = pd.read_csv(gtf_tsv, sep='\t')
case = case_result_dir.split('/')[-1]
out_columns = ['donor_gene_id', 'acceptor_gene_id', 'Frame', 'donor_transcript_name', 'acceptor_transcript_name',
                   'donor_transcript_priority', 'acceptor_transcript_priority']
# for i in range(0, len(case_repo)):
    #get each case
outfile = os.path.join(fusion_dir, case + '.fullResult.csv')
frame_result = os.path.join(fusion_dir, case + '_frame_prediction.csv')
# get fusion prediction
print 'predicting fusion frame'
fusion_predict = pd.read_csv(outfile, sep=',')
pair_result = pd.DataFrame()
for index, item in fusion_predict.iterrows():
    dd = FrameDetection(item['donor_gene_id'], item['acceptor_gene_id'], item['LeftBreak'],
                   item['RightBreak'], gtf)
    pair_result = pair_result.append(dd)
pair_result['Frame'] = pair_result.apply(lambda x: is_shifting(x['distance_startCodon_to_break_x'],
                                   x['distance_startCodon_to_break_y']), axis=1)
# rename some columns
old_columns = ['gene_id_x', 'gene_name_x', 'transcript_id_x', 'distance_startCodon_to_break_x', 'gene_id_y',
               'gene_name_y', 'transcript_id_y', 'distance_startCodon_to_break_y']
new_columns = ['donor_gene_id', 'donor_gene_name', 'donor_transcript_id', 'donor_distance_startCodon_to_break', 'acceptor_gene_id',
               'acceptor_gene_name', 'acceptor_transcript_id', 'acceptor_distance_startCodon_to_break']
pair_result[new_columns] = pair_result[old_columns]
pair_result = pair_result.drop(columns=old_columns)
# annotate each transcript pair
print 'annotating fusion transcripts'
all_genes = pair_result.drop_duplicates(['donor_gene_id', 'acceptor_gene_id'])[['donor_gene_id', 'acceptor_gene_id']]
donor_txs = pd.DataFrame()
acceptor_txs = pd.DataFrame()
for name, row in all_genes.iterrows():
    donor_txs = donor_txs.append(annotate_transcript(row['donor_gene_id'], gtf, 'donor'))
    acceptor_txs = acceptor_txs.append(annotate_transcript(row['acceptor_gene_id'], gtf, 'acceptor'))
merged = pd.merge(pair_result, donor_txs, how='left', on=['donor_transcript_id'])
merged = pd.merge(merged, acceptor_txs, how='left', on=['acceptor_transcript_id'])
merged.to_csv(frame_result, index=None)

# select most prioritized transcript pair
print 'selecting prioritized transcript'
prioritized_transcripts = pd.DataFrame()
Group = merged.groupby(['donor_gene_id', 'acceptor_gene_id'])
for name, group in Group:
    group = group.sort_values(['donor_transcript_name', 'acceptor_transcript_name',
                               'donor_tag', 'acceptor_tag'])
    dg_id = group.donor_gene_id.iloc[0]
    ag_id = group.acceptor_gene_id.iloc[0]
    Frame = group.Frame.iloc[0]
    donor_principal = group['donor_transcript_name'].iloc[0]
    acceptor_principal = group['acceptor_transcript_name'].iloc[0]
    group = group[group['Frame'].isin(['in-frame', 'out-frame'])]
    if len(group) == 0:
        prioritized_transcripts = prioritized_transcripts.append(
        pd.DataFrame(data=[[dg_id, ag_id, Frame, 'N/A', 'N/A', 'N/A', 'N/A']], columns=out_columns))
    else:
        donor_tx_name = group['donor_transcript_name'].iloc[0]
        acceptor_tx_name = group['acceptor_transcript_name'].iloc[0]
        Frame = group.Frame.iloc[0]
        if group.donor_transcript_name.iloc[0] != donor_principal and group.acceptor_transcript_name.iloc[0] == acceptor_principal:
            prioritized_transcripts = prioritized_transcripts.append(pd.DataFrame(data=[[dg_id, ag_id, Frame, donor_tx_name, acceptor_tx_name, 'alternative',
             'principal']], columns=out_columns))
        elif group.donor_transcript_name.iloc[0] == donor_principal and group.acceptor_transcript_name.iloc[0] != acceptor_principal:
            prioritized_transcripts = prioritized_transcripts.append(pd.DataFrame(data=[[dg_id, ag_id, Frame, donor_tx_name, acceptor_tx_name, 'principal',
             'alternative']], columns=out_columns))
        elif group.donor_transcript_name.iloc[0] != donor_principal and group.acceptor_transcript_name.iloc[0] != acceptor_principal:
            prioritized_transcripts = prioritized_transcripts.append(pd.DataFrame(data=[[dg_id, ag_id, Frame, donor_tx_name, acceptor_tx_name, 'alternative',
             'alternative']], columns=out_columns))
        else:
            prioritized_transcripts = prioritized_transcripts.append(pd.DataFrame(data=[[dg_id, ag_id, Frame, donor_tx_name, acceptor_tx_name, 'principal',
             'principal']], columns=out_columns))
fusion_gene_transcript = pd.merge(fusion_predict, prioritized_transcripts, how='inner', on=['donor_gene_id', 'acceptor_gene_id'])
fusion_gene_transcript.to_csv(os.path.join(fusion_dir, case+'.fusion_gene_transcript.csv'), index=None)

print 'Success! fusion analysis done'