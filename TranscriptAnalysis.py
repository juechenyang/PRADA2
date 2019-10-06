"""
created by Juechen Yang at 2/12/19

"""

def annotate_transcript(geneID, gtf, side):
    transcripts = gtf[(gtf.gene_id == geneID) & (gtf.feature == 'transcript')][['transcript_id', 'transcript_type', 'transcript_status', 'transcript_name', 'tag',
       'transcript_support_level']]
    transcripts[[side + '_transcript_id', side + '_transcript_type', side + '_transcript_status', side + '_transcript_name', side + '_tag',
       side + '_transcript_support_level']] = transcripts[['transcript_id', 'transcript_type', 'transcript_status', 'transcript_name', 'tag',
       'transcript_support_level']]
    transcripts = transcripts.drop(columns=['transcript_id', 'transcript_type', 'transcript_status', 'transcript_name', 'tag',
       'transcript_support_level'])
    return transcripts