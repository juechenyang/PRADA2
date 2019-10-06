"""
created by Juechen Yang at 1/8/19

"""
import os
import pandas as pd

def export_to_bed(gtf,intermediate_file_dir, lincRNA):
    if lincRNA:
        lincRNAIDs = pd.read_csv(os.path.join(intermediate_file_dir, 'intersect_total.txt'), names=['ids'], sep='\t')
        exons = gtf[(gtf.feature == 'exon') & (gtf.seqname != 'chrM') & (gtf.gene_type ==
                                                                          'protein_coding') | (gtf.gene_id.isin(lincRNAIDs['ids']))][
            ['seqname', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type']]
    else:
        exons = gtf[(gtf.feature == 'exon') & (gtf.seqname != 'chrM') & (gtf.gene_type ==
                                                                         'protein_coding')][
            ['seqname', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type']]
    exons.start = exons.start - 1
    exons.to_csv(os.path.join(intermediate_file_dir, 'exon.bed'), index=None, header=False, sep='\t')
def check_gtf_tsv(gtf_tsv, annotation_file):
    if not os.path.isfile(gtf_tsv):
        from gtfparse import read_gtf
        parsed_gtf = read_gtf(annotation_file)
        parsed_gtf.to_csv(gtf_tsv, sep='\t', index=None)

