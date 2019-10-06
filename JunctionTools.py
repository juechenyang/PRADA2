"""
created by Juechen Yang at 1/14/19

"""
import numpy as np

def filter_junction(junctionDF):
    # chromosome filter that can keep desired choromosome
    chromosomeSet = ['chr1', 'chr8', 'chr2', 'chr12', 'chrX', 'chr10', 'chr22', 'chr5', \
                     'chr11', 'chr7', 'chr15', 'chr19', 'chr3', 'chr21', \
                     'chr13', 'chr18', 'chr9', 'chr16', 'chr6', 'chr14', 'chr4', \
                     'chr17', 'chr20', 'chrY']
    junctionDF = junctionDF[junctionDF['donor_chr'].isin(chromosomeSet)]
    junctionDF = junctionDF[junctionDF['acceptor_chr'].isin(chromosomeSet)]

    junctionDF[['donor_fb','acceptor_fb']] = junctionDF[['donor_fb','acceptor_fb']].astype(int)
    junctionDF['donor_start'] = np.where((junctionDF['junc_type'] != -1) & (junctionDF['donor_strand'] == '+'),
                                         junctionDF['donor_fb'] - 1, junctionDF['donor_start'])
    junctionDF['donor_start'] = np.where((junctionDF['junc_type'] != -1) & (junctionDF['donor_strand'] == '-'),
                                         junctionDF['donor_fb'] + 1, junctionDF['donor_start'])
    junctionDF['acceptor_start'] = np.where((junctionDF['junc_type'] != -1) & (junctionDF['acceptor_strand'] == '+'),
                                            junctionDF['acceptor_fb'] + 1, junctionDF['acceptor_start'])
    junctionDF['acceptor_start'] = np.where((junctionDF['junc_type'] != -1) & (junctionDF['acceptor_strand'] == '-'),
                                            junctionDF['acceptor_fb'] - 1, junctionDF['acceptor_start'])
    # create end position for each junction
    junctionDF['donor_end'] = junctionDF['donor_start']
    junctionDF['acceptor_end'] = junctionDF['acceptor_start']

    # convert position features to integer type
    junctionDF[['donor_start', 'donor_end', 'junc_type', 'acceptor_start', 'acceptor_end']] = junctionDF[
        ['donor_start', 'donor_end', 'junc_type', 'acceptor_start', 'acceptor_end']].astype(int)

    # make junctions compatible 0-based bed file format
    junctionDF['donor_start'] = junctionDF['donor_start'] - 1
    junctionDF['acceptor_start'] = junctionDF['acceptor_start'] - 1

    # make a unique ID for each junction
    junctionDF.insert(0, 'juncID', range(1, len(junctionDF) + 1))
    return junctionDF
