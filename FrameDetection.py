"""
created by Juechen Yang at 1/30/19

"""
import pandas as pd


# detect if junction break is in exon
def is_in_exon(start, end, break_point):
    if break_point <= end and break_point >= start:
        return True
    else:
        return False
# main function for frame detection
def FrameDetection(dg_id, ag_id, dg_break, ag_break, gtf):
    # get all exons based on gene_id
    donor_exons = gtf[(gtf.feature == 'exon') & (gtf.gene_id == dg_id)]
    acceptor_exons = gtf[(gtf.feature == 'exon') & (gtf.gene_id == ag_id)]
    # find all exons that overlapped with junction break
    donor_exons['cover_break'] = donor_exons.apply(lambda x: is_in_exon(x['start'], x['end'], dg_break), axis=1)
    acceptor_exons['cover_break'] = acceptor_exons.apply(lambda x: is_in_exon(x['start'], x['end'], ag_break), axis=1)
    # find all qualified transcripts
    donor_qf_tx_name = donor_exons[donor_exons.cover_break == True].transcript_id.unique()
    acceptor_qf_tx_name = acceptor_exons[acceptor_exons.cover_break == True].transcript_id.unique()

    donor_qf_tx = gtf[((gtf.feature == 'start_codon') | (gtf.feature == 'exon') | (gtf.feature == 'stop_codon')) &
                      (gtf.transcript_id.isin(donor_qf_tx_name))][
        ['feature', 'gene_id', 'transcript_id', 'gene_name', 'strand', 'start', 'end']]
    acceptor_qf_tx = gtf[((gtf.feature == 'start_codon') | (gtf.feature == 'exon') | (gtf.feature == 'stop_codon')) &
                         (gtf.transcript_id.isin(acceptor_qf_tx_name))][
        ['feature', 'gene_id', 'transcript_id', 'gene_name', 'strand', 'start', 'end']]

    donor_coding_count = get_conding_count(donor_qf_tx, dg_break, False)
    acceptor_coding_count = get_conding_count(acceptor_qf_tx, ag_break, True)
    donor_coding_count['merge_id'] = 0
    acceptor_coding_count['merge_id'] = 0
    merged = pd.merge(donor_coding_count, acceptor_coding_count, how='outer', on=['merge_id'])
    merged = merged.drop(columns=['merge_id'])
    return merged

# get distance by giving start and end pos for both entity and region
def getDistance(e_start, e_end, r_start, r_end):

    if e_start >= r_start and e_end <= r_end:
        distance = e_end - e_start + 1
    elif e_start < r_start and e_end >= r_start and r_end <= r_end:
        distance = e_end - r_start + 1
    elif e_start >= r_start and e_start <= r_end and e_end > r_end:
        distance = r_end - e_start + 1
    elif e_start < r_start and e_end > r_end:
        distance = r_end - r_start + 1
    else:
        distance = 0
    return distance
'''
compute the number of nucleotide between start codon and junction break, 
if either start or stop codon not found then distance was given as -1 which
indicates N/A for computation. if junction break was not fall into region 
between start and stop codon then distance was given as -2 indicates UTR break

'''
def get_conding_count(qf_tx, threshold, ac_side):
    gene_id = qf_tx.gene_id.iloc[0]
    gene_name = qf_tx.gene_name.iloc[0]
    qf_tx_group = qf_tx.groupby('transcript_id')
    countDF = pd.DataFrame()
    in_frame = list()
    for name, group in qf_tx_group:
        transcript_id = group.transcript_id.iloc[0]
        if ac_side and group.strand.unique()[0] == '+':
            threshold = threshold - 1
        elif ac_side and group.strand.unique()[0] == '-':
            threshold = threshold + 1
        if group.strand.unique()[0] == '+':
            if not ('start_codon' in group['feature'].unique() and 'stop_codon' in group['feature'].unique()):
                number = -1
                in_frame.append([gene_id, gene_name, transcript_id, number])
                continue
            coding_start = group[group.feature == 'start_codon']['start'].iloc[0]
            coding_end = group[group.feature == 'stop_codon']['start'].iloc[0]
            if not (threshold >= coding_start and threshold <= coding_end):
                number = -2
                in_frame.append([gene_id, gene_name, transcript_id, number])
                continue
            group['distance'] = group.apply(lambda x: getDistance(x['start'], x['end'], coding_start, threshold), axis=1)
            number = group['distance'].sum()
        else:
            if not ('start_codon' in group['feature'].unique() and 'stop_codon' in group['feature'].unique()):
                number = -1
                in_frame.append([gene_id, gene_name, transcript_id, number])
                continue
            coding_start = group[group.feature == 'start_codon']['end'].iloc[0]
            coding_end = group[group.feature == 'stop_codon']['end'].iloc[0]
            if not (threshold >= coding_end and threshold <= coding_start):
                number = -2
                in_frame.append([gene_id, gene_name, transcript_id, number])
                continue
            group['distance'] = group.apply(lambda x: getDistance(x['start'], x['end'], threshold, coding_start), axis=1)
            number = group['distance'].sum()
        in_frame.append([gene_id, gene_name, transcript_id, number])
    countDF = countDF.append(pd.DataFrame(data=in_frame, columns=['gene_id', 'gene_name', 'transcript_id', 'distance_startCodon_to_break']))
    return countDF

# determine fusion transcript type
def is_shifting(donor_distance, acceptor_distance):
    if donor_distance == -1 or acceptor_distance == -1:
        return 'N/A'
    if donor_distance == -2 and acceptor_distance != -2:
        return 'UTR/CDS'
    elif donor_distance != -2 and acceptor_distance == -2:
        return 'CDS/UTR'
    elif donor_distance == -2 and acceptor_distance == -2:
        return 'UTR/UTR'
    elif donor_distance%3 == acceptor_distance%3:
        return 'in-frame'
    elif donor_distance%3 != acceptor_distance%3:
        return 'out-frame'
