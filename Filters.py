"""
created by Juechen Yang at 1/14/19

"""
import numpy as np
def switchDonorAcceptor(donorF, acceptorF, df):
    df[donorF], df[acceptorF], = np.where(df['switchOrNot'] == 'switch',
                                                    [df[acceptorF], df[donorF]],
                                                    [df[donorF], df[acceptorF]])


def filterOrientation(df):

    #create an extra column called flipped_acceptor_strand
    df['flipped_acceptor_strand'] = np.where(df['acceptor_strand'] == '+', '-', '+')

    #create another column
    df['switchOrNot'] = np.where(
        ((df['donor_strand'] != df['donorgene_strand']) & (df['flipped_acceptor_strand'] == df['acceptorgene_strand'])),
        'switch', 'fail')
    df['switchOrNot'] = np.where(
        (df['donor_strand'] == df['donorgene_strand']) & (df['flipped_acceptor_strand'] != df['acceptorgene_strand']),
        'noSwitch', df.switchOrNot)
    crNew = df[df['switchOrNot'] != 'fail']

    #switch columns based "switchOrNot" rules
    switchDonorAcceptor('donor_start', 'acceptor_start', crNew)
    switchDonorAcceptor('donor_end', 'acceptor_end', crNew)
    switchDonorAcceptor('donor_chr', 'acceptor_chr', crNew)
    switchDonorAcceptor('donor_strand', 'acceptor_strand', crNew)
    switchDonorAcceptor('donorgene_chromosome', 'acceptorgene_chromosome', crNew)
    switchDonorAcceptor('donorgene_start', 'acceptorgene_start', crNew)
    switchDonorAcceptor('donorgene_end', 'acceptorgene_end', crNew)
    switchDonorAcceptor('donorgene_strand', 'acceptorgene_strand', crNew)
    switchDonorAcceptor('donor_gene', 'acceptor_gene', crNew)
    switchDonorAcceptor('donor_gene_id', 'acceptor_gene_id', crNew)
    switchDonorAcceptor('donorgene_type', 'acceptorgene_type', crNew)
    switchDonorAcceptor('donor_fb', 'acceptor_fb', crNew)

    return crNew

def sameGeneFilter(result):
    result = result[((result['donor_gene'] != result['acceptor_gene']) & ~(
                result['donor_gene'].isnull() | result['acceptor_gene'].isnull()))]
    result = result.drop_duplicates(result.columns)

    return result

# def readThroughFilter(result):
#     result['abs_gap'] = abs(result.donor_end - result.acceptor_end)
#     result = result[~((result.junc_type == -1) & (result.donorgene_strand == result.acceptorgene_strand) & (result.abs_gap <= 1000000)
#                       & (result.donor_strand == '+') & (result.acceptor_strand == '-')
#                       & (result.donor_end < result.acceptor_end))]
#     result = result[~((result.junc_type == -1) & (result.donorgene_strand == result.acceptorgene_strand) & (result.abs_gap <= 1000000)
#                       & (result.donor_strand == '-') & (result.acceptor_strand == '+')
#                       & (result.donor_end > result.acceptor_end))]
#     return result

def readThroughFilter(result):
    result['abs_gap'] = abs(result.donor_end - result.acceptor_end)
    result = result[~((result.junc_type == -1) & (result.donorgene_strand == result.acceptorgene_strand) & (result.abs_gap <= 1000000)
                      & (result.donor_chr == result.acceptor_chr) & (result.donorgene_strand == '+')
                      & (result.donor_end < result.acceptor_end))]
    result = result[~((result.junc_type == -1) & (result.donorgene_strand == result.acceptorgene_strand) & (result.abs_gap <= 1000000)
                      & (result.donor_chr == result.acceptor_chr) & (result.donorgene_strand == '-')
                      & (result.donor_end > result.acceptor_end))]
    return result

def overlap_exon_filter(dstart, dend, astart, aend, dchr, achr):
    if dchr == achr:
        if dstart < astart and dend >= astart:
            return 'Yes'
        elif dstart >= astart and dstart <= aend:
            return 'Yes'
        else:
            return 'No'
    else:
        return 'No'


