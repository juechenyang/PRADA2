#define IO functions. May add more gradually.
import pybedtools,os
import pandas as pd

def getMatching(table, name, exonBed, intermediate_file_dir):
    table.to_csv(os.path.join(intermediate_file_dir, 'junc.bed'), sep='\t', index=None, header=False)
    bed = pybedtools.example_bedtool(os.path.join(intermediate_file_dir, 'junc.bed'))
    OverlapExon = bed.intersect(exonBed, wb=True)
    OverlapExon.saveas(os.path.join(intermediate_file_dir, 'OverlapExon.tsv'))
    matching = pd.read_csv(os.path.join(intermediate_file_dir, 'OverlapExon.tsv'), sep='\t',
                                  names=[name+'_chr', name+'_start', name+'_end', 'juncID', name+'_strand', name+'_fb',
                                         name+'gene_chromosome', name+'gene_start', name+'gene_end', name+'gene_strand', name+'_gene_id',name+'_gene', name+'gene_type'])
    pybedtools.cleanup(remove_all=True)
    return matching





def getCandidates(df, Identifier, samfile, junctionSpanning_th, discordant_th):
    numOfDiscordant = '#discordant_read_pairs'
    discordantReads = 'discordant_reads'
    numOfJunctionSpanning = '#junctionSpanning_reads'
    junctionSpanningReads = 'junctionSpanning_reads'
    numberConsistent = '#Consistent'
    total_reads_at_junction = 'total_reads_at_junction'

    Group = df.groupby(Identifier)
    modf = list()
    extra_columns = [numOfJunctionSpanning, junctionSpanningReads, total_reads_at_junction, 'LeftBreakCandidates', 'LeftBreak', 'RightBreakCandidates', 'RightBreak',
                     numOfDiscordant, numberConsistent, discordantReads]
    for name, group in Group:
        discordantGroup = group[group.junc_type == -1]
        junctionSpanningGroup = group[group.junc_type != -1]

        if len(junctionSpanningGroup.rd_name.unique()) < junctionSpanning_th:
            continue
        if len(discordantGroup.rd_name.unique()) < discordant_th:
            continue
        temp = list(name)
        #add junction spanning reads
        if junctionSpanningGroup.donorgene_strand.unique()[0] == '-':
            donorBreak = junctionSpanningGroup.drop_duplicates(['rd_name']).donor_fb.mode().min() + 1
        else:
            donorBreak = junctionSpanningGroup.drop_duplicates(['rd_name']).donor_fb.mode().max() - 1
        if junctionSpanningGroup.acceptorgene_strand.unique()[0] == '+':
            acceptorBreak = junctionSpanningGroup.drop_duplicates(['rd_name']).acceptor_fb.mode().min() + 1
        else:
            acceptorBreak = junctionSpanningGroup.drop_duplicates(['rd_name']).acceptor_fb.mode().max() - 1
        allreads = list()
        for read in samfile.fetch(junctionSpanningGroup.donor_chr.iloc[0], donorBreak-1, donorBreak):
            allreads.append(read.qname)
        for read in samfile.fetch(junctionSpanningGroup.acceptor_chr.iloc[0], acceptorBreak-1, acceptorBreak):
            allreads.append(read.qname)
        allreads = pd.Series(data=allreads, index=None)
        junctionSpanningGroup = junctionSpanningGroup[junctionSpanningGroup.rd_name.isin(allreads)]
        temp.append(len(junctionSpanningGroup.rd_name.unique()))
        temp.append(''.join(pd.Series(junctionSpanningGroup['rd_name'].unique()).to_string(index=None).replace("\n", "*").split()))
        temp.append(len(allreads.unique()))
        newGroup = junctionSpanningGroup.drop_duplicates(['rd_name']).groupby(
            'donor_fb').size().reset_index(name='num').sort_values(['num'], ascending=False)
        newGroup['new'] = newGroup.donor_fb.map(str) + '#' + newGroup.num.map(str)
        temp.append(''.join(pd.Series(newGroup.new).to_string(index=None).replace("\n", "***").split()))
        temp.append(donorBreak)
        newGroup = junctionSpanningGroup.drop_duplicates(['rd_name']).groupby(
            'acceptor_fb').size().reset_index(name='num').sort_values(['num'], ascending=False)
        newGroup['new'] = newGroup.acceptor_fb.map(str) + '#' + newGroup.num.map(str)
        temp.append(''.join(pd.Series(newGroup.new).to_string(index=None).replace("\n", "***").split()))
        temp.append(acceptorBreak)

        #add discordant reads
        temp.append(len(discordantGroup.rd_name.unique()))
        dop = discordantGroup.donorgene_strand.unique()[0] == '+'
        dom = discordantGroup.donorgene_strand.unique()[0] == '-'
        aop = discordantGroup.acceptorgene_strand.unique()[0] == '+'
        aom = discordantGroup.acceptorgene_strand.unique()[0] == '-'
        if dop and aop:
            temp.append(len(discordantGroup[(discordantGroup.donor_end <= donorBreak) & (
                        discordantGroup.acceptor_end >= acceptorBreak)].rd_name.unique()))
        elif dop and aom:
            temp.append(len(discordantGroup[(discordantGroup.donor_end <= donorBreak) & (
                        discordantGroup.acceptor_end <= acceptorBreak)].rd_name.unique()))
        elif dom and aop:
            temp.append(len(discordantGroup[(discordantGroup.donor_end >= donorBreak) & (
                        discordantGroup.acceptor_end >= acceptorBreak)].rd_name.unique()))
        elif dom and aom:
            temp.append(len(discordantGroup[(discordantGroup.donor_end >= donorBreak) & (
                        discordantGroup.acceptor_end <= acceptorBreak)].rd_name.unique()))
        temp.append(''.join(
            pd.Series(discordantGroup['rd_name'].unique()).to_string(index=None).replace("\n", "*").split()))
        modf.append(temp)
    result = pd.DataFrame(modf, columns=Identifier + extra_columns)

    return result

# def getRealJSR(chr, start, end, samfile, junctionSpanningRds):
#     jsp_total = junctionSpanningRds.split('*')
#     allreads = set()
#     for read in samfile.fetch(chr, start-1, end):
#         allreads.add(read.qname)
#     overlapped = allreads.intersection(set(jsp_total))
#     return len(overlapped), '*'.join(list(overlapped)), len(allreads)