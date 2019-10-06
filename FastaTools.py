"""
created by Juechen Yang at 12/13/18

"""
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline

def get_max_tx_seq(gene_id, seqDB, tx, exs, fasta_repo):
    desired_txs = tx[tx.gene_id == gene_id]
    desired_txs = desired_txs.sort_values('length', ascending=False)
    tx_id = desired_txs['transcript_id'].iloc[0]
    gene = desired_txs['gene_name'].iloc[0]
    desired_exs = exs[exs.transcript_id == tx_id]
    def getExonSeqs(ex_chr, strand, start, end):
        if strand == '+':
            return str(Seq(str(seqDB.get(ex_chr)[start: end])))
        else:
            return str(Seq(str(seqDB.get(ex_chr)[start: end])).reverse_complement())
    desired_exs['seqs'] = desired_exs.apply(lambda x: getExonSeqs(x['seqname'], x['strand'], x['start']-1, x['end']), axis=1)
    tx_seq = desired_exs.seqs.str.cat(sep='')
    length = desired_exs.length.sum()
    tx_seq_obj = Seq(tx_seq)
    tx_fasta = os.path.join(fasta_repo, gene + '_' + tx_id + '.fasta')
    f = open(tx_fasta, 'w')
    tx_seq_record = SeqRecord(tx_seq_obj, id=gene_id + ' | ' + tx_id + ' | ' + str(length))
    SeqIO.write(tx_seq_record, f, 'fasta')
    f.close()
    return tx_fasta

def parsingFasta(fasta_file):
    record = SeqIO.parse(fasta_file, "fasta")
    seqDB = dict()
    for item in record:
        seqDB[item.id] = item.seq
    return seqDB
def blastFastas(fasta1, fasta2):
    output = NcbiblastnCommandline(query=fasta1, subject=fasta2, outfmt=6, task='blastn')()[0]
    try:
        evalue = output.split('\n')[0].split('\t')[10]
        return evalue
    except IndexError:
        return 100
