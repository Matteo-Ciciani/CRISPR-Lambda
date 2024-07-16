import os
import pandas as pd
import numpy as np
from Bio import Seq, SeqIO
import itertools
from collections import Counter
import multiprocessing as mp
import pysam
from glob import glob
import bz2

# READ DATA

print('Reading data')

datadir = 'pileups_indels/'

ampl_name = 'GFP'

indel_data_files = glob(os.path.join(datadir, 'GFP*_indel_reads.tsv.bz2'))
indel_data_files.sort()
indel_data_files = indel_data_files[::-1]

indel_data = {f.split('/')[-1].replace('_indel_reads.tsv.bz2', ''):pd.read_csv(bz2.open(f, 'rt'), sep='\t', index_col=0)
              for f in indel_data_files}

indel_data = {f:pd.concat([indel_data[f]['0'].apply(eval), indel_data[f]['1'].apply(eval)], axis=1) for f in indel_data}
indel_data = {f:indel_data[f].assign(UMI=indel_data[f]['0'].apply(lambda x: x[0])) for f in indel_data}

# FORMAT DATA

print('Formatting data')

ncores = 32
number_UMIs = {f:pd.Series(Counter(indel_data[f]['UMI'])) for f in indel_data}
n_UMIs = {f:len(number_UMIs[f]) for f in indel_data}
nchunks = {f:ncores*100 if f=='GFP-1-B_S3' else ncores for f in indel_data}
chunks = {f:[[i*(n_UMIs[f]//nchunks[f]), (i+1)*(n_UMIs[f]//nchunks[f])] for i in range(nchunks[f])] for f in indel_data}
for f in chunks:
    chunks[f][-1][-1] = n_UMIs[f]
UMI_chunks = {f:{UMI:n for n,chunk in enumerate(chunks[f]) for UMI in number_UMIs[f].index[chunk[0]:chunk[1]]} for f in indel_data}

indel_data_chunks_df = {f:{i:{} for i in range(len(chunks[f]))} for f in indel_data}
for f in indel_data:
    for i in indel_data[f].index:
        UMI = indel_data[f].loc[i, 'UMI']
        chunk = UMI_chunks[f][UMI]
        if UMI not in indel_data_chunks_df[f][chunk]:
            indel_data_chunks_df[f][chunk][UMI] = []
        indel_data_chunks_df[f][chunk][UMI].append(indel_data[f].loc[i])
        
# READ AMPLICONS

basedir = './'

spike = str([rec.seq for rec in SeqIO.parse(os.path.join(basedir, 'amplicons/spike.fasta'), 'fasta')][0])
EGFP = str([rec.seq for rec in SeqIO.parse(os.path.join(basedir, 'amplicons/EGFP.fasta'), 'fasta')][0])

EGFP_F_S = 'CGAGCTGAAGGGCATCGAC'
EGFP_R_S = 'CTCGTCCATGCCGAGAGTGA'
EGFP_amplicon_S = EGFP[EGFP.find(EGFP_F_S):EGFP.find(str(Seq.Seq(EGFP_R_S).reverse_complement()))+len(EGFP_R_S)]

spike_F_S = 'TGCCATCGGCAAGATTCAAG'
spike_R_S = 'TGCCACAAAAGTCGACCCG'
spike_amplicon_S = spike[spike.find(spike_F_S):spike.find(str(Seq.Seq(spike_R_S).reverse_complement()))+len(spike_R_S)]

ampls = {'GFP':EGFP, 'spike':spike}

ampl_name = 'spike' if list(indel_data_chunks_df.keys())[0].startswith('SPIKE') else 'GFP'

guide_seq = {'GFP':{'B':'caacgagaagcgcgatcaca', '1':'gagcTaagaccccaacgaga'}, 'spike':{'1':'CGACCCACCGGAAGCAGAAG', 'B':'CAAATTGATCGCCTGATAAC'}}
guide_pos = {ampl:{g:ampls[ampl].find(guide_seq[ampl][g].upper()) for g in guide_seq[ampl]} for ampl in guide_seq}
guides_by_sample = {
    'GFP-EMPTY_S1':['1'],
    'GFP-1_S2':['1'],
    'GFP-1-B_S3':['1', 'B'],
    'SPIKE-RT-EMPTY_S9':['1', 'B'],
    'SPIKE-RT-1-B_S10':['1', 'B'],
    'SPIKE-EVO-EMPTY_S11':['1', 'B'],
    'SPIKE-EVO-1-B_S12':['1', 'B']
}
nick_sites = {f:{g:guide_pos[ampl_name][g]+3 if g!='B' else guide_pos[ampl_name][g]+17 for g in guides_by_sample[f]}
              for f in indel_data}

# DEFINE FUNCTIONS

def expand_cigar(cigar):
    ec = ''
    for tup in cigar:
        if tup[0] == 0: c = 'M'
        elif tup[0] == 1: c = 'I'
        elif tup[0] == 2: c = 'D'
        elif tup[0] == 4: c = 'S'
        else: c = ''
        ec += c*tup[1]
    return ec

def process_pair_indels(read_in, read_del):
    # make read pair result: for each region around nicks, get # indels and length
    pair_in = {}
    pair_del = {}
    all_del = []
    
    # insertions
    for read_n in ['0', '1']:
        for insertion in read_in[read_n]:
            if insertion[0] not in pair_in or pair_in[insertion[0]]<insertion[1]:
                pair_in[insertion[0]] = insertion[1]
    
    # deletion
    for read_n in ['0', '1']:
        for deletion in read_del[read_n]:
            for i in range(deletion[0], deletion[0]+deletion[1]):
                all_del.append(i)
    deleted_bases = np.unique(all_del)
    deleted_bases.sort()
    if len(deleted_bases)>0:
        pair_del = {}
        start = deleted_bases[0]
        length = 1
        for base in deleted_bases[1:]:
            if base == start+length:
                length += 1
            else:
                pair_del[start] = length
                start = base
                length = 1
        pair_del[start] = length
        
    return pair_in, pair_del

def process_read_pair(pair_data):
    read_n = []
    read_cov = {'0':[], '1':[]}
    read_in = {'0':[], '1':[]}
    read_del = {'0':[], '1':[]}

    for n_r in ['0', '1']:
        refpos = pair_data[n_r][1] # pos wrt reference
        cigar = pair_data[n_r][2]
        seq = pair_data[n_r][3]
        qual = pair_data[n_r][4]
        seqpos = 0 # pos wrt the read
        qual_33 = [ord(c)-33 for c in qual]
        is_read1 = pair_data[n_r][5]
        read_n.append('R1' if is_read1 else 'R2')

        # get expanded cigar
        expanded_cigar = expand_cigar(cigar)

        iread = 0 # iterates on read
        iref = refpos # iterates on ref

        indel_len = 0 # count indel len
        indel_start = 0 # track where the indels started

        in_insertion = False
        in_deletion = False

        for i in range(len(expanded_cigar)):
            if expanded_cigar[i] == 'M' or expanded_cigar[i] == 'D':
                read_cov[n_r].append(iref) # add the position in the reference to the bases covered by the read

            # if already in indel, continue expanding it until it ends
            if expanded_cigar[i] == 'D':
                if in_deletion:
                    indel_len += 1
                else:
                    if in_insertion:
                        read_in[n_r].append([indel_start, indel_len])
                        in_insertion = False
                    indel_len = 1 # count indel len
                    in_deletion = True
                    indel_start = iref # track where the indels started

            if expanded_cigar[i] == 'I':
                if in_insertion:
                    indel_len += 1
                else:
                    if in_deletion:
                        read_del[n_r].append([indel_start, indel_len])
                        in_deletion = False
                    indel_len = 1 # count indel len
                    in_insertion = True
                    indel_start = iref # track where the indels started                             

            # increase position on read and/or reference
            if expanded_cigar[i] == 'S' or expanded_cigar[i] == 'I' or expanded_cigar[i] == 'M':
                iread += 1
            if expanded_cigar[i] == 'D' or expanded_cigar[i] == 'M':
                iref += 1
        # save last indel if present
        if in_insertion:
            read_in[n_r].append([indel_start, indel_len])
        if in_deletion:
            read_del[n_r].append([indel_start, indel_len])
    return read_n, read_cov, read_in, read_del

def collapse_UMI_insertions(UMI_ins):
    UMI_ins = pd.DataFrame(UMI_ins).transpose()
    UMI_ins = UMI_ins.sort_index()
    indel_window = 5
    indexes_in_window = [UMI_ins.index[0]]
    collapsed_indels = {}

    for i in UMI_ins.index[1:]:
        # i is close enough to other in window and all indels are mutually exclusive
        if i-indexes_in_window[-1]<=indel_window and all(UMI_ins.loc[indexes_in_window+[i]].notna().sum(axis=0)<=1):
            indexes_in_window.append(i)
        else:
            # check that indel window is supported by more than half the reads
            if UMI_ins.loc[indexes_in_window].notna().sum().sum()>0.5*UMI_ins.shape[1]:
                length = int(UMI_ins.loc[indexes_in_window].mean().mode()[0])
                best_pos = UMI_ins.loc[indexes_in_window].notna().sum(axis=1).idxmax()
                collapsed_indels[best_pos] = length
            indexes_in_window = [i]
    # append last window
    if UMI_ins.loc[indexes_in_window].notna().sum().sum()>0.5*UMI_ins.shape[1]:
        length = int(UMI_ins.loc[indexes_in_window].mean().mode()[0])
        best_pos = UMI_ins.loc[indexes_in_window].notna().sum(axis=1).idxmax()
        collapsed_indels[best_pos] = length
    return pd.Series(collapsed_indels, dtype='int')

def collapse_UMI_deletions(UMI_deletions):
    del_pos_counter = pd.Series(Counter([i for pair in UMI_deletions for pos in pair for i in range(
        pos,pos+pair[pos])]))
    del_pos_counter = del_pos_counter[del_pos_counter>0.5*len(UMI_deletions)]
    consensus_del = {}
    if len(del_pos_counter)>0:
        start = del_pos_counter.index[0]
        length = 1
        for i in del_pos_counter.index[1:]:
            if i==start+length:
                length +=1
            else:
                consensus_del[start] = length
                start = i
                length = 1
        consensus_del[start] = length
    return pd.Series(consensus_del, dtype='int')

def process_UMI(UMI_insertions, UMI_deletions):
    # run a scanning window ~2 nt along the reads, allowing for indels prositioned slightly differently
    # there's an argument to treat deletions differently, considering the overlap based on the deletion length
    collapsed_in = collapse_UMI_insertions(UMI_insertions) if any([len(x)>0 for x in UMI_insertions]) else pd.Series(dtype='float64')
    collapsed_del = collapse_UMI_deletions(UMI_deletions) if any([len(x)>0 for x in UMI_deletions]) else pd.Series(dtype='float64')
    return collapsed_in, collapsed_del

def count_indels(chunk_data, ampl, UMI_depth=5):
    insertions_UMI = {i:{} for i in range(len(ampl))} # positions should be shifted by -0.5 when plotted
    deletions_UMI = {i:{} for i in range(len(ampl))}
    n_indel_reads_UMI = {'Insertions':0, 'Deletions':0}
    
    for UMI in chunk_data:
        if len(chunk_data[UMI])>=UMI_depth:
            pair_in_list = []
            pair_del_list = []
            for pair_data in chunk_data[UMI]:
                read_n, read_cov, read_in, read_del = process_read_pair(pair_data)
                pair_in, pair_del = process_pair_indels(read_in, read_del)
                pair_in_list.append(pair_in)
                pair_del_list.append(pair_del)
            UMI_in, UMI_del = process_UMI(pair_in_list, pair_del_list)

            # add indels to results
            for i in UMI_in.index:
                if UMI not in insertions_UMI[i]:
                    insertions_UMI[i][UMI] = int(UMI_in.loc[i])
            
            for i in UMI_del.index:
                if UMI not in deletions_UMI[i]:
                    deletions_UMI[i][UMI] = int(UMI_del.loc[i])
                    
            if len(UMI_in)>0:
                n_indel_reads_UMI['Insertions'] +=1
            if len(UMI_del)>0:
                n_indel_reads_UMI['Deletions'] +=1
            
    return insertions_UMI, deletions_UMI, n_indel_reads_UMI

def make_log_result(res, nchunks, f):
    def log_result(retval):
        res[f].append(retval)
        print('{}: {}/{} chunks processed...'.format(f, len(res[f]), nchunks[f]))
    return log_result

# RUN ANALYSIS

ncores = 32
pool = mp.Pool(ncores)

results = {f:[] for f in indel_data_chunks_df}
print('Running analysis...')
for f in indel_data_chunks_df:
    for chunk in indel_data_chunks_df[f]:
        pool.apply_async(count_indels, args=[indel_data_chunks_df[f][chunk], ampls[ampl_name]],
                         callback=make_log_result(results, nchunks, f))

pool.close()
pool.join()

# PROCESS RESULTS

overall_insertions = {f:{} for f in results}
overall_deletions = {f:{} for f in results}
overall_counts = {f:{'Insertions':0, 'Deletions':0} for f in results}
for f in results:
    for res in results[f]:
        for pos in res[0]: # insertions
            if pos not in overall_insertions[f]:
                overall_insertions[f][pos] = {}
            for UMI in res[0][pos]:
                overall_insertions[f][pos][UMI] = res[0][pos][UMI]
        for pos in res[1]:
            if pos not in overall_deletions[f]:
                overall_deletions[f][pos] = {}
            for UMI in res[1][pos]:
                overall_deletions[f][pos][UMI] = res[1][pos][UMI]
        for indel in res[2]:
            overall_counts[f][indel] += res[2][indel]
            
overall_counts = pd.DataFrame(overall_counts).transpose()
overall_counts.to_csv(os.path.join(datadir, ampl_name + '_indel_read_counts_UMI_5.tsv'), sep='\t')

for f in overall_insertions:
    with bz2.open(os.path.join(datadir, f+'_insertions_by_UMI_5.tsv.bz2'), 'wt') as fout:
        pd.Series(overall_insertions[f]).to_csv(fout, sep='\t')
        
for f in overall_deletions:
    with bz2.open(os.path.join(datadir, f+'_deletions_by_UMI_5.tsv.bz2'), 'wt') as fout:
        pd.Series(overall_deletions[f]).to_csv(fout, sep='\t')