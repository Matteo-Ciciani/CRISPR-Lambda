import os
import pandas as pd
import numpy as np
from Bio import Seq, SeqIO
import itertools
from collections import Counter
import multiprocessing as mp
import pysam
from functools import reduce
from glob import glob
import bz2

# READ DATA

print('Reading data')

datadir = 'pileups_indels/'

indel_data_files = glob(os.path.join(datadir, 'SPIKE*_indel_reads.tsv.bz2'))
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

def make_pileup(chunk_data, ampl, UMI_depth=5):
    pileup_data_UMI = {i:[] for i in range(len(ampl))}
    UMI_counter = {}
    
    for UMI in chunk_data:
        single_UMI_pileup = {i:[] for i in range(len(ampl))}
        if len(chunk_data[UMI])>=UMI_depth:
            for pair_data in chunk_data[UMI]:
                reads_cov = [{}, {}]
                has_indel = False
                is_clipped = False
                UMI = pair_data[0][0]
                if UMI not in UMI_counter:
                    UMI_counter[UMI] = 0
                UMI_counter[UMI] += 1
                
                for n_r in [0,1]:
                    refpos = pair_data[n_r][1] # pos wrt reference
                    seq = pair_data[n_r][3]
                    cigar = pair_data[n_r][2]
                    qual = pair_data[n_r][4]
                    seqpos = 0 # pos wrt the read
                    qual_33 = [ord(c)-33 for c in qual]
                    is_read1 = pair_data[n_r][5]

                    # get expanded cigar
                    expanded_cigar = expand_cigar(cigar)

                    iread = 0 # iterates on read
                    iref = refpos # iterates on ref

                    for i in range(len(expanded_cigar)):
                        if expanded_cigar[i] == 'M':
                            reads_cov[n_r][iref] = (seq[iread], qual_33[iread], iread)

                        # increase position on read and/or reference
                        if expanded_cigar[i] == 'S' or expanded_cigar[i] == 'I' or expanded_cigar[i] == 'M':
                            iread += 1
                        if expanded_cigar[i] == 'D' or expanded_cigar[i] == 'M':
                            iref += 1
                
                # make pileup and check for sequencing errors
                # p = position on the reference
                common_pos = [p for p in reads_cov[0] if p in reads_cov[1]]
                # for all position not in common, just add to pileup and to error frequency
                for read in [0,1]:
                    for p in reads_cov[read]:
                        if p not in common_pos:
                            if reads_cov[read][p][1]>28 and reads_cov[read][p][0]!='N':
                                single_UMI_pileup[p].append(reads_cov[read][p][0])

                for p in common_pos:
                    # error correction on overlapping part of the read
                    if reads_cov[0][p][0]==reads_cov[1][p][0] and reads_cov[0][p][0]!='N': 
                        if min(reads_cov[0][p][1], reads_cov[1][p][1]) > 28:
                            single_UMI_pileup[p].append(reads_cov[0][p][0])
            
            # process UMI pileup
            single_UMI_pileup_consensus = {}
            for p in single_UMI_pileup:
                if len(single_UMI_pileup[p])>0.5*UMI_counter[UMI] and len(single_UMI_pileup[p])>=UMI_depth:
                    UMI_bases = Counter(single_UMI_pileup[p])
                    max_base = UMI_bases.most_common()[0][0]
                    if UMI_bases[max_base]/len(single_UMI_pileup[p])>0.5:
                        pileup_data_UMI[p].append(max_base)
    
    # process pileup_data_UMI
    pileup_data_UMI = {i:Counter(pileup_data_UMI[i]) for i in pileup_data_UMI}
    pileup_data_UMI = pd.DataFrame(pileup_data_UMI).transpose()[['A', 'C', 'G', 'T']].fillna(0)
    return pileup_data_UMI

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
        pool.apply_async(make_pileup, args=[indel_data_chunks_df[f][chunk], ampls[ampl_name]],
                         callback=make_log_result(results, nchunks, f))

pool.close()
pool.join()

# PROCESS RESULTS ____________________________________________________________________________________________________________________________

print('Processing results...')
pileup_counter_UMI = {f:reduce(lambda x, y: x.add(y, fill_value=0), [r for r in results[f] if r is not None]) for f in results}

def porcess_pileup(pileup, ampls, ampl_name):
    pileup = pileup.assign(**{'chr':[ampl_name]*len(ampls[ampl_name]), 'pos':[i+1 for i in range(len(ampls[ampl_name]))],
                              'ref':list(ampls[ampl_name]), 'cov':pileup.sum(axis=1).values})
    af_list = []
    for i in pileup.index:
        a_count = 0
        for b in ['A', 'C', 'G', 'T']:
            if b!=pileup.loc[i, 'ref']:
                a_count += pileup.loc[i, b]
        af_list.append(a_count/pileup.loc[i, 'cov'] if pileup.loc[i, 'cov']!=0 else 0.0)
    pileup = pileup.assign(af=af_list)
    return pileup[['chr','pos','ref','A','C','G','T','af','cov']]

pileup_counter_UMI = {f:porcess_pileup(pileup_counter_UMI[f], ampls, ampl_name) for f in pileup_counter_UMI}
for f in pileup_counter_UMI:
    pileup_counter_UMI[f].to_csv(os.path.join(basedir, 'pileups_indels/{}_pileup_indel_reads.tsv'.format(f)), sep='\t')