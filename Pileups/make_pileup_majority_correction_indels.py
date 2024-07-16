import os
import sys
import pandas as pd
import numpy as np
from Bio import Seq, SeqIO
import itertools
from collections import Counter
from Bio import Align
import multiprocessing as mp
import pysam
import gzip
import seaborn as sns
from scipy import stats
import math
import logomaker
from scipy.spatial.distance import hamming
from glob import glob
from functools import reduce
import bz2

basedir = './'

# AMPLICONS __________________________________________________________________________________________________________________________________

spike = str([rec.seq for rec in SeqIO.parse(os.path.join(basedir, 'amplicons/spike.fasta'), 'fasta')][0])
EGFP = str([rec.seq for rec in SeqIO.parse(os.path.join(basedir, 'amplicons/EGFP.fasta'), 'fasta')][0])

EGFP_F_S = 'CGAGCTGAAGGGCATCGAC'
EGFP_R_S = 'CTCGTCCATGCCGAGAGTGA'
EGFP_amplicon_S = EGFP[EGFP.find(EGFP_F_S):EGFP.find(str(Seq.Seq(EGFP_R_S).reverse_complement()))+len(EGFP_R_S)]

spike_F_S = 'TGCCATCGGCAAGATTCAAG'
spike_R_S = 'TGCCACAAAAGTCGACCCG'
spike_amplicon_S = spike[spike.find(spike_F_S):spike.find(str(Seq.Seq(spike_R_S).reverse_complement()))+len(spike_R_S)]

ampls = {'GFP':EGFP, 'spike':spike}

# READ BAM ___________________________________________________________________________________________________________________________________

print('Reading data...')
bam_file_l1 = sys.argv[1]
bam_file_l2 = sys.argv[2]

name = bam_file_l1.replace('_L001.filtered.bam', '')

ampl_name = 'spike' if bam_file_l1.startswith('SPIKE') else 'GFP'

fh = pysam.AlignmentFile(os.path.join(basedir, 'bwa/filtered', bam_file_l1), 'rb')
bam_data_l1 = [read for read in fh.fetch(until_eof=True)]
fh.close()

fh = pysam.AlignmentFile(os.path.join(basedir, 'bwa/filtered', bam_file_l2), 'rb')
bam_data_l2 = [read for read in fh.fetch(until_eof=True)]
fh.close()

## NICK SITES ________________________________________________________________________________________________________________________________

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
nick_sites = {g:guide_pos[ampl_name][g]+3 if g!='B' else guide_pos[ampl_name][g]+17 for g in guides_by_sample[name]}
    
# UMI STATS __________________________________________________________________________________________________________________________________

UMI_list_1 = [r.qname.split(' ')[0].split('UMI_')[1] for r in bam_data_l1]
UMI_list_2 = [r.qname.split(' ')[0].split('UMI_')[1] for r in bam_data_l2]

UMI_list_1 = pd.Series(Counter(UMI_list_1))
UMI_list_2 = pd.Series(Counter(UMI_list_2))

number_UMIs = UMI_list_1.add(UMI_list_2, fill_value=0)
number_UMIs.to_csv(os.path.join(basedir, 'UMI_stats', name+'_UMI_counts.tsv'), sep='\t')
number_UMIs = number_UMIs.sample(frac=1) # shuffle UMIs to get more balanced chunks

# MAKE PILEUPS FUNCTIONS _____________________________________________________________________________________________________________________

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

def make_pileups(UMI_data, ampl, nick_sites):

    pileup_data = {i:[] for i in range(len(ampl))} # raw pileup
    pileup_data_UMI = {i:{} for i in range(len(ampl))}
    errors_by_position = {'R1':{i:[] for i in range(251)}, 'R2':{i:[] for i in range(251)}}
    errors_by_position_with_correction = {'R1':{i:[] for i in range(251)}, 'R2':{i:[] for i in range(251)}}
    failed_read_names = []
    exception_messages = []
    indel_reads_by_UMI = {}
    UMI_counter = {}
    clipped_reads = []
    
    nick_range = 20

    for qname in UMI_data:
        try:
            pair_data = UMI_data[qname]
            if len(pair_data)==2:
                reads_cov = [{}, {}]
                read_n = []
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
                    read_n.append('R1' if is_read1 else 'R2')

                    # get expanded cigar
                    expanded_cigar = expand_cigar(cigar)

                    iread = 0 # iterates on read
                    iref = refpos # iterates on ref

                    for i in range(len(expanded_cigar)):
                        if expanded_cigar[i] == 'M':
                            reads_cov[n_r][iref] = (seq[iread], qual_33[iread], iread)
                            
                        if iread > 20 and iread < (215 if is_read1 else 230):
                            if expanded_cigar[i] == 'D' or expanded_cigar[i] == 'I':
                                if any([abs(nick_sites[nick]-iref)<nick_range for nick in nick_sites]):
                                    has_indel = True
                                    break
                            # elif expanded_cigar[i] == 'S':
                            #     if any([abs(nick_sites[nick]-iref)<nick_range for nick in nick_sites]):
                            #         is_clipped = True
                            #         break

                        # increase position on read and/or reference
                        if expanded_cigar[i] == 'S' or expanded_cigar[i] == 'I' or expanded_cigar[i] == 'M':
                            iread += 1
                        if expanded_cigar[i] == 'D' or expanded_cigar[i] == 'M':
                            iref += 1
                            
                if has_indel:
                    if UMI not in indel_reads_by_UMI:
                        indel_reads_by_UMI[UMI] = []
                    indel_reads_by_UMI[UMI].append(pair_data)
                elif is_clipped:
                    clipped_reads.append(pair_data)
                else:
                    # make pileup and check for sequencing errors
                    # p = position on the reference
                    common_pos = [p for p in reads_cov[0] if p in reads_cov[1]]
                    # for all position not in common, just add to pileup and to error frequency
                    for read in [0,1]:
                        for p in reads_cov[read]:
                            if p not in common_pos:
                                if reads_cov[read][p][1]>28 and reads_cov[read][p][0]!='N':
                                    pileup_data[p].append(reads_cov[read][p][0])
                                    # check if error
                                    errors_by_position[read_n[read]][reads_cov[read][p][2]].append(reads_cov[read][p][0]!=ampl[p])
                                    errors_by_position_with_correction[read_n[read]][reads_cov[read][p][2]].append(
                                        reads_cov[read][p][0]!=ampl[p])
                                    # add to UMI pileup positions without overlap
                                    if UMI not in pileup_data_UMI[p]:
                                        pileup_data_UMI[p][UMI] = []
                                    pileup_data_UMI[p][UMI].append(reads_cov[read][p][0])

                    for p in common_pos:
                        # error correction on overlapping part of the read
                        if reads_cov[0][p][0]==reads_cov[1][p][0] and reads_cov[0][p][0]!='N': 
                            if min(reads_cov[0][p][1], reads_cov[1][p][1]) > 28:
                                pileup_data[p].append(reads_cov[0][p][0])
                                errors_by_position_with_correction[read_n[0]][reads_cov[0][p][2]].append(reads_cov[0][p][0]!=ampl[p])
                                errors_by_position_with_correction[read_n[1]][reads_cov[1][p][2]].append(reads_cov[1][p][0]!=ampl[p])
                                # add to UMI pileup positions with overlap
                                if UMI not in pileup_data_UMI[p]:
                                    pileup_data_UMI[p][UMI] = []
                                pileup_data_UMI[p][UMI].append(reads_cov[0][p][0])
                        # check if error w/o correction
                        errors_by_position[read_n[0]][reads_cov[0][p][2]].append(reads_cov[0][p][0]!=ampl[p])
                        errors_by_position[read_n[1]][reads_cov[1][p][2]].append(reads_cov[1][p][0]!=ampl[p])
        except Exception as e:
            failed_read_names.append(qname)
            exception_messages.append(str(e))
    
    try:
        # make raw base counts
        pileup_counter = {i:Counter(pileup_data[i]) for i in pileup_data}
        pileup_counter = pd.DataFrame(pileup_counter).transpose()[['A', 'C', 'G', 'T']].fillna(0)
    except Exception as e:
        exception_messages.append(str(e))
        pileup_counter = None
    
    try:
        # make error counts
        error_counter = {R:pd.Series({i:sum(errors_by_position[R][i]) for i in errors_by_position[R]}) for R in errors_by_position}
        error_coverage = {R:pd.Series({i:len(errors_by_position[R][i]) for i in errors_by_position[R]}) for R in errors_by_position}

        # make error counts with correction
        error_counter_with_correction = {R:pd.Series({i:sum(errors_by_position_with_correction[R][i]) for i in errors_by_position_with_correction[R]}
                                                    ) for R in errors_by_position_with_correction}
        error_coverage_with_correction = {R:pd.Series({i:len(errors_by_position_with_correction[R][i]) for i in errors_by_position_with_correction[R]}
                                                     ) for R in errors_by_position_with_correction}
    except Exception as e:
        exception_messages.append(str(e))
        error_counter, error_coverage, error_counter_with_correction, error_coverage_with_correction = None, None, None, None
        
    try:
        indel_reads = []
        indel_UMIs = []
        for UMI in indel_reads_by_UMI:
            if 'N' not in UMI:
                if len(indel_reads_by_UMI[UMI])/UMI_counter[UMI]>0.5:
                    indel_UMIs.append(UMI)
                    for read_pair in indel_reads_by_UMI[UMI]:
                        indel_reads.append(read_pair)
        indel_UMIs = set(indel_UMIs)
    except Exception as e:
        exception_messages.append(str(e))
        indel_reads = []
    
    try:
        # make UMI counts
        pileup_UMI_2 = {i:[] for i in pileup_data_UMI}
        pileup_UMI_3 = {i:[] for i in pileup_data_UMI}
        pileup_UMI_4 = {i:[] for i in pileup_data_UMI}
        pileup_UMI_5 = {i:[] for i in pileup_data_UMI}
        for i in pileup_data_UMI:
            for UMI in pileup_data_UMI[i]:
                if 'N' not in UMI and UMI not in indel_UMIs:
                    UMI_read_bases = pileup_data_UMI[i][UMI]
                    if len(UMI_read_bases)>1:
                        UMI_bases = Counter(UMI_read_bases)
                        max_base = UMI_bases.most_common()[0][0]
                        if UMI_bases[max_base]/len(UMI_read_bases)>0.5:
                            pileup_UMI_2[i].append(max_base)
                            if len(UMI_read_bases)>2:
                                pileup_UMI_3[i].append(max_base)
                                if len(UMI_read_bases)>3:
                                    pileup_UMI_4[i].append(max_base)
                                    if len(UMI_read_bases)>4:
                                        pileup_UMI_5[i].append(max_base)

        
        pileup_counter_UMI_5 = {i:Counter(pileup_UMI_5[i]) for i in pileup_UMI_5}
        pileup_counter_UMI_5 = pd.DataFrame(pileup_counter_UMI_5).transpose()[['A', 'C', 'G', 'T']].fillna(0)
    except Exception as e:
        exception_messages.append(str(e))
        pileup_counter_UMI_5 = None
    
    return [pileup_counter, error_counter, error_coverage, error_counter_with_correction, error_coverage_with_correction,
             failed_read_names, exception_messages,
            pileup_counter_UMI_5, indel_reads, clipped_reads]

# log function
def make_log_result(res, n_chunks):
    def log_result(retval):
        res.append(retval)
        print('{}/{} chunks processed...'.format(len(res), n_chunks))
    return log_result

# RUN ANALYSIS________________________________________________________________________________________________________________________________

ncores = 32
pool = mp.Pool(ncores)

results = []

n_UMIs = len(number_UMIs)
nchunks = ncores*100
chunks = [[i*(n_UMIs//nchunks), (i+1)*(n_UMIs//nchunks)] for i in range(nchunks)]
chunks[-1][-1] = n_UMIs
UMI_chunks = {UMI:n for n,chunk in enumerate(chunks) for UMI in number_UMIs.index[chunk[0]:chunk[1]]}

bam_dict = {i:{} for i in range(len(chunks))}
for r in bam_data_l1:
    UMI = r.qname.split(' ')[0].split('UMI_')[1]
    chunk = UMI_chunks[UMI]
    if r.qname not in bam_dict[chunk]:
        bam_dict[chunk][r.qname] = []
    bam_dict[chunk][r.qname].append([UMI, r.pos, r.cigar, r.seq, r.qual, r.is_read1])
    
for r in bam_data_l2:
    UMI = r.qname.split(' ')[0].split('UMI_')[1]
    chunk = UMI_chunks[UMI]
    if r.qname not in bam_dict[chunk]:
        bam_dict[chunk][r.qname] = []
    bam_dict[chunk][r.qname].append([UMI, r.pos, r.cigar, r.seq, r.qual, r.is_read1])

print('Running analysis...')
for chunk in bam_dict:
    pool.apply_async(make_pileups, args=[bam_dict[chunk], ampls[ampl_name], nick_sites], callback=make_log_result(results, len(chunks)))

pool.close()
pool.join()

# PROCESS RESULTS ____________________________________________________________________________________________________________________________

print('Processing results...')
failed_reads = [pd.Series(r[7]) for r in results if len(r[7])!=0]
if len(failed_reads)!=0:
    failed_reads = pd.concat(failed_reads)
    failed_reads.to_csv(os.path.join(basedir, 'pileups_indels/{}_failed_reads.tsv'.format(name)), sep='\t')

error_messages = [pd.Series(r[8]) for r in results if len(r[8])!=0]
if len(error_messages)!=0:
    error_messages = pd.concat(error_messages)
    error_messages.to_csv(os.path.join(basedir, 'pileups_indels/{}_error_messages.tsv'.format(name)), sep='\t')

pileup_counter_total = reduce(lambda x, y: x.add(y, fill_value=0), [r[0] for r in results if r is not None])
error_counter_total = {R:reduce(lambda x, y: x.add(y, fill_value=0), [r[1][R] for r in results if r is not None]) for R in ['R1', 'R2']}
error_coverage_total = {R:reduce(lambda x, y: x.add(y, fill_value=0), [r[2][R] for r in results if r is not None]) for R in ['R1', 'R2']}
error_counter_with_correction_total = {R:reduce(lambda x, y: x.add(y, fill_value=0), [r[3][R] for r in results if r is not None]) for R in ['R1', 'R2']}
error_coverage_with_correction_total = {R:reduce(lambda x, y: x.add(y, fill_value=0), [r[4][R] for r in results if r is not None]) for R in ['R1', 'R2']}
pileup_counter_UMI_5_total = reduce(lambda x, y: x.add(y, fill_value=0), [r[5] for r in results if r is not None])

error_rate = {R:error_counter_total[R].divide(error_coverage_total[R]).fillna(0) for R in ['R1', 'R2']}
error_rate_df = pd.DataFrame([error_counter_total['R1'], error_coverage_total['R1'], error_rate['R1'],
                              error_counter_total['R2'], error_coverage_total['R2'], error_rate['R2']],
                            index = ['# errors R1', 'cov R1', 'error rate R1', '# errors R2', 'cov R2', 'error rate R2']).transpose()
error_rate_df.to_csv(os.path.join(basedir, 'pileups_indels/{}_error_rate.tsv'.format(name)), sep='\t')

error_rate_with_correction = {R:error_counter_with_correction_total[R].divide(error_coverage_with_correction_total[R]).fillna(0) for R in [
    'R1', 'R2']}
error_rate_with_correction_df = pd.DataFrame([error_counter_with_correction_total['R1'], error_coverage_with_correction_total['R1'],
                                              error_rate_with_correction['R1'], error_counter_with_correction_total['R2'],
                                              error_coverage_with_correction_total['R2'], error_rate_with_correction['R2']],
                                index = ['# errors R1', 'cov R1', 'error rate R1', '# errors R2', 'cov R2', 'error rate R2']).transpose()
error_rate_with_correction_df.to_csv(os.path.join(basedir, 'pileups_indels/{}_error_rate_with_correction.tsv'.format(name)), sep='\t')

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

pileup_counter_total = porcess_pileup(pileup_counter_total, ampls, ampl_name)
pileup_counter_total.to_csv(os.path.join(basedir, 'pileups_indels/{}_pileup.tsv'.format(name)), sep='\t')

pileup_counter_UMI_5_total = porcess_pileup(pileup_counter_UMI_5_total, ampls, ampl_name)
pileup_counter_UMI_5_total.to_csv(os.path.join(basedir, 'pileups_indels/{}_pileup_UMI_5.tsv'.format(name)), sep='\t')

indel_reads_total = pd.DataFrame(reduce(lambda x,y: x + y, [r[11] for r in results]))
clipped_reads_total = pd.DataFrame(reduce(lambda x,y: x + y, [r[12] for r in results]))
with bz2.open(os.path.join(basedir, 'pileups_indels/{}_indel_reads.tsv.bz2'.format(name)), 'wt') as fout:
    indel_reads_total.to_csv(fout, sep='\t')
with bz2.open(os.path.join(basedir, 'pileups_indels/{}_clipped_reads.tsv.bz2'.format(name)), 'wt') as fout:
    clipped_reads_total.to_csv(fout, sep='\t')

print('Done!')