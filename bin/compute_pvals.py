#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys,itertools
import mmh3
import argparse
from scipy import stats
import scipy
import pickle
import nltk

### inputs:
##### either infile (explicit stratified anchors.txt file), or inpath as
###### /oak/stanford/groups/horence/kaitlin/results_nomad/bulk_RNAseq/paper/AA_antibody_secreting_cells/abundant_stratified_anchors/ and slurmID (can be renamed) ranging from 0 to 63 indicating which file to analyze
##### samplesheet: e.g. /oak/stanford/groups/horence/kaitlin/running_nomad/bulk_RNAseq/samplesheets/samplesheet_AA_antibody_secreting_cells.csv
### outputs:
##### args.outfldr+'/pvals_{}.csv'

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( ## if infile exists, use this
        "--infile",
        type=str,
        default=""
    )
    parser.add_argument(
        "--kmer_size",
        type=int
    )
    parser.add_argument( ### takes as input filepath to samplesheet
        "--samplesheet",
        type=str,
        default=""
    )
    parser.add_argument( ### number of hashes to use
        "--K",
        type=int,
        default=10
    )
    parser.add_argument( ### number of random cj to use
        "--L",
        type=int,
        default=50
    )
    parser.add_argument(
        "--outfile_scores",
        type=str,
        default=""
    )
    args = parser.parse_args()
    return args


### method for taking in dij column and generating Sj
def convertDijToSj(df,idx):
    stridx = str(idx)
    dijcol = 'dij_'+stridx
    prodcol = 'sum_'+stridx
    mucol = 'mu_'+stridx
    sjcol = 'sj_'+stridx
    sampSumcol = 'sampSum_'+stridx

    df[prodcol] = df.counts*df[dijcol] ## for each row entry, compute its count * d_i
    df[sampSumcol] = df.groupby(['anchor','sample'])[prodcol].transform('sum') ## for each sample, compute its sum of D_jk
    ### sampSumCol is essentially "mu of this anchor sample"
    df[mucol] = df.groupby(['anchor'])[prodcol].transform('sum')/df.M ## comute mean across this anchor of D_jk
    df[sjcol]=(df[sampSumcol] - df.nj*df[mucol]) / (df.nj**(.5)) ## compute overall sj
    return df


def hamming(x,y): ### works in conjunction with dMax
    return sum(c1 != c2 for c1, c2 in zip(x,y))


def normalizevec(x): ### shift and scale vector so that it's in [-1,1]
    return 2*(x-x.min()) / (x.max()-x.min())-1

def main():
    args = get_args()
    K = args.K
    L = args.L
    dMax=5 ### get rid of soon.

    ### load samplesheet c_j's
    useHandCjs = False
    samplesheet = args.samplesheet
    if samplesheet != "":
        with open(samplesheet,'r') as f:
            cols = f.readline().split(',')
        if len(cols)==1: ### if len(cols) is 1, then only samplesheet name, no ids
            print("Only 1 samplesheet column, using random cjs")
        elif len(cols)>2:
            print("Improperly formatted samplesheet")
        else:
            useHandCjs = True
            sheetdf = pd.read_csv(samplesheet,names=['fname','handCjs'])
            sheetdf['sample'] = (sheetdf.fname
                            .str.rsplit('/',1,expand=True)[1]
                            .str.split('.',1,expand=True)[0])
            sheetdf['handCjs'] = normalizevec(sheetdf['handCjs'])
            sheetdf = sheetdf.drop(columns='fname')
            sheetdf['handCjs'] = normalizevec(sheetdf.handCjs) ### normalize in case time-series
            print("Successfully loaded custom cj")


    df = pd.read_csv(args.infile, delim_whitespace=True, names=['counts','seq','sample'])
    if useHandCjs:
        df = pd.merge(df,sheetdf)
    else:
        df['handCjs']=1
    print('done reading')

    ### split seq into anchor and target
    ### for some reason some files have duplicates, so drop those
    df['anchor'] = df.seq.str[:args.kmer_size]
    df['target'] = df.seq.str[args.kmer_size:]
    df = df.drop(columns='seq').drop_duplicates()

    ### add in informational columns for each anchor for filtering purposes
    ### this is faster than filter command
    df['anch_uniqTargs'] = df.groupby('anchor').target.transform('nunique')
    df['anch_samples']= df.groupby('anchor')['sample'].transform('nunique')
    df['anchSample_cts'] = df.groupby(['anchor','sample']).counts.transform('sum')
    df = df[(df.anch_uniqTargs>1) & (df.anch_samples > 1) & (df.anchSample_cts>5)]
    df['anch_cts'] = df.groupby('anchor').counts.transform('sum') ## number of reads per anchor
    df = df[df.anch_cts > 30]
    ### above line filters, can change parameters
    print('done filtering')

    ### if no anchors left after filtering, exit
    if df.empty:
        return

    df['nj'] = df.groupby(['anchor','sample']).counts.transform('sum')
    df = df.drop(columns='anchSample_cts') ### this is essentially "old" n_j
    df['M'] = df.groupby(['anchor']).counts.transform('sum')


    #### add in handcrafted dij of distance to most frequent target
    ### find most abundant target for every anchor
    df['targ_cts']=df.groupby(['anchor','target']).counts.transform('sum')
    tmpDf = df.sort_values('targ_cts', ascending=False).drop_duplicates(['anchor'])[['anchor','target']].rename(columns={'target':'maxTarget'}) #### keeps only the first anchor occurence
    mergedDf = pd.merge(df,tmpDf)

    ### compute distances from each target to max target and propagate
    ### first call is to generate hand mu (no limit)
    mergedDf['dij_0'] = np.vectorize(lambda x,y : hamming(x,y))(mergedDf.target,mergedDf.maxTarget)
    df = convertDijToSj(mergedDf,0)
    df['mu_ham'] = df['mu_0']
    
    ### generate levenstein distance based mu
    mergedDf['dij_0'] = np.vectorize(lambda x,y : nltk.edit_distance(x,y))(mergedDf.target,mergedDf.maxTarget)
    df = convertDijToSj(mergedDf,0)
    df['mu_lev'] = df['mu_0']

    ### overwrite dij_0 and associated columns with handcrafted
    mergedDf['dij_0'] = np.vectorize(lambda x,y : min(hamming(x,y),dMax)/dMax)(mergedDf.target,mergedDf.maxTarget)
    df = mergedDf.drop(columns=['targ_cts','maxTarget'])
    df = convertDijToSj(df,0)

    #### hash based dij, randomly assign each target to 0 or 1
    for k in range(1,K):
        df['dij_'+str(k)] = (df['target'].apply(lambda x : (mmh3.hash(x,seed=k)>0) ))*1.0

        df = convertDijToSj(df,k)

    df = df.copy() ## get rid of memory issues

    print("starting p value computation")
    #### start efficient computation of p values ####
    sjcols = ['sj_'+str(t)for t in range(K)]
    sampSumcols = ['sampSum_'+str(t)for t in range(K)]

    dftmp = (df[['anchor','sample','nj'] + sjcols+sampSumcols].drop_duplicates().pivot(index='anchor',values=['nj']+ sjcols+sampSumcols,columns='sample'))

    ### define dimensions
    A = len(dftmp)
    p = df['sample'].nunique()

    ### get sample names and samplesheet cjs in order
    sampleNames = dftmp.columns.get_level_values(1)[:p].values
    sheetCj = np.ones(p)
    if useHandCjs:
        sheetCj = sheetdf.set_index('sample').T[sampleNames].to_numpy().flatten()

    njMat = np.zeros((A,p))
    njMat = dftmp['nj'].fillna(0).to_numpy()
    njMat = njMat.T #### p x A
    Marr =njMat.sum(axis=0)

    sjMat = np.zeros((A,p,K)) ### A x p x K
    for i,col in enumerate(sjcols):
        sjMat[:,:,i]=dftmp[col].fillna(0).to_numpy()
    sjMat = np.swapaxes(sjMat,0,2)### K x p x A

    np.random.seed(0)
    cjs = np.random.choice([-1,1],size=(L+1,p))
    cjs[0,:]=sheetCj ### set first row to be samplesheet cjs

    zerodCjs = (cjs[:,:,None]*(njMat>0)[None,:,:])
    ### zero out Cjs that correspond to njs of 0 to increase power
    ### L x p x A

    testStats = np.abs(cjs@sjMat)
    ### K x L x A
    ### sj are 0 when nj is 0, so don't need to use masked cj

    with np.errstate(divide='ignore', invalid='ignore'): ### discard divide by 0 warnings
        ### with zeroed cjs
        num = (zerodCjs**2).sum(axis=1)*Marr ### L x A
        denom=(cjs@np.sqrt(njMat))**2 ### can use normal cjs, since 0 pattern matches
        a = 1.0/(np.sqrt(num/denom)+1.0)
        a[np.isnan(a)] = 0
        ### L x A

        ### compute p values
        term1HC = 2*np.exp(- 2*(1-a)[None,:,:]**2*testStats**2
                       /(zerodCjs**2).sum(axis=1)[None,:,:])
        term2HC = 2*np.exp(-2*a[None,:,:]**2*Marr*testStats**2
                           /(cjs@np.sqrt(njMat))**2 )
        term2HC[:,a==0]=0
        term1HC = np.nan_to_num(term1HC,nan=1)
        # term1HC = np.maximum(term1HC,1) ### if user defined cjs sum to 0, yields nans
        pv = term1HC+term2HC
        ### K x L x A

    pv_hand_sheetCjs = np.minimum(pv[0,0,:],1)
    pv_hash_sheetCjs = np.minimum((K-1)*pv[1:,0,:].min(axis=0),1)
    pv_hand = np.minimum(L*pv[0,1:,:].min(axis=0),1) ### axis 1 (L) becomes axis 0 after slicing
    pv_hash = np.minimum(L*(K-1)*pv[1:,1:,:].min(axis=(0,1)),1)

    ### get cjOpt and best hash
    B = np.swapaxes(pv[1:,1:,:],0,2)
    min_idx = B.reshape(B.shape[0],-1).argmin(1)
    minpos_vect = np.column_stack(np.unravel_index(min_idx, B[0,:,:].shape))+1
    ### A x 2 array
    ### minpos_vect[:,0] ranges 1 to L
    ### minpos_vect[:,1] ranges 1 to K
    cjOpts = zerodCjs[minpos_vect[:,0],:,list(range(A))]
    hashOpts = minpos_vect[:,1]

    ###compute effect size
    sampSumMat = np.zeros((A,p,K)) ### A x p x K
    for i,col in enumerate(sampSumcols):
        sampSumMat[:,:,i]=dftmp[col].fillna(0).to_numpy()
    sampSumMat = np.swapaxes(sampSumMat,0,2)### K x p x A

    tmpMat = sampSumMat[minpos_vect[:,1],:,list(range(A))] * cjOpts
    effectSize = np.abs((tmpMat*(cjOpts>0)).sum(axis=1)/np.maximum(1,(njMat.T*(cjOpts>0)).sum(axis=1))
                  + (tmpMat*(cjOpts<0)).sum(axis=1)/np.maximum(1,(njMat.T*(cjOpts<0)).sum(axis=1))
                 )

    outdf = pd.DataFrame({'anchor':dftmp.index.to_list(),'pv_hash':pv_hash, 'pv_hand':pv_hand,
             'pv_hand_sheetCjs':pv_hand_sheetCjs, 'pv_hash_sheetCjs':pv_hash_sheetCjs,
             'effect_size':effectSize,'optHash':hashOpts, 'M':Marr})

    entDf = (df.groupby(['anchor','target']).counts.sum()
             .groupby('anchor').apply(lambda x : scipy.stats.entropy(x,base=2))
             .reset_index().rename({'counts':'entropy'},axis=1))

    outdf = outdf.merge(entDf)
    outdf = outdf.merge(df[['anchor','mu_hand']].drop_duplicates())

    if not useHandCjs:
        outdf = outdf.drop(columns=['pv_hand_sheetCjs','pv_hash_sheetCjs'])

    outdf= outdf.join(pd.DataFrame(cjOpts,columns=['cj_rand_opt_'+samp for samp in sampleNames])) ## add in cj rand opt
    outdf.sort_values('pv_hash',inplace=True)


    print('writing')

    outdf.to_csv(args.outfile_scores, sep='\t', index=False)

    print("done")

main()
