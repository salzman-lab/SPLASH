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
import time #### just for testing

### inputs:
#####  stratified anchors.txt file
##### samplesheet: e.g. /oak/stanford/groups/horence/kaitlin/running_nomad/bulk_RNAseq/samplesheets/samplesheet_AA_antibody_secreting_cells.csv
### outputs:
##### args.outfldr+'/pvals_{}.csv'

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( ### input abundant_stratified_{}.txt.gz file 
        "--infile",
        type=str
    )
    parser.add_argument( ### output file path, required
        "--outfile_scores",
        type=str
    )
    parser.add_argument( ### samplesheet file, if samplesheet-based cj are to be used
        "--samplesheet",
        type=str,
        default=""
    )
    parser.add_argument( ### length of k-mer used
        "--kmer_size",
        type=int,
        default=27
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
    parser.add_argument( ### remove all anchors with y or fewer unique targets
        "--anchor_unique_targets_threshold",
        type=int,
        default=1
    )
    parser.add_argument( ### remove all anchors with y or fewer total counts 
        "--anchor_count_threshold",
        type=int,
        default=30
    )
    parser.add_argument( ### remove all anchors that only appear in y or fewer samples
        "--anchor_samples_threshold",
        type=int,
        default=1
    )
    parser.add_argument( ### remove all (anchor,sample) pairs with y or fewer counts
        "--anchor_sample_counts_threshold",
        type=int,
        default=1
    )
    parser.add_argument( ### number of anchors batched together in p-value computation
        "--anchor_batch_size",
        type=int,
        default=1000
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


def hamming(x,y):
    return sum(c1 != c2 for c1, c2 in zip(x,y))


def normalizevec(x): ### shift and scale vector so that it's in [-1,1]
    return 2*(x-x.min()) / (x.max()-x.min())-1

def main():
    args = get_args()
    K = args.K
    L = args.L

    print("running compute_pvals script on {}".format(args.infile))

    ### load samplesheet c_j's
    useSheetCjs = False
    samplesheet = args.samplesheet
    if samplesheet != "":
        with open(samplesheet,'r') as f:
            cols = f.readline().split(',')
        if len(cols)==1: ### if len(cols) is 1, then only samplesheet name, no ids
            print("Only 1 samplesheet column, using random cjs")
        elif len(cols)>2:
            print("Improperly formatted samplesheet")
        else:
            useSheetCjs = True
            sheetdf = pd.read_csv(samplesheet,names=['fname','sheetCj'])
            sheetdf['sample'] = (sheetdf.fname
                            .str.rsplit('/',1,expand=True)[1]
                            .str.split('.',1,expand=True)[0])
            sheetdf['sheetCj'] = normalizevec(sheetdf['sheetCj'])
            sheetdf = sheetdf.drop(columns='fname')
            sheetdf['sheetCj'] = normalizevec(sheetdf.sheetCj) ### normalize in case time-series
            print("Successfully loaded custom cj")


    df = pd.read_csv(args.infile, delim_whitespace=True, names=['counts','seq','sample'])
    if useSheetCjs:
        df = pd.merge(df,sheetdf)
    else:
        df['sheetCj']=1
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
    df = df[
        (df.anch_uniqTargs > args.anchor_unique_targets_threshold) &
        (df.anch_samples > args.anchor_samples_threshold) &
        (df.anchSample_cts > args.anchor_sample_counts_threshold)
    ]
    df['anch_cts'] = df.groupby('anchor').counts.transform('sum') ## number of reads per anchor
    df = df[df.anch_cts > args.anchor_count_threshold]

    print('done filtering')

    ### if no anchors left after filtering, exit
    if df.empty:
        print("no anchors left in dataframe, exiting")
        return

    df['nj'] = df.groupby(['anchor','sample']).counts.transform('sum')
    df = df.drop(columns='anchSample_cts') ### this is essentially "old" n_j
    df['M'] = df.groupby(['anchor']).counts.transform('sum')
    df['number_nonzero_samples'] = df.groupby('anchor')['sample'].nunique() ## count number of samples that this anchor appears in
    
    print('starting hamming and levenshtein computation')
    #### add in handcrafted dij of distance to most frequent target
    ### find most abundant target for every anchor
    df['targ_cts']=df.groupby(['anchor','target']).counts.transform('sum')
    tmpDf = df.sort_values('targ_cts', ascending=False).drop_duplicates(['anchor'])[['anchor','target']].rename(columns={'target':'maxTarget'}) #### keeps only the first anchor occurence
    mergedDf = pd.merge(df,tmpDf)

    ### compute distances from each target to max target and propagate
    ### first call is to generate hand mu (no limit)
    mergedDf['dij_0'] = np.vectorize(lambda x,y : hamming(x,y))(mergedDf.target,mergedDf.maxTarget)
    df = convertDijToSj(mergedDf,0)
    df['mean_target_hamming_distance'] = df['mu_0']

    ### generate levenstein distance based mu
    # mergedDf['dij_0'] = np.vectorize(lambda x,y : nltk.edit_distance(x,y))(mergedDf.target,mergedDf.maxTarget)
    mergedDf['dij_0'] = np.vectorize(lambda x,y : hamming(x,y))(mergedDf.target,mergedDf.maxTarget) ######### DELETE THIS######################
    df = convertDijToSj(mergedDf,0)
    df['mean_target_levenshtein_distance'] = df['mu_0']

    ### overwrite dij_0 and associated columns with handcrafted
    # mergedDf['dij_0'] = np.vectorize(lambda x,y : hamming(x,y))(mergedDf.target,mergedDf.maxTarget) ##### to do, delete?
    # df = convertDijToSj(df,0)

    df = mergedDf.drop(columns=['targ_cts','maxTarget'])


    #### hash based dij, randomly assign each target to 0 or 1
    for k in range(1,K):
        df['dij_'+str(k)] = (df['target'].apply(lambda x : (mmh3.hash(x,seed=k)>0) ))*1.0

        df = convertDijToSj(df,k)

    df = df.copy() ## get rid of memory issues

    print("starting p value computation on {} anchors".format(df.anchor.nunique()))
    startTime = time.time()
    #### start efficient computation of p values ####
    sjcols = ['sj_'+str(t)for t in range(K)]
    sampSumcols = ['sampSum_'+str(t)for t in range(K)]

    df_pivoted_full = (df[['anchor','sample','nj'] + sjcols+sampSumcols].drop_duplicates().pivot(index='anchor',values=['nj']+ sjcols+sampSumcols,columns='sample'))

    ##### Following step can be batched across anchors to reduce memory inefficiency
    ##### i.e., consider only 10k anchors at once

    ### define dimensions
    Atotal = len(df_pivoted_full)
    p = df['sample'].nunique()

    ### get sample names and samplesheet cjs in order
    sampleNames = df_pivoted_full.columns.get_level_values(1)[:p].values
    sheetCj = np.ones(p)
    if useSheetCjs:
        sheetCj = sheetdf.set_index('sample').T[sampleNames].to_numpy().flatten()

    

    #### batched computation start
    anchor_batch_size = args.anchor_batch_size

    pval_random = np.ones(Atotal)
    pval_samplesheet = np.ones(Atotal)
    effect_size_random = np.zeros(Atotal)
    effect_size_samplesheet = np.zeros(Atotal) 
    optHash = np.zeros(Atotal)
    fullMarr = np.zeros(Atotal)
    cjOptMat = np.zeros((Atotal,p))


    #### can change to something like the below?
    # anchor_batch_size = int(10**8 / len(sampleNames) / numRandomCj)

    for i in range(int(Atotal//anchor_batch_size)+1):
        ### operate on dftmp[i*anchor_batch_size : (i+1)*anchor_batch_size]
        # A = len(dftmp[i*anchor_batch_size : (i+1)*anchor_batch_size])
        idx_start = i*anchor_batch_size
        idx_end = min((i+1)*anchor_batch_size,len(df_pivoted_full))
        A = idx_end-idx_start ### number of anchors used here
        # print(dftmp.iloc[0:5]['nj'])
        # continue

        dftmp = df_pivoted_full.iloc[i*anchor_batch_size : (i+1)*anchor_batch_size]
        print(A,len(dftmp),idx_start,idx_end,len(df_pivoted_full))
        assert A==len(dftmp)

        njMat = np.zeros((A,p))
        njMat = dftmp['nj'].fillna(0).to_numpy()
        njMat = njMat.T #### p x Atotal
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

        ###compute effect size for random cj
        sampSumMat = np.zeros((A,p,K)) ### A x p x K
        for i,col in enumerate(sampSumcols):
            sampSumMat[:,:,i]=dftmp[col].fillna(0).to_numpy()
        sampSumMat = np.swapaxes(sampSumMat,0,2)### K x p x A

        tmpMat = sampSumMat[minpos_vect[:,1],:,list(range(A))] * cjOpts
        effectSize = np.abs((tmpMat*(cjOpts>0)).sum(axis=1)/np.maximum(1,(njMat.T*(cjOpts>0)).sum(axis=1))
                    + (tmpMat*(cjOpts<0)).sum(axis=1)/np.maximum(1,(njMat.T*(cjOpts<0)).sum(axis=1))
                    )
        ### where all samples have same cluster id, set effect size to 0
        sameLocs = np.abs(((njMat.T>0)*(cjOpts)).sum(axis=1)/(njMat.T>0).sum(axis=1))==1
        effectSize[sameLocs]=0

        ### comput sheetCj effect size
        ### preserves sign, can be between -1 and 1
        ### take fresh argmin to find minimizing hash
        tmpMat = sampSumMat[pv[1:,0,:].argmin(axis=0)+1,:,list(range(A))] * sheetCj ### A x p
        effectSize_sheet = ((tmpMat*(sheetCj>0)).sum(axis=1)/np.maximum(1,(njMat.T*(sheetCj>0)).sum(axis=1))
                    + (tmpMat*(sheetCj<0)).sum(axis=1)/np.maximum(1,(njMat.T*(sheetCj<0)).sum(axis=1)))
        ### where all samples have same cluster id, set effect size to 0
        sameLocs = np.abs(((njMat.T>0)*(sheetCj)).sum(axis=1)/(njMat.T>0).sum(axis=1))==1
        effectSize_sheet[sameLocs]=0

        pval_random[idx_start:idx_end] = pv_hash
        pval_samplesheet[idx_start:idx_end] = pv_hash_sheetCjs
        effect_size_random[idx_start:idx_end] = effectSize
        effect_size_samplesheet[idx_start:idx_end] = effectSize_sheet
        optHash[idx_start:idx_end] = hashOpts
        fullMarr[idx_start:idx_end] = Marr
        cjOptMat[idx_start:idx_end] = cjOpts




    print('finished with p-value computation, {:.1F} sec'.format(time.time()-startTime))

    
    outdf = pd.DataFrame({'anchor':df_pivoted_full.index.to_list(),'pval_random':pval_random, 'pval_samplesheet':pval_samplesheet,
             'effect_size_random':effect_size_random,'effect_size_samplesheet':effect_size_samplesheet,'optHash':optHash, 'num_observations':fullMarr})

    entDf = (df.groupby(['anchor','target']).counts.sum()
             .groupby('anchor').apply(lambda x : scipy.stats.entropy(x,base=2))
             .reset_index().rename({'counts':'target_entropy'},axis=1))

    outdf = outdf.merge(entDf)
    outdf = outdf.merge(df[['anchor','mean_target_hamming_distance','mean_target_levenshtein_distance','number_nonzero_samples']].drop_duplicates())

    outdf= outdf.join(pd.DataFrame(cjOptMat,columns=['cj_rand_opt_'+samp for samp in sampleNames])) ## add in cj rand opt

    if useSheetCjs:
        outdf['pval_aggregated'] = 2*np.minimum(outdf.pval_random,outdf.pval_samplesheet)
        outdf.sort_values('pval_aggregated',inplace=True)
    else:
        outdf = outdf.drop(columns=['pval_samplesheet','effect_size_samplesheet'])
        outdf.sort_values('pv_hash',inplace=True)

    print('writing')
    outdf.to_csv(args.outfile_scores, sep='\t', index=False)
    print("done")

main()
