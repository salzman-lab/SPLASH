#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys,itertools
import mmh3
import argparse
from scipy import stats
import pickle



def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--infile",
        type=str
    )
    parser.add_argument(
        "--kmer_size",
        type=int
    )
    parser.add_argument(
        "--outfile",
        type=str
    )
    parser.add_argument( ### takes as input filepath to samplesheet
        "--samplesheetIDs",
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


### input: dft, dataframe resuling from df.groupby('anchor') and iterating over anchors
### outputs: updated dataframe dft, (vector of optimal cjs used, handcjs, index of hash function used)
def computePVals(dft,L=10,K=20,useHandCjs=False):

    njs = dft.nj.to_numpy()
    sampSumMat = dft[['sampSum_'+str(k) for k in range(K)]].to_numpy() ### p x K matrix
    sjMat = dft[['sj_'+str(k) for k in range(K)]].to_numpy() ### p x K matrix
    M=njs.sum()
    p = len(njs)
    cjs = np.random.choice([-1,1],size=(L+1,p)) ### number of samples p, last row is for hand cj
    if useHandCjs:
        cjs[-1] = dft.handCjs.to_numpy()

    ### compute our test statistics
    testStat = np.abs(cjs @ sjMat)
    b= np.sqrt(M*(cjs**2).sum(axis=1)/ (cjs*np.sqrt(njs)).sum(axis=1)**2)
    a = 1.0/(b+1)
    a[np.isnan(a)] = 0 #### cancel out warning, in case denom of b is 0

    ### compute two terms of p value
    term1HC = 2*np.exp(- 2*(1-a)[:,None]**2*testStat**2/ (cjs**2).sum(axis=1)[:,None] )
    term2HC = 2*np.exp( - 2*a[:,None]**2*M*testStat**2 / (cjs*np.sqrt(njs)).sum(axis=1)[:,None]**2 )
    term2HC[a==0]=0 #### cancel out warning
    pv = term1HC+term2HC  ##### L+1 x K matrix

    dft['pv_hand'] = min(L*pv[:-1,0].min(),1)
    dft['pv_hash'] = min(L*(K-1)*pv[:-1,1:].min(),1)
#     print(useHandCjs)
    if useHandCjs:
        dft['pv_hand_sheetCjs'] = min(1,L*pv[-1,0])
        dft['pv_hash_sheetCjs'] = min(1,L*pv[-1,1:].min())

    ### idx1 range to L, idx2 ranges to K
    idx1,idx2 = np.unravel_index(np.argmin(pv[:-1,1:]),pv[:-1,1:].shape)
    idx2+=1 ### idx2 should range from 1 to K, not 0 to K-1 (since handcrafted is 0)
    cjOpt = cjs[idx1]
    sampSumOpt = sampSumMat[:,idx2]
    posSet = cjOpt>0

    dft['hash_max'] = idx2

    ### compute mean of + cluster, mean of - cluster, and take the difference (+ because negative cj have a - in front)
    dft['effect_size'] = np.abs(
        cjOpt[posSet]@ sampSumOpt[posSet]/max(1,njs[posSet].sum()) + cjOpt[~posSet]@ sampSumOpt[~posSet]/max(1,njs[~posSet].sum()))
    return dft, (cjOpt,cjs[-1], idx2)


def hamming(x,y): ### works in conjunction with dMax
    return sum(c1 != c2 for c1, c2 in zip(x,y))
def normalizevec(x): ### shift and scale vector so that it's in [-1,1]
    return 2*(x-x.min()) / (x.max()-x.min())-1

def main():
    args = get_args()

    ### load samplesheet c_j's
    useHandCjs = False
    samplesheet = args.samplesheetIDs
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


    K=20    ### number of random hashes to use (actually, use K-1, 0 is for handcrafted)
    L = 20  #### number of random cj to use
    dMax = 5 #### can be modified in input
    #### can also change hamming to levenshtein / other

    df = pd.read_csv(args.infile, delim_whitespace=True, names=['counts','seq','sample'])
    if useHandCjs:
        df = pd.merge(df,sheetdf)
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
    mergedDf['dij_0'] = np.vectorize(lambda x,y : min(hamming(x,y),dMax)/dMax)(mergedDf.target,mergedDf.maxTarget)
    df = mergedDf.drop(columns=['targ_cts','maxTarget'])
    df = convertDijToSj(df,0)


    #### hash based dij, randomly assign each target to 0 or 1
    for k in range(1,K):
        df['dij_'+str(k)] = (df['target'].apply(lambda x : (mmh3.hash(x,seed=k)>0) ))*1.0

        df = convertDijToSj(df,k)


    ### drop columns to remove duplicates (primarily dij,target, and counts), to allow for single row per anchor sample pair
    anchorDF = (df.loc[:,((~df.columns.str.startswith('dij_'))
                                 & (~df.columns.str.startswith('sum_'))
#                                  & (~df.columns.str.startswith('mu_'))
                                )]
                       .drop(columns=['anch_cts','anch_samples','anch_uniqTargs','counts','target'])).drop_duplicates()


    np.random.seed(0) ### for reproducibility, random c_j's
    arr = []
    cjDict = {}
    ### currently just doing a for loop over anchors, since we need cj output
    ###   (can't easily do .apply() )
    for (anch,dft) in anchorDF.groupby('anchor'):
        tmp,extraInfo = computePVals(dft,L,K,useHandCjs)
        arr.append(tmp)
        cjDict[anch] = extraInfo

    dfall = pd.concat(arr)


    ### add in extra columns to output df
    dfall["mu_hand"] = dfall["mu_0"] ### mean of handcrafted
    ### entropy of target vec averaged across samples
    entDf = (df.groupby(['anchor','target']).counts.sum()
         .groupby('anchor').apply(lambda x : stats.entropy(x,base=2))
         .reset_index().rename({'counts':'entropy'},axis=1))
    dfall = dfall.merge(entDf)

    #### construct pvals df with anchor + new statistics
    dfpvals = dfall[['anchor']+dfall.columns[-8:].to_list()].drop_duplicates()


    print('writing')

    (dfpvals
     .sort_values('pv_hash')
     .to_csv(args.outfile, sep='\t', index=False))
    # dfall.to_csv(args.outfile+'extra_info.csv', sep='\t', index=False)

    # with open(args.outfile+'extra_anchor_dict.pkl','wb') as f:
    #     pickle.dump(cjDict,f)

    print("done")

main()
