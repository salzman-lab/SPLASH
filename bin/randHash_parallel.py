#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys,itertools
import mmh3
import argparse
from scipy import stats


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


def main():
    args = get_args()

    K=20    ### number of random hashes to use (actually, use K-1, 0 is for handcrafted)
    L = 10  #### number of random cj to use
    dMax = 5 #### can be modified in input

    df = pd.read_csv(args.infile, delim_whitespace=True, names=['counts','seq','sample'])
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

    ### pivot the datafame to compute aggregate statistics
    dfNew = df[['anchor','sample','nj']+['sj_'+str(k) for k in range(K)] + ['sampSum_'+str(k) for k in range(K)]]
    dfNew = dfNew.drop_duplicates()
    #### pass in sampSum to help calculate effect size
    #### note that columns are reordered alphabetically; watch out
    pvtdf = pd.pivot_table(dfNew,index=['anchor','sample'], values=['nj']+['sj_'+str(k) for k in range(K)] + ['sampSum_'+str(k) for k in range(K)])


    mat_iter = pvtdf.groupby('anchor').apply(pd.DataFrame.to_numpy) ### convert all matrices (for each anchor) to numpy

    #### only store 1 value per anchor
    pvHandcrafted = np.zeros(len(mat_iter))
    pvHash = np.zeros(len(mat_iter))
    eSizeArr = np.zeros(len(mat_iter))

    print('computing p values')
    ### iterate over all anchors
    for i,x in enumerate(mat_iter):
        njs = x[:,0] ## unpack data stored in x
        sampSumMat = x[:,1:K+1] ### p x K matrix
        sjMat = x[:,K+1:] ### p x K matrix
        M=njs.sum()
        p = x.shape[0]
        cjs = np.random.choice([-1,1],size=(L,p)) ### number of samples p

        ### compute our test statistics
        testStat = np.abs(cjs @ sjMat)
        b= np.sqrt(M*(cjs**2).sum(axis=1)/ (cjs*np.sqrt(njs)).sum(axis=1)**2)
        a = 1.0/(b+1)
        a[np.isnan(a)] = 0 #### cancel out warning, in case denom of b is 0

        ### compute two terms of p value
        term1HC = 2*np.exp(- 2*(1-a)[:,None]**2*testStat**2/ (cjs**2).sum(axis=1)[:,None] )
        term2HC = 2*np.exp( - 2*a[:,None]**2*M*testStat**2 / (cjs*np.sqrt(njs)).sum(axis=1)[:,None]**2 )
        term2HC[a==0]=0 #### cancel out warning
        pv = term1HC+term2HC  ##### L x K matrix

        pvHandcrafted[i] = pv[:,0].min()
        pvHash[i] = pv[:,1:].min()

        ### idx1 range to L, idx2 ranges to K
        idx1,idx2 = np.unravel_index(np.argmin(pv),pv.shape)
        cjOpt = cjs[idx1]
        sampSumOpt = sampSumMat[:,idx2]
        posSet = cjOpt>0

        ### compute mean of + cluster, mean of - cluster, and take the difference (+ because negative cj have a - in front)
        eSizeArr[i] = np.abs(
            cjOpt[posSet]@ sampSumOpt[posSet]/max(1,njs[posSet].sum()) + cjOpt[~posSet]@ sampSumOpt[~posSet]/max(1,njs[~posSet].sum()) )

    ### construct dfres, dataframe of corrected p values
    anchs = pvtdf.index.get_level_values(0).unique().to_list() ### get list of anchors names
    dct = {'anchor':anchs}
    dct['pv_hand']=np.minimum(1,L*pvHandcrafted)
    dct['pv_hash']=np.minimum(1,L*(K-1)*pvHash)
    dct['effectSize'] = eSizeArr
    dfres = pd.DataFrame(dct)


    dfall = pd.DataFrame()

    ### should never reach this, due to earlier catch
    if not df.empty and not dfres.empty:
        dfall = pd.merge(df,dfres)

    ### add in extra columns to output df
    dfall["mu_hand"] = dfall["mu_0"] ### mean of handcrafted
    ### entropy of target vec averaged across samples
    entDf = (df.groupby(['anchor','target']).counts.sum()
         .groupby('anchor').apply(lambda x : stats.entropy(x,base=2))
         .reset_index().rename({'counts':'entropy'},axis=1))
    dfall = dfall.merge(entDf)


    print('writing')

    (dfall[['anchor','pv_hand','pv_hash','effectSize','mu_hand','entropy']]
     .drop_duplicates().sort_values('pv_hash')
     .to_csv(args.outfile, sep='\t', index=False))
    dfall.to_csv(args.outfile+'extra_info.csv', sep='\t', index=False)



main()
