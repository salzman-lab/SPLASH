#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys,itertools
import mmh3
import argparse

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

    df[prodcol] = df.counts*df[dijcol]
    df[mucol] = df.groupby(['anchor'])[prodcol].transform('sum')/df.M
    df[sjcol]=(df.groupby(['anchor','sample'])[prodcol].transform('sum') - df.nj*df[mucol]) / (df.nj**(.5))
    return df


def main():
    args = get_args()

    K=20    ### number of random hashes to use
    L = 10  #### number of random cj to use

    df = pd.read_csv(args.infile, delim_whitespace=True, names=['counts','seq','sample'])

    ### split seq into anchor and target
    ### for some reason some files have duplicates, so drop those
    df['anchor'] = df.seq.str[:args.kmer_size]
    df['target'] = df.seq.str[args.kmer_size:]
    df = df.drop(columns='seq').drop_duplicates()

    ### add in informational columns for each anchor for filtering purposes
    ### this is faster than filter command
    df['anch_uniqTargs'] = df.groupby('anchor').target.transform('nunique')
    df['anch_cts'] = df.groupby('anchor').counts.transform('sum')
    df['anch_samples']= df.groupby('anchor')['sample'].transform('nunique')
    df['anchSample_cts'] = df.groupby(['anchor','sample']).counts.transform('sum')
    df = df[(df.anch_uniqTargs>1) & (df.anch_cts > 10) & (df.anch_samples > 1) & (df.anchSample_cts>4)]
    ### above line filters, can change parameters

    df['nj'] = df.groupby(['anchor','sample']).counts.transform('sum')
    df = df.drop(columns='anchSample_cts') ### this is essentially "old" n_j
    df['M'] = df.groupby(['anchor']).counts.transform('sum')

    #### hash based dij, randomly assign each target to 0 or 1
    for k in range(K):
        df['dij_'+str(k)] = (df['target'].apply(lambda x : (mmh3.hash(x,seed=k)>0) ))*1.0

        df = convertDijToSj(df,k)

    ### pivot the datafame to compute aggregate statistics
    dfNew = df[['anchor','sample','nj','M']+['sj_'+str(k) for k in range(K)]]
    dfNew = dfNew.drop_duplicates()

    # pvtdf = pd.pivot_table(dfNew,index=['anchor'], columns=['sample'],values=['nj']+['sj'+str(k) for k in range(K)])
    pvtdf = pd.pivot_table(dfNew,index=['anchor','sample'], values=['nj']+['sj_'+str(k) for k in range(K)])

    anchs = pvtdf.index.get_level_values(0).unique().to_list()
    mat_iter = pvtdf.groupby('anchor').apply(pd.DataFrame.to_numpy)

    outMat = np.zeros((len(mat_iter),K))
    for i,x in enumerate(mat_iter):
        njs = x[:,0]
        M=njs.sum()
        p = x.shape[0]
        cjs = np.random.choice([-1,1],size=(L,p))

        testStat = np.abs(cjs @ x[:,1:])
        b= np.sqrt(M*(cjs**2).sum(axis=1)/ (cjs*np.sqrt(njs)).sum(axis=1)**2)
        a = 1.0/(b+1)
        a[np.isnan(a)] = 0 #### cancel out warning

        term1HC = 2*np.exp(- 2*(1-a)[:,None]**2*testStat**2/ (cjs**2).sum(axis=1)[:,None] )
        term2HC = 2*np.exp( - 2*a[:,None]**2*M*testStat**2 / (cjs*np.sqrt(njs)).sum(axis=1)[:,None]**2 )
        term2HC[a==0]=0 #### cancel out warning
        pv = term1HC+term2HC  ##### L x K matrix

        pvals = pv.min(axis=0) ### uncorrected

        outMat[i] = pvals

    ### construct df of corrected p values
    dct = {'anchor':anchs}
    for k in range(K):
        dct['pv'+str(k)]=np.minimum(1, L*K*outMat[:,k])
    dfres = pd.DataFrame(dct)

    dfres['pv']=dfres.drop(columns='anchor').min(axis=1)
    dfres['pv_Rand'] = dfres[['pv'+str(k) for k in range(1,K)]].min(axis=1)

    dfall = pd.DataFrame()

    if not df.empty and not dfres.empty:
        dfall = pd.merge(df,dfres)

    print('writing')

    dfres[['anchor','pv_Rand']].sort_values('pv_Rand').to_csv(args.outfile, sep='\t', index=False)
    # dfres[['anchor','pv_Rand']].sort_values('pv_Rand').anchor.to_csv('{}/anchors_{}.txt'.format(outfldr,slurmID),sep='\t',index=False)


main()
