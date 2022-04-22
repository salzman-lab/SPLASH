#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--use_std",
        type=str
    )
    parser.add_argument(
        "--anchor_sample_scores",
        type=str
    )
    parser.add_argument(
        "--anchor_target_counts",
        type=str
    )
    parser.add_argument(
        "--samplesheet",
        type=str
    )
    parser.add_argument(
        "--kmer_size",
        type=int
    )
    parser.add_argument(
        "--outfile_norm_scores",
        type=str
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # read in files
    counts = pd.read_csv(args.anchor_target_counts, sep='\t')
    anchor_sample_scores = pd.read_csv(args.anchor_sample_scores, sep='\t')

    # read in samples from samplesheet
    # [fastq_file, optional group_id]
    sample_list = pd.read_csv(
        args.samplesheet,
        header=None,
        sep='\t'
    )

    # redefine bool
    if args.use_std == "true":
        use_std = True
    elif args.use_std == "false":
        use_std = False
    if sample_list.shape[1] == 1:
        use_std = True

    # get list of samples from fastq_files
    # if fastq_file = "file1.fastq.gz", sample = "file1"
    samples = (
        sample_list
        .iloc[:,0]
        .apply(
            lambda x:
            os.path.basename(x).split('.')[0]
        )
        .tolist()
    )

    # make group_ids dict, default to 1 if no group_ids provided
    group_ids_dict = {}
    if sample_list.shape[1] == 1:
        for i in range(0, len(samples)):
            group_ids_dict[samples[i]] = 1
    else:
        group_ids = sample_list.iloc[:,1].tolist()
        for i in range(0, len(samples)):
            group_ids_dict[samples[i]] = group_ids[i]

    # define terms
    n_j = (
        counts                          # file of anchor-target-sample counts
        .drop('distance', axis=1)       # df of anchor-targets x sample counts
        .groupby('anchor')              # over all anchors
        .sum()                          # get the per-sample sums
    )

    sqrt_n_j = (
        n_j                             # per-sample sums
        .apply(np.sqrt)                 # square-root the per anchor-sample sums
    )

    anchor_sample_scores = (
        anchor_sample_scores                    # file of anchor_sample_scores per anchor and per anchor-sample
        .set_index('anchor')                    # df of anchor x S_j
        .drop('anchor_sample_std', axis=1)      # drop anchor_score column
    )

    if use_std:
        norm_summary_score = (
            anchor_sample_scores                # df of anchor x S_j
            .mul(pd.Series(group_ids_dict))     # multiply by group_ids
            .mul(sqrt_n_j, axis=1)              # multiply by n_j (square-root of per anchor-sample sums)
            .std(axis=1)                        # std over samples, to get a per-anchor score
        )
    else:
        norm_summary_score = (
            anchor_sample_scores                # df of anchor x S_j
            .mul(pd.Series(group_ids_dict))     # multiply by group_ids
            .mul(sqrt_n_j, axis=1)              # multiply by n_j (square-root of per anchor-sample sums)
            .sum(axis=1)                        # sum over samples, to get a per-anchor score
        )

    sum_sqrt_n_j_c_j = (
        sqrt_n_j                                # square-root of per anchor-sample sums
        .mul(pd.Series(group_ids_dict))         # multiply by group_ids
        .sum(axis=1)                            # sum over samples, to get a per-anchor score
    )

    # intialise the scores table
    scores_table = n_j.copy()
    scores_table['total_anchor_counts'] = scores_table.sum(axis=1)
    scores_table = (
        pd.merge(
            scores_table,
            pd.DataFrame(norm_summary_score, columns=['norm_summary_score']),
            on='anchor'
        )
    )

    # initialise columns
    scores_table['first_moment'] = None
    scores_table['second_moment'] = None
    scores_table['third_moment'] = None
    scores_table['fourth_moment'] = None
    scores_table['fourth_central_moment'] = None
    scores_table['expectation'] = None
    scores_table['variance_distance'] = None

    #### Tavor added
    # scores_table['p-value_cj_v1'] = None
    # scores_table['p-value_cj_vGamma'] = None
    # scores_table['p-value_cj_H'] = None
    # scores_table['p-value_cj_Hv2'] = None
    scores_table['l1_norm_Sj'] = None
    scores_table['l2_norm_Sj'] = None
    scores_table['p-value_cj'] = None
    scores_table['p-value_noCj'] = None
    test_stat_cj = ( ### regenerate in case use_std is True
            anchor_sample_scores                # df of anchor x S_j
            .mul(pd.Series(group_ids_dict))     # multiply by group_ids
            .mul(sqrt_n_j, axis=1)              # multiply by n_j (square-root of per anchor-sample sums)
            .sum(axis=1)                        # sum over samples, to get a per-anchor score
        )

    #### due to sqrt(n_j) factor off;
    sj_vals = (
        anchor_sample_scores                # df of anchor x S_j
                .mul(sqrt_n_j, axis=1)              # multiply by n_j (square-root of per anchor-sample sums)
    )

    np.random.seed(0) #### for random choices of cj, fix randomness

    # for each anchor, get expectation and variance_distance
    for anchor, df in counts.groupby('anchor'):

        # subset and rename distance as i
        distances = (
            df
            .drop(['anchor', 'target'], axis=1)
            .rename(columns={'distance' : 'i'})
        )

        # initialise
        p = pd.DataFrame()
        p['i'] = distances['i'].unique()
        p['counts'] = None

        for i, df in distances.groupby('i'):
            # get total counts for targets with that distance
            count = (
                df
                .drop('i', axis=1)
                .values
                .sum()
            )
            # set the count
            p.loc[p['i']==i, 'counts'] = count

        # get fraction of times each distance occurs as p_hat
        p['p_hat'] = p['counts']/p['counts'].sum()

        # define moments
        first_moment = (p['p_hat'] * p['i']).sum()
        second_moment = (p['p_hat'] * p['i']**2).sum()
        third_moment = (p['p_hat'] * p['i']**3).sum()
        fourth_moment = (p['p_hat'] * p['i']**4).sum()

        mu = (p['i'] * p['p_hat']).sum()
        fourth_central_moment = (p['i'] * (p['i'] - mu)**4).sum()

        # define variance_distance
        variance_distance = (p['p_hat'] * p['i']**2).sum() - ((p['p_hat'] * p['i']).sum()) ** 2

        # add these values to the score table
        scores_table.loc[anchor, 'variance_distance'] = variance_distance
        scores_table.loc[anchor, 'first_moment'] = first_moment
        scores_table.loc[anchor, 'second_moment'] = second_moment
        scores_table.loc[anchor, 'third_moment'] = third_moment
        scores_table.loc[anchor, 'fourth_moment'] = fourth_moment
        scores_table.loc[anchor, 'fourth_central_moment'] = fourth_central_moment



        ### Tavor p-value addition
        varEst = second_moment-first_moment**2 ### ignoring normalizing issue
        cjs = pd.Series(group_ids_dict).to_numpy()
        njs = n_j.loc[anchor].to_numpy() 
        # sjs = anchor_sample_scores.loc[anchor].to_numpy()[:-1] ### drop the "num_anchor_sample_scores" col, not normalized by sqrt(nj)
        # s_val = norm_summary_score.loc[anchor] ### value of test statistic
        s_val = test_stat_cj.loc[anchor]

        M = njs.sum() ### total number of observations for this anchor
        D = 8 #### maximum distance value,  ###### todo add as arg
        eps = np.finfo(np.float32).eps ### to prevent division by 0
        nj_invSqrt = 1.0/np.sqrt(njs+eps)
        nj_invSqrt[njs==0]=0 ## set to 0


        ### Bernstein-based concentration, p-value_cj_v1
        term1 = 2*np.exp(-s_val**2 / 
            (16*varEst*(cjs**2).sum() + 8/3*s_val*D*(np.abs(cjs)*nj_invSqrt).max() ))
        term2 = 2*np.exp(-M**2*s_val**2 / 
            (16*M*varEst*(cjs*np.sqrt(njs)).sum()**2 + 4/3*s_val*D* np.abs( (cjs*np.sqrt(njs)).sum() )))
        term3 = np.exp(-(M-1)/(4*D**4))
        pvalv1 = min(term1+term2+term3,1)
        # scores_table.loc[anchor, 'p-value_cj_v1'] = pvalv1


        ### Bernsteion, coarsely optimizing over gamma, pvalvGamma
        a = 8*M*(cjs**2).sum()
        b = 8*M*s_val*D/3*(np.abs(cjs)*nj_invSqrt).max()
        c = -2*D**4*s_val**2
        gamma=(-b+np.sqrt(b**2-4*a*c))/(2*a)

        term1g = 2*np.exp(-s_val**2 / 
            (8*(varEst+gamma)*(cjs**2).sum() + 8/3*s_val*D*(np.abs(cjs)*nj_invSqrt).max() ))
        term2g = 2*np.exp(-M**2*s_val**2 / 
            (8*M*(varEst+gamma)*(cjs*np.sqrt(njs)).sum()**2 + 4/3*s_val*D* np.abs( (cjs*np.sqrt(njs)).sum() )))
        term3g = np.exp(-(M-1)*gamma**2/(2*D**4*(varEst+gamma)))
        pvalvGamma = min(term1g+term2g+term3g,1)
        # scores_table.loc[anchor, 'p-value_cj_vGamma'] = pvalvGamma


        ### just using Hoeffdings directly with boundedness, p-value_cj_H
        term1v2 = 2*np.exp(-s_val**2 / 
            (8*D**2*(cjs**2).sum()))
        term2v2 = 2*np.exp(-M*s_val**2 / 
            (8*D**2*(cjs*np.sqrt(njs)).sum()**2 ))
        pvalH=min(term1v2+term2v2,1)
        # scores_table.loc[anchor, 'p-value_cj_H'] = pvalH


        ### optimized Hoeffding's bound, p-value_cj_Hv2
        b= np.sqrt(M*(cjs**2).sum()/ (cjs*np.sqrt(njs)).sum()**2)
        a = 1.0/(b+1)

        term1Ha = 2*np.exp(- 2*(1-a)**2*s_val**2/ (D**2 * (cjs**2).sum()))
        term2Ha = 2*np.exp( - 2*a**2*M*s_val**2 / (D**2 *(cjs*np.sqrt(njs)).sum()**2 ))
        pvalHv2 = min(term1Ha + term2Ha,1)
        # scores_table.loc[anchor, 'p-value_cj_Hv2'] = pvalHv2


        ### set overall p-value to be minimum of all candidate pvalues
        scores_table.loc[anchor, 'p-value_cj']=min([pvalv1,pvalH,pvalvGamma,pvalHv2])


        #### no c_j
        numCj = 5
        cjArr = np.zeros((numCj,len(njs)))
        for l in range(numCj):
            cjArr[l] = np.random.choice([-1,1],size=len(njs))

        sjs = sj_vals.loc[anchor].to_numpy()[:-1] ### pull sj values
        nanLocs = np.isnan(sjs) ### remove nans, zero out both
        sjs[nanLocs]=0
        cjArr[:,nanLocs] = 0

        testEvals = cjArr@sjs #### look at the numCj linear statistics
        testStat = np.abs(testEvals).max()

        pvalArr = np.zeros(numCj)
        for l in range(numCj):
            cjs = cjArr[l]

            #### use simple Hoeffding's v2, union bound
            b= np.sqrt(M*(cjs**2).sum()/ (cjs*np.sqrt(njs)).sum()**2)
            a = 1.0/(b+1)

            term1HC = 2*np.exp(- 2*(1-a)**2*testStat**2/ (D**2 * (cjs**2).sum()))
            term2HC = 2*np.exp( - 2*a**2*M*testStat**2 / (D**2 *(cjs*np.sqrt(njs)).sum()**2 ))

            pvalArr[l] = min(term1HC + term2HC,1)

        pvalHC = min(pvalArr.sum(),1)
        scores_table.loc[anchor, 'p-value_noCj'] = pvalHC

        #### adding in l1 and l2 (centered)
        sjs = sjs[sjs!=0]
        scores_table.loc[anchor, 'l1_norm_Sj'] = np.linalg.norm(sjs,ord=1)
        scores_table.loc[anchor, 'l2_norm_Sj'] = np.linalg.norm(sjs,ord=2)


    # output scores table
    (
        scores_table
        .reset_index()
        .to_csv(args.outfile_norm_scores, sep='\t', index=False)
    )


main()
