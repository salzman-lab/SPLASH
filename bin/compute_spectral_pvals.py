import numpy as np
import pandas as pd
from tqdm import tqdm,tqdm_notebook, tqdm_pandas
import os
import glob
import pickle
import argparse
import scipy
import scipy.stats
import statsmodels.api as sm
from pathlib import Path 

from stats_utils import *

##### Additional details / theory will appear in an upcoming submission


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
    parser.add_argument( ### samplesheet file, if samplesheet-based cj (metadata) are to be used
        "--samplesheet",
        type=str,
        default=""
    )
    parser.add_argument( ### columns to output. Options are default, metadata, full
        "--output_verbosity",
        type=str,
        default="default"
    )
    parser.add_argument( ### optional file indicating subset of anchors to be used
        "--anchor_list",
        type=str,
        default=""
    )
    parser.add_argument( ### optional file indicating subset of anchors to be used
        "--kmer_size",
        type=int,
        default=27
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
    ############### This option is currently not supported
    parser.add_argument( ### flag, whether to save cj or not
        "--save_c_f",
        type=bool,
        default=False
        ### if save_c_f flag, then can read in optimizing c and f as below:
        #### with open(args.outfile_scores+'/spectral_cj.npy','rb') as f:
        ####     a = np.load(f)
        #### with open(args.outfile_scores+'/spectral_f.pkl', 'rb') as handle:
        ####     b = pickle.load(handle)
    )
    parser.add_argument( ### fraction of data partitioned for training
        "--train_fraction",
        type=float,
        default=.25
    )
    parser.add_argument( ### fraction of data partitioned for training
        "--num_rand_cf",
        type=int,
        default=50
    )
    args = parser.parse_args()
    return args




#### main function
def main():
    args = get_args()

    if args.output_verbosity not in ['default','metadata','full']:
        print('invalid option for output_verbosity')
        return

    if args.save_c_f:
        print('save_c_f not supported')
        return

    ### read in anchor list file
    if len(args.anchor_list)>0:
        print("using passed in anchor list")
        anchLst = pd.read_csv(args.anchor_list,names=['anchor']).anchor.to_list()
        print(len(anchLst), "anchors")
    else:
        print('using all anchors')
        anchLst = []

    print('constructing counts dataframe')
    countsDf = constructCountsDf(args,anchLst)

    if anchLst == []:
        anchLst = countsDf.anchor.unique()
        print('generated all anchors, ', len(anchLst))

    print("parsing samplesheet")
    useSheetCj, samplesheetDf = parseSamplesheet(args.samplesheet)


    anchLst = set(anchLst)
    nuniqueAnchors = countsDf.anchor.nunique()
    anchsUsed = np.ones(nuniqueAnchors,dtype='bool')
    resultsDf = pd.DataFrame()

    print("Starting loop over anchors")
    for anch_idx,(anch,anch_table) in tqdm(enumerate(countsDf.groupby('anchor')), total = nuniqueAnchors):
        if anch not in anchLst:
            anchsUsed[anch_idx]=False
            continue

        ### Row to be added to dataframe
        newRow = {'anchor':anch}

        ### Get the relevant data from the table, pivot it
        anch_pivot_table = (anch_table.drop(columns='anchor')
                            .pivot(index=['target'], columns='sample', values='counts')
                            .fillna(0)) ### index is targets, for levenshtein / etc computation
        
        ### this is the contingency table to operate on
        anch_contingency_table = anch_pivot_table.to_numpy()


        ### compute asymptotically valid comparison tests
        if args.output_verbosity=="experimental":
            newRow['pval_chi2'] = computeChi2Test(anch_contingency_table)
            newRow['pval_lrt'] = computeLRT_Test(anch_contingency_table)


        #### split data into train and test portions
        np.random.seed(0) #### to make it deterministic 
        X = splitCountsColwise(anch_contingency_table,args.train_fraction)
        Xtrain = X
        Xtest = anch_contingency_table-Xtrain


        ### compute simple c,f from spectral approach (correspondence analysis style)
        ###   and compute nomad_simpleSVD_pv
        cOpt,fOpt = get_spectral_cf_svd(Xtrain,anch_contingency_table.shape)
        newRow['pval_SVD_corrAnalysis'] = testPval(Xtest,cOpt,fOpt)
        newRow['pval_asymp_SVD_corrAnalysis'] = computeAsympNOMAD(Xtest,cOpt,fOpt)
        newRow['effect_size_cts_SVD'] = effectSize_cts(Xtest,cOpt,fOpt)


        ### compute pvalsRandOpt
        cOpt,fOpt,_ = generateRandOptcf(Xtrain,Xtest.shape)
        newRow['pval_rand_init_EM']=testPval(Xtest,cOpt,fOpt)    


        ### compute nomad's base p-value
        nomadpvminarr = np.zeros(args.num_rand_cf)
        nomadasympArr = np.zeros(args.num_rand_cf)
        randCs = np.random.choice([-1,1],size=(args.num_rand_cf,len(cOpt)))
        randFs = np.random.choice([0,1],size=(args.num_rand_cf, len(fOpt)))
        for k in range(args.num_rand_cf):
            nomadpvminarr[k] = testPval(anch_contingency_table,randCs[k], randFs[k])
            nomadasympArr[k] = computeAsympNOMAD(anch_contingency_table,cOpt,fOpt)
        newRow['pval_base'] = min(1,args.num_rand_cf*nomadpvminarr.min())
        newRow['pval_asymp_base'] = min(1,args.num_rand_cf*nomadasympArr.min())

        ### compute effect size for 
        minimizerIdx = nomadpvminarr.argmin()
        newRow['effect_size_base'] = effectSize_bin(anch_contingency_table,randCs[minimizerIdx],randFs[minimizerIdx])


        if args.output_verbosity=='experimental':
            ### compute for continuous c,f
            nomadasympArr = np.zeros(args.num_rand_cf)
            nomadctsArr = np.zeros(args.num_rand_cf)
            randCs = np.random.uniform(low=-1,high=1,size=(args.num_rand_cf,len(cOpt)))
            randFs = np.random.uniform(size=(args.num_rand_cf, len(fOpt)))
            for k in range(args.num_rand_cf):
                nomadctsArr[k] = testPval(anch_contingency_table,randCs[k], randFs[k])
                nomadasympArr[k] = computeAsympNOMAD(anch_contingency_table,randCs[k],randFs[k])
            
            newRow['pval_asymp_cts_base'] = min(1,args.num_rand_cf*nomadasympArr.min())
            newRow['pval_cts_base'] = min(1,args.num_rand_cf*nomadctsArr.min())

            ### compute pvalsSpectral
            cOpt,fOpt,_ = generateSpectralOptcf(Xtrain,Xtest.shape)
            newRow['pval_spectral_EM']=testPval(Xtest,cOpt,fOpt)


        ##### hasn't been thoroughly tested, but seems to be working as expected
        if args.output_verbosity == 'metadata' or (useSheetCj and args.output_verbosity=='experimental'): ### not fully tested, use with caution
            sheetCj = samplesheetDf[anch_pivot_table.columns].to_numpy().flatten()

            cOpt,fOpt,_ = generateSignedSheetCjOptcf(Xtrain,sheetCj,Xtest.shape)
            newRow['pval_metadata_EM']=testPval(Xtest,cOpt,fOpt)
            
            cOpt,fOpt = generateSheetCjOptcf(Xtrain,sheetCj,Xtest.shape)
            newRow['pval_metadata_optF']=testPval(Xtest,cOpt,fOpt)

            nomadasympArr = np.zeros(args.num_rand_cf)
            nomadpvArr = np.zeros(args.num_rand_cf)
            randFs = np.random.choice([0,1], size=(args.num_rand_cf, len(fOpt)))
            for k in range(args.num_rand_cf):
                nomadpvArr[k] = testPval(anch_contingency_table,sheetCj, randFs[k])
                nomadasympArr[k] = computeAsympNOMAD(anch_contingency_table,sheetCj,randFs[k])

            newRow['pval_metadata_asymp_base'] = min(1,args.num_rand_cf*nomadasympArr.min())
            newRow['pval_metadata_base'] = min(1,args.num_rand_cf*nomadpvArr.min())

             ### compute effect size for base nomad with sheetCj
            minimizerIdx = nomadpvminarr.argmin()
            newRow['effect_size_metadata_base'] = effectSize_bin(anch_contingency_table,sheetCj,randFs[minimizerIdx])


        ### compute additional quantities (e.g. M, number of unique targets, etc)
        rowMetadata = computeBaseQuantities(anch_contingency_table)
        newRow = newRow | rowMetadata
        newRow['mean_target_levenshtein_distance'] = computeAverageDist(anch_pivot_table,nltk.edit_distance)
        newRow['mean_target_hamming_distance'] = computeAverageDist(anch_pivot_table,hamming)

        resultsDf = resultsDf.append(newRow,ignore_index=True)

    outdf = resultsDf

    if args.output_verbosity != 'experimental':
        outdf.drop(columns=outdf.columns[outdf.columns.str.contains('asymp')],inplace=True)        

    outdf = outdf.sort_values('pval_SVD_corrAnalysis')

    filepath = Path(args.outfile_scores)  
    filepath.parent.mkdir(parents=True, exist_ok=True)

    outdf.to_csv(filepath, sep='\t', index=False)



    # if args.save_c_f:
    #     if not useSheetCj:
    #         cjArr = cjArr[:,:4]
    #     cjArr = cjArr[anchsUsed]
    #     with open(args.outfile_scores[:-4]+'_spectral_cj.npy', 'wb') as f:
    #         np.save(f,cjArr)

    #     with open(args.outfile_scores[:-4]+'_spectral_f.pkl', 'wb') as handle:
    #         pickle.dump(fArr, handle, protocol=pickle.HIGHEST_PROTOCOL)

        #### to be read in as below
        # with open(args.outfile_scores[:-4]+'_spectral_cj.npy','rb') as f:
        #     a = np.load(f)
        # with open(args.outfile_scores[:-4]+'_spectral_f.pkl', 'rb') as handle:
        #     b = pickle.load(handle)


print('starting spectral p value computation')
main()