import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import itertools
import glob
import argparse
import tqdm




def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( ### anchors_pvals.tsv file
        "--anchor_pvals",
        type=str,
        default=""
    )
    parser.add_argument( ### take in additional_summary.tsv file
        "--pvals_all_file",
        type=str,
        default=""
    )
    parser.add_argument( ### folder to fractions.tab files
        "--consensus_fldr",
        type=str,
        default="/oak/stanford/groups/horence/kaitlin/results_nomad/bulk_RNAseq/paper/HLCA_SS2_P2_macrophage_capillary/consensus_anchors/"
    )
    parser.add_argument( ### folder to samplesheet, to get samplesheet cjs
        "--samplesheet",
        type=str,
        default="/oak/stanford/groups/horence/kaitlin/running_nomad/bulk_RNAseq/samplesheets/samplesheet_HLCA_SS2_P2_capillary.csv"
    )
    parser.add_argument( ### folder to save figures to
        "--outfldr",
        type=str,
        default="/oak/stanford/groups/horence/tavorb/HLCA_SS2_P2_capillary/figs/"
    )
    parser.add_argument( ### anchor target counts file
        "--anch_cts_file",
        type=str,
        default="/oak/stanford/groups/horence/kaitlin/results_nomad/bulk_RNAseq/paper/HLCA_SS2_P2_capillary/anchor_top250_targets_counts.tsv"
    )
    parser.add_argument( ### file of 1 column, no header, 1 anchor per line
        "--anch_list",
        type=str,
        default=""
    )
    parser.add_argument( ### if no ctsDf is provided, construct manually
        "--strat_fldr", 
        type=str,
        default="" #/oak/stanford/groups/horence/kaitlin/results_nomad/bulk_RNAseq/paper/HLCA_SS2_P2_capillary/abundant_stratified_anchors/
    )
    
    args = parser.parse_args()
    return args

def normalizevec(x): ### shift and scale vector so that it's in [-1,1]
    return 2*(x-x.min()) / (x.max()-x.min())-1

#### For pretty plotting of targets, to convert basepair to int (color)
def bpToInt(x):
    if x=='A':
        return 0
    if x=='T':
        return 1
    if x=='C':
        return 2
    if x=='G':
        return 3
    return 4 #### for nan


### Inputs:
#### dfpvals: primary output csv file of previous step, currently additional_summary.tsv
#### dfcts: sample x target counts matrix
#### anch: anchor to plot
#### cjSheet: cj values from samplesheet
#### savePath: folder to save output figure
def plotContingency(dfpvals,dfcts,anch,cjSheet,saveFldr,consensusFldr='',plotConsensus=False):
    
    plotBySheet = ~np.all(cjSheet==1)
    
    ### load in list of samples where this anchor is present
    sampleNames = (dfpvals.columns[dfpvals.columns.str.startswith('cj_rand')]
              .str.lstrip('cj_rand_opt_').to_list())
    
    if anch not in dfcts.index:
        print(anch, " not in dfcts")
        return
    
    tmp = dfcts.loc[anch,sampleNames]
    tbl = tmp.fillna(0).to_numpy().T
    
    ### need to aggregate rows (different sample / anchors)
    pvRow = dfpvals[dfpvals.anchor==anch].mode(axis=0).iloc[0].copy()
    
    ### read in cjOpt
    cjOpt = pvRow[dfpvals.columns.str.startswith('cj_rand')].to_numpy().flatten().astype('float')
    
    ### if it's facing the wrong way, flip it
    if cjOpt @ cjSheet < 0 :
        cjOpt = -cjOpt
        
    ### construct some condition to filter out samples
#     relevantSamples = tbl.sum(axis=1)>0 ### one possible condition
    relevantSamples = cjOpt!=0
    tbl=tbl[relevantSamples,:]
    cjOpt=cjOpt[relevantSamples]
    cjSheet=cjSheet[relevantSamples]

    ### order targets by abundance
    targOrdering = np.argsort(-tbl.sum(axis=0))
    tbl = tbl[:,targOrdering]
    thresh = max(-np.sort(-tbl.sum(axis=0))[min(10,len(targOrdering)-1)],5)
    relevantCols = tbl.sum(axis=0)>=thresh
    
    
    ### sort samples first by cj, then by abundance of first target
    if plotBySheet: ###sort samples by cjSheet
        toSort = list(zip(cjSheet,-tbl[:,0]/np.maximum(tbl.sum(axis=1),1)))
    else: ### sort samples by cjRand
        toSort = list(zip(cjOpt,-tbl[:,0]/np.maximum(tbl.sum(axis=1),1))) 

    sampOrdering = [x[1] for x in sorted((e,i) for i,e in enumerate(toSort))]
    tbl=tbl[sampOrdering,:]
    cjOpt=cjOpt[sampOrdering]
    cjSheet=cjSheet[sampOrdering]
    
    if plotBySheet:
        n1 = (cjSheet>0).sum()
        n2 = (cjSheet<0).sum()
    else:
        n1 = (cjOpt>0).sum()
        n2 = (cjOpt<0).sum()
        
    if 'pv_hand_sheetCjs' not in dfpvals.columns:
        pvRow.loc['pv_hand_sheetCjs']=np.nan
        pvRow.loc['pv_hash_sheetCjs']=np.nan
        pvRow.loc['effect_size_sheetCjs']=np.nan
    
    if plotConsensus:
        fig, axs = plt.subplots(nrows=2,ncols=3)
        fig.set_size_inches(15, 10)
    else:
        fig, axs = plt.subplots(nrows=2,ncols=2)
        fig.set_size_inches(10, 10)
    fig.set_facecolor('white')
    
    (fig.suptitle('''Normalized target distribution for each sample for anchor {}
    Local gene {}, consensus gene {}
    Corrected (both rand and sheet c)hash-based p value = {:.2E}
    (rand f, rand c) p value = {:.2E}, (hand f, rand c) p value = {:.2E}
    (rand f, hand c) p value = {:.2E}, (hand f, hand c) p value = {:.2E}
    Cosine similarity between random c and samplesheet c = {:.2F}
    Rand effect size = {:.2F}, SheetCj effect size = {:.2F}
    Mu lev = {:.2F}, Mu ham = {:.2F}, total number of reads used for inference M = {}
    Number of +1 samples = {}, number of -1 samples = {}'''
              .format(anch,
                      pvRow.anchor_local_gene,
                      pvRow.consensus_gene_mode,
                      pvRow['pv_hash_both_corrected'],
                     pvRow['pv_hash'],
                     pvRow['pv_hand'],
                     pvRow['pv_hash_sheetCjs'],
                     pvRow['pv_hand_sheetCjs'],
                     cjOpt@cjSheet/len(cjSheet),
                      pvRow['effect_size_randCjs'],
                     pvRow['effect_size_sheetCjs'],
                      pvRow['mu_lev'],
                      pvRow['mu_ham'],
                     int(pvRow.M),
                     n1,
                     n2),
                 y=1.05))

    fig.subplots_adjust(hspace=.3)


    ### construct target matrix for visual printing
    targs = dfcts.loc[anch].index.values[targOrdering][relevantCols]
    targMat = np.zeros((27,len(targs)))
    for i,targ in enumerate(targs):
        targMat[:,i] = [bpToInt(x) for x in targ]


    ### plot contingency table
    plt.sca(axs[0,1])
    plt.imshow(tbl[:,relevantCols]/tbl.sum(axis=1)[:,None],aspect='auto')
    plt.xticks([])
    plt.yticks(ticks=range(tbl.shape[0]),labels=tbl.sum(axis=1).astype('int'))
    plt.ylabel(r'$n_j$ for sample j')
    plt.xlabel('Target (value is # occurences)')
    plt.xticks(ticks=range(relevantCols.sum()),labels=tbl[:,relevantCols].sum(axis=0).astype('int'),rotation=-80)
    plt.colorbar()

    ### plot target matrix
    plt.sca(axs[1,1])
    plt.imshow(targMat,aspect='auto')
    plt.xticks(ticks=range(len(targs)),labels=targs,rotation=-80, ha="left",rotation_mode="anchor")
    plt.ylabel('basepair')
    plt.colorbar()
    
    ### plot best random cj found
    plt.sca(axs[plotBySheet+0,0]) ### 1 if plotBySheet = 1
    plt.imshow(cjOpt.reshape(-1,1),vmin=-1,vmax=1)
    plt.xticks([])
    plt.ylabel('cj found by random')
    plt.colorbar(location='left',pad=.4)
    
    ### plot samplesheet cj
    plt.sca(axs[1-plotBySheet,0]) ### 0 if plotBySheet=0
    plt.imshow(cjSheet.reshape(-1,1),vmin=-1,vmax=1)
    plt.xticks([])
    plt.ylabel('Samplesheet cj')
    plt.colorbar(location='left',pad=.4)

    if plotConsensus:
        sampleNames = np.array(sampleNames)[relevantSamples][sampOrdering]
        df = pd.DataFrame(data={'sample':sampleNames,'cjRand':cjOpt,'cjSheet':cjSheet})
        df = df.set_index('sample')

        for samp in sampleNames:
            frac_fname = consensusFldr+'/{}_down_fractions.tab'.format(samp)
            with open(frac_fname,'r') as f:
                lines = f.readlines()
                colcount = max([len(l.strip().split('\t')) for l in lines])
            colcount-=2 ### anchor and consensus
            dfc = pd.read_csv(frac_fname
                              ,sep='\t',names=['anchor','consensus']+['frac_'+str(x) for x in range(colcount)])
            if anch not in dfc.anchor.values:
                continue
            row=dfc.set_index('anchor').loc[anch]

            df.loc[samp,['consensus']+['frac_'+str(x) for x in range(colcount)]] = row

        df.reset_index().to_csv(
            saveFldr+'/{}_consensus.csv'.format(anch),
                                sep='\t',index=False)

        ### plotting
        maxLen = len(df.columns)-3 ### cjRand,cjSheet,consensus (sample is index)
        consensi = df.loc[sampleNames,'consensus'].to_list()
        consensusMat = -1*np.ones((len(consensi),maxLen)) ### negative 1 is nan
        for i,cons in enumerate(consensi):
            consensusMat[i,:len(cons)] = [bpToInt(x) for x in cons]

        fracs = df.loc[sampleNames,['frac_'+str(x) for x in range(maxLen)]]
        fracMat = fracs.fillna(0).to_numpy()

        ### plot consensi
        plt.sca(axs[0,2])
        plt.imshow(consensusMat,aspect='auto')
        plt.ylabel('sample')
        plt.xlabel('basepair')
        plt.ylabel('sample')
        plt.colorbar()

        ### plot fracs
        plt.sca(axs[1,2])
        plt.imshow(fracMat,aspect='auto')
        plt.ylabel('sample')
        plt.xlabel('basepair')
        plt.colorbar()
    
    
    ### outside loop, save
    plotSaveName = '/anch_{}_gene_{}_consensus_{}_n1_{}_n2_{}_n_{}_pv_{:.1E}_eSize_{:.2F}_muLev_{:.2F}.pdf'.format(
                            anch,
        pvRow.anchor_local_gene,
        pvRow.consensus_gene_mode,
        int(n1//5)*5, int(n2//5)*5, int(min(n1,n2)//5)*5,
    pvRow.pv_hash_both_corrected,
    max(pvRow['effect_size_randCjs'],pvRow['effect_size_sheetCjs']),
    pvRow.mu_lev)
    plt.savefig(saveFldr+plotSaveName,bbox_inches='tight')
    #### name contains n1,n2 floored to nearest 5, and p is minimum of these two
    
    plt.close()
        
    
    
#### currently merged into plotContingency     
# def saveConsensuses(dfpvals,anch,cjSheet,consensusFldr,saveFldr,
#                     cjOptOrd,cjSheetOrd,sampleNamesOrd):
# #     sampleNames = (dfpvals.columns[dfpvals.columns.str.startswith('cj_rand')]
# #               .str.lstrip('cj_rand_opt_').to_list())
#     pvRow = dfpvals[dfpvals.anchor==anch].mode(axis=0).iloc[0].copy()

# #     cjOpt = pvRow[dfpvals.columns.str.startswith('cj_rand')].to_numpy().flatten().astype('float')

# #     if cjOpt @ cjSheet < 0:
# #         cjOpt = -cjOpt
    
# #     sampleNames = np.array(sampleNames)[sampOrderingIdx]
# #     cjOpt = cjOpt[sampOrderingIdx]
# #     cjSheet = cjSheet[sampOrderingIdx]
#     sampleNames = sampleNamesOrd
#     cjOpt = cjOptOrd
#     cjSheet = cjSheetOrd

#     df = pd.DataFrame(data={'sample':sampleNames,'cjRand':cjOpt,'cjSheet':cjSheet})
#     df = df.set_index('sample')

#     for samp in sampleNames:
#         frac_fname = consensusFldr+'/{}_down_fractions.tab'.format(samp)
#         with open(frac_fname,'r') as f:
#             lines = f.readlines()
#             colcount = max([len(l.strip().split('\t')) for l in lines])
#         colcount-=2 ### anchor and consensus
#         dfc = pd.read_csv(frac_fname
#                           ,sep='\t',names=['anchor','consensus']+['frac_'+str(x) for x in range(colcount)])
#         if anch not in dfc.anchor.values:
#             continue
#         row=dfc.set_index('anchor').loc[anch]

#         df.loc[samp,['consensus']+['frac_'+str(x) for x in range(colcount)]] = row

#     df.reset_index().to_csv(
#         saveFldr+'/{}_consensus.csv'.format(anch),
#                             sep='\t',index=False)
    
    
#     ### plotting
#     maxLen = len(df.columns)-3 ### cjRand,cjSheet,consensus (sample is index)
#     consensi = df.loc[sampleNames,'consensus'].to_list()
#     consensusMat = -1*np.ones((len(consensi),maxLen)) ### negative 1 is nan
#     for i,cons in enumerate(consensi):
#         consensusMat[i,:len(cons)] = [bpToInt(x) for x in cons]
        
        
#     fracs = df.loc[sampleNames,['frac_'+str(x) for x in range(maxLen)]]
#     fracMat = fracs.fillna(0).to_numpy()
    
    
#     ### figure
#     fig, axs = plt.subplots(nrows=2,ncols=2)
#     fig.set_facecolor('white')
#     fig.set_size_inches(10, 10)
#     fig.suptitle('Normalized target distribution for each sample for anchor {}'.format(anch))
    
#     ### plot consensi
#     plt.sca(axs[0,1])
#     plt.imshow(consensusMat,aspect='auto')
#     plt.ylabel('sample')
# #     plt.yticks(ticks=range(len(sampleNames)),labels=sampleNames)
#     plt.xlabel('basepair')
#     plt.ylabel('sample')
#     plt.colorbar()
    
#     ### plot fracs
#     plt.sca(axs[1,1])
#     plt.imshow(fracMat,aspect='auto')
#     plt.ylabel('sample')
# #     plt.yticks(ticks=range(len(sampleNames)),labels=sampleNames)
#     plt.xlabel('basepair')
#     plt.colorbar()
    
#     plotBySheet = ~np.all(cjSheet==1)
#     ### plot best random cj found
#     plt.sca(axs[plotBySheet+0,0]) ### 1 if plotBySheet = 1
#     plt.imshow(cjOpt.reshape(-1,1),vmin=-1,vmax=1)
#     plt.xticks([])
#     plt.ylabel('cj found by random')
#     plt.colorbar(location='left',pad=.4)
    
#     ### plot samplesheet cj
#     plt.sca(axs[1-plotBySheet,0]) ### 0 if plotBySheet=0
#     plt.imshow(cjSheet.reshape(-1,1),vmin=-1,vmax=1)
#     plt.xticks([])
#     plt.ylabel('Samplesheet cj')
#     plt.colorbar(location='left',pad=.4)
    
#     plotSaveName = '/anch_{}_consensus_plot.pdf'.format(anch)
#     plt.savefig(saveFldr+plotSaveName,bbox_inches='tight')
#     plt.close()
    
    
### construct anchor target counts matrix from raw abundant_stratified files
def constructCountsDf(strat_fldr,anchLst):
    kmer_size = 27
    dfs = []
    paths = glob.glob(strat_fldr+"/abundant_stratified_*.txt.gz")
    print('reading in abudant_stratified files')
    for path in tqdm.tqdm(paths,total=len(paths)):
        df = pd.read_csv(path.strip(), delim_whitespace=True, names=['counts','seq','sample'])
        df['anchor'] = df.seq.str[:kmer_size]
        df['target'] = df.seq.str[kmer_size:]
        df = df.drop(columns='seq').drop_duplicates()
        df = df[df.anchor.isin(anchLst)]
        dfs.append(df)

    print('concat and pivoting')
    ctsDf = pd.concat(dfs)
    ctsDf = (ctsDf
        .pivot(index=['anchor', 'target'], columns='sample', values='counts')
        .reset_index()
        .fillna(0))
    return ctsDf
    
    
    
def main():
    print('starting')
    args = get_args()
    
    if not os.path.exists(args.pvals_all_file):
        print("no additional_summary tsv file")
        return
    
    if not os.path.exists(args.outfldr):
        os.makedirs(args.outfldr)
    
    ### count number of lines to check if chunking is needed
    numLines = int(os.popen('wc -l '+args.pvals_all_file).read().split()[0])
    print("summary file has {:.2E} lines".format(numLines))
    lineThresh = 2000000
    
    if numLines<lineThresh:
        print("File is small enough, reading p value file directly")
        dfpvals = pd.read_csv(args.pvals_all_file,sep='\t')
    else:
        print("numLines>lineThresh, chunking p value file")
        chunksize=10**6
        dfs=[]

        
        with pd.read_csv(args.pvals_all_file, sep='\t',chunksize=chunksize) as reader:
            ### stream over chunks of size chunksize
            for df in tqdm.tqdm(reader,total=numLines//chunksize+1):
                ### apply some filtering to make the dataframe smaller
                ### ***note that this is prefiltering for anchList***
                df = df[(~df.consensus_gene_mode.isna())]
                df = df[(df.consensus_gene_mode.apply(lambda x : len(x.split(' ')))==1) &
                       (df.consensus_gene_mode.apply(lambda x : len(x.split(',')))==1) &
                       (df.consensus_gene_mode != '[]') ]
                df = df[(df['target_top_ann']=='target_hits_grch38_1kgmaj') & 
                        (df['anchor_top_ann']=='anchor_hits_grch38_1kgmaj') & 
                        (df['anchor_num_ann']==1) & (df['target_num_ann']==1) &
                        (df['target_top_ann_hit'] == df['anchor_top_ann_hit'])]
                df = df.drop_duplicates()
                dfs.append(df)
                
        ### concatenate all the streamed chunks
        dfpvals = pd.concat(dfs)

        ### write to file, in case we want to check
        dfpvals.to_csv(args.outfldr+'/dfpvals_short.csv',sep='\t',index=False)
    
    dfpvanch = pd.read_csv(args.anchor_pvals,sep='\t')
    
    dfpvals = dfpvals.merge(dfpvanch) ### I don't think this should be losing anything?
    
    sampleNames = (dfpvals.columns[dfpvals.columns.str.startswith('cj_rand')]
              .str.lstrip('cj_rand_opt_').to_list())
    
    print('reading cjs')
    sheetCj = np.ones(len(sampleNames))
    if args.samplesheet!='':
        with open(args.samplesheet,'r') as f:
            cols = f.readline().split(',')
        if len(cols)==1: ### if len(cols) is 1, then only samplesheet name, no ids
            print("Only 1 samplesheet column, using random cjs")
        elif len(cols)>2:
            print("Improperly formatted samplesheet")
        else:
            

            sheetdf = pd.read_csv(args.samplesheet,names=['fname','handCjs'])
            sheetdf['sample'] = (sheetdf.fname
                            .str.rsplit('/',1,expand=True)[1]
                            .str.split('.',1,expand=True)[0])
            sheetdf = sheetdf.drop(columns='fname')
            sheetdf['handCjs'] = normalizevec(sheetdf.handCjs)

            sheetCj = sheetdf.set_index('sample').T[sampleNames].to_numpy().flatten()
    
    
    print('determining anchor list')
    
    if args.anch_list =='':
        print('no anchor list passed in, generating handcrafted one')
        df= dfpvals
        if np.all(sheetCj==1):
            df=df[df.effect_size_randCjs  >.3]
        else:
            df=df[df.effect_size_sheetCjs  >.3]
            
        ### *** note that prefiltering is applied on read in for chunking*** ###
        df = df[(~df.consensus_gene_mode.isna())]
        df = df[(df.consensus_gene_mode.apply(lambda x : len(x.split(' ')))==1) &
               (df.consensus_gene_mode.apply(lambda x : len(x.split(',')))==1) &
               (df.consensus_gene_mode != '[]') ]
        
        df = df[(df['M']>100) &
                (df.mu_lev>5)
               ]
        df = df[
            (df['target_top_ann']=='target_hits_grch38_1kgmaj') & 
            (df['anchor_top_ann']=='anchor_hits_grch38_1kgmaj') & 
            (df['anchor_num_ann']==1) & (df['target_num_ann']==1) &
            (df['target_top_ann_hit'] == df['anchor_top_ann_hit'])]
        df = df.sort_values('pv_hash_both_corrected')
        anchLst = df[['anchor']].drop_duplicates().anchor.to_list()
    else:
        print("using passed in anchor list")
        anchLst = pd.read_csv(args.anch_list,names=['anchor']).anchor.to_list()

    
    if True or args.anch_cts_file == '' or not os.path.exists(args.anch_cts_file):
#         print('no anch cts file, generating from scratch')
        print('hardwired to regenerate anchor target counts file from scratch')
        ctsDf = (constructCountsDf(args.strat_fldr,
                                   anchLst
                                   )
                 .set_index(['anchor','target']))
        for samp in sampleNames:
            if samp not in ctsDf.columns:
                ctsDf[samp] = 0
        ctsDf.reset_index().to_csv(args.outfldr+'/ctsDfmanual.csv',
                            sep='\t',index=False)
    else:
        print('reading in anchor target counts df')
        ctsDf = pd.read_csv(args.anch_cts_file, sep='\t').set_index(['anchor','target'])

        
    print('plotting')
    for anch in tqdm.tqdm(anchLst,total=len(anchLst)):
        if anch not in dfpvals.anchor.values:
            print("anch not in dfpvals")
            continue
        plotContingency(dfpvals,ctsDf,anch,sheetCj,args.outfldr,
                        args.consensus_fldr,True)
#         cjOptOrd,cjSheetOrd,sampleNamesOrd = plotContingency(dfpvals,ctsDf,anch,sheetCj,args.outfldr)
#         saveConsensuses(dfpvals,anch,sheetCj,args.consensus_fldr,args.outfldr, cjOptOrd,cjSheetOrd,sampleNamesOrd)
main()
            