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
    
    args = parser.parse_args()
    return args


#### For pretty plotting of targets, to convert basepair to int (color)
def bpToInt(x):
    if x=='A':
        return 0
    if x=='T':
        return 1
    if x=='C':
        return 2
    return 3


### Inputs:
#### dfpvals: primary output csv file of previous step, currently additional_summary.tsv
#### dfcts: sample x target counts matrix
#### anch: anchor to plot
#### cjSheet: cj values from samplesheet
#### savePath: folder to save output figure
def plotContingency(dfpvals,dfcts,anch,cjSheet,saveFldr = ''):
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
    pvRow = dfpvals.loc[anch].mode().iloc[0].copy()
    

    cjOpt = pvRow[dfpvals.columns.str.startswith('cj_rand')].to_numpy().flatten().astype('float')
    
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
    

    n1 = (cjSheet>0).sum()
    n2 = (cjSheet<0).sum()

    if 'pv_hand_sheetCjs' not in dfpvals.columns:
        pvRow.loc['pv_hand_sheetCjs']=np.nan
        pvRow.loc['pv_hash_sheetCjs']=np.nan
        pvRow.loc['effect_size_sheetCjs']=np.nan
    fig, axs = plt.subplots(nrows=2,ncols=2)
    fig.set_facecolor('white')
    fig.set_size_inches(10, 10)
    (fig.suptitle('''Normalized target distribution for each sample for anchor {}
    Local gene {}, consensus gene {}, corrected hash-based p value = {:.2E}
    (rand f, rand c) p value = {:.2E}, (hand f, rand c) p value = {:.2E}
    (rand f, hand c) p value = {:.2E}, (hand f, hand c) p value = {:.2E}
    Cosine similarity between random c and samplesheet c = {:.2F}
    Rand effect size = {:.2F}, SheetCj effect size = {:.2F}
    Mu lev = {:.2F}, Mu ham = {:.2F}, total number of reads used for inference M = {}
    Number of +1 samples = {}, number of -1 samples = {}'''
              .format(anch,
                      pvRow.anchor_local_gene,
                      pvRow.consensus_gene_mode,
                      pvRow['pv_hash_corrected'],
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

    ### save or show
    if saveFldr== '':
        plt.show()
    else:
        plotSaveName = '/anch_{}_gene_{}_consensus_{}_n1_{}_n2_{}_p_{}.pdf'.format(
                                anch,
            pvRow.anchor_local_gene,
            pvRow.consensus_gene_mode,
            int(n1//5)*5, int(n2//5)*5, int(min(n1,n2)//5)*5)
        plt.savefig(saveFldr+plotSaveName,bbox_inches='tight')
        #### name contains n1,n2 floored to nearest 5, and p is minimum of these two
    plt.close()
        
        
def saveConsensuses(dfpvals,anch,cjSheet,consensusFldr,saveFldr):
    sampleNames = (dfpvals.columns[dfpvals.columns.str.startswith('cj_rand')]
              .str.lstrip('cj_rand_opt_').to_list())
    pvRow = dfpvals.loc[anch].mode().iloc[0].copy()

    cjOpt = pvRow[dfpvals.columns.str.startswith('cj_rand')].to_numpy().flatten().astype('float')
    
    if cjOpt @ cjSheet < 0 :
        cjOpt = -cjOpt
    
    df = pd.DataFrame(data={'sample':sampleNames,'cjRand':cjOpt,'cjSheet':cjSheet})
    df = df.set_index('sample')

    for samp in sampleNames:
        dfc = pd.read_csv(consensusFldr+'/{}_down_fractions.tab'.format(samp)
                          ,sep='\t',usecols=range(52),names=['anchor','consensus']+['frac_'+str(x) for x in range(50)]) ###### assuming 1000 is an upper bound
        dfc['consensus'] = dfc.consensus.str[:50]
        if anch not in dfc.anchor.values:
            continue
        row=dfc.set_index('anchor').loc[anch][:52]

        df.loc[samp,['consensus']+['frac_'+str(x) for x in range(50)]] = row

        
#     print(df)
    df.reset_index().to_csv(saveFldr+'/{}_consensus.csv'.format(anch),
                            sep='\t',index=False)

def main():
    print('starting')
    args = get_args()
    
    if not os.path.exists(args.outfldr):
        os.makedirs(args.outfldr)
    
    print("reading p value sheet")
    dfpvals = pd.read_csv(args.pvals_all_file,nrows = 1000000,sep='\t')
    
    sampleNames = (dfpvals.columns[dfpvals.columns.str.startswith('cj_rand')]
              .str.lstrip('cj_rand_opt_').to_list())
    
    print('reading cjs')
    sheetCj = np.ones(len(sampleNames))
    if args.samplesheet!='': #### add in
        with open(args.samplesheet,'r') as f:
            cols = f.readline().split(',')
        if len(cols)==1: ### if len(cols) is 1, then only samplesheet name, no ids
            print("Only 1 samplesheet column, using random cjs")
        elif len(cols)>2:
            print("Improperly formatted samplesheet")
        else:
            def normalizevec(x): ### shift and scale vector so that it's in [-1,1]
                return 2*(x-x.min()) / (x.max()-x.min())-1

            sheetdf = pd.read_csv(args.samplesheet,names=['fname','handCjs'])
            sheetdf['sample'] = (sheetdf.fname
                            .str.rsplit('/',1,expand=True)[1]
                            .str.split('.',1,expand=True)[0])
            sheetdf = sheetdf.drop(columns='fname')
            sheetdf['handCjs'] = normalizevec(sheetdf.handCjs)

            sheetCj = sheetdf.set_index('sample').T[sampleNames].to_numpy().flatten()
    
    print('reading in anchor target counts df')
    ctsDf = pd.read_csv(args.anch_cts_file, sep='\t').set_index(['anchor','target'])
    
    
    print('plotting')

        
    
    if args.anch_list =='': 
        print('no anchor list passed in, generating handcrafted one')
        df= dfpvals.copy()
        if np.all(sheetCj==1):
            df=df[df.pv_hash < .01]
        else:
            df = df[df.pv_hash_sheetCjs < .01]
        df = df[(df['mu_lev']>5) 
                & (df['M']>300)]
        anchLst = (df[
            (df['target_top_ann']=='target_hits_grch38_1kgmaj') & 
            (df['anchor_top_ann']=='anchor_hits_grch38_1kgmaj') & 
            (df['anchor_num_ann']==1) & (df['target_num_ann']==1) &
            (df['target_top_ann_hit'] == df['anchor_top_ann_hit'])]
                   [['anchor']].drop_duplicates().anchor.to_list())
    else:
        print("using passed in anchor list")
        anchLst = pd.read_csv(args.anch_list,names=['anchor']).anchor.to_list()

    dfpvals = dfpvals.set_index('anchor')
    for anch in tqdm.tqdm(anchLst,total=len(anchLst)):
        if anch not in dfpvals.index:
            print("anch not in dfpvals")
            continue
        plotContingency(dfpvals,ctsDf,anch,sheetCj,args.outfldr)
        saveConsensuses(dfpvals,anch,sheetCj,args.consensus_fldr,args.outfldr)
main()
            