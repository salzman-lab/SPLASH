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
    parser.add_argument(
        "--pvals_all_file",
        type=str,
        default="/oak/stanford/groups/horence/tavorb/HLCA_SS2_P2_capillary/pvals_stratified/pvals_all_ann.csv"
    )
    parser.add_argument(
        "--ann_file_consensus",
        type=str,
        default="/oak/stanford/groups/horence/kaitlin/results_nomad/bulk_RNAseq/paper/HLCA_SS2_P2_capillary/consensus_anchors/splicing_annotations/consensus_called_exons.tsv"
    )
    parser.add_argument(
        "--samplesheet",
        type=str,
        default="/oak/stanford/groups/horence/kaitlin/running_nomad/bulk_RNAseq/samplesheets/samplesheet_HLCA_SS2_P2_capillary.csv"
    )
    parser.add_argument(
        "--outfldr",
        type=str,
        default="/oak/stanford/groups/horence/tavorb/HLCA_SS2_P2_capillary/figs/"
    )
    parser.add_argument(
        "--anch_cts_file",
        type=str,
        default="/oak/stanford/groups/horence/kaitlin/results_nomad/bulk_RNAseq/paper/HLCA_SS2_P2_capillary/anchor_top250_targets_counts.tsv"
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
#### dfpvals: primary output csv file of previous step, currently anchors_pvals.tsv
#### dfcts: sample x target counts matrix
#### anch: anchor to plot
#### cjSheet: cj values from samplesheet
#### savePath: location to save output figure
def plotContingency(dfpvals,dfcts,anch,cjSheet=-1,savePath = ''):

    ### load in list of samples where this anchor is present
    sampleNames = (dfpvals.columns[dfpvals.columns.str.startswith('cj_rand')]
              .str.lstrip('cj_rand_opt_').to_list())
    
    if anch not in dfcts.index:
        print(anch, " not in dfcts")
        return
    
    tmp = dfcts.loc[anch,sampleNames]
    tbl = tmp.fillna(0).to_numpy().T
    


    cjOpt = dfpvals.loc[anch
                        ,dfpvals.columns.str.startswith('cj_rand')].to_numpy().flatten().astype('float')
    cjHand = cjSheet
#     print(cjSheet)
    if cjOpt @ cjHand < 0 :
        cjOpt = -cjOpt


    ### order targets by abundance
    targOrdering = np.argsort(-tbl.sum(axis=0))
    tbl = tbl[:,targOrdering]
    
    
    ### sort samples first by cj, then by abundance of first target
    toSort = list(zip(cjOpt,-tbl[:,0]/np.maximum(tbl.sum(axis=1),1)))
    sampOrdering = [x[1] for x in sorted((e,i) for i,e in enumerate(toSort))]
    tbl=tbl[sampOrdering,:]
    cjOpt=cjOpt[sampOrdering]
    cjHand=cjHand[sampOrdering]
    
    ### construct some condition to filter out rows
#     relevantSamples = tbl.sum(axis=1)>0 ### one possible condition
    relevantSamples = cjOpt!=0
    tbl=tbl[relevantSamples,:]
    cjOpt=cjOpt[relevantSamples]
    cjHand=cjHand[relevantSamples]

    dfrow = dfpvals.loc[anch]
    if 'pv_hand_sheetCjs' not in dfpvals.columns:
        dfrow.loc['pv_hand_sheetCjs']=np.nan
        dfrow.loc['pv_hash_sheetCjs']=np.nan
    fig, axs = plt.subplots(nrows=2,ncols=2)
    fig.set_facecolor('white')
    fig.set_size_inches(10, 10)
    (fig.suptitle('''Normalized target distribution for each sample for anchor {}
    Local gene {}, consensus gene {}, corrected hash-based p value = {:.2E}
    handcrafted di p value ={:.2E}, hash-based p value = {:.2E}
    hand cs and di p value = {:.2E}, hand cs and hash di p value = {:.2E}
    Cosine similarity between random c and samplesheet c = {:.2F}
    Effect size = {:.2F}, mu hand = {:.2F}, total number of reads M = {}'''
              .format(anch,
                      dfrow.agg_local_gene,
                      dfrow.agg_consensus_gene,
                      dfrow['pv_hash_corrected'],
                     dfrow['pv_hand'],
                     dfrow['pv_hash'],
                     dfrow['pv_hand_sheetCjs'],
                     dfrow['pv_hash_sheetCjs'],
                     cjOpt@cjHand/len(cjHand),
                      dfrow['effect_size'],
                     dfrow['mu_hand'],
                     int(dfrow.M)),
                 y=1.05))

    fig.subplots_adjust(hspace=.3)

    ### maybe should change this to be percentile based
    vals = tbl.sum(axis=0)
    relevantCols = tbl.sum(axis=0)>vals[min(10,len(vals)-1)]


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
    plt.sca(axs[0,0])
    plt.imshow(cjOpt.reshape(-1,1),vmin=-1,vmax=1)
    plt.xticks([])
    plt.ylabel('cj found by random')
    plt.colorbar(location='left',pad=.4)
    
    ### plot samplesheet cj
    plt.sca(axs[1,0])
    plt.imshow(cjHand.reshape(-1,1),vmin=-1,vmax=1)
    plt.xticks([])
    plt.ylabel('hand defined cj')
    plt.colorbar(location='left',pad=.4)

    ### save or show
    if savePath== '':
        plt.show()
    else:
        plt.savefig(savePath,bbox_inches='tight')

    plt.close()
        
        
        
def generateFiles(args):
    print('generating files')
    kmer_size = 27
    
    dfpvals = pd.read_csv(args.pvals_all_file,sep='\t')
    
    sampleNames = (dfpvals.columns[dfpvals.columns.str.startswith('cj_rand')]
              .str.lstrip('cj_rand_opt_').to_list())
    
    print('reading cjs')
    sheetCj = np.ones(len(sampleNames))
    if args.samplesheet!='': #### add in
        def normalizevec(x): ### shift and scale vector so that it's in [-1,1]
            return 2*(x-x.min()) / (x.max()-x.min())-1

        sheetdf = pd.read_csv(args.samplesheet,names=['fname','handCjs'])
        sheetdf['sample'] = (sheetdf.fname
                        .str.rsplit('/',1,expand=True)[1]
                        .str.split('.',1,expand=True)[0])
        sheetdf = sheetdf.drop(columns='fname')
        # sheetdf.loc[0,'handCjs']=10
        sheetdf['handCjs'] = normalizevec(sheetdf.handCjs)
        
        sheetCj = sheetdf.set_index('sample').T[sampleNames].to_numpy().flatten()
    
    print('reading annotations', len(dfpvals))
#     ### merge in annotations for local_gene
#     ann_genome = pd.read_csv(args.ann_file,sep='\t')
#     dfpvals = dfpvals.merge(ann_genome)
        
#     print('merged genome ann', len(ann_genome),len(dfpvals))
    consensusdf = pd.read_csv(args.ann_file_consensus, sep='\t'
                                , usecols=['anchor','anchor_local_gene','consensus_gene'])
    print('read in consensus, now filtering', len(consensusdf))
    consensusdf = consensusdf.drop_duplicates().fillna('')
    consensusdf['agg_consensus_gene'] = consensusdf.groupby('anchor').consensus_gene.transform(lambda x: ','.join(x))
    consensusdf['agg_local_gene'] = consensusdf.groupby('anchor').consensus_gene.transform(lambda x: ','.join(x))
    consensusdf = consensusdf[['anchor','agg_local_gene','agg_consensus_gene']].drop_duplicates()
    consensusdf.agg_consensus_gene = consensusdf.agg_consensus_gene.apply(lambda x: ','.join(list(set(x.split(',')))[:5]))
    consensusdf.agg_local_gene = consensusdf.agg_local_gene.apply(lambda x: ','.join(list(set(x.split(',')))[:5]))

    
    dfpvals = dfpvals.merge(consensusdf,how='left').fillna('')
    print('merged in consensus', len(dfpvals))
    
    if not os.path.exists(args.outfldr):
        os.mkdir(args.outfldr)
        
    print('saving df')
    dfpvals.to_csv(args.outfldr + '/pvals_full_ann_plot.csv', sep='\t', index=False)
#     ctsDf.to_csv(args.outfldr + '/ctsDf_plot.csv', sep='\t')
        
    
    
    

def main():
    print('starting')
    args = get_args()
    
    
    generateFiles(args)
    
#     return
    print('rereading files')
    dfpvals = pd.read_csv(args.outfldr + '/pvals_full_ann_plot.csv', sep='\t')
        
    ctsDf = pd.read_csv(args.anch_cts_file, sep='\t').set_index(['anchor','target'])
    
    sampleNames = (dfpvals.columns[dfpvals.columns.str.startswith('cj_rand')]
              .str.lstrip('cj_rand_opt_').to_list())
    
    sheetCj = np.ones(len(sampleNames))
    if args.samplesheet!='':
        def normalizevec(x): ### shift and scale vector so that it's in [-1,1]
            return 2*(x-x.min()) / (x.max()-x.min())-1

        sheetdf = pd.read_csv(args.samplesheet,names=['fname','handCjs'])
        sheetdf['sample'] = (sheetdf.fname
                        .str.rsplit('/',1,expand=True)[1]
                        .str.split('.',1,expand=True)[0])
        sheetdf = sheetdf.drop(columns='fname')
        # sheetdf.loc[0,'handCjs']=10
        sheetdf['handCjs'] = normalizevec(sheetdf.handCjs)
        
        sheetCj = sheetdf.set_index('sample').T[sampleNames].to_numpy().flatten()
    
    
    print('plotting')
    dfpvals = dfpvals.set_index('anchor')
    anchLst = dfpvals[(dfpvals.effect_size>.5)
             & (dfpvals.mu_hand>.4)].index.values[:500]
    for anch in tqdm.tqdm(anchLst,total=len(anchLst)):
            plotContingency(dfpvals,ctsDf,anch,sheetCj,
                            args.outfldr+'/anch_{}_gene_{}_consensus_{}.pdf'.format(
                                anch,dfpvals.loc[anch].agg_local_gene,dfpvals.loc[anch].agg_consensus_gene))

main()
            