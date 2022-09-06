import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import itertools
import glob
import argparse
import tqdm
import matplotlib as mpl
from matplotlib import colors
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator
import scipy 
import scipy.stats

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( ### anchors_pvals.tsv file
        "--ds_name",
        type=str,
        default=""
    )
    parser.add_argument( ### anchors_pvals.tsv file
        "--anchor_pvals",
        type=str,
        default=""
    )
    parser.add_argument(
        "--kmer_size",
        type=int
    )
    parser.add_argument( ### take in additional_summary.tsv file
        "--pvals_all_file",
        type=str,
        default=""
    )
    parser.add_argument( ### folder to fractions.tab files
        "--consensus_fldr",
        type=str,
        default="" #/oak/stanford/groups/horence/kaitlin/results_nomad/bulk_RNAseq/paper/HLCA_SS2_P2_macrophage_capillary/consensus_anchors/
    )
    parser.add_argument( ### folder to samplesheet, to get samplesheet cjs
        "--samplesheet",
        type=str,
        default="" #/oak/stanford/groups/horence/kaitlin/running_nomad/bulk_RNAseq/samplesheets/samplesheet_HLCA_SS2_P2_capillary.csv
    )
    parser.add_argument( ### folder to save figures to
        "--outfldr",
        type=str,
        default="" #/oak/stanford/groups/horence/tavorb/HLCA_SS2_P2_capillary/figs/
    )
    parser.add_argument( ### anchor target counts file
        "--anch_cts_file",
        type=str,
        default="" #/oak/stanford/groups/horence/kaitlin/results_nomad/bulk_RNAseq/paper/HLCA_SS2_P2_capillary/anchor_top250_targets_counts.tsv
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
    parser.add_argument( ### genome_annotations_anchor.tsv file
        "--gen_ann",
        type=str,
        default=""
    )
    parser.add_argument( ### LEGACY. gtf file, not used
        "--gtf_file",
        type=str,
        default='/oak/stanford/groups/horence/tavorb/stringstats_tavor/gencode_basic_ann_.gtf.gz'
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

### For shortening strings to printable size
def shorten(x,l=15):
    if x!=x:
        return x
    if type(x) is str:
        return x[:l]
    return x

### Inputs:
#### dfpvals: primary output csv file of previous step, currently additional_summary.tsv
#### dfcts: sample x target counts matrix
#### anch: anchor to plot
#### cjSheet: cj values from samplesheet
#### savePath: folder to save output figure
def plotContingency(dfpvals,dfcts,anch,cjSheet,args,plotConsensus=False):
    
    saveFldr = args.outfldr
    consensusFldr = args.consensus_fldr
    
    
    plotBySheet = ~np.all(cjSheet==1)
    
    ### load in list of samples where this anchor is present
    sampleNames = (dfpvals.columns[dfpvals.columns.str.startswith('cj_rand')]
              .str.lstrip('cj_rand_opt_').to_list())
    
    if anch not in dfcts.index:
        print(anch, " not in dfcts")
        return
    
    
    tmp = dfcts.loc[anch,sampleNames]
    tbl = tmp.fillna(0).to_numpy().T
    
#     targsPrint = dfcts.loc[anch].index.values
#     print(anch,targsPrint,tbl.sum(axis=0))
    
    
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
    sortedTargCounts = -np.sort(-tbl.sum(axis=0))
    thresh = max(sortedTargCounts[min(10,len(targOrdering)-1)],5,sortedTargCounts[0]/10)
    thresh = min(thresh,sortedTargCounts[1]) ### need at least 2 targets
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
        fig, axs = plt.subplots(nrows=3,ncols=2,gridspec_kw={'height_ratios': [1, 5,5]})
        fig.set_size_inches(10, 10)
#         fig, axs = plt.subplots(nrows=3,ncols=2)
    else:
        print('here be dragons')
        return
        fig, axs = plt.subplots(nrows=2,ncols=2)
        fig.set_size_inches(10, 10)
        
    fig.set_facecolor('white')
    
    (fig.suptitle('''Dataset {}, anchor {}
    Local gene {}, consensus gene {}, transcript gene {}
    Corrected (both rand and sheet c)hash-based p value = {:.2E}
    (rand f, rand c) p value = {:.2E}, (hand f, rand c) p value = {:.2E}
    (rand f, hand c) p value = {:.2E}, (hand f, hand c) p value = {:.2E}
    Cosine similarity between random c and samplesheet c = {:.2F}
    Rand effect size = {:.2F}, SheetCj effect size = {:.2F}
    Mu lev = {:.2F}, Mu ham = {:.2F}, total number of reads used for inference M = {}
    Number of +1 samples = {}, number of -1 samples = {}'''
              .format(args.ds_name,
                  anch,
                      pvRow.anchor_local_gene,
                      pvRow.consensus_gene_mode,
                      pvRow.transcript_gene,
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

        
    ### coloring    
    dna_colors = ['#ffffff','#8dd3c7','#ffffb3','#fb8072','#80b1d3'] #["N","A","T","C","G"] 
    clust_colors = ['#e78ac3','#2c7bb6'] ##[-1,1]
    frac_cm = plt.cm.viridis_r
    p_aspect='auto'

    ### plot contingency table

    plt.sca(axs[1,0])
    plt.imshow( (tbl[:,relevantCols]/tbl.sum(axis=1)[:,None]).T,aspect=p_aspect,cmap=frac_cm, interpolation='nearest')
    plt.yticks([])
    plt.xticks(ticks=range(tbl.shape[0]),labels=tbl.sum(axis=1).astype('int'),rotation=-80)
    plt.xlabel(r'$n_j$ for sample j')
#     plt.ylabel('Target (value is # occurences)')
    plt.ylabel('Target')
    plt.yticks(ticks=range(relevantCols.sum()),labels=tbl[:,relevantCols].sum(axis=0).astype('int'))
    plt.colorbar(orientation="horizontal",pad=.3)
    
    plottedMat = tbl[:,relevantCols].T
    pd.DataFrame({"sampleName":np.array(sampleNames)[relevantSamples][sampOrdering],
                  "cj":cjSheet,
                  "cts1":plottedMat[0,:],
                  "cts2":plottedMat[1,:],
                 "nj":tbl.sum(axis=1)}).to_csv(
        saveFldr+"/anch_{}_cts.tsv".format(anch) ,sep='\t',index=False)
    
    
    ### plot target matrix with fancy colors
    plt.sca(axs[1,1])
    
    col_dict={0:dna_colors[1],
      1:dna_colors[2],
      2:dna_colors[3],
      3:dna_colors[4]}
    labels = np.array(["A","T","C","G"])
    len_lab = len(labels)
    cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

    norm_bins = np.sort([*col_dict.keys()]) + 0.5
    norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)

    norm = mpl.colors.BoundaryNorm(norm_bins, len_lab, clip=True)
    fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
    diff = norm_bins[1:] - norm_bins[:-1]
    tickz = norm_bins[:-1] + diff / 2
    
    
    im = plt.imshow(targMat.T,aspect=p_aspect, cmap=cm,norm=norm, interpolation='nearest')
    plt.gca().yaxis.tick_right()
    plt.yticks(ticks=range(len(targs)),labels=targs)
    plt.xlabel('basepair')
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    
    
    cb = plt.colorbar(im,orientation="horizontal",pad=0.3) ### just to get cb.ax
    cb = plt.colorbar(im, format=fmt, ticks=tickz,
                      orientation="horizontal", 
                     cax=cb.ax)

    ### samplesheet colorscheme
    ### breaks down for more than 2 clusters
    ##
    col_dict={
        -1:clust_colors[0],
        1:clust_colors[1]}
    labels = np.array([-1,1])
    len_lab = len(labels)
    cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

    norm_bins = np.sort([*col_dict.keys()]) + 0.5
    norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)

    norm = mpl.colors.BoundaryNorm(norm_bins, len_lab, clip=True)
    fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
    diff = norm_bins[1:] - norm_bins[:-1]
    tickz = norm_bins[:-1] + diff / 2
    
    
    ### plot best random cj found
    plt.sca(axs[0,0+plotBySheet]) ### 1 if plotBySheet = 1
    im = plt.imshow(cjOpt.reshape(-1,1).T,aspect=p_aspect,
                    cmap=cm,norm=norm,
               interpolation='nearest')
    plt.yticks([])
    plt.xlabel('cj found by random')
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    if not plotBySheet:
        imcb = im
    
    ### plot samplesheet cj
    plt.sca(axs[0,1-plotBySheet]) ### 0 if plotBySheet=0
    im = plt.imshow(cjSheet.reshape(-1,1).T,aspect=p_aspect,
                    cmap=cm,norm=norm,
               interpolation='nearest')
    plt.yticks([])
    plt.xlabel('Samplesheet cj')
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    if plotBySheet:
        imcb = im

        
    ### joint colorbar
    plt.sca(axs[0,1])
    
    box = plt.gca().get_position()
    cbarHeight = box.width/2
    axColor = plt.axes([box.x0*1.05 + box.width * 1.05, box.y0+box.height/2-cbarHeight/2, 0.01, cbarHeight]) ## left, bottom, width, height 
    cb = plt.colorbar(imcb, cax = axColor, orientation="vertical")
    cb = plt.colorbar(im, format=fmt, ticks=tickz,
                      orientation="vertical", 
                     cax=cb.ax)
    
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
        plt.sca(axs[2,0])
        col_dict={
            -1:dna_colors[0],
            0:dna_colors[1],
          1:dna_colors[2],
          2:dna_colors[3],
          3:dna_colors[4]}
        labels = np.array(["N","A","T","C","G"])
        len_lab = len(labels)
        cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

        norm_bins = np.sort([*col_dict.keys()]) + 0.5
        norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)

        norm = mpl.colors.BoundaryNorm(norm_bins, len_lab, clip=True)
        fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
        diff = norm_bins[1:] - norm_bins[:-1]
        tickz = norm_bins[:-1] + diff / 2
        
        
        im = plt.imshow(consensusMat.T,aspect='auto',
                        cmap=cm,norm=norm,
                        interpolation='nearest')
        plt.xlabel('sample')
        plt.ylabel('basepair')
        plt.xlabel('sample')
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))

        cb = plt.colorbar(im,orientation="horizontal",pad=0.2) ### just to get cb.ax
        cb = plt.colorbar(im, format=fmt, ticks=tickz,
                          orientation="horizontal", 
                         cax=cb.ax)

        ### plot fracs
        plt.sca(axs[2,1])
        plt.imshow(fracMat.T,aspect='auto',
                   cmap=frac_cm,
                   interpolation='nearest')
        plt.xlabel('sample')
        plt.ylabel('basepair')
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.colorbar(orientation="horizontal",pad=.2)

    
    ### outside loop, save
    plotSaveName = '/{}_anch_{}_gene_{}_cons_{}_tran_{}_n_{}_pv_{:.0E}_eSize_{:.2F}_muLev_{:.1F}_ent_{:.1F}.pdf'.format(
                            args.ds_name,
                            anch,
        shorten(pvRow.anchor_local_gene,20),
        shorten(pvRow.consensus_gene_mode,50),
        pvRow.transcript_gene,
        int((n1+n2)//5)*5,
    pvRow.pv_hash_both_corrected,
    max(pvRow['effect_size_randCjs'],np.abs(pvRow['effect_size_sheetCjs'])),
    pvRow.mu_lev,
    pvRow.entropy)
    plt.savefig(saveFldr+plotSaveName,bbox_inches='tight')
    plt.savefig(saveFldr+plotSaveName[:-3]+"jpg",bbox_inches='tight')
    #### name contains n1,n2 floored to nearest 5, and p is minimum of these two
    
    plt.close()
    
    
    #### plot contingency tables
    df = pd.read_csv(saveFldr+"/anch_{}_cts.tsv".format(anch) ,sep='\t')
    df['normed'] = df.cts1/df.nj
    
    plt.figure()
    l1 = len(df[df.cj==1].normed)
    l2 = len(df[df.cj==-1].normed)
    bins=np.arange(-1E-6,1+1E-6,.1+1E-7)

    plt.hist(df[df.cj==1].normed,
         bins=bins,
         weights=np.ones(l1),
#          label=r'$c_j$=1',
        label='capillary cell',
        color=clust_colors[1],alpha=.4, edgecolor='black',hatch='\\')
    plt.hist(df[df.cj==-1].normed,
         bins=bins,
         weights=np.ones(l2),
#          label=r'$c_j$=-1',
         label='macrophage',
        color=clust_colors[0],alpha=.4, edgecolor='black',hatch='/')
    plt.title('{}_{}'.format(args.ds_name,anch))
    plt.xlabel('Fraction occurence target 1',fontsize=15)
    plt.ylabel('Number of samples',fontsize=15)
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend()
    plt.savefig(saveFldr+"/{}_predHist.jpg".format(anch),bbox_inches='tight')
    plt.close()
    

    ### empirical probability 
    p = df.cts1.sum()/df.nj.sum()
    nmin = int(df.nj.min())
    nmax = int(df.nj.max())
    nRange = range(nmin,nmax+1)
    upperThresh = [scipy.stats.binom.ppf(.99,n,p)/n for n in nRange]
    lowerThresh=[scipy.stats.binom.ppf(.01,n,p)/n for n in nRange]

    binArr = np.zeros(len(df))
    for i,row in df.iterrows():
        nval = row.nj
        n_ind = int(nval - nmin)
        binArr[i] = (row.normed > upperThresh[n_ind]) or (row.normed < lowerThresh[n_ind])

    pval = scipy.stats.binomtest(int(binArr.sum()),len(binArr),.02).pvalue
    
    
    plt.figure(figsize=(5,5))
    ### with no samplesheet:
#     plt.scatter(df.nj,df.normed
#                 ,label='cell')
    plt.scatter(df[df.cj==1].nj,df[df.cj==1].normed
                ,label='capillary cell',c=clust_colors[1])
    plt.scatter(df[df.cj==-1].nj,df[df.cj==-1].normed
                ,label='macrophage',c=clust_colors[0])
    plt.plot(nRange,upperThresh,'--r',label='.99 quantile')
    plt.plot(nRange,lowerThresh,'--g',label='.01 quantile')
    # plt.plot(nRange,median,'-.',label='median')
    plt.axhline(y=p,ls='-.',label='mean')
    plt.legend()
    plt.xlabel(r'Number of occurences in sample ($n_j$)',fontsize=15)
    plt.ylabel('Fraction occurence target 1',fontsize=15)
    plt.title('{}_{}, \n {}/{} hits yields p value of {:.2E}'.format(
        args.ds_name,anch,
        int(binArr.sum()),len(binArr),pval))
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.savefig(saveFldr+"/{}_pred_nj.jpg".format(anch),bbox_inches='tight')
    plt.close()
    
    
### construct anchor target counts matrix from raw abundant_stratified files
def constructCountsDf(strat_fldr,anchLst):
    kmer_size = args.kmer_size
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
    
    if args.anch_list !='':
        print("using passed in anchor list")
        anchLst = pd.read_csv(args.anch_list,names=['anchor']).anchor.to_list()
        if len(anchLst)==0:
            print('no anchors in list')
            return
        print(len(anchLst), "anchors")

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
#                 df = df[(df['target_top_ann']=='target_hits_grch38_1kgmaj') & 
#                         (df['anchor_top_ann']=='anchor_hits_grch38_1kgmaj') & 
#                         (df['anchor_num_ann']==1) & (df['target_num_ann']==1) &
#                         (df['target_top_ann_hit'] == df['anchor_top_ann_hit'])]
                df = df.drop_duplicates()
                if args.anch_list != '':
                    df = df[df.anchor.isin(anchLst)]
                dfs.append(df)
                
        ### concatenate all the streamed chunks
        dfpvals = pd.concat(dfs)

        ### write to file, in case we want to check
        dfpvals.to_csv(args.outfldr+'/dfpvals_short.csv',sep='\t',index=False)
    
    dfpvanch = pd.read_csv(args.anchor_pvals,sep='\t') ### need this for rand cj
    
    dfpvals = dfpvals.merge(dfpvanch[['anchor'] +
                                     dfpvanch.columns[dfpvanch.columns.str.startswith('cj_rand_opt')]
                                     .to_list()])
    
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
    
    
    
    print('reading in genome_ann')
    dfgen = pd.read_csv(args.gen_ann, sep='\t')
    dfgen['transcript_first'] = dfgen.local_transcript.str.split('|').str[0]
    dfgen['transcript_gene'] = dfgen.end_to_end_transcript.str.split('|').str[5] ### new
    


    ### merge transcript_gene into dfpvals
    dfpvals = dfpvals.merge(dfgen[['anchor', 'transcript_gene']], how='left')

    
    
    ### determining anchor list
    if args.anch_list =='':
        print('no anchor list passed in, generating handcrafted one')
        ### *** note that prefiltering is applied on read in for chunking*** ###
        
        
        df= dfpvals
#         if np.all(sheetCj==1):
#             df=df[df.effect_size_randCjs  >.3]
#         else:
#             df=df[df.effect_size_sheetCjs  >.3]
            
#         df = df[(~df.consensus_gene_mode.isna())]
#         df = df[(df.consensus_gene_mode.apply(lambda x : len(x.split(' ')))==1) &
#                (df.consensus_gene_mode.apply(lambda x : len(x.split(',')))==1) &
#                (df.consensus_gene_mode != '[]') ]
        
        df = df[(df['M']>100) #&
#                 (df.mu_lev>5)
               ]
        df = df[
            (df['target_top_ann']=='target_hits_grch38_1kgmaj') & 
            (df['anchor_top_ann']=='anchor_hits_grch38_1kgmaj') & 
            (df['anchor_num_ann']==1) & (df['target_num_ann']==1)]
        
        df = df.sort_values('pv_hash_both_corrected')
        anchLst = df[['anchor']].drop_duplicates().anchor.to_list()
        print("Length of original anchor list is {}, reduced to 300".format(len(anchLst)))
        anchLst = anchLst[:500]
        

        

    ctsDfExists=False
    needHandCtsDf = True
#     if False and args.anch_cts_file != '' and os.path.exists(args.anch_cts_file):
#         print("-----These files need to be regenerated, don't use---")
#         print('reading in anchor target counts df')
#         ctsDf = pd.read_csv(args.anch_cts_file, sep='\t')
#         anchIn = pd.Series(anchLst).isin(ctsDf.anchor)
#         print("fraction of anchors in ctsDf = {:.2F}".format(anchIn.mean()))
#         needHandCtsDf = ~anchIn.all()
#         ctsDfExists=True
#         ctsDf = ctsDf.set_index(['anchor','target'])

    if needHandCtsDf:
        print('regenerating anchor target counts file from scratch')
        ctsDfRaw = (constructCountsDf(args.strat_fldr,anchLst)
                 .set_index(['anchor','target']))
        for samp in sampleNames:
            if samp not in ctsDfRaw.columns:
                ctsDfRaw[samp] = 0
        ctsDfRaw.reset_index().to_csv(args.outfldr+'/ctsDfmanual.csv',
                            sep='\t',index=False)
        if ctsDfExists ==False:
            ctsDf = ctsDfRaw
    
        
    skippedAnchs = []
        
    print('plotting')
    for anch in tqdm.tqdm(anchLst,total=len(anchLst)):
        if anch not in dfpvals.anchor.values:
            skippedAnchs.append(anch)
            print("anch {} not in dfpvals".format(anch))
            continue
            
        ctsDfToUse = ctsDf
        if anch not in ctsDf.index:
            print(anch, " not in dfcts, using raw counts")
            ctsDfToUse = ctsDfRaw
                        
        plotContingency(dfpvals,ctsDfToUse,anch,sheetCj,args,True)
        
    
    ### if something went awry, save list of skipped anchors
    if len(skippedAnchs)>0:
        pd.Series(skippedAnchs).to_csv(args.outfldr+'/skippedAnchs.csv',
                            sep='\t',index=False)
    
### run main
main()
            
