#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import itertools
import glob
import argparse
import matplotlib as mpl
from matplotlib import colors
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator
import scipy
import scipy.stats
import time
import seaborn as sns


### command line inputs
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset",
        type=str,
        help='optional argument, for title and file naming',
        default='dataset'
    )
    parser.add_argument(
        "--anchor_pvals",
        type=str,
        help='required: anchors_pvals.tsv file'
    )
    parser.add_argument(
        "--anchor_Cjs",
        type=str,
        help='required: anchors_pvals.tsv file'
    )
    parser.add_argument(
        "--num_heatmap_anchors",
        type=int,
        default=100
    )
    parser.add_argument(
        "--samplesheet",
        type=str,
        help='NOMAD samplesheet file, if metadata is available',
        default=''
    )
    parser.add_argument(
        "--additional_summary",
        type=str,
        help='optional: additional_summary.tsv file',
        default=''
    )
    parser.add_argument(
        "--heatmap_anchor_list",
        type=str,
        help='optional: file of 1 column, no header, 1 anchor per line, plot these specific anchors',
        default=''
    )
    parser.add_argument(
        "--genome_annotations_anchors",
        type=str,
        help='optional: genome_annotations_anchor.tsv file',
        default=''
    )
    parser.add_argument(
        "--outfile_contingency_table",
        type=str
    )
    parser.add_argument(
        "--outfile_skipped_anchors",
        type=str
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
###   dfpvals: primary output csv file of previous step, currently additional_summary.tsv
###   dfcts: sample x target counts matrix
###   anch: anchor to plot
###   cjSheet: cj values from samplesheet
###   args: passed in argparser
def plotContingency(dfpvals, sampleNames, dfcts, anch, cjSheet, useSheetCjs, args,plotConsensus=True):
    start = time.time()

    ### check if samplesheet metadata passed in (cjSheet set to all 1s by default)
    plotBySheet = ~np.all(cjSheet==1)

    if anch not in dfcts.index:
        print(anch, " not in dfcts")
        return

    tmp = dfcts.loc[anch,sampleNames]
    tbl = tmp.fillna(0).to_numpy().T

    ### need to aggregate rows (different sample / anchors)
    pvRow = dfpvals[dfpvals.anchor==anch].mode(axis=0).iloc[0].copy()

    ### read in cjOpt
    cjOpt = pvRow[sampleNames].to_numpy().flatten().astype('float')

    ### if it's facing the wrong way, flip it (for visual purposes)
    if cjOpt @ cjSheet < 0 :
        cjOpt = -cjOpt

    ### construct some condition to filter out samples
    relevantSamples = cjOpt!=0 ### for larger datasets
    tbl=tbl[relevantSamples,:]
    cjOpt=cjOpt[relevantSamples]
    cjSheet=cjSheet[relevantSamples]

    ### order targets by abundance
    targOrdering = np.argsort(-tbl.sum(axis=0))
    tbl = tbl[:,targOrdering]
    sortedTargCounts = -np.sort(-tbl.sum(axis=0))
    # thresh = 2 ### for small datasets
    thresh = max(sortedTargCounts[min(10,len(targOrdering)-1)],5,sortedTargCounts[0]/10) ## for large datasets
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

    if not useSheetCjs:
        pvRow.loc['pval_samplesheet']=np.nan
        pvRow.loc['effect_size_samplesheet']=np.nan

    ### currently always true, but in case it should be changed for later
    if plotConsensus:
        fig, axs = plt.subplots(nrows=3,ncols=2,gridspec_kw={'height_ratios': [1, 5,5]})
        fig.set_size_inches(10, 10)
    else:
        print('here be dragons')
        return
        fig, axs = plt.subplots(nrows=2,ncols=2)
        fig.set_size_inches(10, 10)

    fig.set_facecolor('white')

    if useSheetCjs:
        (fig.suptitle('''Dataset {}, anchor {}
        Gene {}, consensus gene {}, transcript gene {}
        Corrected (both rand and sheet c)hash-based p value = {:.2E}
        (rand f, rand c) p value = {:.2E}, (rand f, sheet c) p value = {:.2E}
        Cosine similarity between random c and samplesheet c = {:.2F}
        Random c effect size = {:.2F}, Sheet c effect size = {:.2F}
        Mu_lev = {:.2F}, total number of observations M = {}
        Number of +1 samples = {}, number of -1 samples = {}'''
                .format(
                        args.dataset,
                        anch,
                        pvRow.gene,
                        pvRow.consensus_gene_mode,
                        pvRow.transcript,
                        pvRow['pval_aggregated_corrected'],
                        pvRow['pval_random'],
                        pvRow['pval_samplesheet'],
                        cjOpt@cjSheet/len(cjSheet),
                        pvRow['effect_size_random'],
                        pvRow['effect_size_samplesheet'],
                        pvRow['mean_target_levenshtein_distance'],
                        int(pvRow.num_observations),
                        n1,
                        n2
                    ),y=1.05))
    else:
        (fig.suptitle('''Dataset {}, anchor {}
        Gene {}, consensus gene {}, transcript gene {}
        (rand f, rand c) p value = {:.2E}, Corrected (rand f, rand c) p value = {:.2E}
        Cosine similarity between random c and samplesheet c = {:.2F}
        Random c effect size = {:.2F}
        Mu_lev = {:.2F}, total number of observations M = {}
        Number of +1 samples = {}, number of -1 samples = {}'''
                .format(
                        args.dataset,
                        anch,
                        pvRow.gene,
                        pvRow.consensus_gene_mode,
                        pvRow.transcript,
                        pvRow['pval_random'],
                        pvRow['pval_random_corrected'],
                        cjOpt@cjSheet/len(cjSheet),
                        pvRow['effect_size_random'],
                        pvRow['mean_target_levenshtein_distance'],
                        int(pvRow.num_observations),
                        n1,
                        n2
                    ),y=1.05))
    fig.subplots_adjust(hspace=.3)


    ### construct target matrix for visual printing
    targs = dfcts.loc[anch].index.values[targOrdering][relevantCols]
    targMat = np.zeros((27,len(targs)))
    for i,targ in enumerate(targs):
        targMat[:,i] = [bpToInt(x) for x in targ]


    ### coloring
    dna_colors = ['#ffffff','#8dd3c7','#ffffb3','#fb8072','#80b1d3'] #["N","A","T","C","G"]
    clust_colors= ['#ffbc54','#C86CE4'] ## [-1,1] orange, purple
    frac_cm = plt.cm.viridis_r
    p_aspect='auto'

    ### plot contingency table
    plt.sca(axs[1,0])
    plt.imshow( (tbl[:,relevantCols]/tbl.sum(axis=1)[:,None]).T,aspect=p_aspect,cmap=frac_cm, interpolation='nearest')
    plt.yticks([])
    plt.xticks(ticks=range(tbl.shape[0]),labels=tbl.sum(axis=1).astype('int'),rotation=-80)
    plt.xlabel(r'$n_j$ for sample j')
    plt.ylabel('Target')
    plt.yticks(ticks=range(relevantCols.sum()),labels=tbl[:,relevantCols].sum(axis=0).astype('int'))
    plt.colorbar(orientation="horizontal",pad=.3)

    plottedMat = tbl[:,relevantCols].T
    pd.DataFrame(
        {"sampleName":np.array(sampleNames)[relevantSamples][sampOrdering],
        "cj":cjSheet,
        "cts1":plottedMat[0,:],
        "cts2":plottedMat[1,:],
        "nj":tbl.sum(axis=1)}
    ).to_csv(f"anchor_{anch}_counts.tsv", sep='\t', index=False)


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

    ### currently always true
    if plotConsensus:
        sampleNames = np.array(sampleNames)[relevantSamples][sampOrdering]
        df = pd.DataFrame(data={'sample':sampleNames,'cjRand':cjOpt,'cjSheet':cjSheet})
        df = df.set_index('sample')

        for samp in sampleNames:
            frac_fname = f'{samp}_down_fractions.tab'
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

        (
            df
            .reset_index()
            .to_csv(
                f'{anch}_consensus.tsv',sep='\t',index=False
            )
        )

        ### plotting
        maxLen = len(df.columns)-3 ### cjRand,cjSheet,consensus (sample is index)
        consensi = df.loc[sampleNames,'consensus'].to_list()
        consensusMat = -1*np.ones((len(consensi),maxLen)) ### negative 1 is nan
        for i,cons in enumerate(consensi):
            consensusMat[i,:len(cons)] = [bpToInt(x) for x in cons]


        fracs = df.loc[sampleNames,['frac_'+str(x) for x in range(maxLen)]]
        fracMat = fracs.fillna(0).to_numpy()

        ### plot consensuses
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

    if useSheetCjs:
        pval_col = 'pval_aggregated_corrected'
    else:
        pval_col = 'pval_random_corrected'

    ### outside loop, save
    plotSaveName = '{}_anchor_{}_gene_{}_cons_{}_transcript_{}_n_{}_pv_{:.0E}_eSize_{:.2F}_muLev_{:.1F}_ent_{:.1F}.pdf'.format(
                            args.dataset,
                            anch,
        shorten(pvRow.gene,20),
        shorten(pvRow.consensus_gene_mode,50),
        pvRow.transcript,
        int((n1+n2)//5)*5,
    pvRow[pval_col],
    max(pvRow['effect_size_random'],np.abs(pvRow['effect_size_samplesheet'])),
    pvRow.mean_target_levenshtein_distance,
    pvRow.target_entropy)
    plt.savefig(plotSaveName,bbox_inches='tight') ### plot as pdf
    plt.savefig(plotSaveName[:-3]+"png",bbox_inches='tight') ### plot as png
    plt.close()


    #### plot contingency tables
    df = pd.read_csv(f"anchor_{anch}_counts.tsv" ,sep='\t')
    df['normed'] = df.cts1/df.nj

    label1=r'$c_j$=1'
    label2=r'$c_j$=-1'
#     if plotBySheet:
#         label1='capillary cell'
#         label2='macrophage'

    plt.figure()
    l1 = len(df[df.cj==1].normed)
    l2 = len(df[df.cj==-1].normed)
    bins=np.arange(-1E-6,1+1E-6,.1+1E-7)

    plt.hist(df[df.cj==1].normed,
         bins=bins,
         weights=np.ones(l1),
         label=label1,
        color=clust_colors[1],alpha=.4, edgecolor='black',hatch='\\')
    plt.hist(df[df.cj==-1].normed,
         bins=bins,
         weights=np.ones(l2),
         label=label2,
        color=clust_colors[0],alpha=.4, edgecolor='black',hatch='//')
    plt.title('{}_{}'.format(args.dataset,anch))
    plt.xlabel('Fraction occurence target 1',fontsize=15)
    plt.ylabel('Number of samples',fontsize=15)
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend()
    plt.savefig(f"{anch}_sampleTargFracHist.png", bbox_inches='tight')
    plt.close()


    ### empirical probability
    p = df.cts1.sum()/df.nj.sum()
    nmin = int(df.nj.min())
    nmax = int(df.nj.max())
    nRange = range(nmin,nmax+1)
    upperThresh = [scipy.stats.binom.ppf(.99,n,p)/n for n in nRange]
    lowerThresh=[scipy.stats.binom.ppf(.01,n,p)/n for n in nRange]

    binArr = np.zeros(len(df))
    pMax = 0
    numFlips = 0
    parr = []
    for i,row in df.iterrows():
        nval = row.nj
        n_ind = int(nval - nmin)
        binArr[i] = (row.normed > upperThresh[n_ind]) or (row.normed < lowerThresh[n_ind])

        pTrue = scipy.stats.binom.cdf(lowerThresh[n_ind]*nval-1,n=nval,p=p) + (1-scipy.stats.binom.cdf(upperThresh[n_ind]*nval,n=nval,p=p))
        pMax = max(pMax,pTrue)
        if pTrue>0:
            numFlips+=1



    ## 2 sided binomial test for heads probability .02 (0->.01,.99->1)
    pval = scipy.stats.binomtest(int(binArr.sum()),len(binArr),.02).pvalue

    ## 2 sided binomial test for heads probability pMax, with only number of feasible flips
    ### ideally, we would want a poisson binomial test, but those are hard to find
    newpv = scipy.stats.binomtest(int(binArr.sum()),numFlips,pMax).pvalue


    if plotBySheet and df.cj.nunique()==2:
        df['ctype_plot'] = 'c_j=+1' ####kait
        df.loc[df.cj==-1, 'ctype_plot'] = 'c_j=-1' ####kait

        g=sns.JointGrid(x=df.nj,y=df.normed,hue=df.ctype_plot, height=5, palette=clust_colors)
        g.plot_joint(sns.scatterplot)
        g.plot_marginals(sns.histplot,bins=10,multiple='layer',hue=df.cj,alpha=.4)
        hatches = ["////", "\\\\"]
        # Loop over the bars
        for bars, hatch in zip(g.ax_marg_y.containers, hatches):
            # Set a different hatch for each group of bars
            for bar in bars:
                bar.set_hatch(hatch)

    else:    ### with no samplesheet:
        g=sns.JointGrid(x=df.nj,y=df.normed,height=5, palette=['gray'])
        g.plot_joint(sns.scatterplot,label='sample',color='gray')
        g.plot_marginals(sns.histplot,bins=15,color=clust_colors[1],alpha=.4)

    g.ax_marg_x.remove()
    g.ax_joint.plot(nRange,upperThresh,ls='--',c='grey',label='.99 quantile')
    g.ax_joint.plot(nRange,lowerThresh,ls='--',c='grey',label='.01 quantile')
    g.ax_joint.legend(bbox_to_anchor=(1.5,1.5), borderaxespad=0)
    g.ax_joint.set_xlabel(r'Number of anchor observations in sample ($n_j$)',fontsize=15)
    g.ax_joint.set_ylim([-.1,1.1])
    g.ax_joint.set_ylabel('Fraction occurence target 1',fontsize=15)
    g.ax_joint.set_title('{} : {}, \n {}/{} hits with p<{:.3F} yields {:.2E}'.format(
        args.dataset,anch,
        int(binArr.sum()), numFlips, pMax, newpv))
    g.ax_joint.yaxis.set_major_locator(MaxNLocator(integer=True))



    plt.savefig(f"{anch}_njBinomScatter.png", bbox_inches='tight')
    plt.close()


### construct anchor target counts matrix from raw abundant_stratified files
def constructCountsDf(anchLst):
    kmer_size = 27
    dfs = []
    paths = glob.glob("abundant_stratified_*.txt.gz")
    print('reading in abudant_stratified files')
    for path in paths:
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


def parseSamplesheet(args):
    ### read in samplsheet, and output samplesheetCj (properly ordered)
    try:
        sheetdf = pd.read_csv(args.samplesheet,names=['fname','sheetCjs'])
    except:
        sheetdf = pd.read_csv(args.samplesheet,names=['fname'])

    sheetdf['sample'] = (
        sheetdf['fname']
        .str.rsplit('/',1,expand=True)[1]
        .str.split('.',1,expand=True)[0]
    )

    sampleNames = sheetdf['sample'].tolist()

    if len(sheetdf.columns)==2:
        useSheetCjs = True

        sheetdf = sheetdf.drop(columns='fname')
        sheetdf['sheetCjs'] = normalizevec(sheetdf['sheetCjs'])
        sheetCj = (
            sheetdf
            .set_index('sample')
            .T[sampleNames]
            .to_numpy()
            .flatten()
        )

    else:
        useSheetCjs = False
        sheetCj = np.ones(len(sampleNames))

    return sampleNames, sheetCj, useSheetCjs


def main():
    print('starting')
    args = get_args()

    ## If we aren't plotting any anchors, break
    if args.num_heatmap_anchors == 0:
        return

    ### read in passed in anchor list
    if args.heatmap_anchor_list!='':
        print("using passed in anchor list")
        anchLst = pd.read_csv(args.heatmap_anchor_list,names=['anchor']).anchor.to_list()
        if len(anchLst)==0:
            print('no anchors in list')
            return
        print(len(anchLst), "anchors")

    ## merge in anchor statistics with anchor cjs
    anchor_pvals = pd.read_csv(args.anchor_pvals,sep='\t')
    anchor_Cjs = pd.read_csv(args.anchor_Cjs,sep='\t')

    # sampleNames = anchor_Cjs.drop('anchor', axis=1).columns.tolist()
    dfpvals = pd.merge(anchor_pvals, anchor_Cjs, on='anchor')

    # Parsing samplesheet for sampleNames and sheetCjs (if they exist)
    sampleNames, sheetCj, useSheetCjs = parseSamplesheet(args)

    ### if we are using additional_summary.tsv file
    if args.additional_summary != '':
        if not os.path.exists(args.additional_summary):
            print("no additional_summary tsv file")
            return

        ### count number of lines to check if chunking is needed
        numLines = int(os.popen('wc -l '+args.additional_summary).read().split()[0])
        print("summary file has {:.2E} lines".format(numLines))
        lineThresh = 2000000

        if numLines<lineThresh:
            print("File is small enough, reading p value file directly")
            dfconcat = pd.read_csv(args.additional_summary,sep='\t')
        else:
            print("numLines>lineThresh, chunking p value file")
            chunksize=10**6
            dfs=[]

            with pd.read_csv(args.additional_summary, sep='\t',chunksize=chunksize) as reader:
                ### stream over chunks of size chunksize
                for df in reader:
                    df = df.drop_duplicates()

                    #### can discard aanchors that will not be plotted
                    if args.heatmap_anchor_list != '':
                        df = df[df.anchor.isin(anchLst)]
                    dfs.append(df)

            ### concatenate all the streamed chunks
            dfconcat = pd.concat(dfs)

        dfpvals = dfconcat.merge(dfpvals[['anchor'] + sampleNames])

    else:
        ### set gene annotations to be NA
        dfpvals['gene']='NA'
        dfpvals['consensus_gene_mode']='NA'

    if args.genome_annotations_anchors != '' and os.path.exists(args.genome_annotations_anchors):
        print('reading in genome_ann')
        dfgen = pd.read_csv(args.genome_annotations_anchors, sep='\t')
        dfgen['transcript'] = dfgen.transcript.astype(str).str.split('|').str[0]

        ### merge transcript into dfpvals
        dfpvals = dfpvals.merge(dfgen[['anchor', 'transcript']], how='left')
    else:
        print('No annotation file given, setting transcript gene to NA')
        dfpvals['transcript']='NA'

    ### determining anchor list
    ### can use filters to speed up, if not all anchors are desired
    if args.heatmap_anchor_list =='':
        print('no anchor list passed in, using top {} most significant anchors'.format(args.num_heatmap_anchors))
        df= dfpvals

        if useSheetCjs:
            df = df.sort_values('pval_aggregated_corrected')
        else:
            df = df.sort_values('pval_random_corrected')
        anchLst = df[['anchor']].drop_duplicates().anchor.to_list()[:args.num_heatmap_anchors]
        print("Length of original anchor list is {}, reduced to {}".format(df.anchor.nunique(),len(anchLst)))

    #### generate contingency tables from abundant_stratified_anchors files
    print('regenerating anchor target counts file from scratch')
    ctsDf = (constructCountsDf(anchLst).set_index(['anchor','target']))
    for samp in sampleNames:
        if samp not in ctsDf.columns:
            ctsDf[samp] = 0
    ctsDf.reset_index().to_csv(args.outfile_contingency_table,sep='\t',index=False)

    ### loop over anchors and plot
    skippedAnchs = []
    print('plotting')
    for anch in anchLst:
        if anch=='anchor':
            continue
        if anch not in dfpvals.anchor.values:
            skippedAnchs.append(anch)
            print("anch {} not in dfpvals".format(anch))
            continue

        plotContingency(dfpvals, sampleNames, ctsDf, anch, sheetCj, useSheetCjs, args, True)

    ### if something went awry, save list of skipped anchors
    if len(skippedAnchs)>0:
        pd.Series(skippedAnchs).to_csv(
            args.outfile_skipped_anchors,sep='\t',index=False)


### run main
print("TESTINGTESTING")
print('running main')
main()
