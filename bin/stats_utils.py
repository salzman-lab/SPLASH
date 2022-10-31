from hashlib import new
import numpy as np
import pandas as pd
from tqdm import tqdm,tqdm_notebook, tqdm_pandas
import scipy
import scipy.stats
import statsmodels.api as sm
from pathlib import Path 
import nltk


#### Usage:
## This file provides statistically functionality for NOMAD.
## Theoretical analysis and additional details to appear in an upcoming preprint
## There are 4 sets of functions in this file: 
### 1. General utility files: samplesheet parsing, abundant stratified parsing, splitting counts matrices, etc.
### 2. c,f construction: These take as input a contingency table, and output c and f that are predicted to perform well on them
### 3. p-value and effect size computation: NOMAD allows for its default p-value computation and asymptotically valid p-value computation. Also included are chi-squared and LRT
### 4. Other data: computing average levenshtein distance, target entropy, etc on dataframes.


###### 1: General utility files

### read in samplsheet, and output samplesheetCj dataframe
def parseSamplesheet(samplesheet):
    sheetCj = 0
    success = False
    if samplesheet!='':
        with open(samplesheet,'r') as f:
            cols = f.readline().split(',')
        if len(cols)==1: ### if len(cols) is 1, then only samplesheet name, no ids
            print("Only 1 samplesheet column, using random cjs")
        elif len(cols)>2:
            print("Improperly formatted samplesheet")
        else:

            sheetdf = pd.read_csv(samplesheet,names=['fname','sheetCjs'])
            #### this approach needs to be changed for 10X
            sheetdf['sample'] = (sheetdf.fname
                            .str.rsplit('/',1,expand=True)[1]
                            .str.split('.',1,expand=True)[0])
            sheetdf = sheetdf.drop(columns='fname')
            sheetdf['sheetCjs'] = normalizevec(sheetdf.sheetCjs,-1,1)

            sheetCj = sheetdf.set_index('sample').T
            print('using samplesheet metadata: not fully tested, be warned')
            success=True
    return success, sheetCj

### construct contingency tables, only for anchors in anchLst
def constructCountsDf(args,anchLst):
    df = pd.read_csv(args.infile, delim_whitespace=True, names=['counts','seq','sample'])
    print('done reading')

    ### split seq into anchor and target
    ### for some reason some files have duplicates, so drop those
    df['anchor'] = df.seq.str[:args.kmer_size]
    df['target'] = df.seq.str[args.kmer_size:]
    df = df.drop(columns='seq').drop_duplicates()
    if len(anchLst)>0:
        df = df[df.anchor.isin(anchLst)]

    ### add in informational columns for each anchor for filtering purposes
    ### this is faster than filter command
    df['anch_uniqTargs'] = df.groupby('anchor').target.transform('nunique')
    df['anch_samples']= df.groupby('anchor')['sample'].transform('nunique')
    df['anchSample_cts'] = df.groupby(['anchor','sample']).counts.transform('sum')
    df['anch_cts'] = df.groupby('anchor').counts.transform('sum') ## number of reads per anchor

    df = df[
        (df.anch_cts > args.anchor_count_threshold) &
        (df.anch_uniqTargs > args.anchor_unique_targets_threshold) &
        (df.anch_samples > args.anchor_samples_threshold) &
        (df.anchSample_cts > args.anchor_sample_counts_threshold)
    ]

    #### would like to do this, but table gets too large, so iterating instead
    # print('pivoting')
    # ctsDf = (df
    #     .pivot(index=['anchor', 'target'], columns='sample', values='counts')
    #     .reset_index()
    #     .fillna(0))
    return df


### split contingency table into train and test data
def splitCounts(mat,downSampleFrac = .5): #slight modification of https://stackoverflow.com/questions/11818215/subsample-a-matrix-python
    keys, counts = zip(*[
    ((i,j), mat[i,j])
        for i in range(mat.shape[0])
        for j in range(mat.shape[1])
        if mat[i,j] > 0
    ])
    # Make the cumulative counts array
    counts = np.array(counts, dtype=np.int64)
    sum_counts = np.cumsum(counts)

    # Decide how many counts to include in the sample
    frac_select = downSampleFrac
    count_select = int(sum_counts[-1] * frac_select)

    # Choose unique counts
    ind_select = sorted(np.random.choice(range(sum_counts[-1]), count_select,replace=False))

    # A vector to hold the new counts
    out_counts = np.zeros(counts.shape, dtype=np.int64)

    # Perform basically the merge step of merge-sort, finding where
    # the counts land in the cumulative array
    i = 0
    j = 0
    while i<len(sum_counts) and j<len(ind_select):
        if ind_select[j] < sum_counts[i]:
            j += 1
            out_counts[i] += 1
        else:
            i += 1

    # Rebuild the matrix using the `keys` list from before
    out_mat = np.zeros(mat.shape, dtype=np.int64)
    for i in range(len(out_counts)):
        out_mat[keys[i]] = out_counts[i]
        
    return out_mat

### simple wrapper of the above to downsample matrix columnwise, as opposed to overall
def splitCountsColwise(mat,downSampleFrac = .5): #slight modification of https://stackoverflow.com/questions/11818215/subsample-a-matrix-python
    out_mat = np.zeros_like(mat)
    for j in range(mat.shape[1]):
        if mat[:,j].sum()>0:
            out_mat[:,j] = splitCounts(np.reshape(mat[:,j],(mat.shape[0],-1)),downSampleFrac).flatten()
        
    return out_mat

### shift and scale an input vector to be in the range [minval, maxval]
def normalizevec(x,minval=0,maxval=1):
    if x.max()==x.min():
        return np.zeros_like(x)
    x = np.array(x)
    x01= (x-x.min())/(x.max()-x.min())
    return x01*(maxval-minval)+minval




######### 2: c,f generation

#### get spectral c,f, from simple SVD (correspondence analysis style)
def get_spectral_cf_svd(X,tblShape=-1):

    ### for backwards comaptability 
    if tblShape==-1:
        tblShape=X.shape

    relevantTargs = X.sum(axis=1)>0
    relevantSamples = X.sum(axis=0)>0
    
    r,c=X.shape

    if relevantTargs.sum()<2 or relevantSamples.sum()<2:
        return np.zeros(c),np.zeros(r)

    X = X[np.ix_(relevantTargs,relevantSamples)]

    E = np.outer(X.sum(axis=1),X.sum(axis=0))/X.sum()
    svdmat = (X-E)@np.diag(1.0/X.sum(axis=0))
    u,s,vh=np.linalg.svd(svdmat,full_matrices=False) ### very important to set full_matrices to false, else memory errors

    cRaw = vh[0,:]
    fRaw = u[:,0]

    cGuess = cRaw
    fGuess = normalizevec(fRaw,0,1)


    fElong = .5*np.ones(tblShape[0])
    fElong[relevantTargs] = fGuess
    fOpt = fElong
    
    cElong = np.zeros(tblShape[1])
    cElong[np.arange(tblShape[1])[relevantSamples]]=cGuess ### fancy indexing
    cOpt = cElong
    
    return cOpt,fOpt


### starting at c, run alternating maximization
def altMaximize(X,c):
    #### if clustering put all in same cluster, perturb
    if np.all(c==c[0]):
        c[0] = -1*c[0]
        
    nj = X.sum(axis=0)
    njinvSqrt = 1.0/np.maximum(1,np.sqrt(nj)) ### avoid divide by 0 errors
    njinvSqrt[nj==0]=0
    
    Xtild = (X - 1.0/X.sum()*np.outer(X@np.ones(X.shape[1]), X.T@np.ones(X.shape[0]))) @ np.diag(njinvSqrt)
    
    Sold = 0
    i=0
    while True:
        ### find optimal f for fixed c
        f = np.sign(Xtild @ c)
        f1 = (f+1)/2 ### to rescale f to be [0,1] valued
        f2 = (1-f)/2
        f = f1
        if np.abs(f2@Xtild@c) > np.abs(f1@Xtild@c):
            f = f2
        
        ### find optimal c for fixed f
        c = Xtild.T @ f
        if np.linalg.norm(c)>0:
            c /= np.linalg.norm(c)
        
        ### compute objective value, if fixed, stop
        S = f @ Xtild @ c
        if S==Sold: ### will terminate once fOpt is fixed over 2 iterations
            break
        Sold = S
        i+=1
        if i>50:
            c = np.zeros_like(c)
            f=np.zeros_like(f)
            S=0
            break
    return c,f,np.abs(S)


### find locally-optimal (unconstrained) c and f from spectral c initialization
def generateSpectralOptcf(X,tblShape=-1):
    np.random.seed(0)### for filling in f

    ### for backwards comaptability 
    if tblShape==-1:
        tblShape=X.shape
    
    relevantTargs = X.sum(axis=1)>0
    relevantSamples = X.sum(axis=0)>0

    r,c=X.shape
    
    if relevantTargs.sum()<2 or relevantSamples.sum()<2:
        return np.zeros(c),np.zeros(r)

    X = X[np.ix_(relevantTargs,relevantSamples)]

    
    ###### Spectral initialization of c
    ### pairwise similarity matrix
    A = X.T @ X
    A = A-np.diag(np.diag(A))
    
    ### remove isolated elements
    sampleMask2= (A.sum(axis=0)>0)
    
    if sampleMask2.sum()==0: ### graph is not connected (I think)
        sampleMask2 = np.ones(X.shape[1],dtype='bool')
        c = np.random.choice([-1,1],size=X.shape[1])
    else:
        A = A[np.ix_(sampleMask2,sampleMask2)]
        X = X[:,sampleMask2]

        ## construct diagonals for Laplacian
        D = np.diag(A.sum(axis=0))

        ### spectrally normalized Laplacian
        Dnorm = np.diag(A.sum(axis=0)**(-1/2))
        L = np.eye(A.shape[0]) - Dnorm@A@Dnorm

        ### normal
        # L=D-A

        ### potentially could merge results from 2 clusterings?

        ### compute fiedler vector
        eigvals,eigvecs =np.linalg.eig(L)
        c = np.sign(normalizevec(np.real(eigvecs[:,np.argsort(np.abs(eigvals))[1]])))
        ## maybe more fancy checking needed if graph isn't connected
    
    #### if clustering put all samples in same cluster, shouldn't happen
    if np.all(c==c[0]):
        c[0] = -1*c[0]
    
    c,fOpt,S=altMaximize(X,c)
    
    ## extend to targets and samples that didn't occur in training data
    fElong = np.random.choice([0,1],size=tblShape[0])
    fElong[relevantTargs] = fOpt
    fOpt = fElong
    
    cElong = np.zeros(tblShape[1])
    cElong[np.arange(tblShape[1])[relevantSamples][sampleMask2]]=c ### fancy indexing
    cOpt = cElong
    
    return np.nan_to_num(cOpt,0),np.nan_to_num(normalizevec(fOpt),0)


### find locally-optimal (unconstrained) c and f from random initialization
### v2
def generate_alt_max_cf(X,tblShape=-1, randSeed=0,numRandInits=10):
    np.random.seed(randSeed) ### random initialization and extension
    
    ### for backwards comaptability 
    if tblShape==-1:
        tblShape=X.shape

    relevantTargs = X.sum(axis=1)>0
    relevantSamples = X.sum(axis=0)>0

    nrows,ncols=X.shape
    
    if relevantTargs.sum()<2 or relevantSamples.sum()<2:
        return np.zeros(ncols),np.zeros(nrows)

    X = X[np.ix_(relevantTargs,relevantSamples)]
    
    Sbase=0
    fMax=0
    cMax=0
    for _ in range(numRandInits):
        c = np.random.choice([-1,1],size=X.shape[1])
        c,f,S = altMaximize(X,c)
        if S > Sbase:
            fMax = f
            cMax = c
            Sbase = S
            

    ## extend to targets and samples that didn't occur previously
    fElong = np.random.choice([0,1],size=nrows)
    fElong[relevantTargs] = fMax
    fOpt = fElong
    
    cElong = np.zeros(ncols)
    cElong[np.arange(ncols)[relevantSamples]]=cMax ### fancy indexing
    cOpt = cElong
    
    return cOpt,fOpt

### find locally-optimal c and f from random initialization
#### constrained to c having the same sign as sheetCj
def generateSignedSheetCjOptcf(X,sheetCj,tblShape = -1):
    np.random.seed(0)### for filling in f
    
    ### for backwards comaptability 
    if tblShape==-1:
        tblShape=X.shape


    relevantTargs = X.sum(axis=1)>0
    relevantSamples = X.sum(axis=0)>0

    X = X[np.ix_(relevantTargs,relevantSamples)]

    cSign = sheetCj[relevantSamples]
    c = cSign

    nj = X.sum(axis=0)
    njinvSqrt = 1.0/np.maximum(1,np.sqrt(nj)) ### avoid divide by 0 errors
    njinvSqrt[nj==0]=0
    
    fOpt = np.sign((X-X@np.outer(np.ones(X.shape[1]),nj)/np.maximum(1,X.sum()))@(c*njinvSqrt))
    fOpt = (fOpt+1)/2 ### to rescale f to be [0,1] valued
    
    Sold = 0
    i=0
    while True:
        
        ### find optimal f for fixed c
        fOpt = np.sign((X-X@np.outer(np.ones(X.shape[1]),nj)/np.maximum(1,X.sum()))@(c*njinvSqrt))
        fOpt = (fOpt+1)/2 ### to rescale f to be [0,1] valued
        
        ### find optimal c for fixed f
        Sj = np.multiply(fOpt @ (X-X@np.outer(np.ones(X.shape[1]),nj)/X.sum()),njinvSqrt)
        ### construct plus and minus variants of Sj
        cplus = Sj * ((Sj*cSign) >0)
        cplus /= np.maximum(1,np.linalg.norm(cplus,2))
        c=cplus
        Splus = fOpt @ (X-X@np.outer(np.ones(X.shape[1]),nj)/X.sum())@(c*njinvSqrt)/np.maximum(1,np.linalg.norm(c,2))
        
        cminus = Sj * ((Sj*cSign) <0)
        cminus /= np.maximum(1,np.linalg.norm(cminus,2))
        c=cminus
        Sminus = fOpt @ (X-X@np.outer(np.ones(X.shape[1]),nj)/X.sum())@(c*njinvSqrt)/np.maximum(1,np.linalg.norm(c,2))

        if Splus >= -1*Sminus:
            c = cplus
            S = Splus
        else:
            c = cminus
            S = Sminus

        if S==Sold: ### will terminate once fOpt is fixed over 2 iterations
            break
        Sold = S
        i+=1
        if i>50:
            c = np.zeros_like(c)
            fOpt=np.zeros_like(fOpt)
            break
            
    ## extend to targets and samples that didn't occur previously
    fElong = np.random.choice([0,1],size=tblShape[0])
    fElong[relevantTargs] = fOpt
    fOpt = fElong
    
    cElong = np.zeros(tblShape[1])
    cElong[np.arange(tblShape[1])[relevantSamples]]=c ### fancy indexing
    cOpt = cElong
    return cOpt,fOpt


### find optimal f for given input sheetCj
def generateSheetCjOptcf(X,sheetCj,tblShape=-1):
    np.random.seed(0) ### for filling in f

    ### for backwards comaptability 
    if tblShape==-1:
        tblShape=X.shape
    
    relevantTargs = X.sum(axis=1)>0
    X = X[relevantTargs]

    ## set cj as sheetCj
    c = sheetCj

    nj = X.sum(axis=0)
    njinvSqrt = 1.0/np.maximum(1,np.sqrt(nj)) ### avoid divide by 0 errors
    njinvSqrt[nj==0]=0
    
    ## compute opt f
    fOpt = np.sign((X-X@np.outer(np.ones(X.shape[1]),nj)/np.maximum(1,X.sum()))@(c*njinvSqrt))
    fOpt = (fOpt+1)/2 ### to rescale f to be [0,1] valued
    
    ## extend to targets and samples that didn't occur previously
    fElong = np.random.choice([0,1],size=tblShape[0])
    fElong[relevantTargs] = fOpt
    fOpt = fElong
    
    cOpt = c ### sheetCj
    return cOpt,fOpt



######## 3. p-value computation

### asymptotic pvalue
### NOMAD is asymptotically normal with variance upper bounded by totalVar under the null
###   allowing is to provide asymptotically valid p-values
def computeAsympNOMAD(X,cOpt,fOpt):
    if (cOpt==0).all():
        return 1
    
    cOpt = np.nan_to_num(cOpt,0)
    fOpt = np.nan_to_num(fOpt,0)
    if not (fOpt.max()<=1) and (fOpt.min()>=0):
        fOpt = normalizevec(np.nan_to_num(fOpt,0))
    
    ### only requirement for valid p-value
    assert((fOpt.max()<=1) and (fOpt.min()>=0))
    nj = X.sum(axis=0)
    njinvSqrt = 1.0/np.maximum(1,np.sqrt(nj))
    njinvSqrt[nj==0]=0
    
    ### compute p value
    S = fOpt @ (X-X@np.outer(np.ones(X.shape[1]),nj)/X.sum())@(cOpt*njinvSqrt)

    if S<1E-10: ### prevent edge issues
        return 1

    M=X.sum()
    
    muhat = (fOpt@X).sum()/M
    
    varF = (fOpt-muhat)**2 @ X.sum(axis=1)/X.sum()
    totalVar = varF * (np.linalg.norm(cOpt)**2 - (cOpt@np.sqrt(nj))**2/M)
    pval = 2*np.exp(-S**2/(2*totalVar ))
                
    return min(np.nan_to_num(pval,1),1)


### test p-value for fixed c and f on contingency table
def testPval(X,cOpt,fOpt):
    
    cOpt = np.nan_to_num(cOpt,0)
    fOpt = np.nan_to_num(fOpt,0)

    if np.all(cOpt == cOpt[0]) or np.all(fOpt==fOpt[0]):
        return 1

    if not (fOpt.max()<=1) and (fOpt.min()>=0):
        fOpt = normalizevec(fOpt,0,1)
    
    ### only requirement for valid p-value
    assert((fOpt.max()<=1) and (fOpt.min()>=0))
    nj = X.sum(axis=0)
    njinvSqrt = 1.0/np.maximum(1,np.sqrt(nj))
    njinvSqrt[nj==0]=0
    
    ### compute test statistic
    S = fOpt @ (X-X@np.outer(np.ones(X.shape[1]),nj)/X.sum())@(cOpt*njinvSqrt)

    M=X.sum()
    
    denom = (np.linalg.norm(cOpt)**2 - (cOpt@np.sqrt(nj))**2/M)
    pval = 2*np.exp(-2*S**2/denom)
                
    # if np.isnan(pval):
    #     print(S,denom,cOpt, fOpt@X/ nj)

    return min(np.nan_to_num(pval,1),1)

### compute chi2 pvalue
def computeChi2Test(X):
    if len(X.shape)==1:
        return 1
    X = X[X.sum(axis=1)>0]
    X = X[:,X.sum(axis=0)>0]
    _,pv,_,_= scipy.stats.contingency.chi2_contingency(X)
    return pv

### compute LRT against full null
def computeLRT_Test(X):
    if len(X.shape)==1:
        return 1
    X = X[X.sum(axis=1)>0]
    X = X[:,X.sum(axis=0)>0]
    _,pv,_,_ = scipy.stats.contingency.chi2_contingency(X, lambda_="log-likelihood")
    return pv

###### effect size computation

### standard effect size definition, for +/-1 (or, binarized) c
def effectSize_bin(X,c,f):
    if (c>0).sum()==0 or (c<0).sum()==0:
        return 0

    return np.abs(f@X@(c>0) / (X@(c>0)).sum() - f@X@(c<0) / (X@(c<0)).sum())


### new effect size definition, for continuous c
def effectSize_cts(X,c,f):
    if (X@np.abs(c)).sum() ==0:
        return 0

    return np.abs(f@X@c / (X@np.abs(c)).sum())



######## 4: table metadata computation

#### compute mean entropy difference to detect V(D)J type events
def simple_ent_dif(X):
    X = X[:,X.sum(axis=0)>0] ### avoid 0 count samples
    normalizedMat = X/X.sum(axis=0)
    agg_ent = scipy.stats.entropy(normalizedMat.mean(axis=1),base=2)
    meanEnt = np.mean([scipy.stats.entropy(normalizedMat[:,j],base=2) for j in range(X.shape[1])])
    return agg_ent-meanEnt

### compute hamming distance between two vectors
def hamming(x,y):
    return sum(c1 != c2 for c1, c2 in zip(x,y))

### takes as input dataframe where index is sequences (e.g. targets, compactors)
### this is slow for levenshtein: can be sped up by approximating the solution by ignoring infrequent targets
def computeAverageDist(df,distfn):
    cont_table = df.to_numpy()
    targCts = cont_table.sum(axis=1)
    targMax = targCts.argmax()
    targLst = df.index.to_numpy()
    distances = np.vectorize(lambda x,y : distfn(x,y))(targLst,targLst[targMax])
    return np.dot(distances,targCts)/targCts.sum()

### compute many base quantities
def computeBaseQuantities(X):
    newRow = {}

    newRow['M'] = X.sum()
    newRow['anch_uniqTargs']=(X.sum(axis=1)>0).sum()
    newRow['number_nonzero_samples'] = (X.sum(axis=0)>0).sum()
    newRow['target_entropy'] = scipy.stats.entropy(X.sum(axis=1),base=2)
    newRow['entropy_difference'] = simple_ent_dif(X)

    return newRow




