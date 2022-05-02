#!/usr/bin/env Rscript

library(data.table)
library(spgs)
library(stringdist)
print('yes!')
return_d<-function(targets, dmethod, MAXD){

    # inputs a list of strings in sorted order
    n.tot=length(targets)
    ds= c(0,rep(MAXD,n.tot-1))

    # run through all needed to be computed
    for (i in 2:(n.tot)){
        for (j in 1:(i-1) ){
            ## define two target strings
            si=strsplit(targets[i], "")[[1]] #string i
            sj=strsplit(targets[j], "")[[1]] #string j

            ## If HAMMING
            if (dmethod=="hamming"){
                ds[i] = min(MAXD ,min(ds[i], sum( si != sj )))
            }
            ### If Jaccard
            if (dmethod =="jaccard"){
                ds[j] = min(MAXD, min(ds[i], stringdist(si, sj,method="jaccard",q=6)))
            }
        }
    }

    ## RETURNS
    ds
}

args <- commandArgs(TRUE)
infile <- args[1]
samplesheet <- args[2]
distance_type <- args[3]
max_targets <- as.integer(args[4])
max_dist <- as.integer(args[5])
bonfer <- as.integer(args[6])
pval_threshold <- as.numeric(args[7])
outfile_scores <- args[8]
outfile_anchors <- args[9]


samplesheet <- fread(samplesheet, header=F)

if (ncol(samplesheet) == 1) {
    ## if there is only one column of files, use_c is set to false
    use_c = FALSE

} else if (ncol(samplesheet) == 2) {
    ## if there is a second column, use_c is set to true
    use_c = TRUE

    colnames(samplesheet) <- c('sample', 'cs')
    ## create sample column for joining later
    samplesheet[, sample := sub('.fastq.gz', '', basename(sample))]
}

m = fread(infile, fill=T, header=F)
names(m) = c("counts", "seq", "sample")

m[ ,anchor := substr(seq,1,27)]
m[ ,target := substr(seq,28,54)]
m[ ,n.targ := length(unique(target)),by=anchor]
m = m[n.targ>1] # JUST COMPUTE FOR MORE THAN 1 TARGET
m[ ,target.count := sum(counts), by=anchor] # gives the count per target for each anchor

i = 0
for (anch in unique(m$anchor)){ ## INEFFICIENT can be done with lapply
    i = i + 1

    ll=length(unique(m[anchor== anch][order(-counts)]$target)) # TOTAL NUM TARGETS

    ## GENERATES LIST OF TARGETS ONLY FOR THE TOP max_targets
    tlist= unique(m[anchor== anch][order(-counts)]$target) [1: min(ll,max_targets)]

    #RETURN HAMMING "hamming" or JACCARD "jaccard" FROM THE unique list
    target.d = return_d(tlist, distance_type, max_dist)

    #CREATE A NEW DATAFRAME AND ADD ANCHOR ID
    into= data.table(cbind (tlist,as.numeric(target.d))); names(into) = c("target","target.d")
    into[,anchor:=anch]

    ### LOGIC TO EITHER GROW A DATAFRAME OR INITIALIZE IT FOR LATER MERGE
    if ( sum(ls() %like% "tomerge.distances")>0){ tomerge.distances=rbind(tomerge.distances,into)}
    if ( sum(ls() %like% "tomerge.distances")==0){ tomerge.distances= into}  ## keeps a
}

######################################################
## CREATE THE KEY TO MERGE TARGET DISTANCES ON ON
tomerge.distances[,seq:=paste(tomerge.distances$anchor,tomerge.distances$target,sep="")]

### MERGES ON DISTANCE  JUST ADD THE 2 fields IN THE LIST, MATCHED ON SEQ
new = merge(m,tomerge.distances[,list(seq,target.d)],by = "seq")

#################### JUST COMPUTING THE LINEAR STATISTICS -- PLS CHECK-- BASED ON VECTOR OF CS BEING SET TO 1 --NEEDS ADJUSTMENT PER FORMULAS
new[ ,n_j:=sum(counts),by=list(anchor,sample)] # counts up all the targets for this anchor,sample
new[ ,mu:=sum(counts* as.numeric(target.d))/sum(counts), by="anchor"]
new[ ,target.d :=  as.numeric(target.d)]
#### CHECK-- eg i'm not accounting for missing values
new[ ,score_per_sample:=sum((counts*(target.d-mu))/sqrt(n_j) ), by=list(sample,anchor) ] # if this is the way you are computing it
new[ ,num.sample:= length(unique(sample)), by=anchor]
#  M is total counts for the anchor
new[ ,M:=sum(counts),by=anchor]

## create a new list of cs by sample by anchor has all the data
compute.a = unique( new[,list(anchor, sample, M, score_per_sample, n_j, mu)])
## define new c per anchor then mergein


if (use_c == FALSE){

    pv=matrix(nrow=dim(compute.a)[1],ncol=bonfer)

    for (j in 1:bonfer){ # 5 projections of cs
        cs=runif(dim(compute.a)[1],-1,1) # generates newcs in the next line by first removing current cs
        compute.a[,cs:=NULL]; compute.a=cbind(compute.a,cs)
        compute.a[,anchor_score :=  sum(cs*(score_per_sample)),by=anchor] #
        ## CONSTANTS NEEDED FOR BOTH TERMS
        compute.a[,a:= 1/(1+sqrt(M*sum(cs^2) /(sum(cs*sqrt(n_j))^2 ))) ,by=anchor]
        compute.a[,sumcsq:=sum(cs^2),by=anchor]
        compute.a[,sumc_sqrtnj := sum( cs*sqrt(n_j)),by=anchor]
        ## since anchorscore is squared, no need for abs value
        compute.a[, pterm1:= 2*(exp(- (2*(1-a)^2*anchor_score^2/(max_dist^2*sumcsq))))]
        compute.a[, pterm2:= 2*(exp(- ((2*a^2*anchor_score^2*M)/(max_dist^2 *sumc_sqrtnj^2)))) ]

        pv[,j]= bonfer*(compute.a$pterm1+compute.a$pterm2)
    }

    ## take the min across the 5 already BH-corrected values
    bf.cor.p=apply(pv, 1, min)
    compute.a=cbind(compute.a, bf.cor.p)

} else {

    compute.a=merge(compute.a, samplesheet, by="sample",all.x=T)# OUTER JOIN, assign 0 to cs 0
    compute.a=compute.a[is.na(cs),cs:=0]

    compute.a[,anchor_score :=  sum(cs*(score_per_sample)),by=anchor] #
    ## CONSTANTS NEEDED FOR BOTH TERMS
    compute.a[,a:= 1/(1+sqrt(M*sum(cs^2) /(sum(cs*sqrt(n_j))^2 ))) ,by=anchor]
    compute.a[,sumcsq:=sum(cs^2),by=anchor]
    compute.a[,sumc_sqrtnj := sum( cs*sqrt(n_j)),by=anchor]
    ## since anchorscore is squared, no need for abs value
    compute.a[, pterm1:= 2*(exp(- (2*(1-a)^2*anchor_score^2/(max_dist^2*sumcsq))))]
    compute.a[, pterm2:= 2*(exp(- ((2*a^2*anchor_score^2*M)/(max_dist^2 *sumc_sqrtnj^2)))) ]
    compute.a[,bf.cor.p:=pterm1+pterm2] ## NO BF CORRECTION
}


# anchors = compute.a[pval<pval_threshold]
anchors = head(compute.a[order(compute.a$bf.cor.p, decreasing=T), "anchor"], 1000)

write.table(compute.a, outfile_scores, col.names=T, row.names=F, quote=F, sep='\t')
write.table(anchors, outfile_anchors, col.names=F, row.names=F, quote=F, sep='\t')
