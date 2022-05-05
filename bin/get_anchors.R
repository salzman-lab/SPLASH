#!/usr/bin/env Rscript

library(data.table)
library(spgs)
library(stringdist)
library(R.utils)


args <- commandArgs(TRUE)
infile <- args[1]
samplesheet <- args[2]
distance_type <- args[3]
max_targets <- as.integer(args[4])
max_distance <- as.integer(args[5])
bonfer <- as.integer(args[6])
pval_threshold <- as.numeric(args[7])
run_type <- args[8]
outfile_scores <- args[9]
outfile_anchors <- args[10]

# infile="/oak/stanford/groups/horence/kaitlin/results_stringstats/SS2_test/merged_abundant_anchor_counts.tsv.gz"
# samplesheet="/scratch/groups/horence/kaitlin/running_stringstats/samplesheets/samplesheet_SS2_test.csv"
# distance_type="lev"
# max_targets=4
# max_distance=10
# bonfer=10
# pval_threshold=0.1
# outfile_scores="scores_merged_abundant_anchor_counts.tsv"
# outfile_anchors="anchors_merged_abundant_anchor_counts.tsv"
# chunk_size=5


## Function to return the minimum distances in an ordered list of strings
get_distances <- function(targets, distance_type, max_distance, chunk_size){
    n.tot = length(targets)
    ds= c(0, rep(max_distance, n.tot-1))
    for (i in 2:(n.tot)){ ## run through all needed to be computed
        j=1
        if (distance_type == "hamming"){
            ds[i] = min(
                max_distance,
                min (ds[i], floor(stringdist(  targets[i],targets[j],method="hamming") / chunk_size))
            )
        }
        if (distance_type == "lev"){
                ds[i] = min(
                max_distance,
                min(ds[i], stringdist(targets[i], targets[j], method="lv"))
            )
        }
    }
    ## Return distances
    ds
}

## Function to return ____
get_cs <- function(use_c, in_matrix, dimensions, samplesheet) {
    if (use_c == FALSE) {
        ## repeat this dimensions times or take the first dimensions vectors of U; if cs are userdef, dimensions=1
        nsubs = min(dim(in_matrix)[1], 20000)
        myinds =sample(dim(in_matrix)[1], nsubs)
        inp = in_matrix[order(-anchor.var)][myinds]
        print("starting reshape")
        x = reshape(in_matrix[myinds, list(anchor, sample, score_per_sample)], timevar="anchor", idvar="sample", direction="wide")

        ## remove the first row whish is samples, remove missing values, take the svd
        print ("STARTING SVD")
        sample = x[,1]
        x=x[,-1]
        x[is.na(x)] = 0

        ## IDEALLY WILL ENFORCE SOMETHING such as CS are +/- a constant or 0.
        xsvd = svd(x)$u/abs(svd(x)$u)
        dimensions = min(dimensions,dim(xsvd)[2])

        cs = xsvd[,1:dimensions]
        out = cbind(sample,cs)

    } else {

        out = samplesheet

    }

    ## return matrix with column ofcs
    out
}


chunk_size = 4

## Function to return the minimum distances in an ordered list of strings
get_distances <- function(targets, distance_type, max_distance, chunk_size){
    n.tot = length(targets)
    ds= c(0, rep(max_distance, n.tot-1))
    for (i in 2:(n.tot)){ ## run through all needed to be computed
        j=1
        if (distance_type == "hamming"){
            ds[i] = min(max_distance, min ( ds[i], floor ( stringdist(  targets[i],targets[j],method="hamming") / chunk_size)))

        }
        if (distance_type == "lev"){
	           ds[i] = min(
                max_distance,
                min(ds[i], stringdist(targets[i], targets[j], method="lv"))
            )
        }
    }
    ## Return distances
    ds
}

## Function to return ____
get_cs <- function(use_c, in_matrix, dimensions, samplesheet) {
    if (use_c == FALSE) {
        ## repeat this dimensions times or take the first dimensions vectors of U; if cs are userdef, dimensions=1
        nsubs = min(dim(in_matrix)[1], 20000)
        myinds =sample(dim(in_matrix)[1], nsubs)
        inp = in_matrix[order(-anchor.var)][myinds]
        print("starting reshape")
        x = reshape(in_matrix[myinds, list(anchor, sample, score_per_sample)], timevar="anchor", idvar="sample", direction="wide")

        ## remove the first row whish is samples, remove missing values, take the svd
        print ("STARTING SVD")
        sample = x[,1]
        x=x[,-1]
        x[is.na(x)] = 0

        ## IDEALLY WILL ENFORCE SOMETHING such as CS are +/- a constant or 0.
        xsvd = svd(x)$u/abs(svd(x)$u)
        dimensions = min(dimensions,dim(xsvd)[2])

        cs = xsvd[,1:dimensions]
        out = cbind(sample,cs)

    } else {
      out = samplesheet#  merge(in_matrix, samplesheet, by="sample", all.x=T)

    }

    ## return matrix with column ofcs
    out
}


chunk_size = 4

samplesheet <- fread(samplesheet, header=F)

if (ncol(samplesheet) == 1) {
    ## if there is only one column of files, use_c is set to false
    use_c = FALSE
    samplesheet = data.table()

} else if (ncol(samplesheet) == 2) {
    ## if there is a second column, use_c is set to true
    use_c = TRUE

    colnames(samplesheet) <- c('sample', 'cs')
    ## create sample column for joining later
    samplesheet[, sample := sub('.fastq.gz', '', basename(sample))]
}

m = fread(infile, fill=T, header=F)
print("READ IN")
names(m) = c("counts", "seq", "sample")

m[ ,anchor := substr(seq,1,27)]
m[ ,target := substr(seq,28,54)]

## Only continue computation for anchors with more than one target
m[ ,n.targ := length(unique(target)),by=anchor]
m = m[n.targ>1]
m[,M:=sum(counts),by=anchor]
# if impatient m=m[M>100]
## gives the count per target for each anchor
m[ ,target.count := sum(counts), by=anchor]
if (dim(m)[1]>0){
    ## gives the count per target for each anchor
    m[ ,target.count := sum(counts), by=anchor]
    ## get per target counts
    m[ ,pertarget.counts := sum(counts), by=list(target,anchor)]
    print("STARTING GET DISTANCES");
    for (anch in unique(m$anchor)){
        ## gives the total number of targets
        ll = length(unique(m[anchor== anch][order(-counts)]$target))

        ## generates list of top targets for targets
	    ## changed
        tlist= unique(m[anchor== anch][order(-pertarget.counts)]$target)[1: min(ll, max_targets)]
        ## gets list of distances from unique list of targets

        target.d = get_distances(tlist, distance_type, max_distance,  chunk_size)
 if (anch %like% "GGG"){print(anch)}
        ## CREATE A NEW DATAFRAME AND ADD ANCHOR ID
        into = data.table(cbind (tlist,as.numeric(target.d)))
        names(into) = c("target", "target.d")
        into[,anchor := anch]

    ## LOGIC TO EITHER GROW A DATAFRAME m with its TARGET DISTANCES OR INITIALIZE IT FOR LATER MERGE

    into[ ,seq := paste(into$anchor,into$target,sep="")]
    if ( sum(ls() %like% "tomerge.distances")>0){ tomerge.distances=rbind(tomerge.distances,into)}
    if ( sum(ls() %like% "tomerge.distances")==0){ tomerge.distances= into}
    }


    print("ENDING GET DISTANCES")
    ## CREATE THE KEY TO MERGE TARGET DISTANCES ON ON

    ## MERGES ON DISTANCE  JUST ADD THE 2 fields IN THE LIST, MATCHED ON SEQ
    new = merge(m, tomerge.distances[,list(seq,target.d)], by = "seq")

    ## Computing linear statistics
    ## counts up all the targets for this anchor,sample
    new[ ,n_j :=sum(counts), by=list(anchor, sample)]
    new[ ,mu := sum(counts * as.numeric(target.d)) / sum(counts), by="anchor"]
    new[ ,anchor.var := sum(counts * (as.numeric(target.d)^2))/sum(counts)-mu^2, by="anchor"]
    new[ ,target.d :=  as.numeric(target.d)]

    ## compute the per sample score, and number of samples
    new[ ,score_per_sample := sum((counts*(target.d-mu))/sqrt(n_j)) , by=list(sample,anchor) ]
    new[ ,num.sample:= length(unique(sample)), by=anchor]

    ## all is the variable with all the statistical information
    if(sum(ls() %like% "all")>0){ all=rbind(all,new)}
    if(sum(ls() %like% "all")==0){ all= new}

    all = unique(all)
    ## reduces the matrix to unique anchor, samples with associated

    compute.a = unique( all[ ,list(anchor,sample, M, mu, score_per_sample, num.sample, n_j, anchor.var)])
    print("moving to get cs")
}

## add column of cs
c.mx = get_cs(  use_c=0, compute.a, bonfer, samplesheet)

## in case we ask for more tests than dimensions of the mx, set ot the min
bonfer = min(bonfer, (dim(c.mx)[2]-1))

## matrix of pvalues
pv = matrix(nrow=dim(compute.a)[1], ncol=bonfer)

## start bonferroni correction
for (j in 1:bonfer){ # bonfer is number of projections of cs
    print(paste("starting ",j,"th bonfer"))
    ## merge in new cs
    newc = data.table(data.frame(c.mx)[,c(1,(j+1))])
    ## creates a temp matrix of sample, col1, and Jth c vector
    names(newc)[1]="sample"
    names(newc)[2]="cs"

    ## if a col of cs exists, remove
    if (sum(names(compute.a) == "cs")>0){compute.a[,cs:=NULL]}
    ## added all.x=T
    compute.a = merge(unique(compute.a), unique(newc),all.x=T, by="sample") ## MERGES THIS C into the
    # 2 added lines in case cs are not presen
    compute.a[is.na(cs),cs:=0]
    compute.a[cs==0, n_j:=0]
    compute.a[ ,M:=sum(n_j), by=anchor] ## new
    compute.a[ ,num.sample := length(unique(sample)), by=anchor]
    compute.a[ ,anchor_score := sum(cs*(score_per_sample)),by=anchor]

    ## constants needed for pvalue calculation
    compute.a[ ,a:= 1/(1+sqrt(M*sum(cs^2) /(sum(cs*sqrt(n_j))^2 ))), by=anchor]
    compute.a[ ,sumcsq:=sum(cs^2), by=anchor]
    compute.a[ ,sumc_sqrtnj := sum(cs*sqrt(n_j)), by=anchor]

    ## since anchorscore is squared, no need for abs value
    compute.a[, pterm1:= 2*(exp(-(2*(1-a)^2*anchor_score^2/(max_distance^2*sumcsq))))]
    compute.a[, pterm2:= 2*(exp(-((2*a^2*anchor_score^2*M)/(max_distance^2 *sumc_sqrtnj^2))))]

    pv[ ,j] = bonfer * (compute.a$pterm1+compute.a$pterm2)

}

## take the min across the 5 already BH-corrected values
bf.cor.p = apply(pv,1,min)
compute.a = cbind(compute.a, bf.cor.p)

## L1 is ballpark calculation for reference for now, likely will be removed.
compute.a[ ,l1 := sum(abs(score_per_sample)), by=anchor]

## add significance
compute.a[ ,l1 := sum(abs(score_per_sample)), by=anchor]

## units of variance * num samples -- like "sds"
compute.a[ ,l1.sdlike.units := l1 / (anchor.var*num.sample)]
compute.a = (compute.a[anchor.var>.5][mu>1][order(-l1.sdlike.units)])

## write out anchor scores
summary.file = unique(compute.a[order(-l1.sdlike.units)])
write.table(summary.file, outfile_scores, col.names=T, row.names=F, quote=F, sep='\t')

print(head(summary.file))

## write out significant anchors
anchors = summary.file[bf.cor.p < pval_threshold]
anchors = head(unique(summary.file[order(bf.cor.p, decreasing=T), "anchor"]), 1000)

write.table(anchors, outfile_anchors, col.names=F, row.names=F, quote=F, sep='\t')
write.table(file=paste("cmx_",outfile_anchors,sep=""), c.mx, col.names=F, row.names=F, quote=F, sep='\t')

write.table(file=paste("cmx_",outfile_anchors,sep=""), c.mx, col.names=F, row.names=F, quote=F, sep='\t')
