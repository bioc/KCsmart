setClass("probeAnnotation", representation(chromosome="character", maploc="integer", name="character"))
setClass("KcghData", representation(probeAnnotation="probeAnnotation", data="matrix"))
setClass("KcghDataSplit", representation(probeAnnotation="probeAnnotation", pos="matrix", neg="matrix"))
setClass("KcghDataSum", representation(probeAnnotation="probeAnnotation", pos="vector", neg="vector", nrSamples="integer"))
setClass("KcghDataMirror", representation(mirrorLocs="list"), contains="KcghDataSum")
setClass("KCData", representation(data="list"))
setClass("samplePointMatrix", representation(totalLength="integer", maxy="numeric", miny="numeric", sampleDensity="integer", sigma="integer", mirrorLocs="list", probeAnnotation="probeAnnotation"), contains="KCData")
setClass("sigSegments", representation(gains="list", losses="list", sigma="integer", sigLevels="list"))
setClass("sigRegions", representation(sigma="integer"), contains="KCData")
setClass("scaleSpace", contains="KCData")

setClass("spmCollection", representation(data="matrix", annotation="probeAnnotation", cl="vector", mirrorLocs="list", sampleDensity="integer", sigma="integer"))
setClass("snrResult", representation(snrValues="numeric", fudge="numeric", permutations="matrix"))
setClass("compKc", representation(spmCollection="spmCollection", method="character", siggenesResult="SAM", snrResult="snrResult"))
setClass("compKcSigRegions", representation(regionTable="data.frame", method="character", cutoff="numeric", fdr="numeric"))

setMethod("initialize", "KcghData", function(.Object, cghData){
    .Object@probeAnnotation <- new("probeAnnotation", cghData$chrom, cghData$maploc, row.names(cghData))
    .Object@data <- as.matrix(cghData[,3:ncol(cghData)])
    .Object
})

setMethod("initialize", "probeAnnotation", function(.Object, chromosome=NULL, maploc=NULL, name=NULL){
    .Object@chromosome <- as.character(chromosome)
    .Object@maploc <- as.integer(maploc)
    .Object@name <- as.character(name)
    chromNA <- which(is.na(.Object@chromosome))
    maplocNA <- which(is.na(.Object@maploc))
    if(length(chromNA) > 0){
	stop(paste('There is no chromosome annotation for probe',chromNA, ', please remove probe from dataset or annotate'))
    }
    if(length(maplocNA) > 0){
	stop(paste('There is no maploc annotation for probe', maplocNA, ', please remove probe from dataset or annotate'))
    }
    if(length(.Object@maploc) != length(.Object@chromosome)){
	stop('The number of chromosome and maploc annotations are not equal')
    }
    .Object
})

setMethod("initialize", signature("KcghDataSplit"), function(.Object, KcghData){
    if(!is(KcghData,"KcghData")){stop("Need KcghData object as input")}

    pos <- KcghData@data
    neg <- pos
    
    pos.mask <- pos < 0

    pos[pos.mask] <- 0
    neg[!pos.mask] <- 0
    
    .Object@probeAnnotation <- KcghData@probeAnnotation
      
    .Object@pos <- pos
    .Object@neg <- neg
    .Object
})


setMethod("initialize", signature("KcghDataSum"), function(.Object, KcghDataSplit=NULL){
    if(is.null(KcghDataSplit)){
        .Object@probeAnnotation <- new("probeAnnotation")
        return(.Object)
    }
    if(!is(KcghDataSplit,"KcghDataSplit")){stop("Need KcghDataSplit object as input")}

    .Object@probeAnnotation <- KcghDataSplit@probeAnnotation
    
    .Object@pos <- rowSums(KcghDataSplit@pos)
    .Object@neg <- rowSums(KcghDataSplit@neg)
	
    .Object@nrSamples <- ncol(KcghDataSplit@pos)
    .Object
})

setMethod("initialize", signature("KcghDataMirror"), function(.Object, mirrorLocs, KcghDataSum=NULL){
    if(is.null(KcghDataSum)){
        .Object@probeAnnotation <- new("probeAnnotation")
        return(.Object)
    }
    if(!is(KcghDataSum,"KcghDataSum")){stop("Need KcghDataSum object as input")}

    .Object@probeAnnotation <- KcghDataSum@probeAnnotation
    
    .Object@pos <- KcghDataSum@pos
    .Object@neg <- KcghDataSum@neg
    .Object@nrSamples <- KcghDataSum@nrSamples
    .Object@mirrorLocs <- mirrorLocs

    .Object
})



setMethod("initialize", signature("compKc"), function(.Object, spmCollection, method=c("siggenes", "perm"), comparedata) {
	
	method <- match.arg(method)

	if(!is(spmCollection,"spmCollection")){stop("Need spmCollection object as input")}

	.Object@spmCollection <- spmCollection
	
	if(is(comparedata,"snrResult"))
		.Object@snrResult <- comparedata
	else if (is(comparedata,"SAM"))
		.Object@siggenesResult <- comparedata
	else
		stop("comparedata needs to be either SAM or snrResult")

	
	.Object@method <- method
	.Object	
})


