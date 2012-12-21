setMethod("sort", signature("KcghData"), function(x, decreasing=FALSE, ...){
    dataOrdering <- order(x@probeAnnotation@chromosome, x@probeAnnotation@maploc)
    x@probeAnnotation <- x@probeAnnotation[dataOrdering]
    x@data <- x@data[dataOrdering,]
    x
})

setMethod("sort", signature("KcghDataSum"), function(x, decreasing=FALSE, ...){
    dataOrdering <- order(x@probeAnnotation@chromosome, x@probeAnnotation@maploc)
    x@probeAnnotation <- x@probeAnnotation[dataOrdering]
    x@pos <- x@pos[dataOrdering]
    x@neg <- x@neg[dataOrdering]
    x
})
