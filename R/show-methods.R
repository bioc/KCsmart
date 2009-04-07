setMethod("show", signature("samplePointMatrix"), function(object){
	cat("Sample point matrix\n")
	cat(paste('Contains', length(object@data),"chromosomes\n"))
	cat(paste("The chromosome start, end and centromere locations of", attr(object@mirrorLocs, 'Version'),"were used\n"))
	cat(paste("The original data contained", length(object@probeAnnotation@maploc), "probes\n\n"))
	cat("The following parameters have been used:\n")
	cat(paste("Sigma:", object@sigma), "\n")
	cat(paste("Sample density:", object@sampleDensity, "\n\n"))
	cat("Use 'plot' to display the results visually\n")
})


setMethod("show", signature("sigSegments"), function(object){
	cat("Significantly gained and lost segments\n")
	cat(paste("Cutoff value was", object@sigLevels$pos, "for gains and", object@sigLevels$neg, "for losses\n"))
	cat(paste("A total of", length(object@gains),"gains and", length(object@losses),"losses were detected\n\n"))
	cat("Use 'write.table' to write the file to disk\n")
})

setMethod("show", signature("spmCollection"), function(object){
	cat("Sample point matrix collection\n")
	cat(paste("Contains", ncol(object@data), "samples\n"))
	cat(paste("The default class vector is:", object@cl, "\n"))
	cat("The following parameters have been used:\n")
	cat(paste("The chromosome start, end and centromere locations of", attr(object@mirrorLocs, 'Version'),"were used\n"))
	cat(paste("Sigma:", object@sigma), "\n")
	cat(paste("Sample density:", object@sampleDensity, "\n\n"))
})

setMethod("show", signature("compKc"), function(object){
	cat("Comparison of ", length(object@spmCollection@cl), 
	   " samples (",sum(object@spmCollection@cl==0)," vs. ",sum(object@spmCollection@cl==1),")\n", sep="")
	cat("Using", ifelse(object@method=="siggenes","siggenes", "SNR and sample permutations"), "to find significant regions\n")
})

setMethod("show", signature("compKcSigRegions"), function(object){
	cat("Significantly different regions using",object@method, "\n")
	
	cat(ifelse(object@method=="siggenes", "delta", "snr cutoff"), "value was", object@cutoff, ", FDR =",object@fdr, "\n")

	print(object@regionTable)
	cat("Use 'write.table' to write the file to disk\n")
})