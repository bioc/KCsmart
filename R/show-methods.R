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