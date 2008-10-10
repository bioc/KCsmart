setMethod("unlist", signature("KCData"), functionunlistSpm <- function(x, recursive="missing", use.names ="missing"){
	spmUnlisted <- unlist(x@data)
	spmPos <- spmUnlisted[grep("pos",names(spmUnlisted))]
	spmNeg <- spmUnlisted[grep("neg",names(spmUnlisted))]
	
	return(list(pos=spmPos,neg=spmNeg))
})
