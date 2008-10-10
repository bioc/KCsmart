setMethod("write.table", signature(x="sigSegments"), function(x, file="", append = FALSE, quote = 7, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                 col.names=c("Status", "Chromosome", "Start", "End", "Average KC score", "Mode KC score", "Probes"), qmethod = c("escape", "double")){
	#first convert the object to a data.frame, separately for gains and for losses
	gains <- unlist(lapply(x@gains, function(y){return(list(y$chromosome,as.numeric(y$start),as.numeric(y$end),as.numeric(y$avgy),as.numeric(y$modey),paste(y$probenames,collapse=", ")))}))
	losses <- unlist(lapply(x@losses, function(y){return(list(y$chromosome,as.numeric(y$start),as.numeric(y$end),as.numeric(y$avgy),as.numeric(y$modey),paste(y$probenames,collapse=", ")))}))
	if(is.null(gains)){gains=NA}
	if(is.null(losses)){losses=NA}
	
	gainsMatrix <- matrix(gains, ncol=6, byrow=TRUE)
	lossesMatrix <- matrix(losses, ncol=6, byrow=TRUE)
	
	writeMatrix <- rbind(cbind('G', gainsMatrix), cbind('L', lossesMatrix))
	
	write.table(writeMatrix, file=file, append=append, quote=quote, sep=sep, eol=eol, na=na, dec=dec, row.names=row.names,
                 col.names=c("Status", "Chromosome", "Start", "End", "Average KC score", "Mode KC score", "Probes"), qmethod=qmethod)
	cat(paste("Output written to file", file, "\n"))
})
