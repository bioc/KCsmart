#input functions

.checkMirrorLocs <- function(mirrorLocs, data){
	chroms <- unique(data@probeAnnotation@chromosome)
	nrChroms <- length(chroms)
	mirrored <- matrix(rep(0,(3*length(mirrorLocs))),ncol=3,nrow=length(mirrorLocs))
	chromNames <- attr(mirrorLocs, 'chromNames')
	
	for (i in chroms){
		chromIndex  <- which(chromNames == i)
		if(length(chromIndex ) == 0){
			chromIndex  = as.numeric(i)
			warning(paste('Could not find chromosome',i,'in mirror locations object, using chromosome #', i, 'as an alternative'))
		}
		for(j in 1:length(mirrorLocs[[chromIndex ]])){
			max.maploc <- max(data@probeAnnotation@maploc[data@probeAnnotation@chromosome == i])
			if(mirrorLocs[[chromIndex ]][j] == 0){
				#print(paste('Mirroring at start of chromosome',i))
				mirrored[chromIndex ,1] <- mirrored[chromIndex ,1] + 1
			}
			if(mirrorLocs[[chromIndex ]][j] < max.maploc & mirrorLocs[[chromIndex]][j] > 0){
				#print(paste('Mirroring around position',mirrorLocs[[chromIndex]][j],'on chromosome',i))
				mirrored[chromIndex ,2] <- mirrored[chromIndex ,2] + 1
				
				#fix mirrorLocs ..
				if(j == length(mirrorLocs[[chromIndex ]])){
					warning(paste('Probe annotation AFTER chromosome ', i, ' end, adjusting mirror locations ..', sep=""), immediate.=TRUE)
					mirrorLocs[[chromIndex ]][j] <- max.maploc + 1
				}
			}
			if(mirrorLocs[[chromIndex ]][j] > max.maploc){
				#print(paste('Mirroring at end of chromosome',i))
				mirrored[chromIndex ,3] <- mirrored[chromIndex ,3] + 1
			}
		
		}
	}
	
	if(sum(mirrored[,3]) < nrChroms){
		print('Warning, not all chromosomes are mirrored at end!')
	}
	else{
		print('Mirror locations looking fine')
	}
	
	return(mirrorLocs)
}


.mirrorData <- function(data, mirrorLocs, mirrorLength=1500000){

       if(!is(data, "KcghDataSum")){stop("Need KcghDataSum object as input")}

	cc <- new("KcghDataSum")
	data <- new("KcghDataMirror", mirrorLocs, data)
	data <- sort(data)
	
	dataChromosomes <- unique(data@probeAnnotation@chromosome)
	chromNames <- attr(mirrorLocs, 'chromNames')
	
	for (i in dataChromosomes){
		#data should be ordered!
		
		current.data.index <- which(data@probeAnnotation@chromosome == i)
		cc@probeAnnotation <- data@probeAnnotation[current.data.index]

		cc@pos <- data@pos[current.data.index]
		cc@neg <- data@neg[current.data.index]
		
		first.probe.index <- current.data.index[1]
		new.first.probe.name <- paste(data@probeAnnotation@name[first.probe.index], 'first', i)
		data@probeAnnotation@name[first.probe.index] <- new.first.probe.name
		
		last.probe.index <- current.data.index[length(current.data.index)]
		new.last.probe.name <- paste(data@probeAnnotation@name[last.probe.index], 'last', i)
		data@probeAnnotation@name[last.probe.index] <- new.last.probe.name
		
		
		chromIndex  <- which(chromNames == i)
		if(length(chromIndex ) == 0){
			chromIndex  = as.numeric(i)
			warning(paste('Could not find chromosome',i,'in mirror locations object, using chromosome #', i, 'as an alternative'))
			data@probeAnnotation@chromosome[data@probeAnnotation@chromosome == i] <-  chromNames[chromIndex]
		}
		
		for(j in 1:length(mirrorLocs[[chromIndex ]])){
		
			#print(paste('At:', i, j))
			cm <- mirrorLocs[[chromIndex ]][j]
			
			#reset mirrorLength
			mirrorLength.current <- mirrorLength
			
			if(sum(cc@probeAnnotation@maploc < cm)){
				last.probe.before.index <- max(which(cc@probeAnnotation@maploc < cm))
				last.probe.before <- cc@probeAnnotation@maploc[last.probe.before.index]
			}
			else{
				#print('No probes before')
				last.probe.before.index <- 0
				last.probe.before <- 0
			}
			
			if(last.probe.before.index == length(cc@pos)){
				#print('No probes after')
				first.probe.after <- last.probe.before
			}
			else{
				first.probe.after <- cc@probeAnnotation@maploc[last.probe.before.index + 1]
				
				#if mirror probes around centromere would influence probes on other side of centromere, adjust mirrorLength
				if((cm > 0) & ((last.probe.before + mirrorLength.current) > first.probe.after)){
					mirrorLength.current <- first.probe.after - last.probe.before
					#print(paste("Yup, at: ", i , ',', j))
					#print(paste('Changed to', mirrorLength.current))
				}
			}
			
			#print(paste('before:', last.probe.before, 'after:', first.probe.after))
			
			#first calculate position of to be mirrored probes
			left.mirror.probes <- which((cc@probeAnnotation@maploc < last.probe.before) & (cc@probeAnnotation@maploc > (last.probe.before - mirrorLength.current)))
			right.mirror.probes <- which((cc@probeAnnotation@maploc > first.probe.after) & (cc@probeAnnotation@maploc < (first.probe.after + mirrorLength.current)))

			#now add maplocs to data maploc

			#adjust maploc positions
			maplocs.left <- last.probe.before + (last.probe.before - cc@probeAnnotation@maploc[left.mirror.probes])
			maplocs.right <- first.probe.after - (cc@probeAnnotation@maploc[right.mirror.probes] - first.probe.after)
			maplocs <- c(maplocs.left, maplocs.right)
			
			data@probeAnnotation@chromosome <- c(data@probeAnnotation@chromosome, rep(chromNames[chromIndex], length(maplocs)))
			data@probeAnnotation@maploc <- c(data@probeAnnotation@maploc, as.integer(maplocs))
			
			virtualNames <- vector()
			if(length(left.mirror.probes) > 0){
				virtualNames <- c(virtualNames, paste(cc@probeAnnotation@name[left.mirror.probes], 'virtual left',i))
			}
			if(length(right.mirror.probes) > 0){
				virtualNames <- c(virtualNames, paste(cc@probeAnnotation@name[right.mirror.probes], 'virtual right',i))
			}

			data@probeAnnotation@name <- c(data@probeAnnotation@name, virtualNames)
			
			data@pos <- c(data@pos, cc@pos[left.mirror.probes], cc@pos[right.mirror.probes])
			data@neg <- c(data@neg, cc@neg[left.mirror.probes], cc@neg[right.mirror.probes])
		}
	}

	data <- sort(data)
	data@mirrorLocs <- mirrorLocs

	return(data)
}

.spm2spmc <- function(spm){
  if(!is(spm, "samplePointMatrix")){stop("Need samplePointMatrix object as input")}
    
	sampleDensity <- spm@sampleDensity
	mirrorLocs <- spm@mirrorLocs

	chromosomes <- names(spm@data)
	chromNames <- attr(mirrorLocs, 'chromNames')
	chromosomesOrdered <- chromNames[chromNames %in% chromosomes]
	chromosomesOrdered <- c(chromosomesOrdered, chromosomes[!(chromosomes %in% chromNames)])
	chromosomes <- chromosomesOrdered

	spmCollection <- new("spmCollection")
	spmCollection@mirrorLocs <- mirrorLocs
	spmCollection@sigma <- spm@sigma
	spmCollection@sampleDensity <- spm@sampleDensity
	spmCollection@data <- matrix(ncol=1, nrow=0)
	tmpAnnotation <- vector(mode="list")
	
	for(i in chromosomes){	
		#browser()	
		spmCollection@data <- rbind(spmCollection@data, matrix(as.numeric(spm[[i]]$pos) + as.numeric(spm[[i]]$neg)))
		
		nrSamplePoints <- length(spm[[i]]$pos)
		samplePointPositions <- seq(1, (nrSamplePoints*sampleDensity), by=sampleDensity)
		tmpAnnotation$chrom <- c(tmpAnnotation$chrom, rep(i, nrSamplePoints))
		tmpAnnotation$maploc <- c(tmpAnnotation$maploc, as.numeric(samplePointPositions))
	}

	spmCollection@annotation <- new("probeAnnotation", tmpAnnotation$chrom, tmpAnnotation$maploc)
	
	
	return(spmCollection)
}

.add2spmc <- function(spmc, spm){
	chromosomes <- unique(spmc@annotation@chromosome)
	tempData <- matrix(ncol=1, nrow=0)
	
	for(i in chromosomes){		
		tempData <- rbind(tempData, matrix(as.numeric(spm[[i]]$pos) + as.numeric(spm[[i]]$neg)))
	}
	
	spmc@data <- cbind(spmc@data, tempData)
	
	return(spmc)
}


.convertCGHbase <-function(cghBase){
	#convert cghRaw to internal data representation
	#get midpositions
	maploc <- apply(cbind(CGHbase::bpstart(cghBase), CGHbase::bpend(cghBase)),1,function(x)return(round(mean(x))))
	tData <- data.frame(chrom=CGHbase::chromosomes(cghBase), maploc=maploc, data=CGHbase::copynumber(cghBase))
	row.names(tData) <- Biobase::featureNames(cghBase) 
	KCGHdata <- new("KcghData", tData)
	return(KCGHdata)	
}


#wrapper functions
calcSpm <- function(data, mirrorLocs, sigma=1000000, sampleDensity=50000, maxmem=1000, verbose=T){

	#do checks
	#is it CGHbase data?
	if(is(data, "cghRaw")){
		data <- .convertCGHbase(data)
	}

	if(!is(data, "KcghData")){
		data <- new("KcghData", data)
	}
	
	mirrorLocs <- KCsmart:::.checkMirrorLocs(mirrorLocs, data)
	mirrorLength <- sigma * 4
	
	print("Splitting data ..")
	data <- new("KcghDataSplit", data)
	
	print("Summing data ..")
	data <- new("KcghDataSum", data)
	
	print("Mirroring data ..")
	data <- KCsmart:::.mirrorData(data, mirrorLocs, mirrorLength)

	print("Calculating sample point matrix ..")
	spm <- KCsmart:::.samplePointMatrix(data, sampleDensity=sampleDensity, sigma=sigma, maxmem=maxmem, verbose=verbose)
	rm(data)

	print('Done')
	return(spm)
}

calcSpmCollection <- function(data, mirrorLocs, cl=NULL, data2=NULL, sigma=1000000, sampleDensity=50000, maxmem=1000, verbose=F, doChecks=T){

	#do checks
	#is it CGHbase data?
	if(doChecks){
	
		if(is.null(data2) & is.null(cl)){
			stop('Please provide a default class vector or a second data set')
		}
	
		if(is(data, "cghRaw")){
			data <- KCsmart:::.convertCGHbase(data)
		}
	
		if(!is(data, "KcghData")){
			data <- new("KcghData", data)
		}
		
		#checks for 2nd data set
		if(!is.null(data2)){
			if(is(data2, "cghRaw")){
				data2 <- KCsmart:::.convertCGHbase(data2)
			}
	
			if(!is(data2, "KcghData")){
				data2 <- new("KcghData", data2)
			}
			if(!all.equal(data@probeAnnotation, data2@probeAnnotation)) {stop("The data has different annotations, unable to continue")}
		}
		
		
		mirrorLocs <- KCsmart:::.checkMirrorLocs(mirrorLocs, data)
	}
	
	nrSamples <- ncol(data@data)
	mirrorLength <- sigma * 4
	
	#print("Splitting data ..")
	data <- new("KcghDataSplit", data)
	
	originalData <- data

	for(i in 1:nrSamples){
		cat(paste("Processing sample", i, "/", nrSamples, "\r"))
		#single sample construction
		
		#copy probeAnnotations
		singleSample <- originalData
		singleSample@pos <- originalData@pos[,i, drop=F]
		singleSample@neg <- originalData@neg[,i, drop=F]
		singleSample <- new("KcghDataSum", singleSample)  
		
		data <- singleSample
		
		#cat(paste("Mirroring data ..","\r"))
		data <- KCsmart:::.mirrorData(data, mirrorLocs, mirrorLength)
	
		#cat(paste("Calculating sample point matrix ..","\r"))
		spm <- KCsmart:::.samplePointMatrix(data, sampleDensity=sampleDensity, sigma=sigma, maxmem=maxmem, verbose=verbose)
		rm(data)
		
		#cat(paste("Converting to spmc ..","\r"))
		if(!exists("spmc")){
			spmc <- KCsmart:::.spm2spmc(spm)
		}
		else{
			spmc <- KCsmart:::.add2spmc(spmc, spm)
		}
	}
	
	cat("\n")
	
	if(!is.null(cl)){
		if(sum(cl == 0 | cl == 1) != nrSamples) {stop('Invalid class vector given')}
		else{
			spmc@cl <- cl
		}
	}
	
	if(!is.null(data2)){
		spmc2 <- calcSpmCollection(data2, spmc@mirrorLocs, cl=cl, data2=NULL, sigma=sigma, sampleDensity=sampleDensity, maxmem=maxmem, verbose=verbose, doChecks=F)
		nrSamplesFirstClass <- ncol(spmc@data)
		spmc@data <- cbind(spmc@data, spmc2@data)
		rm(spmc2)
		spmc@cl <- c(rep(0, nrSamplesFirstClass), rep(1, (ncol(spmc@data) - nrSamplesFirstClass)))
	}
	
	
	return(spmc)
}


#sig level
findSigLevelTrad <- function(data, observedSpm, n=1, p=0.05, maxmem=1000){
	#method 1 (as published)
	if(!is(data, "KcghData")){data <- new("KcghData", data)}

	print(paste('Calculating alpha = ',p,'significance cut-off'))
	
	sampleDensity <- observedSpm@sampleDensity
	sigma <- observedSpm@sigma
	
	peaks <- vector(mode="list")
	
	dataSplit <- new("KcghDataSplit", data)
	
	spmUnlisted <- unlist(observedSpm)
		
	posPeaks <- KCsmart:::.findPeaks(spmUnlisted$pos)
	negPeaks <- KCsmart:::.findPeaks(spmUnlisted$neg, mode='neg')
	
	print(paste('Found ', length(posPeaks), ' pos peaks and ', length(negPeaks), ' neg peaks in observed sample point matrix'))
	
	p.corrected.pos <- p / length(posPeaks)
	p.corrected.neg <- p / length(negPeaks)
	
	mirrorLength <- sigma * 4
	mirrorLocs <- observedSpm@mirrorLocs
	
	print('Starting permutations ..')
	
	for(i in 1:n){
		cat(paste("\rAt iteration", i,'of',n))
		
		spmPermuted <- KCsmart:::.permutedSpm(dataSplit, mirrorLocs, mirrorLength, sampleDensity=sampleDensity, sigma=sigma, maxmem=maxmem)
		spmPermuted <- unlist(spmPermuted)

		peaks$pos <- c(peaks$pos, KCsmart:::.findPeaks(spmPermuted$pos))
		peaks$neg <- c(peaks$neg, KCsmart:::.findPeaks(spmPermuted$neg, mode='neg'))
	}
	
	
	cutoff.pos <- min(sort(peaks$pos, decreasing=TRUE)[1:(p.corrected.pos*length(peaks$pos))])
	cutoff.neg <- max(sort(peaks$neg)[1:(p.corrected.neg*length(peaks$neg))])

	print("\n")	
	return(list(pos=cutoff.pos, neg=cutoff.neg))
}


findSigLevelFdr <- function(data, observedSpm, n=1, fdrTarget=0.05, maxmem=1000){
		#FDR approach
		
		if(!is(data, "KcghData")){data <- new("KcghData", data)}
		print(paste('Calculating fdr = ',fdrTarget,'significance cut-off'))
		
		sampleDensity <- observedSpm@sampleDensity
		sigma <- observedSpm@sigma
		mirrorLength <- sigma * 4
		mirrorLocs <- observedSpm@mirrorLocs
		
		totalSpmPermuted <- vector(mode="list")
		observedSpmUnlisted <- unlist(observedSpm)
		dataSplit <- new("KcghDataSplit", data)
		
		for(i in 1:n){
			cat(paste("\rAt iteration", i,'of',n))
		
			spmPermuted <- KCsmart:::.permutedSpm(dataSplit, mirrorLocs, mirrorLength, sampleDensity=sampleDensity, sigma=sigma, maxmem=maxmem)
			spmPermuted <- unlist(spmPermuted)
			totalSpmPermuted$pos <- c(totalSpmPermuted$pos, spmPermuted$pos)
			totalSpmPermuted$neg <- c(totalSpmPermuted$neg, spmPermuted$neg)
		}
		
		finalCutoffPos <- KCsmart:::.findCutoffByFdr(observedSpmUnlisted$pos, totalSpmPermuted$pos, fdr=fdrTarget)
		finalCutoffNeg <- -KCsmart:::.findCutoffByFdr(-observedSpmUnlisted$neg, -totalSpmPermuted$neg, fdr=fdrTarget)

		print("\n")	
		return(list(pos=finalCutoffPos,neg=finalCutoffNeg))
}


.findPeaks_old <- function(data, mode='pos'){
	data <- data[!is.na(data)]
	dataOriginal <- data
	
	data <- c(data[-1], 0) - data
	data <- apply(as.matrix(data), 1,function(x){if(x>0){return(1)};if(x<0){return(-1)}else{return(0)}})
	nonFlats <- which(data!=0)
	data <- data[nonFlats]
	data <- data - c(data[-1], 0)
	if(mode=='pos'){
		peaksIndices <- which(data==2) + 1
	}
	else{
		peaksIndices <- which(data==-2) + 1
	}
	
	peaks <- dataOriginal[nonFlats][peaksIndices]

	return(peaks)
}

.findPeaks <- function(data, mode='pos'){
    dir <- ifelse(mode=="pos", -2, 2)
    data <- data[!is.na(data)]
    dataOriginal <- data[-1]
   
    data <- diff.default(data)
    data[data > 0] <- 1
    data[data < 0] <- -1
    nonFlats <- data != 0

    data2 <- diff.default(data[nonFlats])
   
    dataOriginal[nonFlats][data2 == dir]
}


.findCutoffByFdr <- function(observed, permuted, fdr=0.05, precision=0.01){
	#righ-tailed!!
	t.fdr <- 0
	tCutoff <- max(observed, na.rm=TRUE)
	tCutoffOld <- 0
	maxTFdr <- fdr + (precision*fdr)
	minTFdr <-  fdr - (precision*fdr)
	
	multiplier <- length(permuted) / length(observed)
	
	#WARNING: this is 1-tailed (right)!
	
	while(!((t.fdr > minTFdr) & (t.fdr < maxTFdr))){
		if(length(tCutoffOld) > 1){
			tCutoffOldComp <- tCutoffOld[(length(tCutoffOld) - 1)]
		}
		else{
			tCutoffOldComp <- 0
		}
		
		if(t.fdr > fdr){
			tCutoff <- tCutoff + 0.5*abs(tCutoffOldComp - tCutoff)
		}
		else{
			tCutoff <- tCutoff - 0.5*abs(tCutoffOldComp - tCutoff)
		}
		
		tCutoffOld <- c(tCutoffOld, tCutoff)
		
		fp <- sum(permuted >= tCutoff, na.rm=TRUE)
		called <- sum(observed >= tCutoff, na.rm=TRUE)
		
		if(called == 0){print('Uhoh, 0 called!');return()}
		
		t.fdr <-  (fp/multiplier)  / called
		
		if(length(tCutoffOld) > 50){print('Uhoh, maximum number of iterations reached!'); return()}
	}

	return(tCutoff)
}


.permutedSpm <- function(data, mirrorLocs, mirrorLength=sigma*4, sampleDensity=50000, sigma=1000000, maxmem=1000){
	#takes sample data and returns sample point matrix of permuted probe label data
	if(!is(data, "KcghDataSplit")){stop("Need KcghDataSplit as input")}
	
	#permute data
	data@pos <- apply(data@pos, 2, function(x){return(sample(x))})
	data@neg <- apply(data@neg, 2, function(x){return(sample(x))})

	dataSum <- new("KcghDataSum", data)
	dataMirrored <- KCsmart:::.mirrorData(dataSum, mirrorLocs, mirrorLength)

	#get sample point matrix
	spmPermuted <- KCsmart:::.samplePointMatrix(dataMirrored, sampleDensity, sigma, maxmem=maxmem, verbose=FALSE)
	
	return(spmPermuted)
}


.samplePointMatrix <- function(data, sampleDensity=50000,sigma=1000000, mirrorLocs=data@mirrorLocs, maxmem=1000, verbose=TRUE){
        if(!is(data, "KcghDataSum")){stop("Need KcghDataSum object as input")}
	
	cl <- unlist(lapply(mirrorLocs, max))
	attr(cl, 'chromNames') <- attr(mirrorLocs, 'chromNames')
	
	spm <- new("samplePointMatrix")
	sigma.4 <- sigma * 4
	nrsamples <- data@nrSamples
	dataChromosomes <- unique(data@probeAnnotation@chromosome)
	chromNames <- attr(cl, 'chromNames')
	
	#prepare scaffold matrix -> maximal matrix
	sample.points <- seq.int(0,max(cl), by=sampleDensity)
	max.nrprobes <- max(aggregate(data@probeAnnotation@chromosome, by=list(chromosome=data@probeAnnotation@chromosome), length)[,2])
	
	#the maximum number of probes that are allowed to be processed simultaneously, this directly affects memory usage
	#use lower values for less memory usage
	max.nrprobes <- maxmem
	
	sp <- rep.int(sample.points, max.nrprobes)
	scaffold.matrix <- matrix(sp, nrow=max.nrprobes, byrow=TRUE)

	for(i in dataChromosomes){
		if(verbose){cat(paste("\nProcessing chromosome", i,"\n"))}
		
		chromIndex  <- which(chromNames == i)
		if(length(chromIndex ) == 0){
			chromIndex  = as.numeric(i)
			warning(paste('Could not find chromosome',i,'in mirror locations object'))
		}
		
		all.cprobes <- which(data@probeAnnotation@chromosome == i)
		
		#set max chunk size
		chunk.size <- max.nrprobes
		chunks <- c(seq(0,length(all.cprobes), by=chunk.size), length(all.cprobes))
		
		#find virtual probes
		real.first.probe <- grep('first', data@probeAnnotation@name[all.cprobes])
		real.last.probe <- grep('last', data@probeAnnotation@name[all.cprobes])
		virtual.probes <- grep('virtual', data@probeAnnotation@name[all.cprobes])
		virtual.probes <- virtual.probes[((virtual.probes > real.first.probe) & (virtual.probes < real.last.probe))]
		
		c.scaffold.matrix <- scaffold.matrix[1, 1:(cl[chromIndex ]/sampleDensity)]

		#identify virtual probes
		virtual <- c.scaffold.matrix < data@probeAnnotation@maploc[all.cprobes][real.first.probe]
		virtual <- virtual | (c.scaffold.matrix > (data@probeAnnotation@maploc[all.cprobes][real.last.probe]))
		if(length(virtual.probes > 0)){
			virtual <- virtual | ((c.scaffold.matrix > data@probeAnnotation@maploc[all.cprobes][(min(virtual.probes) - 1)]) & (c.scaffold.matrix < (data@probeAnnotation@maploc[all.cprobes][(max(virtual.probes) + 1)])))
		}
		
		stored.cspm <- vector(mode="list")
		stored.w.regression <- vector(length=(cl[chromIndex ]/sampleDensity))
		
		for(current.chunk in 1:(length(chunks) - 1)){
		
			if(verbose){cat(paste("\rProcessing chunk", current.chunk, 'of', (length(chunks) - 1), '(',(chunks[current.chunk] + 1),'-', chunks[current.chunk + 1],')'))}
			cprobes <- all.cprobes[(chunks[current.chunk]+1):chunks[current.chunk + 1]]
			nrprobes <- length(cprobes)
		
			cspm <- scaffold.matrix[1:nrprobes, 1:(cl[chromIndex ]/sampleDensity), drop=FALSE]

			probe.locs <- data@probeAnnotation@maploc[cprobes]
			cspm <- abs(cspm - probe.locs)
			
			cspm[cspm > sigma.4] <- NA
			cspm[,virtual] <- NA
			
			cspm <- dnorm(cspm, mean=0, sd=sigma)
			cspm[is.na(cspm)] <- 0
			cspm.pos <- cspm * data@pos[cprobes]
			cspm.neg <- cspm * data@neg[cprobes]
			
			w.regression <- colSums(cspm, na.rm=TRUE) * nrsamples
			
			stored.w.regression <- rbind(stored.w.regression, w.regression)
			rm(cspm)
			rm(w.regression)
			
			csspm.pos <- colSums(cspm.pos, na.rm=TRUE)
			csspm.neg <- colSums(cspm.neg, na.rm=TRUE)
			rm(cspm.pos)
			rm(cspm.neg)
			
			stored.cspm$pos <- rbind(stored.cspm$pos, csspm.pos)
			stored.cspm$neg <- rbind(stored.cspm$neg, csspm.neg)
			
			rm(csspm.pos)
			rm(csspm.neg)
		}
		
		stored.w.regression <- colSums(stored.w.regression, na.rm=TRUE)
		stored.cspm$pos <- colSums(stored.cspm$pos, na.rm=TRUE) / stored.w.regression
		stored.cspm$neg <- colSums(stored.cspm$neg, na.rm=TRUE) / stored.w.regression

		spm[[i]] <- list(pos = stored.cspm$pos, neg = stored.cspm$neg)
		
		attr(spm[[i]], 'chromosome') <- chromIndex 
		rm(stored.cspm)
		
		#gc()
	}
	
	total <- 0
	maxy <- 0
	miny <- 0

	for(i in 1:length(spm@data)){
		total <- total + length(spm[[i]]$pos)
		if(max(spm[[i]]$pos, na.rm=TRUE) > maxy){
			maxy <- max(spm[[i]]$pos, na.rm=TRUE)
		}
		if(min(spm[[i]]$neg, na.rm=TRUE) < miny){
			miny <- min(spm[[i]]$neg, na.rm=TRUE)
		}
	}
	
	spm@totalLength <- as.integer(total)
	spm@maxy <- maxy
	spm@miny <- miny
	spm@sampleDensity <- as.integer(sampleDensity)
	spm@sigma <- as.integer(sigma)
	spm@mirrorLocs <- data@mirrorLocs

	#attach probeAnnotation to spm, less the virtual probes
	virtualProbes <- grep("virtual", data@probeAnnotation@name)
	spm@probeAnnotation <- data@probeAnnotation[-virtualProbes]
	
	if(verbose){cat("\n\n")}
	
	return(spm)
}


getSigSegments <- function(spm, sigLevels, chromosomes=NULL){
    if(!is(spm, "samplePointMatrix")){stop("Need samplePointMatrix object as input")}
    
	sigSegments <- new("sigSegments")
	
	sampleDensity <- spm@sampleDensity
	mirrorLocs <- spm@mirrorLocs
	
	if(is.null(chromosomes)){
		chromosomes <- names(spm@data)
		chromNames <- attr(mirrorLocs, 'chromNames')
		chromosomesOrdered <- chromNames[chromNames %in% chromosomes]
		chromosomesOrdered <- c(chromosomesOrdered, chromosomes[!(chromosomes %in% chromNames)])
		chromosomes <- chromosomesOrdered
	}
	else{
		chromosomes <- as.character(chromosomes)
		spmLengths <- unlist(lapply(spm@data, function(x){return(length(x$pos))}))
		total <- sum(spmLengths[chromosomes])
	}
	
	for(i in chromosomes){
		cprobes <- (spm@probeAnnotation@chromosome == i)
		
		#gains
		sigGains <- which(spm[[i]]$pos >= sigLevels$pos)
		
		if(length(sigGains) > 0){
		
		#now get chromosome region
		t <- diff(sigGains)
		t2 <- t - c(2,t[-c(length(t))])
		t2[1] <- 1
		t2[length(t2)] <- 1
		t3 <- which(t2 != 0)
		t <- t3
		
		for(j in seq(1,(length(t) - 1), by=2)){
			currentSegment <- vector(mode="list")
			startPosition <- t[j]
			end.position <- t[j+1]
			currentSegment$chromosome <- i
			currentSegment$x <- sigGains[startPosition:end.position]
			currentSegment$y <- signif(spm[[i]]$pos[sigGains[startPosition:end.position]],4)
			currentSegment$avgy <- signif(mean(currentSegment$y, na.rm=TRUE), 4)
			currentSegment$modey <- signif(max(currentSegment$y, na.rm=TRUE), 4)
			currentSegment$start <- (sigGains[startPosition]) * sampleDensity
			currentSegment$end <- (sigGains[end.position]) * sampleDensity
			currentSegment$probes <- which((spm@probeAnnotation@maploc[cprobes] >= currentSegment$start) & (spm@probeAnnotation@maploc[cprobes] <= currentSegment$end))
			currentSegment$probenames <- spm@probeAnnotation@name[cprobes][currentSegment$probes]
			sigSegments@gains <- c(sigSegments@gains, list(currentSegment))
		}
		}
		
		
		#losses
		sigLosses <- which(spm[[i]]$neg <= sigLevels$neg)
		
		if(length(sigLosses) > 0){
		
		#now get chromosome region
		t <- diff(sigLosses)
		t2 <- t - c(2,t[-c(length(t))])
		t2[1] <- 1
		t2[length(t2)] <- 1
		t3 <- which(t2 != 0)
		t <- t3
		
		#sigSegments[[i]]$negsegments <- vector(mode='list', length=length(t)/2)

		for(j in seq(1,(length(t) - 1), by=2)){
			currentSegment <- vector(mode="list")
			startPosition <- t[j]
			end.position <- t[j+1]
			currentSegment$chromosome <- i
			currentSegment$x <- sigLosses[startPosition:end.position]
			currentSegment$y <- signif(spm[[i]]$neg[sigLosses[startPosition:end.position]], 4)
			currentSegment$avgy <- signif(mean(currentSegment$y, na.rm=TRUE), 4)
			currentSegment$modey <- signif(min(currentSegment$y, na.rm=TRUE), 4)
			currentSegment$start <- (sigLosses[startPosition]) * sampleDensity
			currentSegment$end <- (sigLosses[end.position]) * sampleDensity
			currentSegment$probes <- which((spm@probeAnnotation@maploc[cprobes] >= currentSegment$start) & (spm@probeAnnotation@maploc[cprobes] <= currentSegment$end))
			currentSegment$probenames <- spm@probeAnnotation@name[cprobes][currentSegment$probes]
			
			sigSegments@losses <- c(sigSegments@losses, list(currentSegment))
		}
		}
	}
	
	sigSegments@sigma <- spm@sigma
	sigSegments@sigLevels <- sigLevels
	return(sigSegments)
}


.getSigRegions <- function(spm, sigLevels, chromosomes=NULL){
	sigRegions <- new("sigRegions")
	mirrorLocs <- spm@mirrorLocs
	
	if(is.null(chromosomes)){
		chromosomes <- names(spm@data)
		chromNames <- attr(mirrorLocs, 'chromNames')
		chromosomesOrdered <- chromNames[chromNames %in% chromosomes]
		chromosomesOrdered <- c(chromosomesOrdered, chromosomes[!(chromosomes %in% chromNames)])
		chromosomes <- chromosomesOrdered
	}
	else{
		chromosomes <- as.character(chromosomes)
		spmLengths <- unlist(lapply(spm@data, function(x){return(length(x$pos))}))
		total <- sum(spmLengths[chromosomes])
	}
	
	for(i in chromosomes){
		sigRegions[[i]]$posx <- which(spm[[i]]$pos >= sigLevels$pos)
		sigRegions[[i]]$posy <- spm[[i]]$pos[sigRegions[[i]]$posx]
		sigRegions[[i]]$negx <- which(spm[[i]]$neg <= sigLevels$neg)
		sigRegions[[i]]$negy <- spm[[i]]$neg[sigRegions[[i]]$negx]
		chromIndex  <- which(attr(mirrorLocs, 'chromNames') == i)
		if(length(chromIndex ) == 0){
			chromIndex  = as.numeric(i)
			warning(paste('Could not find chromosome',i,'in mirror locations object'))
		}
		attr(sigRegions[[i]], 'chromosome') <- chromIndex 
	}
	
	sigRegions@sigma <- spm@sigma
	return(sigRegions)
}

#plot functions
plotScaleSpace <- function(spms, sigLevels, chromosomes=NULL, type='b'){
	mirrorLocs <- spms[[1]]@mirrorLocs
	scaleSpace <- new("scaleSpace")
	
	for(i in 1:length(spms)){
		#since only calculating the gains or the losses (instead of both) in the .getSigRegions function hardly makes any difference at all (performance-wise)
		#the 'type' distinction is not propagated into that function
		scaleSpace@data[[i]] <- KCsmart:::.getSigRegions(spms[[i]], sigLevels[[i]], chromosomes)
	}
	
	#plot sig regions, pass any spm as argument
	plot(scaleSpace, spm=spms[[1]], type=type)
}


idPoints <- function(spm, mode='pos', dev=2, chromosomes=NULL){
	options(locatorBell = FALSE)
	if(!is.null(chromosomes)){
		chromosomes <- as.character(chromosomes)
	}
	else{
		chromosomes <- names(spm@data)
		mirrorLocs <- spm@mirrorLocs
		chromNames <- attr(mirrorLocs, 'chromNames')
		chromosomesOrdered <- chromNames[chromNames %in% chromosomes]
		chromosomesOrdered <- c(chromosomesOrdered, chromosomes[!(chromosomes %in% chromNames)])
		chromosomes <- chromosomesOrdered
	}
	
	spmLengths <- unlist(lapply(spm@data, function(x){return(length(x$pos))}))
	spmLengths <- spmLengths[chromosomes]
	colin <- cumsum(spmLengths)
	
	yvalues <- vector()
	for(i in chromosomes){
		if(mode=='pos'){
			yvalues <- c(yvalues, signif(spm[[i]]$pos, 3))
		}
		else{
			yvalues <- c(yvalues, signif(spm[[i]]$neg, 3))
		}
	}
	totalLength <- length(yvalues)
	
	dev.set(dev)
	identifiedPoints <- identify(seq(1,totalLength), yvalues, labels=yvalues)
	if(length(identifiedPoints) > 0){
		selectedy <- yvalues[identifiedPoints]
		selectedPoints <- data.frame(colin=identifiedPoints, KCscore=selectedy)
		row.names(selectedPoints) <- NULL
		#get chromosome + chromosome position
		chromFind <- apply(matrix(selectedPoints$colin),1, function(x){return(abs(x-colin))})
		chromFind <- apply(matrix(chromFind, ncol=length(selectedPoints$colin)), 2, which.min)
		
		chromDistance <- selectedPoints$colin - colin[chromFind]
		chromFind[chromDistance > 0] <- chromFind[chromDistance > 0] + 1
		selectedPoints$chromosome <- chromosomes[chromFind]
		selectedPoints$chromPosition <- (selectedPoints$colin - c(0,colin)[chromFind]) * (spm@sampleDensity / 1000)
		selectedPoints <- selectedPoints[,c('KCscore', 'chromosome','chromPosition','colin')]
		return(selectedPoints)
	}
}

compareSpmCollection <- function(spmCollection, nperms=20, method=c("siggenes", "perm"), siggenes.args=NULL, altcl=NULL) {

	method <- match.arg(method)
	stopifnot(is(spmCollection,"spmCollection"))
	if(!is.null(altcl)){
		if((sum(altcl == 0 | altcl == 1)) != length(spmCollection@cl)){stop('Invalid class vector given')}
		else{
			spmCollection@cl <- altcl
		}
	}

	res <- switch(method,
		siggenes = KCsmart:::.comparativeKcSiggenes(spmCollection@data, spmCollection@cl, nperms=nperms, siggenes.args),
		perm = KCsmart:::.comparativeKcPerms(spmCollection@data, spmCollection@cl , nperms=nperms))

	new("compKc", spmCollection, method, res)
}

getSigRegionsCompKC <- function(compKc, fdr=.01, maxRegionGap=10) {

	stopifnot(is(compKc,"compKc"))

	#prepare the result stuff
	siglist <- list()
	
	if(compKc@method == "siggenes") {
		#get all dvalues from sam
		a <- compKc@siggenesResult
		siglist$testval <- a@d
		#find delat, suppress output
		sink(tempfile())
		deltas <- findDelta(a,fdr=fdr)
		sink()
		siglist$cutoff <- ifelse(is.matrix(deltas), deltas[2,1], deltas[1])
		siggeneslist <- multtest:::summary(a, siglist$cutoff)
		siglist$issignificant <- 1:nrow(compKc@spmCollection@data) %in% siggeneslist@mat.sig$Row
	}
	else {
		siglist$cutoff <- KCsmart:::.findfdrcutoff(compKc@snrResult@permutations, compKc@snrResult@snrValues, fdr)
		siglist$testval <- compKc@snrResult@snrValues
		siglist$issignificant <- abs(siglist$testval) > siglist$cutoff
	}
	
	#replace NA in sig vector with FALSE
	siglist$issignificant[is.na(siglist$issignificant)] <- FALSE
	#determine regions
	regions <- KCsmart:::.getRegions(siglist$issignificant, maxRegionGap)

	#add annotation to regiontable
	regions$startchrom <- compKc@spmCollection@annotation@chromosome[regions$startrow]
	regions$endchrom <- compKc@spmCollection@annotation@chromosome[regions$endrow]

	#split regions on chromsome border
	if(any(regions$startchrom != regions$endchrom)) {
		tosplit <- subset(regions, regions$startchrom != regions$endchrom)
		notsplit <- subset(regions, regions$startchrom == regions$endchrom)
	
		for(i in 1:nrow(tosplit)) {
			s <- tosplit[i,]
			done <- F
			repeat {
				chendrow <- max(which(compKc@spmCollection@annotation@chromosome == s$startchrom))
				notsplit <- rbind(notsplit, c(s$startrow, chendrow, s$startchrom, s$startchrom))
				if(done) break;

				s$startchrom <- compKc@spmCollection@annotation@chromosome[chendrow+1]
				s$startrow <- chendrow+1

				if(s$startchrom == s$endchrom) done <- T
			}
		}
		#sort table
		regions <- notsplit[order(notsplit$startrow),]
		
	}
	startposition <- compKc@spmCollection@annotation@maploc[regions$startrow]
	endposition <- compKc@spmCollection@annotation@maploc[regions$endrow]

	new("compKcSigRegions", 
	   regionTable=data.frame(startrow=regions$startrow, endrow=regions$endrow, chromosome=regions$startchrom, startposition, endposition), 
	   method=compKc@method, fdr=fdr, cutoff=siglist$cutoff)
}

.comparativeKcSiggenes <- function(data, cl, nperms, args=NULL) {
	require(siggenes)

	callargs <- c(list(data, cl, B=nperms), args)
	do.call(sam, callargs)
}

.comparativeKcPerms <- function(data, cl, nperms) {

	snrresults <- KCsmart:::.snr(data[,cl==0], data[,cl==1])

	#prepare the permutations
	perms <- KCsmart:::.makePermutations(sum(cl==0), sum(cl==1), nperms)

	m <- matrix(NA, nrow=nrow(data), ncol=length(perms))
	for(p in 1:length(perms)) {
		m[,p] <- KCsmart:::.snr(data[,perms[[p]]$a], data[,perms[[p]]$b], s0=snrresults$s0)
	}
	
	new("snrResult", snrValues=snrresults$snr, fudge=snrresults$s0, permutations=m)
}

.makePermutations <- function(sizeA, sizeB, nperms=2000) {
	#How many perms are there?
	totalperms <- choose(sizeA+sizeB, sizeA)
	n <- sizeA+sizeB
	a <- 1:sizeA
	b <- (sizeA+1):n

	#FIXME: if totalperms !>> nperms we should select perms from allperms to avoid duplicates
	if(totalperms > nperms) {
		return(lapply(1:nperms, function(x){z <- sample(n); return(list(a=z[a], b=z[b]))}))
	}
	else {
		warning("There are only ", totalperms," permutations. Returning all")
		#TODO: implemtent all perms

	}
}

.snr <- function(a, b, s0=NULL) {
	n1 <- ncol(a)
	n2 <- ncol(b)
	
	n1row <- rep(n1, nrow(a)) - rowSums(is.na(a))
	n2row <- rep(n2, nrow(a)) - rowSums(is.na(b))

	m1 <- rowMeans(a,na.rm=T)
	m2 <- rowMeans(b,na.rm=T)

	var1 <- KCsmart:::.varr(a, meanx=m1, n1row)
	var2 <- KCsmart:::.varr(b, meanx=m2, n2row)
	
	sigma12 <- sqrt( (var1/n1row) +  (var2/n2row ) )

	if(is.null(s0)){
		s0 <- quantile(sigma12, .95, na.rm=T)
		return(list(snr=(m1-m2)/(sigma12 + s0), s0=s0))
	} else {
		return(	(m1-m2)/(sigma12 + s0))
	}
}

.varr <- function(x, meanx, nna=NULL) {
	n <- ncol(x)
	p <- nrow(x)
	Y <-matrix(1,nrow=1,ncol=n)
	#calc nna vector if not provided
	if(is.null(nna)) nna <- rep(n, p) - rowSums(is.na(x))
	nnam <- (1/(nna-1))
	xdif <- x - (meanx %*% Y)
	ans <- rowSums(xdif^2, na.rm=T) * nnam
	ans[nnam==1] <- NA
	drop(ans)
}

.findfdrcutoff <- function(permdata, realdata, fdr=.01, tails=c("both", "up", "down")) {
	tails <- match.arg(tails)

	if(tails == "down") {
		realdata <- -realdata
		permdata <- -permdata
	}

	if(tails == "both") {
		permdata <- abs(permdata)
	}

	#sort by default removes NA values
	sr <- switch(tails, 
	   both=sort(abs(realdata), decreasing=T),
	   sort(realdata, decreasing=T))

	l <- length(sr)
	#pick the 50% point to test for fdr
	p <- floor(l/2)

	cat("Start fdr test value =", sr[p])
	tfdr <- mean(colSums(permdata > sr[p], na.rm=T))/p
	cat(sr[p]," fdr =", tfdr,"\n")
	stepfraction <- 0.5
	stepsize <- ceiling(p*stepfraction)
	up <- 0
	while(stepsize > 1) {
		#up or down?
		newup <- ifelse(( tfdr - fdr) > 0,0,1)
		if(newup != up) stepfraction <- stepfraction/2
		up <- newup
		stepsize <- ceiling(p*stepfraction)

		p <- ifelse(up == 0,p-stepsize, p+stepsize)
		tfdr <- mean(colSums(abs(permdata) > sr[p], na.rm=T))/p
		cat("p=", p, "step=", stepsize, ", ", sr[p], " fdr =", tfdr,"\n")
	}
	return (sr[p])
}


.getRegions <- function(v, allowedGapsize=10) {
		ss <- c(0,v) - c(v,0)
		starts <- which(ss == -1)
		ends <- which(ss == 1)-1
		if(length(starts) == 0) {
			return(data.frame())
		}
		
		if(length(starts) == 1) {
			return(data.frame(startrow=starts, endrow=ends))
		}
		gapsizes <- c(0, starts[-1] - ends[-length(ends)] -1)

		joinedstarts <- numeric()
		joinedends <- numeric()
		pos <- 1
		ingap <- 0
		
		if(allowedGapsize > 0) {
			for(i in 1:(length(starts)-1)) {
				if(ingap == 0) joinedstarts[pos] <- starts[i]

				if(gapsizes[i+1] <= allowedGapsize) {
					joinedends[pos] <- ends[i+1]
					ingap <- 1
				}
				else {
					joinedends[pos] <- ends[i]
					ingap <- 0
					pos <- pos+1
					if(i == (length(starts)-1)) {
						joinedstarts[pos] <- starts[i+1]
						joinedends[pos] <- ends[i+1]
					}
				}
			}
			return(data.frame(startrow=joinedstarts, endrow=joinedends)) 
		}
	return(data.frame(startrow=starts, endrow=ends)) 
}


