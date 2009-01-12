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


.findPeaks <- function(data, mode='pos'){
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
		
		gc()
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
