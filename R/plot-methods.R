setMethod("plot", signature(x="scaleSpace", y="missing"), function(x, y, spm, type='b', ...){
	
	
	nrScales <- length(x)
	mirrorLocs <- spm@mirrorLocs
	
	chromosomes <- names(x[[1]]@data)
	spmLengths <- unlist(lapply(spm@data, function(x){return(length(x$pos))}))
	chrom.indices <- unlist(lapply(x[[1]]@data, function(x){return(attr(x, 'chromosome'))}))

	total <- sum(spmLengths[chromosomes])
	
	if(type=='b' | type=='g'){
		plot(0,0,xlim=c(0, total), ylim=c(-0.1,nrScales),col='white', xaxt="n", yaxt="n", main='Scale space gains', xlab='Genomic position (in mb)', ylab='Scale space',...)
		gain.dev <- dev.cur()
	}
	if(type=='b' | type=='l'){
		if(type=='b'){x11()}
		plot(0,0,xlim=c(0, total), ylim=c(-0.1,nrScales),col='white', xaxt="n", yaxt="n", main='Scale space losses',xlab='Genomic position (in mb)',ylab='Scale space', ...)
		loss.dev <- dev.cur()
	}
	
	scale.spaces <- vector(length=nrScales)
	heatcolors <- rev(heat.colors(32))
	
	for(i in 1:nrScales){
		xOffset <- 0
		unlisted <- unlist(x[[i]]@data)
		posy <- unlisted[grep("posy", names(unlisted))]
		negy <- unlisted[grep("negy", names(unlisted))]
		if(length(posy>0)){
			maxy <- max(posy, na.rm=TRUE)
		}
		else{
			maxy <- 1
		}
		if(length(negy)>0){
			miny <- min(negy, na.rm=TRUE)
		}
		else{
			miny <- -1
		}
		
		for(j in chromosomes){	
			
			if(type=='b' | type=='g'){
				dev.set(gain.dev)
				
				poscolors <- (x[[i]][[j]]$posy / maxy) * 32
				poscolors <- heatcolors[poscolors]
				
				xloc <- x[[i]][[j]]$posx + xOffset
				yloc <- rep((i-1), length(x[[i]][[j]]$posx))
				if(length(xloc) > 0){
					segments(xloc,yloc ,xloc, (yloc+0.98), col=poscolors)
				}
			}
			
			if(type=='b' | type=='l'){
				dev.set(loss.dev)
				
				negcolors <- (x[[i]][[j]]$negy / miny) * 32
				negcolors <- heatcolors[negcolors]
				
				xloc <- x[[i]][[j]]$negx + xOffset
				yloc <- rep((i-1), length(x[[i]][[j]]$negx))
				if(length(xloc) > 0){
					segments(xloc,yloc ,xloc, (yloc+0.98), col=negcolors)
				}
			}
			
			xOffset <- xOffset + length(spm[[j]]$pos)
			
			scale.spaces[i] <- paste(x[[i]]@sigma / 1000000, 'Mb')
		}
	}
	
	
	colin <- cumsum(unlist(lapply(mirrorLocs[chrom.indices], max)))
	#get 2nd position out of mirrorLocs, this either the centromere or the end position
	#add these values to the cumsum and use the results to plot centromeres (or end positions that get overwritten by the chromosome borders)
	centromeres <- c(0, colin) + c(unlist(lapply(mirrorLocs[chrom.indices], function(x){return(x[2])})), 0)
	#centromeres <- centromers[-24]
	totalbp <- spm@sampleDensity * as.numeric(total)
	
	axisBy <- 5000
	sampleDensity <- spm@sampleDensity
	
	if(type=='b' | type=='g'){
		dev.set(gain.dev)
		abline(v=centromeres/sampleDensity, col='lightblue', lty=2)
		abline(v=c(0, colin/sampleDensity), col='darkblue')
		text((colin/sampleDensity) ,-0.1, labels=chromosomes, pos=2, cex=0.6)
		axis(1,seq(0,totalbp, by=(axisBy * spm@sampleDensity))/1000000, at=seq(0,total,by=axisBy))
		axis(2,scale.spaces, at=seq(0.5,nrScales,by=1), las=1)
	}
	
	if(type=='b' | type=='l'){
		dev.set(loss.dev)
		abline(v=centromeres/sampleDensity, col='lightblue', lty=2)
		abline(v=colin/spm@sampleDensity, col='darkblue')
		text((colin/sampleDensity) ,-0.1, labels=chromosomes, pos=2, cex=0.6)
		axis(1,seq(0,totalbp, by=(axisBy * spm@sampleDensity))/1000000, at=seq(0,total,by=axisBy))
		axis(2,scale.spaces, at=seq(0.5,nrScales,by=1), las=1)
	}
})

setMethod("plot", signature(x="samplePointMatrix", y="missing"), function(x, y, type="b", sigLevels=NULL, chromosomes=NULL, colinAxis=NULL, fillColor=NULL, maploc=NULL, interpolation=1, ...){
	mirrorLocs <- x@mirrorLocs
	
	if(!is.null(fillColor) & is.null(sigLevels)){
		warning('Fill color given but no significance levels, unable to color significant regions')
	}
	
	if(is.null(fillColor) & !is.null(sigLevels)){
		fillColor <- vector(mode='list')
		fillColor$pos='red'
		fillColor$neg='green'
	}
	
	if(!is.list(fillColor)){
		fillColor=NULL
	}

	if(type == 'b'){
		layout(c(1,2))
	}
	
	total <- x@totalLength
	
	if(is.null(chromosomes)){
		chromosomes <- names(x@data)
		chromNames <- attr(mirrorLocs, 'chromNames')
		chromosomesOrdered <- chromNames[chromNames %in% chromosomes]
		chromosomesOrdered <- c(chromosomesOrdered, chromosomes[!(chromosomes %in% chromNames)])
		chromosomes <- chromosomesOrdered
		#set colinear axis if parameter is not set and showing all chromosomes
		if(is.null(colinAxis)){
			colinAxis <- TRUE
		}
	}
	else{
		chromosomes <- as.character(chromosomes)
		spmLengths <- unlist(lapply(x@data, function(x){return(length(x$pos))}))
		total <- sum(spmLengths[chromosomes])
		#set colinear axis to false if parameter is not set and not showing all chromosomes
		if(is.null(colinAxis)){
			colinAxis <- FALSE
		}
	}
	
	totalbp <- as.numeric(total) * x@sampleDensity
	
	#determine how to scale axis
	if(colinAxis){
		axisBy <- 10^(floor(log10(total)))
	}
	else{
		axisBy <- 10^(floor(log10(total/length(chromosomes))))
	}
	
		
	maxy <- x@maxy
	miny <- x@miny
	
	#gains
	
	if(type == 'g' | type == 'b' | type == 1){

	if(type==1){
		posminy <- miny
	}
	else{
		posminy <- 0
	}
	
	plot(0,0,xlim=c(0, total), ylim=c(posminy,maxy),col='white', xaxt="n", main='Gains', xlab='Genomic position (in mb)', ylab='Normalized KC score', ...)
	xOffset <- 0
	abline(v=xOffset,col='darkblue')
	for(i in chromosomes){
		chromosome <- attr(x@data[[i]], 'chromosome')
		#color area under the curve
		if(!is.null(fillColor) & !is.null(sigLevels)){
			sigRegions <- which(x[[i]]$pos > sigLevels$pos)
			if(length(sigRegions) > 1){
				#make separate polygons by inserting 'NA's between segments
				t <- diff(sigRegions)
				endPoints <- sigRegions[c(1,which(t>1)+1, length(sigRegions))]
				sigRegions2 <- rep(NA,length(sigRegions) + length(endPoints))
				fillHeight <- rep(NA,length(sigRegions) + length(endPoints))
				for(k in 1:length(endPoints)){
					currentSigRegion <- which((sigRegions<endPoints[k+1]) & (sigRegions>=endPoints[k]))
					sigRegions2[currentSigRegion + k] <- sigRegions[currentSigRegion]	
					fillHeight[currentSigRegion + k] <- x[[i]]$pos[sigRegions[currentSigRegion]]
					fillHeight[c((currentSigRegion[1] + k),(tail(currentSigRegion, n=1)+k))] <- sigLevels$pos
				}
				sigRegions <- sigRegions2
				polygon(xOffset + seq(1,length(x[[i]]$pos))[sigRegions], fillHeight, col=fillColor$pos, border=NA)
			}
		}
		#to avoid getting really large images the user can set an interpolation
		plottingPoints <- seq(1,length(x[[i]]$pos), by=interpolation)
		lines(xOffset + plottingPoints, x[[i]]$pos[plottingPoints], type="l")
		
		chromosome.length <- 0
		#if centromere is present, plot it
		if(length(mirrorLocs[[chromosome]]) == 3){
			centromereLoc <- xOffset + ((mirrorLocs[[chromosome]][2]/mirrorLocs[[chromosome]][3]) * length(x[[i]]$pos))
			abline(v=centromereLoc, col='lightblue', lty=2)
			chromosome.length=mirrorLocs[[chromosome]][3]
		}
		if(!is.null(maploc) & FALSE){
			if(chromosome.length<1){
				chromosome.length <- mirrorLocs[[chromosome]][2]
			}	
			
			locs <- xOffset + ((maploc[[i]]/chromosome.length) * length(x[[i]]$pos))
			segments(locs,0,locs,1, col="purple")	
		}
		
		text(xOffset,0.001, labels=i, pos=4, cex=0.6)
		if(!colinAxis){
			axis(1,seq(0, length(x[[i]]$pos), by=axisBy) * x@sampleDensity/1000000, at=seq(xOffset, (xOffset + length(x[[i]]$pos)), by=axisBy))
		}
		xOffset <- xOffset + length(x[[i]]$pos)
		abline(v=xOffset,col='darkblue')
	}
	
	if(!is.null(sigLevels)){
		abline(h=sigLevels$pos, col="red", lty=2)	
	}
	
	if(colinAxis){
		axis(1,seq(0, totalbp, by=(axisBy * x@sampleDensity))/1000000, at=seq(0,total,by=axisBy))
	}
	}
	
	#losses
	if(type == 'l' | type == 'b' | type == 1){
	
	#if not in 1 plot, open new device
	if(type != 1){
		plot(0,0,xlim=c(0, total), ylim=c(miny, 0), col='white', xaxt="n",main='Losses', xlab='Genomic position (in mb)', ylab='Normalized KC score', ...)
	}
	xOffset <- 0
	abline(v=xOffset,col='darkblue')
	for(i in chromosomes){
		chromosome <- attr(x[[i]], 'chromosome')
		
		if(!is.null(fillColor) & !is.null(sigLevels)){
			sigRegions <- which(x[[i]]$neg < sigLevels$neg)
			if(length(sigRegions) > 1){
				t <- diff(sigRegions)
				endPoints <- sigRegions[c(1,which(t>1)+1, length(sigRegions))]
				sigRegions2 <- rep(NA,length(sigRegions) + length(endPoints))
				fillHeight <- rep(NA,length(sigRegions) + length(endPoints))
				for(k in 1:length(endPoints)){
					currentSigRegion <- which((sigRegions<endPoints[k+1]) & (sigRegions>=endPoints[k]))
					sigRegions2[currentSigRegion + k] <- sigRegions[currentSigRegion]	
					fillHeight[currentSigRegion + k] <- x[[i]]$neg[sigRegions[currentSigRegion]]
					fillHeight[c((currentSigRegion[1] + k),(tail(currentSigRegion,n=1)+k))] <- sigLevels$neg
				}
				sigRegions <- sigRegions2
				polygon(xOffset + seq(1,length(x[[i]]$neg))[sigRegions], fillHeight, col=fillColor$neg, border=NA)
			}
		}
		
		plottingPoints <- seq(1,length(x[[i]]$neg), by=interpolation)
		lines(xOffset + plottingPoints, x[[i]]$neg[plottingPoints], type="l")
		
		#if centromere is present, plot it
		if(length(mirrorLocs[[chromosome]]) == 3){
			centromereLoc <- xOffset + ((mirrorLocs[[chromosome]][2]/mirrorLocs[[chromosome]][3]) * length(x[[i]]$neg))
			abline(v=centromereLoc, col='lightblue', lty=2)
		}
		if(type != 1){
			text(xOffset,-0.001, labels=i, pos=4, cex=0.6)
		}
		if(!colinAxis){
			axis(1,seq(0, length(x[[i]]$pos), by=axisBy) * x@sampleDensity/1000000, at=seq(xOffset, (xOffset + length(x[[i]]$pos)), by=axisBy))
		}
		xOffset <- xOffset + length(x[[i]]$neg)
		abline(v=xOffset, col='darkblue')
	}
	
	if(!is.null(sigLevels)){
		abline(h=sigLevels$neg, col="red", lty=2)	
	}
	
	if(colinAxis){
		axis(1,seq(0,totalbp, by=(axisBy * x@sampleDensity))/1000000, at=seq(0,total,by=axisBy))
	}
	}
})