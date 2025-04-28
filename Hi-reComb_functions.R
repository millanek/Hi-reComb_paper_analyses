
loadGenomes <- function() {
	fAstCal14ChrSizes <- read.table("genomes/fAstCal14_final.genome")[1:22,]; 
	fAstCal14ChrSizes <- cbind(fAstCal14ChrSizes, cumsum(fAstCal14ChrSizes[,2])); names(fAstCal14ChrSizes) <- c("chr","size","cumSize")

	fNeoMul12ChrSizes <- read.table("genomes/fNeoMul12_final.genome")[1:22,]; 
	fNeoMul12ChrSizes <- cbind(fNeoMul12ChrSizes, cumsum(fNeoMul12ChrSizes[,2])); names(fNeoMul12ChrSizes) <- c("chr","size","cumSize")
	fNeoMul12ChrSizesAll <- read.table("genomes/fNeoMul12_final.genome"); 
	fNeoMul12ChrSizesAll <- cbind(fNeoMul12ChrSizesAll, cumsum(fNeoMul12ChrSizesAll[,2])); names(fNeoMul12ChrSizesAll) <- c("chr","size","cumSize")

	stickGenome <- read.table("genomes/GCF_016920845.1_GAculeatus_UGA_version5_genomic.genome")[1:21,]; 
	stickGenome <- cbind(stickGenome, cumsum(stickGenome[,2])); names(stickGenome) <- c("chr","size","cumSize")
	
	genomes <- list(fAstCal14ChrSizes, fNeoMul12ChrSizes, stickGenome)
	return(genomes)
}

# Rates manually copied from Hi-reComb Simulate output
loadSimulationErrorRates <- function() {
	fpEstimates05pc <-c(0.00619463, 0.00549218, 0.00445164, 0.00495447, 0.00499885, 0.00561732, 0.00569197, 0.00636681, 0.00416328, 0.00560568)
	fpEstimates1pc <- c(0.0100932, 0.0100579, 0.00815174, 0.012543, 0.00638092, 0.0101994, 0.00949528, 0.0092311, 0.0107887, 0.0117928)
	fpEstimates2pc <- c(0.0215488, 0.021312, 0.0195243, 0.0187395, 0.0173437, 0.0195534, 0.0208934, 0.0215381, 0.0205156, 0.0244012)
	fpEstimates5pc <- c(0.0550622,0.0500238, 0.0504625, 0.0594678, 0.0512728, 0.0489242, 0.0516474, 0.0502548, 0.0469761, 0.0456089)
	return(list(fpEstimates05pc, fpEstimates1pc, fpEstimates2pc, fpEstimates5pc))
}

# Rates manually copied from Hi-reComb RecombMap output
loadEmpiricalErrorRates <- function() {
	fpEstimatesAulStu <- c(0.011785, 0.0152884, 0.0417516, 0.0137903, 0.0102613, 0.0107384, 0.0108437, 0.00830416, 0.0136961, 0.0119457, 0.00884437, 0.0176149, 0.00579438, 0.0147112, 0.0135998, 0.0107724, 0.0080795, 0.013282, 0.00764981, 0.00965774, 0.0146057, 0.00993886)
	fpEstimatesAstCal <- c(0.00924034, 0.014941, 0.0440328, 0.035085, 0.0176831, 0.0149784, 0.0124832, 0.0164024, 0.0151912, 0.0132663, 0.020781, 0.0427739, 0.0141437, 0.0163788, 0.013343, 0.0135502, 0.0135397, 0.016274, 0.0173753, 0.00920529, 0.0287199, 0.018175)
	fpEstimatesAstNub <- c(0.00718598,0.00681244,0.017547,0.0095622,0.00671574,0.00701416,0.00799517,0.00610466, 0.00795083, 0.00809538, 0.00794872, 0.00825306, 0.00764594, 0.0082155, 0.00792008, 0.00659299, 0.00710879, 0.00778075, 0.00523764, 0.00683005, 0.00883134, 0.00664865)
	fpEstimatesT2 <- c(0.0145719, 0.0104147, 0.0476935, 0.0200516, 0.00740117, 0.00990681, 0.00911055, 0.0131951, 0.0187807,  0.0111496, 0.00919284, 0.0118015, 0.00582867, 0.00656461, 0.0143198, 0.00809787, 0.0121379, 0.0105535, 0.00971165, 0.0110106)
	fpEstimatesWB1 <- c(0.00701286, 0.0065505, 0.00693688, 0.00642361, 0.00636616, 0.00613475, 0.00586282, 0.00548617, 0.00572161, 0.00581034, 0.00820836, 0.00617624, 0.00650975, 0.0073319, 0.00571139, 0.00665045, 0.00618253, 0.00483517, 0.00387824, 0.00594342, 0.00652188)
	fpEstimatesWB1trio <- c(0.00373433, 0.00319813, 0.00375655, 0.00363486, 0.0051728, 0.00333759, 0.00321091, 0.0026603, 0.00311272, 0.00307625, 0.00502331, 0.00399733, 0.00287288, 0.00300303, 0.00358339, 0.00407407, 0.00321577, 0.00340345, 0.00252107, 0.00398956, 0.00399984)
	allRates <- list(fpEstimatesAulStu, fpEstimatesAstCal, fpEstimatesAstNub, fpEstimatesT2, fpEstimatesWB1, fpEstimatesWB1trio)
	names(allRates) <- c("A. stuartgranti", "A. calliptera", "A. nubila", "N. brichardi", "stickleback", "sticklebackTrio")
	return(allRates)
}

plotTenReconstructions <- function(mapAulStuR, tit="", maxY=3e-7) {
	plot(mapAulStuR[,1], mapAulStuR[,4],type='l',xlab="chr 2 (LS420020) coordinates (Mb)",ylab="r (per bp)", main=paste("Simulate and reconstruct", tit), ylim=c(0, maxY),col= "gray50",lwd=0.3,xaxt='n')
	axis(1,at=c(0,10000000,20000000,30000000),labels=c(0,10,20,30))
	#lines(mapAulStuR[,1], mapAulStuR[,5], col=cols[5-3],lwd=0.5)
	lines(mapAulStuR[,1], mapAulStuR[,5], col="gray50",lwd=0.3)
	#for (i in 5:13) { lines(mapAulStuR[,1], mapAulStuR[,i], col=cols[i-3],lwd=0.5) } 
	for (i in 5:13) { lines(mapAulStuR[,1], mapAulStuR[,i], col="gray50",lwd=0.3) } 
	lines(mapAulStuR[,1], mapAulStuR[,3], col="black",lwd=2)
	legend("top",c("original", "reconstructed"), col=c("black", "gray50"),lwd=c(2,1),lty=1,bty="n")
}

getCorrelationsWithRefMapAtDifferentScales <- function(mapsSimRec) {
	library("zoo")
	corAcrossScales <- cor(mapsSimRec[,3:13],use="pairwise",method="spearman")[1,]
	corAcrossScales <- rbind(corAcrossScales, cor(rollmean(mapsSimRec[,3:13],5),method="spearman",use="pairwise")[1,]) # 10k
	corAcrossScales <- rbind(corAcrossScales, cor(rollmean(mapsSimRec[,3:13],25),method="spearman",use="pairwise")[1,]) # 50k
	corAcrossScales <- rbind(corAcrossScales, cor(rollmean(mapsSimRec[,3:13],50),method="spearman",use="pairwise")[1,]) # 100k
	corAcrossScales <- rbind(corAcrossScales, cor(rollmean(mapsSimRec[,3:13],125),method="spearman",use="pairwise")[1,]) # 250k
	corAcrossScales <- rbind(corAcrossScales, cor(rollmean(mapsSimRec[,3:13],250),method="spearman",use="pairwise")[1,]) # 500k
	corAcrossScales <- rbind(corAcrossScales, cor(rollmean(mapsSimRec[,3:13],500),method="spearman",use="pairwise")[1,]) # 1Mb
	corAcrossScales <- rbind(corAcrossScales, cor(rollmean(mapsSimRec[,3:13],1000),method="spearman",use="pairwise")[1,]) # 2Mb
	corAcrossScales <- rbind(corAcrossScales, cor(rollmean(mapsSimRec[,3:13],1500),method="spearman",use="pairwise")[1,]) # 3Mb
	corAcrossScales <- rbind(corAcrossScales, cor(rollmean(mapsSimRec[,3:13],2000),method="spearman",use="pairwise")[1,]) # 4Mb
	corAcrossScales <- rbind(corAcrossScales, cor(rollmean(mapsSimRec[,3:13],2500),method="spearman",use="pairwise")[1,]) # 5Mb
	corAcrossScales <- corAcrossScales[,-1]
	rownames(corAcrossScales) <- c("2kb","10kb","50kb","100kb","250kb","500kb","1Mb","2Mb","3Mb","4Mb","5Mb")
	colnames(corAcrossScales) <- paste("r",1:10,sep="")
	return(corAcrossScales)
}

plotCorrelations <- function(correlationsList,lg,limYbottom=0.5,limYtop=1) {
	cols <- palette.colors(palette = "Okabe-Ito")[-1]
	corScales <- c(2000,10000,50000,100000,250000,500000,1000000,2000000,3000000,4000000,5000000)
	sdTop <- apply(correlationsList[[1]],1,mean) + apply(correlationsList[[1]],1,sd)*0.5
	sdBottom <- apply(correlationsList[[1]],1,mean) - apply(correlationsList[[1]],1,sd)*0.5
	plot(corScales, apply(correlationsList[[1]],1,mean),type='b',ylim=c(limYbottom, limYtop),ylab="Spearman correlation",xlab="Map scale", col=cols[1],pch=16,log='x',xaxt='n')
	axis(1,at=c(2000,10000,50000,100000,250000,500000,1000000,2000000,5000000),labels=c("2kb","10kb","50kb","100kb","250kb","500kb","1Mb","2Mb","5Mb"))
	for (j in 1:length(corScales)) {
		lines(c(corScales[j],corScales[j]),c(sdTop[j], sdBottom[j]), col=cols[1])
	} 
	
	for (i in 2:length(correlationsList)) {
		sdTop <- apply(correlationsList[[i]],1,mean) + apply(correlationsList[[i]],1,sd)*0.5
		sdBottom <- apply(correlationsList[[i]],1,mean) - apply(correlationsList[[i]],1,sd)*0.5
		lines(corScales, apply(correlationsList[[i]],1,mean),type='b', col=cols[i],pch=16)
		for (j in 1:length(corScales)) {
			lines(c(corScales[j],corScales[j]),c(sdTop[j], sdBottom[j]), col=cols[i])
		} 
	}
	legend("bottomright",lg,col=cols,lty=1,pch=16)
}

BinMean <- function (vec, every, na.rm = FALSE) {
  n <- length(vec)
  x <- .colMeans(vec, every, n %/% every, na.rm)
  r <- n %% every
  if (r) x <- c(x, mean.default(vec[(n - r + 1):n], na.rm = na.rm))
  x
 }
 
plotAndLoadCichlidMap <- function(thisSample="H5a2",m="", doPlot=TRUE) { 
	map_all <- list(22);
	if(doPlot) { pdf(paste0("output_plots/reconstructed_maps/", thisSample, m, "_all.pdf"),width=18,height=40) }
	if(doPlot) { par(mfcol=c(11,2),oma=c(0,0,2,0)) }
	for(i in c(1:22)) {
		map <- read.table(paste0("maps/reconstructed/", thisSample, m, "_chr",i,".txt"),row.names=NULL,header=FALSE)
		map_all[[i]] <- map;
		if(doPlot) { plot(map[,1],map[,3],type='l',xlab=paste0("chr ",i," coordinates"),ylab="rho (per bp)",ylim=c(0,3e-7)) }
	}
	if(doPlot) { dev.off() }
	return(map_all)
}


plotAndLoadGasAcuMap <- function(trio="") { 
	mapGasAcu_all <- list(20);
	pdf(paste0("output_plots/reconstructed_maps/WB1", trio, "_all.pdf"),width=18,height=40)
	par(mfcol=c(10,2),oma=c(0,0,2,0))
	for(i in c(1:18,20,21)) {
	map <- read.table(paste0("maps/reconstructed/WB1", trio, "_chr",i,".txt"),row.names=NULL,header=FALSE)
	if (i <= 18) { mapGasAcu_all[[i]] <- map;} else { mapGasAcu_all[[i-1]] <- map;}
	plot(map[,1],map[,3],type='l',xlab=paste0("chr ",i," coordinates"),ylab="rho (per bp)",ylim=c(0,3e-7))
	}
	dev.off()
	return(mapGasAcu_all)
}

plotAndLoadCichlidBootstrap <- function(thisSample="H5a2") { 
	map_all <- list(22);
	pdf(paste0("output_plots/", thisSample, "_bootstrap_all.pdf"),width=18,height=40)
	par(mfcol=c(11,2),oma=c(0,0,2,0))
	for(i in c(1:22)) {
		mapNAul <- read.table(paste0("maps/reconstructed/wBootstrap/", thisSample, "_chr",i,"_boot.txt"),row.names=NULL,header=FALSE)
		pct95<-apply(mapNAul[,4:dim(mapNAul)[2]],1,quantile, 0.95,na.rm=T)
		pct5<-apply(mapNAul[,4:dim(mapNAul)[2]],1,quantile, 0.05,na.rm=T)
		avgBoot<-apply(mapNAul[,4:dim(mapNAul)[2]],1,median,na.rm=T)
		boot_chr <- cbind(mapNAul[,1:3], avgBoot, pct95, pct5)
		map_all[[i]] <- boot_chr
		plot(mapNAul[,1], avgBoot,type='l',xlab=paste0("chr ",i," coordinates"),ylab="r (per bp)",ylim=c(0,3e-7),col="black",lwd=2)
		lines(mapNAul[,1], pct95, lty=2,lwd=1,col="gray70")
		lines(mapNAul[,1], pct5, lty=2,lwd=1,col="gray70")
	}
	dev.off()
	return(map_all)
}

getMeanRates <- function(mapAulStu, mapAstCal, mapAstNub, mapNeoMul, mapGasAcu, mapGasAcuTrio) {
	meanRatesAulStu <- numeric(0); meanRatesAstCal <- numeric(0); meanRatesAstNub <- numeric(0); meanRatesNeoMul <- numeric(0);
	meanRatesGasAcu <- numeric(0); meanRatesGasAcuTrio <- numeric(0);
for (i in 1:22) { 
	meanRatesAulStu <- c( meanRatesAulStu, mean(mapAulStu[[i]][,3],na.rm=T) ) 
	meanRatesAstCal <- c( meanRatesAstCal, mean(mapAstCal[[i]][,3],na.rm=T) ) 
	meanRatesAstNub <- c( meanRatesAstNub, mean(mapAstNub[[i]][,3],na.rm=T) ) 
	meanRatesNeoMul <- c( meanRatesNeoMul, mean(mapNeoMul[[i]][,3],na.rm=T) ) 
	if (i <= 20) {
		meanRatesGasAcu <- c( meanRatesGasAcu, mean(mapGasAcu[[i]][,3],na.rm=T) )
		meanRatesGasAcuTrio <- c( meanRatesGasAcuTrio, mean(mapGasAcuTrio[[i]][,3],na.rm=T) )
	}
	}
	meanRates <- list(meanRatesAulStu, meanRatesAstCal, meanRatesAstNub, meanRatesNeoMul, meanRatesGasAcu, meanRatesGasAcuTrio)
	names(meanRates) <- c("AulStu","AstCal","AstNub","NeoMul","GasAcu","GasAcuTrio")
	return(meanRates)
}

plotChrSizeVsMeanRateCichlid <- function(genomes, meanRates) {
	par(mfrow=c(2,2))
	lmAulStu <- lm(meanRates$AulStu~genomes$fAstCal[,2])
	plot(genomes$fAstCal[,2], meanRates$AulStu, col=colsSpecies[1],pch=16, xlab="Chromosome size", ylab="Mean r")
	abline(lmAulStu,col=colsSpecies[1])
	plot(genomes$fAstCal[,2], meanRates$AstNub, col=colsSpecies[2],pch=16, xlab="Chromosome size", ylab="Mean r")
	lmAstNub <-lm(meanRates$AstNub~genomes$fAstCal[,2])
	abline(lmAstNub,col=colsSpecies[2])
	plot(genomes$fNeoMul[,2], meanRates$NeoMul, col=colsSpecies[3],pch=16, xlab="Chromosome size", ylab="Mean r")
	lmNeoMul <- lm(meanRates$NeoMul~genomes$fNeoMul[,2])
	abline(lmNeoMul,col=colsSpecies[3])
	plot(genomes$fAstCal[,2], meanRates$AstCal, col=colsSpecies[4],pch=16, xlab="Chromosome size", ylab="Mean r")
	lmAstCalAll <- lm(meanRates$AstCal~genomes$fAstCal[,2])
	abline(lmAstCalAll,col=colsSpecies[4])
	lmAstCalNoBig <- lm(meanRates$AstCal[-c(3,7)]~genomes$fAstCal[-c(3,7),2])
	abline(lmAstCalNoBig,col=colsSpecies[4])
	allLms <- list(lmAulStu, lmAstNub, lmNeoMul, lmAstCalAll, lmAstCalNoBig)
	names(allLms) <- c("AulStu", "AstNub", "NeoMul", "AstCalAll", "AstCalNoBig")
	return(allLms)
}


getMapLengths_Cm <- function(mapAulStu, mapAstCal, mapAstNub, mapNeoMul, mapGasAcu, mapGasAcuTrio,windowSize=2000) {
	lengthsAulStu <- numeric(0); lengthsAstCal <- numeric(0); lengthsAstNub <- numeric(0); lengthsNeoMul <- numeric(0);
	lengthsGasAcu <- numeric(0); lengthsGasAcuTrio <- numeric(0);
for (i in 1:22) { 
	lengthsAulStu <- c( lengthsAulStu, sum(mapAulStu[[i]][,3]*windowSize,na.rm=T)*100 ) 
	lengthsAstCal <- c( lengthsAstCal, sum(mapAstCal[[i]][,3]*windowSize,na.rm=T)*100 ) 
	lengthsAstNub <- c( lengthsAstNub, sum(mapAstNub[[i]][,3]*windowSize,na.rm=T)*100 ) 
	lengthsNeoMul <- c( lengthsNeoMul, sum(mapNeoMul[[i]][,3]*windowSize,na.rm=T)*100 ) 
	if (i <= 20) {
		lengthsGasAcu <- c( lengthsGasAcu, sum(mapGasAcu[[i]][,3]*windowSize,na.rm=T)*100 )
		lengthsGasAcuTrio <- c( lengthsGasAcuTrio, sum(mapGasAcuTrio[[i]][,3]*windowSize,na.rm=T)*100 )
	}
	}
	lengths_Cm <- list(lengthsAulStu, lengthsAstCal, lengthsAstNub, lengthsNeoMul, lengthsGasAcu, lengthsGasAcuTrio)
	names(lengths_Cm) <- c("AulStu","AstCal","AstNub","NeoMul","GasAcu","GasAcuTrio")
	return(lengths_Cm)
}


getMapLengths_Cm_linkageMaps <- function(AlbertsonMap, KocherMap, roestiMap) {
	LGlengthsAlbertson <- aggregate(cM~LG, AlbertsonMap,max); 
	LGlengthsKocher <- aggregate(cM~LG, KocherMap,max); 
	LGlengthsRoesti <- aggregate(cM~chromosome_reassembled, roestiMap,max);
	LGlengths_all <- list(LGlengthsAlbertson[,2], LGlengthsKocher[,2], LGlengthsRoesti[-19,2])
	names(LGlengths_all) <- c("Albertson", "Kocher", "Roesti")
	return(LGlengths_all)
}

plotBootRegion <- function(boot1Chr,boot2Chr,start,end,bootYmin=0,bootYmax=3e-7, bootLty=1,meanAdjust=FALSE, zScore=FALSE) {
	
	if(meanAdjust) {
		boot1Chr$pct95 <- (boot1Chr$pct95 - mean(boot1Chr$avgBoot,na.rm=T))
		boot1Chr$pct5 <- (boot1Chr$pct5 - mean(boot1Chr$avgBoot,na.rm=T))
		boot1Chr$avgBoot <- (boot1Chr$avgBoot - mean(boot1Chr$avgBoot,na.rm=T))
		boot2Chr$pct95 <- (boot2Chr$pct95 - mean(boot2Chr$avgBoot,na.rm=T))
		boot2Chr$pct5 <- (boot2Chr$pct5 - mean(boot2Chr$avgBoot,na.rm=T))
		boot2Chr$avgBoot <- (boot2Chr$avgBoot - mean(boot2Chr$avgBoot,na.rm=T))
	} else { 
		if(zScore) {
			boot1Chr$pct95 <- (boot1Chr$pct95 - mean(boot1Chr$avgBoot,na.rm=T))/sd(boot1Chr$avgBoot,na.rm=T)
			boot1Chr$pct5 <- (boot1Chr$pct5 - mean(boot1Chr$avgBoot,na.rm=T))/sd(boot1Chr$avgBoot,na.rm=T)
			boot1Chr$avgBoot <- (boot1Chr$avgBoot - mean(boot1Chr$avgBoot,na.rm=T))/sd(boot1Chr$avgBoot,na.rm=T)
			boot2Chr$pct95 <- (boot2Chr$pct95 - mean(boot2Chr$avgBoot,na.rm=T))/sd(boot2Chr$avgBoot,na.rm=T)
			boot2Chr$pct5 <- (boot2Chr$pct5 - mean(boot2Chr$avgBoot,na.rm=T))/sd(boot2Chr$avgBoot,na.rm=T)
			boot2Chr$avgBoot <- (boot2Chr$avgBoot - mean(boot2Chr$avgBoot,na.rm=T))/sd(boot2Chr$avgBoot,na.rm=T)
		}
	}
	
	plot(boot1Chr[start:end,1], boot1Chr[start:end,4],type='l',xlab=paste0("Coordinates"),ylab="r (per bp)",ylim=c(bootYmin, bootYmax),col=colsSpecies[1],lwd=1)
	#lines(boot1Chr[start:end,1], boot1Chr$pct95[start:end], lty=bootLty,lwd=1,col=colsSpecies[1])
	#lines(boot1Chr[start:end,1], boot1Chr$pct5[start:end], lty=bootLty,lwd=1,col=colsSpecies[1])
	val <- which(!is.na(boot1Chr[,6]))
	polygon(c(boot1Chr[val,1], rev(boot1Chr[val,1])),c(boot1Chr[val,]$pct95,rev(boot1Chr[val,]$pct5)),col=rgb(col2rgb(colsSpecies[1])[1]/255,col2rgb(colsSpecies[1])[2]/255,col2rgb(colsSpecies[1])[3]/255,alpha=0.5),border=NA)
	lines(boot2Chr[start:end,1], boot2Chr[start:end,4],lwd=1,col=colsSpecies[2])
	#lines(boot2Chr[start:end,1], boot2Chr$pct95[start:end], lty=bootLty,lwd=1,col=colsSpecies[2])
	#lines(boot2Chr[start:end,1], boot2Chr$pct5[start:end], lty= bootLty,lwd=1,col=colsSpecies[2])
	val <- which(!is.na(boot2Chr[,6]))
	polygon(c(boot2Chr[val,1], rev(boot2Chr[val,1])),c(boot2Chr[val,5],rev(boot2Chr[val,6])),col=rgb(col2rgb(colsSpecies[2])[1]/255,col2rgb(colsSpecies[2])[2]/255,col2rgb(colsSpecies[2])[3]/255,alpha=0.5),border=NA)
}

plotBootOneMap <- function(boot1Chr,start,end,bootYmin=0,bootYmax=3e-7, bootLty=1) {
	plot(boot1Chr[start:end,1], boot1Chr[start:end,4],type='l',xlab=paste0("Coordinates"),ylab="r (per bp)",ylim=c(bootYmin, bootYmax),col=colsSpecies[1],lwd=1)
	#lines(boot1Chr[start:end,1], boot1Chr$pct95[start:end], lty=bootLty,lwd=1,col=colsSpecies[1])
	#lines(boot1Chr[start:end,1], boot1Chr$pct5[start:end], lty=bootLty,lwd=1,col=colsSpecies[1])
	val <- which(!is.na(boot1Chr[,6]))
	polygon(c(boot1Chr[val,1], rev(boot1Chr[val,1])),c(boot1Chr[val,]$pct95,rev(boot1Chr[val,]$pct5)),col=rgb(col2rgb(colsSpecies[1])[1]/255,col2rgb(colsSpecies[1])[2]/255,col2rgb(colsSpecies[1])[3]/255,alpha=0.5),border=NA)
}

makeSigDiffRegions <- function(x) {
  diffs <- c(1, diff(x))
  start_indexes <- c(1, which(diffs > 1))
  end_indexes <- c(start_indexes - 1, length(x))
  startEnd <- cbind(x[start_indexes], x[end_indexes], x[end_indexes]-x[start_indexes]+1); colnames(startEnd) <- c("start","end","length")
  return(startEnd)
}

getAllSigDiffRegionsPerChr <- function(boot1Chr,boot2Chr,ma=FALSE) {
	boot1ChrL <- dim(boot1Chr)[1]; boot2ChrL <- dim(boot2Chr)[1]; sharedL <- min(boot1ChrL, boot2ChrL)
	if(ma) {
		boot1Chr$pct95 <- boot1Chr$pct95 - mean(boot1Chr$avgBoot,na.rm=T)
		boot1Chr$pct5 <- boot1Chr$pct5 - mean(boot1Chr$avgBoot,na.rm=T)
		boot2Chr$pct95 <- boot2Chr$pct95 - mean(boot2Chr$avgBoot,na.rm=T)
		boot2Chr$pct5 <- boot2Chr$pct5 - mean(boot2Chr$avgBoot,na.rm=T)
	}
	boot1Lower <- which(boot1Chr[1:sharedL,]$pct95 < boot2Chr[1:sharedL,]$pct5)
	boot1Higher <- which(boot1Chr[1:sharedL,]$pct5 > boot2Chr[1:sharedL,]$pct95)
	sigDiffs <- list(boot1Lower, boot1Higher)
	sigDiffRegions <- lapply(sigDiffs, makeSigDiffRegions)
	sigDiffRegionsDf <- rbind(cbind(sigDiffRegions[[1]],rep(1,dim(sigDiffRegions[[1]])[1])), cbind(sigDiffRegions[[2]],rep(2,dim(sigDiffRegions[[2]])[1]))) 
	sigDiffRegionsDf <- sigDiffRegionsDf[order(sigDiffRegionsDf[,1]),]
	colnames(sigDiffRegionsDf) <- c(colnames(sigDiffRegionsDf)[1:3],"elevated_in")
	return(sigDiffRegionsDf)
}

getAllSigDiffRegions <- function(boot1, boot2,nChr=22,meanAdjust=FALSE) {
	allChr <- list(nChr)
	for (i in 1:nChr) {
		thisChr <- getAllSigDiffRegionsPerChr(boot1[[i]], boot2[[i]], ma=meanAdjust)
		allChr[[i]] <- thisChr
	}
	return(allChr)
}

getCorrelationsOfTwoMapsAtDifferentScales <- function(map1all, map2all,Nchr=c(1:20)) {
	corAcrossScalesAllChr <- matrix(rep(0,length(Nchr)*11),nrow=11,ncol=length(Nchr))
	rownames(corAcrossScalesAllChr) <- c("2kb","10kb","50kb","100kb","250kb","500kb","1Mb","2Mb","3Mb","4Mb","5Mb")
	colnames(corAcrossScalesAllChr) <- 1:length(Nchr)
	for (i in Nchr) {
		print(paste0("chr ", i))
		map1 <- map1all[[i]]; map2 <- map2all[[i]];
		maxL <- min(length(map1[,3]),length(map2[,3]))
		corAcrossScales <- cor(x=map1[1:maxL,3],y=map2[1:maxL,3], use="pairwise",method="spearman")
		corAcrossScales <- c(corAcrossScales, cor(x=rollmean(map1[1:maxL,3],5),y=rollmean(map2[1:maxL,3],5),method="spearman",use="pairwise")) # 10k
		corAcrossScales <- c(corAcrossScales, cor(x=rollmean(map1[1:maxL,3],25),y=rollmean(map2[1:maxL,3],25),method="spearman",use="pairwise")) # 50k
		corAcrossScales <- c(corAcrossScales, cor(x=rollmean(map1[1:maxL,3],50),y=rollmean(map2[1:maxL,3],50),method="spearman",use="pairwise")) # 100k
		corAcrossScales <- c(corAcrossScales, cor(x=rollmean(map1[1:maxL,3],125),y=rollmean(map2[1:maxL,3],125),method="spearman",use="pairwise")) # 250k
		corAcrossScales <- c(corAcrossScales, cor(x=rollmean(map1[1:maxL,3],250),y=rollmean(map2[1:maxL,3],250),method="spearman",use="pairwise")) # 500k
		corAcrossScales <- c(corAcrossScales, cor(x=rollmean(map1[1:maxL,3],500),y=rollmean(map2[1:maxL,3],500),method="spearman",use="pairwise")) # 1Mb
		corAcrossScales <- c(corAcrossScales, cor(x=rollmean(map1[1:maxL,3],1000),y=rollmean(map2[1:maxL,3],1000),method="spearman",use="pairwise")) # 2Mb
		corAcrossScales <- c(corAcrossScales, cor(x=rollmean(map1[1:maxL,3],1500),y=rollmean(map2[1:maxL,3],1500),method="spearman",use="pairwise")) # 3Mb
		corAcrossScales <- c(corAcrossScales, cor(x=rollmean(map1[1:maxL,3],2000),y=rollmean(map2[1:maxL,3],2000),method="spearman",use="pairwise")) # 4Mb
		corAcrossScales <- c(corAcrossScales, cor(x=rollmean(map1[1:maxL,3],2500),y=rollmean(map2[1:maxL,3],2500),method="spearman",use="pairwise")) # 5Mb
		#corAcrossScales <- corAcrossScales[,-1]
		names(corAcrossScales) <- c("2kb","10kb","50kb","100kb","250kb","500kb","1Mb","2Mb","3Mb","4Mb","5Mb")
		corAcrossScalesAllChr[,i] <- corAcrossScales;
		#colnames(corAcrossScales) <- paste("r",1:10,sep="")
	}
	return(corAcrossScalesAllChr)
}

get_GasAcu_WB_LDmap <- function(LDmaps) {
	LDmapList <- list(21)
	for (i in 12:32) {
		LDmaps_chr <- LDmaps[which(LDmaps$chr == paste0("NC_0532",i,".1")),c(2,3,4)]
		LDmapList[[i-11]] <- LDmaps_chr
	}
	return(LDmapList)
}

reformatRoestiMap <- function(roestiMap) {
	mapVec <- numeric(0)
	for (i in 1:21) {
		map_chr <- roestiMap[which(roestiMap$chromosome_reassembled == i),c(5:7)]
		map_chr_r <- ( (map_chr[-1,3] - map_chr[-dim(map_chr)[1],3]) / (map_chr[-1,2] - map_chr[-dim(map_chr)[1],2] + 1) ) / 100
		map_chr_start <- c(0,map_chr[-dim(map_chr)[1],2])
		map_chr_r <- c(NA, map_chr_r)
		map_chr <- cbind(map_chr[,1],map_chr_start,map_chr[,2], map_chr_r)
		colnames(map_chr) <- c("chr","start","end","r")
		if (i == 1) { mapVec <- map_chr; } 
		else { mapVec <- rbind(mapVec, map_chr);}
	}
	return(mapVec)
}


divideRoestiMapByChr <- function(roestiMap) {
	mapList <- list(21)
	for (i in 1:21) {
		map_chr <- roestiMap[which(roestiMap$chr == i),c(2:4)]
		mapList[[i]] <- map_chr
	}
	return(mapList)
}

divideLDmapByChrAstCal <- function(LDmaps) {
	LDmapList <- list(22)
	for (i in 19:40) {
		LDmaps_chr <- LDmaps[which(LDmaps$chr == paste0("LS4200",i,".2")),c(2,3,4)]
		LDmapList[[i-18]] <- LDmaps_chr
	}
	return(LDmapList)
}
