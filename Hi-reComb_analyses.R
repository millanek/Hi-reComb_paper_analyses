#############################################################################
#### Author: Milan Malinsky
#### Analyses of recombination for the paper
#### Malinsky, M. et al. Hi-reComb: constructing recombination maps from bulk gamete Hi-C sequencing. bioRxiv (2025) doi:10.1101/2025.03.06.641907. 
#### Last edit: 27th April 2025

#############################################################################
#### Pre-load functions and data
library("zoo")
setwd("~/Documents/Hi-reComb_analyses/")
source("Hi-reComb_functions.R")
colsSpecies <- palette.colors(palette = "R4")[-c(1,3,5)]
genomes <- loadGenomes(); names(genomes) <- c("fAstCal","fNeoMul","GasAcu")
hetSites <- c(1000705, 811335, 1137937, 1342749, 1344307); names(hetSites) <- c("AstCal","AulStu","GasAcu","NeoMul","AstNub")
# Linkage-based maps for comparison:
AlbertsonMap <- read.table("publishedMaps/AlbertsonMap.csv"); KocherMap <- read.table("publishedMaps/KocherMap.csv")
names(AlbertsonMap) <- c("sc","pos","LG","cM"); names(KocherMap) <- c("sc","pos","LG","cM"); 
roestiMap <- read.table("publishedMaps/roesti_map.txt",row.names=NULL,header=TRUE)
# LD-based maps for comparison:
LDmapStickleback <- read.table("maps/LD/Gacul_pyrho_2Kb_WB.txt",row.names=NULL,header=F)
names(LDmapStickleback) <- c("chr","start","end", "WB_a")

#############################################################################
#### 1) Analyses of Hi-reComb Simulate
#############################################################################
# 1.1) Error rate estimation (Figure 2D)
loadSimulationErrorRates();
pdf("output_plots/simulations_error_rate.pdf",4,9)
boxplot(fpEstimates05pc, fpEstimates1pc, fpEstimates2pc, fpEstimates5pc,ylim=c(0,0.06),names=c("0.5%","1%","2%","5%"),xlab="Simulated error rates",ylab="Inferred error rates")
dev.off()

# 1.2) AulStu chr2 reference; 1pc error, varying effective coverage
# Figure 2B
mapAulStuR_500x <- read.table("maps/simulated/simulatedMaps_241210_Aul_500_1pc_FW_2000.txt",row.names=NULL,header=FALSE)
mapAulStuR_1000x <- read.table("maps/simulated/simulatedMaps_241210_Aul_1000_1pc_FW_2000.txt",row.names=NULL,header=FALSE)
mapAulStuR_2000x <- read.table("maps/simulated/simulatedMaps_241210_Aul_2000_1pc_FW_2000.txt",row.names=NULL,header=FALSE)
mapAulStuR_3000x <- read.table("maps/simulated/simulatedMaps_241210_Aul_3000_1pc_m20_FW_2000.txt",row.names=NULL,header=FALSE)
plotTenReconstructions(mapAulStuR_500x,"1pc 500x")
plotTenReconstructions(mapAulStuR_1000x,"1pc 1000x")
plotTenReconstructions(mapAulStuR_2000x,"1pc 2000x")
pdf("output_plots/simulations_1pc_3000x_reconstruction.pdf",9,5)
plotTenReconstructions(mapAulStuR_3000x,"1pc 3000x") ### And this is Figure 2A
dev.off()
plotTenReconstructions(mapAulStuR_3000x_m25,"1pc 3000x")
cor_500x_Aul <- getCorrelationsWithRefMapAtDifferentScales(mapAulStuR_500x)
cor_1000x_Aul <- getCorrelationsWithRefMapAtDifferentScales(mapAulStuR_1000x)
cor_2000x_Aul <- getCorrelationsWithRefMapAtDifferentScales(mapAulStuR_2000x)
cor_3000x_Aul <- getCorrelationsWithRefMapAtDifferentScales(mapAulStuR_3000x)
pdf("output_plots/simulations_1pc_depthDependence.pdf",9,9)
plotCorrelations(list(cor_3000x_Aul, cor_2000x_Aul, cor_1000x_Aul, cor_500x_Aul),c("3000x","2000x","1000x","500x"))
dev.off()

# 1.3) AulStu chr2 reference; 1000x effective coverage; varying error rates
# Figure 2C
mapAulStuR_05pc <- read.table("maps/simulated/simulatedMaps_241210_Aul_1000_05pc_FW_2000.txt",row.names=NULL,header=FALSE)
mapAulStuR_1pc <- read.table("maps/simulated/simulatedMaps_241210_Aul_1000_1pc_FW_2000.txt",row.names=NULL,header=FALSE)
mapAulStuR_2pc <- read.table("maps/simulated/simulatedMaps_241210_Aul_1000_2pc_FW_2000.txt",row.names=NULL,header=FALSE)
mapAulStuR_5pc <- read.table("maps/simulated/simulatedMaps_241210_Aul_1000_5pc_FW_2000.txt",row.names=NULL,header=FALSE)
plotTenReconstructions(mapAulStuR_05pc,"0.5pc")
plotTenReconstructions(mapAulStuR_1pc,"1pc")
plotTenReconstructions(mapAulStuR_2pc,"2pc")
plotTenReconstructions(mapAulStuR_5pc,"5pc")
cor_c1000_05pc_Aul <- getCorrelationsWithRefMapAtDifferentScales(mapAulStuR_05pc)
cor_c1000_1pc_Aul <- getCorrelationsWithRefMapAtDifferentScales(mapAulStuR_1pc)
cor_c1000_2pc_Aul <- getCorrelationsWithRefMapAtDifferentScales(mapAulStuR_2pc)
cor_c1000_5pc_Aul <- getCorrelationsWithRefMapAtDifferentScales(mapAulStuR_5pc)
pdf("output_plots/simulations_1000x_errorRateDependence.pdf",9,9)
plotCorrelations(list(cor_c1000_05pc_Aul, cor_c1000_1pc_Aul, cor_c1000_2pc_Aul, cor_c1000_5pc_Aul),c("0.5pc","1pc","2pc","5pc"),0.45)
dev.off()

# Looking at mean rate estimation
boxplot(apply(mapAulStuR_3000x[,4:13],2,mean,na.rm=T),ylim=c(min(mean(mapAulStuR_3000x[,3]),apply(mapAulStuR_3000x[,4:13],2,mean,na.rm=T)), max(mean(mapAulStuR_3000x[,3]),apply(mapAulStuR_3000x[,4:13],2,mean,na.rm=T))))
abline(h=mean(mapAulStuR_3000x[,3]),col="red")


#############################################################################
#### 2) Analyses of reconstructed maps
#############################################################################

# Map for Supplementary Figure 4:
mapAulStu2 <- read.table("maps/reconstructed/recombMap241209_H5a2_chr2_FW_2000.txt",row.names=NULL,header=FALSE)
pdf("output_plots/recombMap241209_H5a2_chr2_FW_2000.pdf",width=10,height=6)
plot(mapAulStu2[,1], mapAulStu2[,3],type='l',xlab="chr 2 (LS420020) coordinates",ylab="r (per bp)", main="A. stuartgranti - original", ylim=c(0,3e-7))
points(fAstCal14ChrSizes[2,2],0,pch=16)
dev.off()

# Figure 3A - empirical error rates
empiricalErrorRates <- loadEmpiricalErrorRates()
pdf("output_plots/EmpiricalErrorRates.pdf",height=7,width=4)
boxplot(empiricalErrorRates,las=2,ylab="Estimated error rate (per chromosome)",ylim=c(0,0.05))
dev.off()

# Load (and plot) the reconstructed maps:
mapAulStu <- plotAndLoadCichlidMap(thisSample="H5a2", doPlot=FALSE); 
mapAstCal <- plotAndLoadCichlidMap(thisSample="H3C2", doPlot=FALSE)
mapAstNub <- plotAndLoadCichlidMap(thisSample="AN01", doPlot=FALSE)
mapNeoMul <- plotAndLoadCichlidMap(thisSample="T2-2", doPlot=FALSE)
mapGasAcu <- plotAndLoadGasAcuMap(); mapGasAcuTrio <- plotAndLoadGasAcuMap("_trio")

# Figure 3B - Get map lengths for Hi-Recomb and F2 linkage-based maps
mapLengths_Cm <- getMapLengths_Cm(mapAulStu, mapAstCal, mapAstNub, mapNeoMul, mapGasAcu, mapGasAcuTrio)
mapLengths_Cm_LG <- getMapLengths_Cm_linkageMaps(AlbertsonMap, KocherMap, roestiMap)
totalLengths <- c(sapply(mapLengths_Cm,sum)[c(1,3,4,2,5,6)], sapply(mapLengths_Cm_LG,sum))
pdf("output_plots/TotalMapLengths.pdf",height=6,width=4)
barplot(totalLengths[-6],las=2,ylab="Total map length (cM)",col=c(colsSpecies,"black","black","black"))
dev.off()

# Figure 3C - Get mean rates per chromosome and plot that for cichlids
meanRates <- getMeanRates(mapAulStu, mapAstCal, mapAstNub, mapNeoMul, mapGasAcu, mapGasAcuTrio)
pdf("output_plots/CichlidRatePerChr.pdf",height=8,width=8)
lms_cichlid <- plotChrSizeVsMeanRateCichlid(genomes, meanRates)
dev.off()

pdf("output_plots/SticklebackRatePerChr.pdf",height=4,width=4)
lmGasAcu <- lm(meanRates$GasAcuTrio~genomes$GasAcu[-19,2])
plot(genomes$GasAcu[-19,2], meanRates$GasAcuTrio, col=colsSpecies[5],pch=16, xlab="Chromosome size", ylab="Mean r")
abline(lmGasAcu,col=colsSpecies[5])
dev.off()

# Supplementary Figure 1 - effective depth along the chromosome
depthNeoMul16 <- read.table("maps/reconstructed/wDepth/recombMap_wDepth_T2-2_chr16.txt",header=T)
pdf("output_plots/recombMap_wDepth_T2-2_chr15.pdf",width=10,height=6)
plot(depthNeoMul16$pos, depthNeoMul16$EffectiveDepth,type='l',ylab="Effective coverage",xlab="Coordinates")
abline(h=0.2*mean(depthNeoMul16$EffectiveDepth))
dev.off()

## Figure 4 - Bootstrap
bootAulStu <- plotAndLoadCichlidBootstrap(thisSample="H5a2")
bootAstNub <- plotAndLoadCichlidBootstrap(thisSample="AN01")
bootAstCal <- plotAndLoadCichlidBootstrap(thisSample="H3C2")

# Figure 4A - An example chromosome with bootstrap
pdf("output_plots/AulStuChr4.pdf",9,4)
plotBootOneMap(bootAulStu[[4]],1,length(bootAulStu[[4]][,1]),bootYmax=2e-7)
dev.off()

# Figures 4B and 4C - bootstrap differences between AulStu and AstNub
sigDiffAulStuAstNub <- getAllSigDiffRegions(bootAulStu, bootAstNub)

plotBootRegion(bootAulStu[[4]], bootAstNub[[4]], 1, length(bootAulStu[[4]][,1]), bootYmax=2e-7)
pdf("output_plots/AulStuAstNub1.pdf",5,4)
plotBootRegion(bootAulStu[[4]], bootAstNub[[4]], 3000, 5000, bootYmax=2e-7)
sig <- sigDiffAulStuAstNub[[4]]
segments(x0=sig[,1]*2000,y0=0,x1=sig[,2]*2000,y1=0,col="red",lwd=2)
dev.off()
pdf("output_plots/AulStuAstNub2.pdf",5,4)
plotBootRegion(bootAulStu[[4]], bootAstNub[[4]], 6500, 8000, bootYmax=2e-7)
segments(x0=sig[,1]*2000,y0=0,x1=sig[,2]*2000,y1=0,col="red",lwd=2)
dev.off()

# Statistics regarding delta(r) regions:
sigDiffLengths <- sapply(sigDiffAulStuAstNub,function(x) {sum(x[,3])})
sigDiffLongest <- sapply(sigDiffAulStuAstNub,function(x) {max(x[,3])})
sigDiffOverallProportion <- (sum(sigDiffLengths)*2000)/sum(genomes$fAstCal[1:22,2])
sigDiffProportions <- (sigDiffLengths*2000)/genomes$fAstCal[,2]

# Statistics regarding "mean-normalised" delta(r) regions:
sigDiffAulStuAstNubMa <- getAllSigDiffRegions(bootAulStu, bootAstNub, meanAdjust=TRUE)
sigDiffLengthsMa <- sapply(sigDiffAulStuAstNubMa,function(x) {sum(x[,3])})
sigDiffLongestMa <- sapply(sigDiffAulStuAstNubMa,function(x) {max(x[,3])})
sigDiffOverallProportionMa <- (sum(sigDiffLengthsMa)*2000)/sum(genomes$fAstCal[1:22,2])
sigDiffProportionsMa <- (sigDiffLengthsMa*2000)/genomes$fAstCal[,2]


#############################################################################
#### 3) Comparisons with trio-based phasing and LD-based maps
#############################################################################

# ---------- Figure 5A - error rates ---------
empiricalErrorRates <- loadEmpiricalErrorRates()
pdf("output_plots/EmpiricalErrorRatesStickleback.pdf",height=7,width=4)
boxplot(empiricalErrorRates[5:6],las=2,ylab="Estimated error rate (per chromosome)",ylim=c(0,0.01))
dev.off() 

# ---------- Figure 5B correlations with LD-based maps ------------
LDmapSticklebackList <- get_GasAcu_WB_LDmap(LDmapsStickleback); LDmapSticklebackList <- LDmapSticklebackList[-19] # Exclude the X chromosome
corGasAcu_LD <- getCorrelationsOfTwoMapsAtDifferentScales(mapGasAcu, LDmapSticklebackList)
corGasAcuTrio_LD <- getCorrelationsOfTwoMapsAtDifferentScales(mapGasAcuTrio, LDmapSticklebackList)
pdf("output_plots/SticklebackCorrelations.pdf",height=9,width=9)
plotCorrelations(list(corGasAcu_LD, corGasAcuTrio_LD),c("Hi-C","trio"),limYbottom=0.15,limYtop=0.8)
dev.off()

# -------- Figure 5C number of phased SNPs hapcut2 vs. trio
nSNPsTrio <- c(54314, 46994, 34918, 71454, 30199, 37767, 62321, 40587, 42564, 38101, 33250, 42065, 42181, 33720, 34202, 36055, 37619, 31476, 13119, 53428, 36228)
nSNPsHapcut2 <- c(75386, 64152, 48499, 99880, 41911, 51582, 84728, 55508, 57539, 50697, 46602, 58058, 57594, 46679, 46718, 50998, 50274, 43042, 18731, 69920, 47707)
pdf("output_plots/SNPDensityStickleback.pdf",height=7,width=4)
boxplot(nSNPsHapcut2[-19]/genomes$GasAcu[-19,2], nSNPsTrio[-19]/genomes$GasAcu[-19,2],las=2,ylab="Density of phased SNPs (per chromosome)",ylim=c(0.0015,0.004))
dev.off()


#############################################################################
#### 4) Supplementary
#############################################################################
# 4.1) Proportion unmappable regions

amountUncallable <- c(5990107, 6821876, 27319741, 7081125, 5588721, 7135983, 12191986, 4551199, 7653867, 5828483, 5542726, 5435144, 4182366, 6043555, 6666845, 6428874, 6444046, 5674946, 5177836, 5327234, 7329078, 9325122)
amountUncallableNeoMul <- c(5264138, 6286556, 13214761, 5524453, 5619581, 6957406, 23960890, 4718727, 5786479, 5857871, 6629990, 5600420, 4951214, 8393375, 6762722, 5862192, 5868791, 7589555, 4033579, 5569304, 6904865, 8791352)
amountUncallableGasAcu <- c(1474197, 936097, 1164422, 2282413, 1086870, 1014780, 1901008, 1205936, 1725117, 1869977, 1293103, 1518652, 1273726, 765349, 801949, 1166709, 705242, 985952, 2588096, 1151318, 2567993, 5823286)
pdf("output_plots/proportionUnmappable.pdf",height=9,width=6)
par(mfrow=c(3,1),mar=c(4.1,4.1,3.1,2.1))
plot(amountUncallable/genomes$fAstCal[,2],xlab="Chromosome",ylab=c("Proportion unmappable"),main="fAstCal1.5",xaxt='n',ylim=c(0,0.42),pch=16)
axis(1,1:22,1:22)
plot(amountUncallableNeoMul/genomes$fNeoMul[,2],xlab="Chromosome",ylab=c("Proportion unmappable"),main="fNeoMul1.2",xaxt='n',ylim=c(0,0.42),pch=16)
axis(1,1:22,1:22)
plot(amountUncallableGasAcu[-22]/genomes$GasAcu[,2],xlab="Chromosome",ylab=c("Proportion unmappable"),main="stickleback v5",xaxt='n',ylim=c(0,0.42),pch=16)
axis(1,1:21,1:21)
dev.off()

