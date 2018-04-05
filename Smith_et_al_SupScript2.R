###################################################################################
## Additional File 2
##
## This script allows the reproduction of all the figures that were created using R
## The data can be obtained from:
## http://www.compbio.group.cam.ac.uk/Resources/baloc/
## This script expect the data to be in two folders named Expression and CNV
## It this isn't the case then the script will need to be modified accordingly
###################################################################################

## All libraries are available from Bioconductor
library(beadarray)
library(BeadDataPackR)  
library(EBImage)
source("allFunctions.R")

pdf(file = "allScriptsOutput.pdf")

############################################################
## Figure 1 - BeadArray Structure
###########################################################

############################################################
## Figure 2 - Association between fraction part & intensity
###########################################################

setwd("Expression")
t1 <- read.table("4343238066_A_1.txt", header = T, sep = "\t")
t1[,3] <- floor(t1[,3]) + 0.5
t1[,4] <- floor(t1[,4]) + 0.5
write.table(t1, file = "4343238066_A_1.txt.round", sep = "\t", col.name = TRUE, row.names = FALSE, quote = FALSE)

BLData1 <- readIllumina(arrayNames = "4343238066_A_1", textType = ".txt")
BLData2 <- readIllumina(arrayNames = "4343238066_A_1", textType = ".txt.round")

tmpX <- getArrayData(BLData1, what = "GrnX") - floor(getArrayData(BLData1, what = "GrnX"))
tmpY <- getArrayData(BLData1, what = "GrnY") - floor(getArrayData(BLData1, what = "GrnY"))
inten1 <- getArrayData(BLData1)
inten1[which(is.na(inten1))] = 0
inten2 <- getArrayData(BLData2)
inten2[which(is.na(inten2))] = 0

resids1 <- getArrayData(BLData1, what = "residG")
resids1[which(is.na(resids1))] = 0
resids2 <- getArrayData(BLData2, what = "residG")
resids2[which(is.na(resids2))] = 0

par(mfrow = c(2,2))
tmp <- imagePlotAvg(inten1, tmpX, tmpY)
tmp <- imagePlotAvg(inten2, tmpX, tmpY)
tmp <- imagePlotAvg(resids1, tmpX, tmpY, divisions = 100)
tmp <- imagePlotAvg(resids2, tmpX, tmpY, divisions = 100)
setwd("../")


############################################################
## Figure 3 - Localised mapping errors
###########################################################

setwd("CNV")

temp<-read.delim("4127130020_A_4.txt",header=T,as.is=T)
tifgreen<-readTIFF("4127130020_A_4_Grn.tif")
tifred<-readTIFF("4127130020_A_4_Red.tif")

tempmatch<-sqrt((temp$GrnX-1370.615)^2+(temp$GrnY-2474.508)^2)
subtemp<-temp[which(tempmatch<13),]
tempord<-c(5,1,16,14,3,13,19,4,18,6,9,15,2,8,7,12,17,10,11)

subtemp<-subtemp[tempord,]

exptemp<-temp[!is.na(match(temp$Code,subtemp$Code)),]

exptemp<-exptemp[is.na(match(rownames(exptemp),rownames(subtemp))),]
lettertemp<-LETTERS[match(exptemp$Code,subtemp$Code)]

meds<-unlist(lapply(split(log2(exptemp$Red),lettertemp),median))

par(mar=c(3,4,1,1))

layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),heights=c(2,1))

plotTIFF(tifgreen,1330,1410,2435,2515,high="green",axes=F)
axis(1,at=seq(1330,1410,20)+0.5,labels=seq(1330,1410,20))
axis(2,at=seq(2435,2515,20)+0.5,labels=seq(2435,2515,20))
text(subtemp$GrnX+0.5,subtemp$GrnY+0.5,LETTERS[1:19],col="white")

plotTIFF(tifred,1338,1418,2435,2515,high="red",axes=F)
meanY<-median(temp$RedY-temp$GrnY)
meanX<-median(temp$RedX-temp$GrnX)
tempdist<-sqrt((temp$RedX-temp$GrnX-meanX)^2+(temp$RedY-temp$GrnY-meanY)^2)

axis(1,at=seq(1330+8,1410+8,20)+0.5,labels=seq(1330+8,1410+8,20))
axis(2,at=seq(2435,2515,20)+0.5,labels=seq(2435,2515,20))
text(subtemp$RedX+0.5,subtemp$RedY+0.5,LETTERS[1:19],col="white")

rectcol<-rep("grey50",19)
rectcol[c(1,4,5,6,7,14,15,16,18)]<-"#5D5DB5"

plot((log2(subtemp$Red)-meds),type="n",axes=F,xlab="",ylab="")
rect(1:19-0.4,rep(0,19),1:19+0.4,(log2(subtemp$Red)-meds),col=rectcol)
box()
axis(2,at=c(-4,-2,0,2))
axis(1,at=1:19,labels=LETTERS[1:19])
mtext("departure from\nbead-type median",2,2,cex=0.8)
setwd("../")


############################################################
## Figure 4
###########################################################

setwd("CNV")
t1 <- read.table("4127130020_B_8.txt", sep = "\t", header = T, as.is = T)
locsRed <- readLocs("4127130020_B_8_Red.locs")

combRed <- BeadDataPackR:::combineFiles(t1[,c(1,5:7)], locsRed);

combRed2 <- combRed[order(combRed[,5]),]
seg1 <- combRed2[1:147519,]
seg2 <- combRed2[147520 : (2*147519),]
seg3 <- combRed2[(2 * 147519 + 1) : (3*147519),]
seg4 <- combRed2[(3 * 147519 + 1) : (4*147519),]

good <- rbind(seg1, seg4)
bad <- rbind(seg2, seg3)
bad <- bad[-which(!bad[,1] %in% good[,1]),]
good <- good[-which(!good[,1] %in% bad[,1]),]

m1 <- unlist(lapply(split(log2(good[,2]), good[,1]), mean, na.rm = T))
m2 <- unlist(lapply(split(log2(bad[,2]), bad[,1]), mean, na.rm = T))
remIndex <- which( (m1 == -Inf) | (m2 == -Inf) )
m1 <- m1[ -remIndex ]
m2 <- m2[ -remIndex ]

par(mfrow = c(1,2))
plot(m1, m2, pch = ".", ylab = "Mean Intensity for probes found on segment 2 & 3", xlab = "Mean Intensity for probes found on segments 1 & 4", col = "grey40", main = "Original Coordinates", ylim = c(5,15), xlim = c(6,15))
abline(0, 1)
text(x = 7, y = 5, labels = paste("Correlation", round(cor(m1, m2, method = "spearman"), 4) ))

fixed <- rbind(gridShift(seg2, 6, 0, 333, 443), gridShift(seg3, 6, 0, 333, 443))
good <- rbind(seg1, seg4)
good <- good[-which(!good[,1] %in% fixed[,1]),]
fixed <- fixed[-which(!fixed[,1] %in% good[,1]),]
m3 <- unlist(lapply(split(log2(fixed[,2]), fixed[,1]), mean, na.rm = T))
m1 <- unlist(lapply(split(log2(good[,2]), good[,1]), mean, na.rm = T))
remIndex <- which( (m1 == -Inf) | (m3 == -Inf) )
m1 <- m1[ -remIndex ]
m3 <- m3[ -remIndex ]


plot(m1, m3, pch = ".", ylab = "Mean Intensity for probes found on segment 2 & 3", xlab = "Mean Intensity for probes found on segments 1 & 4", col = "grey40", main = "Remapped Coordinates", ylim = c(5,15), xlim = c(6,15))
abline(0, 1)
text(x = 7, y = 5, labels = paste("Correlation", round(cor(m1, m3, method = "spearman"), 4) ))
setwd("../")

############################################################
## Figure 5
###########################################################
setwd("Expression")
par(mar=c(1,1,3,1))

layout(matrix(c(1,2,3,4,9,5,6,7,8,9),byrow=T,nrow=2),widths=c(3,3,3,3,1))

mycolramp<-c(colorRampPalette(c("yellow","white"),bias=2)(50),colorRampPalette(c("white","blue"),bias=2)(50))

zlim1 <- -1.5
zlim2 <- 1.5

textfilevec<-c("4343238066_D_1.txt", "4343238080_A_1.txt", "4343238080_A_2.txt", "4343238080_B_2.txt", "4343238080_C_1.txt", "4343238080_C_1.txt", "4343238080_E_1.txt", "4343238080_F_2.txt")

mainvec<-c("4343238066_D_1 S 4","4343238080_A_1 S 9","4343238080_A_2 S 5","4343238080_B_2 S 2","4343238080_C_1 S 1","4343238080_C_1 S 1","4343238080_E_1 S 7","4343238080_F_2 S 9")

ylowvec<-c(6125,16205,8100,2150,0,10100,12075,16025)
yhighvec<-c(8100,18000,10075,4150,2150,12075,14050,18000)

for(ind in 1:8){
    textfile<-read.delim(textfilevec[ind])
    tmp<-sapply(split(textfile$Grn,textfile$Code),median)
    textfile<-cbind(textfile,tmp[match(textfile$Code,names(tmp))])
    colnames(textfile)[5]<-"Med"
    textfile<-cbind(textfile,(log2(textfile$Grn)-log2(textfile$Med)))
    colnames(textfile)[6]<-"Res"
    textfile<-textfile[which((textfile$GrnY>ylowvec[ind])&(textfile$GrnY<yhighvec[ind])),]
    textfile$GrnX<-ceiling((textfile$GrnX-min(textfile$GrnX))/20+0.000000001)
    textfile$GrnY<-ceiling((textfile$GrnY-min(textfile$GrnY))/20+0.000000001)
    textfile$Res[textfile$Res==-Inf]<-min(textfile$Res[textfile$Res != -Inf])
    temp<-split(x=textfile$Res,f=list(textfile$GrnX,textfile$GrnY))
    temp2<-(matrix(sapply(temp,mean,na.rm=T),nrow=103))
    temp2[which(temp2>zlim2)]<-zlim2
    temp2[which(temp2<zlim1)]<-zlim1
    image(t(temp2),col=mycolramp,axes=F,xlab="",main=mainvec[ind],zlim=c(zlim1,zlim2))
    newtemp<-temp2
    newtemp[!is.na(temp2)]<-NA
    newtemp[is.na(temp2)]<-1
    image(t(newtemp),col="orange",axes=F,xlab="",main=mainvec[ind],zlim=c(zlim1,zlim2),add=T)
    }
par(mar=c(1,0,3,4))
image(t(matrix(seq(zlim1,zlim2,0.1),ncol=1)),col=mycolramp,axes=F,xlab="",main="",zlim=c(zlim1,zlim2))
axis(4,labels=paste(c("<= ","","","",">= "), seq(zlim1, zlim2, length.out = 5), sep=""), at=seq(0,1,0.25), las=1)
setwd("../")


###########################################
## Figure 6
## Influence of bright beads on their neighbours
###########################################

setwd("Expression")
BLData <- readIllumina(arrayNames = "4343238080_B_2", useImages = FALSE);
tiff <- readTIFF("4343238080_B_2_Grn.tif")
t1 <- matrix(unlist(scan(file = "4343238080_B_2.txt", what = list(0, 0, double(), double()), skip = 1, sep = "\t", quiet = TRUE)), ncol = 4)


neighbours <- generateNeighbours(BLData, array = 1)
index <- 818270

par(mfrow = c(1,2), mar = c(5,4,4,2))
for(i in index) {
  x <- floor(t1[i,3])
  y <- floor(t1[i,4])

  plotTIFF(tiff, x-13, x+10, y-10, y+11, xaxt = "n", yaxt = "n", accountForZero = TRUE, high = rgb(29, 0, 150, maxColorValue = 255), low = "white")

  currentNeighbours <- neighbours[i,neighbours[i,] != 0]
  k = 1

  values <- matrix(ncol = 2, nrow = length(currentNeighbours))
  
  for(j in currentNeighbours) {
    foregroundWeights(tiff, t1[j,3], t1[j,4], lwd = 1, border = "grey50")
    points(floor(t1[j,3]) - 2, t1[j,4], pch = paste(k, sep = ""), col = "black", cex = 1.5)

    ProbeID <- getArrayData(BLData, what = "ProbeID", array = 1)[j]
    values[k, 1] <- getArrayData(BLData, array = 1)[j]
    ProbeIdx <- which(getArrayData(BLData, what = "ProbeID", array = 1) == ProbeID)
    values[k, 2] <- median(getArrayData(BLData, array = 1)[ProbeIdx])
    k = k+1
  }
  #draw the bead nearby with a more even weights matrix
  foregroundWeights(tiff, t1[368926,3], t1[368926,4], lwd = 1, border = "grey50")
  points(floor(t1[368926,3]) - 2, t1[368926,4], pch = "7", col = "black", cex = 1.5)
  #add the bead centres
  points(t1[,3], t1[,4], col = "black", pch = 4, cex = 1)
  points(t1[index,3], t1[index,4], col = "white", pch = 4)
  barplot(values[,1] - values[,2], ylab = "departure from bead-type median", xlab = "bead label", names.arg = c(1:length(currentNeighbours)), col = "#5D5DB5")
}
setwd("../")

###############################################################
## Figure 7
## Dead pixels affecting background
##################################################

setwd("Expression")
BLData1 <- readIllumina(arrayName = "4343238066_A_2")
BLData2 <- readIllumina(arrayName = "4343238066_A_2", backgroundCalc = "median")
BSData1 <- createBeadSummaryData(BLData1, log = TRUE)
BSData2 <- createBeadSummaryData(BLData2, log = TRUE)


t1 <- matrix(unlist(scan(file = "4343238066_A_2.txt", what = list(0, 0, double(), double()), skip = 1, sep = "\t", quiet = TRUE)), ncol = 4)
tiff <- readTIFF("4343238066_A_2_Grn.tif")

par(mfrow = c(1,2))
plotTIFF(tiff, 964-8, 964+8, 6081-8, 6081+8, xaxt = "n", yaxt = "n", high = "darkblue", mid = "white", low = "white", values = FALSE, textCol = "yellow")

idx <- proximalBeads(c(6081, 964), t1)

xs <- (getArrayData(BLData1, what = "GrnX") + getArrayData(BLData1, what = "xOffset"))[idx];
ys <- (getArrayData(BLData1, what = "GrnY") + getArrayData(BLData1, what = "yOffset"))[idx];

points(xs, ys, pch = 4, col = "white", cex = 1.5)
points(xs+0.7, ys, pch = paste(1:7), col = "grey20", cex = 1.5)

tmp <- NULL
for(i in 1:length(idx)) {
    tmp <- c(tmp, getArrayData(BLData1, log = TRUE)[idx[i]] - (exprs(BSData1[paste(ids[i])])));
    tmp <- c(tmp, getArrayData(BLData2, log = TRUE)[idx[i]] - (exprs(BSData2[paste(ids[i])])));
}
barplot(tmp, col = c("#5D5DB5", "purple"), ylab = "departure from summarized log intensity for bead-type", xlab = "bead label", xaxt = "n")
axis(1, at = seq(1, 16.5, 2.5), 1:length(idx), tick = FALSE)
legend(0, -0.5, legend = c("Mean of five lowest", "Median of five lowest"), fill = c("#5D5DB5", "purple"))

setwd("../")

##################################################
## Figure 8 - Impact on biological interpretation
##################################################

setwd("Expression")

## negative controls from the beadarray library
data(ExpressionControlData) 
controls <- ExpressionControlData[["Humanv2"]]
controls <- controls[which(controls[,2] == "negative"),]

arrays <- c("4343238080_B_1", "4343238080_B_2", "4343238080_D_1", "4343238080_D_2", "4343238080_F_1", "4343238080_F_2");

BLData.standard <- readIllumina(arrayNames = arrays);
BLData.removed <- readIllumina(arrayName = arrays, backgroundCalc = "median");

NN <- neighboursMatrixForAll();

weights <- list()
for(i in 1:length(arrays)) {
    weights[[1]] <- identifyAffectedBeads(arrays[i], NN);
    BLData.removed <- setWeights(BLData.removed, wts = weights[1], array = i);
}
## create the summarized data
BSData.standard <- createBeadSummaryData(BLData.standard, log = TRUE);
BSData.removed <- createBeadSummaryData(BLData.removed, log = TRUE);

ids <- c("5900598")

## get the bead-level data from section D2
d2.standard <- BLData.standard[[4]][which(BLData.standard[[4]][,1] == ids[1]),]
d2.removed <- BLData.removed[[4]][which(BLData.removed[[4]][,1] == ids[1]),]

par(mfrow = c(1,3), mar = c(3,4,3,2))
## plot the summarized log intensity for the 6 sections using both analysis methods
plot(c(exprs(BSData.standard)[ids[1],]), pch = 20, col = "blue", ylim = c(5, 7), ylab = "log intensity", xaxt = "n")
axis(1, at = 1:6, label = c("B_1", "B_2", "D_1", "D_2", "F_1", "F_2"))
points(c(exprs(BSData.removed)[ids[1],]), col = "purple", pch = 20)
legend(0.8, 5.18, legend = c("Standard Analysis", "Affected Beads Removed"), fill = c("blue", "purple"))

## examine the bead-level data and summarization cutoffs for the standard analysis
h.nw <- hist(exprs(BSData.standard)[paste(controls[,1]),4], plot = FALSE)
plot(0,0, ylim = c(3.5, max(log2(d2.standard[,2]))), xlim = c(0, 10 + nrow(d2.standard)), xaxs = "i", xaxt = "n", ylab = "log intensity", xlab = "")
points(x = 10 + seq(1,nrow(d2.standard)), y = log2(d2.standard[,2]), pch = 20, col = "blue",  ylab = "Log2 Intensity")
abline(h = exprs(BSData.standard)[ids[1],4],col = "black", lwd = 2 )
abline(h = median(log2(d2.standard[,2])) + 3*mad(log2(d2.standard[,2])), lty = 3 )
abline(h = median(log2(d2.standard[,2])) - 3*mad(log2(d2.standard[,2])), lty = 3 )
rect(rep(0, 14), h.nw$breaks[1:14], 10 * h.nw$counts / max(h.nw$counts), h.nw$breaks[2:15], col = "grey50")

## examine the bead-level data and summarization cutoffs for the two-step analysis
h.w <- hist(exprs(BSData.removed)[paste(controls[,1]),4], plot = FALSE)
plot(0,0, ylim = c(3.5, max(log2(d2.standard[,2]))), xlim = c(0, 10 + nrow(d2.removed)), xaxs = "i", xaxt = "n", ylab = "log intensity", xlab = "")
points(x = 1:nrow(d2.removed), y = log2(d2.removed[,2]), pch = 20, col = "blue",  ylab = "Log2 Intensity")
points(which(d2.removed[,6] == 0), log2(d2.standard[which(d2.removed[,6] == 0),2]), col = "red", pch = 4, cex = 2)
abline(h = exprs(BSData.removed)[ids[1],4],col = "black" , lwd = 2)
abline(h = median(log2(d2.removed[-which(d2.removed[,6] == 0),2])) + 3*mad(log2(d2.removed[-which(d2.removed[,6] == 0),2])), lty = 3)
abline(h = median(log2(d2.removed[-which(d2.removed[,6] == 0),2])) - 3*mad(log2(d2.removed[-which(d2.removed[,6] == 0),2])), lty = 3)
rect(rep(nrow(d2.removed) + 10, 14) - 10 * h.w$counts / max(h.w$counts), h.w$breaks[1:14], rep(nrow(d2.removed) + 10, 14), h.w$breaks[2:15], col = "grey50")

setwd("../");

## ***************************************************************************** ##
## SUPPLEMENTARY FIGURES
## ============================================================================= ##

####################################################
## Additional File 1 - Script containing functions 
####################################################

####################################################
## Additional File 2 - This script
####################################################

####################################################
## Additional File 3 - Pixels used in packground for multiple beads
####################################################

par(mfrow = c(1,1))
setwd("Expression")
mytxt <- matrix(unlist(scan(file = "4343238066_D_1.txt", what = list(0, 0, double(), double()), skip = 1, sep = "\t", quiet = TRUE)), ncol = 4)
mytif<-readTIFF("4343238066_D_1_Grn.tif")

plotTIFF(mytif,1466,1482,293,309, high =  rgb(29, 0, 150, maxColorValue = 255), low = "white")
points(mytxt[,3],mytxt[,4],pch=20,col="grey20", cex = 2)
lines(c(1473.5,1474.5),c(300.5,300.5),col="grey20")
lines(c(1473.5,1474.5),c(301.5,301.5),col="grey20")
lines(c(1473.5,1473.5),c(300.5,301.5),col="grey20")
lines(c(1474.5,1474.5),c(300.5,301.5),col="grey20")
lines(c(1473.5,1474.5),c(300.5,301.5),col="grey20")
lines(c(1473.5,1474.5),c(301.5,300.5),col="grey20")

lines(c(1473.5,1474.5),c(301.5,301.5),col="grey20")
lines(c(1473.5,1474.5),c(302.5,302.5),col="grey20")
lines(c(1473.5,1473.5),c(301.5,302.5),col="grey20")
lines(c(1474.5,1474.5),c(301.5,302.5),col="grey20")
lines(c(1473.5,1474.5),c(301.5,302.5),col="grey20")
lines(c(1473.5,1474.5),c(302.5,301.5),col="grey20")
setwd("../")

####################################################
## Additional File 5 - Flowchart of method
####################################################


####################################################
## Additional File 6 - Bead simulation
####################################################

results1 <- matrix(NA,ncol = 25, nrow = 25)
offsets <- seq(0.5, 1.49, 1/25)
for(j in 1:25) {
  cat( j, "\n")
  for(i in 1:25) {
      tmp1 <- createBead(xOffset = offsets[i], yOffset = offsets[j], gridSize = 13)
      tmp2 <- pixelateBead(tmp1, gridSize = 13)
      results1[j, i] <- matrixForeground(tmp2, xOffset = offsets[i], yOffset = offsets[j])   
  }
}    
imagePlot(log2(results1), high = "blue", low = "white")
results2 <- matrix(NA,ncol = 25, nrow = 25)
offsets <- seq(0.5, 1.49, 1/25)
for(j in 1:25) {
  cat( j, "\n")
  for(i in 1:25) {
      tmp1 <- createBead(xOffset = offsets[i], yOffset = offsets[j], gridSize = 11)
      tmp2 <- pixelateBead(tmp1, gridSize = 11, maxInten = 2^12)
      results2[j, i] <- matrixForeground(tmp2, xOffset = offsets[i], yOffset = offsets[j])   
  }
}    
imagePlot(log2(results2), high = "blue", low = "white")
results3 <- matrix(NA,ncol = 25, nrow = 25)
offsets <- seq(0.5, 1.49, 1/25)
for(j in 1:25) {
  cat( j, "\n")
  for(i in 1:25) {
      tmp1 <- createBead(xOffset = offsets[i], yOffset = offsets[j], gridSize = 9)
      tmp2 <- pixelateBead(tmp1, gridSize = 9, maxInten = 2^10)
      results3[j, i] <- matrixForeground(tmp2, xOffset = offsets[i], yOffset = offsets[j])   
  }
}    
imagePlot(log2(results3), high = "blue", low = "white")
results4 <- matrix(NA,ncol = 25, nrow = 25)
offsets <- seq(0.5, 1.49, 1/25)
for(j in 1:25) {
  cat( j, "\n")
  for(i in 1:25) {
      tmp1 <- createBead(xOffset = offsets[i], yOffset = offsets[j], gridSize = 7)
      tmp2 <- pixelateBead(tmp1, gridSize = 7, maxInten = 2^8)
      results4[j, i] <- matrixForeground(tmp2, xOffset = offsets[i], yOffset = offsets[j])   
  }
}    
imagePlot(log2(results4), high = "blue", low = "white")

####################################################
## Additional File 7 - Association between residual log intensity and departure from grid
####################################################

setwd("Expression")

BLData <- readIllumina(arrayNames = "4343238080_D_1", singleChannel = TRUE, useImages = FALSE)
t1 <- matrix(unlist(scan(file = "4343238080_D_1.txt", what = list(0, 0, double(), double()), skip = 1, sep = "\t", quiet = TRUE)), ncol = 4)
locs <- readLocs("4343238080_D_1_Grn.locs")
d <- deviationFromGrid(t1, locs);
tmp <- split((getArrayData(BLData, what = "residG")), d[[3]] %/% 0.25)

boxplot(tmp, xaxt = "n", col = "grey70", xlab = "Distance from predicted grid position in pixels", ylab = "Residual Intensity")
axis(1, at = 1:length(tmp), labels = seq(0, max(d[[3]]), 0.25))
abline(h = 0, lty = 3)

setwd("../")


####################################################
## Additional File 8 - Coordinate Shift
####################################################

setwd("CNV")
txtfile<-read.delim("4127130188_B_9.txt",header=T,as.is=T)
GTIFF<-readTIFF("4127130188_B_9_Grn.tif")
RTIFF<-readTIFF("4127130188_B_9_Red.tif")

## Divide the chip into sections
#Y should break up nicely into 4 segments

GY<-cut(txtfile$GrnY,4,0:3)
GX<-cut(txtfile$GrnX,10,0:9)

Xshift<-matrix(0,nrow=4,ncol=10)
Yshift<-matrix(0,nrow=4,ncol=10)

for(i in 1:10){
for(j in 1:4){
Xshift[j,i]<-median(txtfile$GrnX[(GX==(i-1))&(GY==(j-1))]-txtfile$RedX[(GX==(i-1))&(GY==(j-1))])
Yshift[j,i]<-median(txtfile$GrnY[(GX==(i-1))&(GY==(j-1))]-txtfile$RedY[(GX==(i-1))&(GY==(j-1))])
}}


par(mfrow=c(1,2))
par(mar=c(5.1,4.1,1.1,2.1))
plot(Xshift[1,],col="#7E2271",pch=16,ylab="shift in X co-ordinates from green to red channel",xlab="Decile of segment",ylim=c(-9,-6.5),main="",axes=F)
axis(1,at=1:10)
axis(2)
box()
points(Xshift[2,],col="#00626D",pch=16)
points(Xshift[3,],col="#CD5806",pch=16)
points(Xshift[4,],col="#EF6ABF",pch=16)
legend("topleft",legend=c("segment1","segment2","segment3","segment4"),fill=c("#7E2271","#00626D","#CD5806","#EF6ABF"))

plot(Yshift[1,],col="#7E2271",pch=16,ylab="shift in Y co-ordinates from green to red channel",xlab="Decile of segment",ylim=c(-1.2,-0.9),main="",axes=F)
axis(1,at=1:10)
axis(2)
box()
points(Yshift[2,],col="#00626D",pch=16)
points(Yshift[3,],col="#CD5806",pch=16)
points(Yshift[4,],col="#EF6ABF",pch=16)
legend("topleft",legend=c("segment1","segment2","segment3","segment4"),fill=c("#7E2271","#00626D","#CD5806","#EF6ABF"))
setwd("../")

####################################################
## Additional File 9 - Mis-registration
####################################################

setwd("CNV")

locsGrn <- readLocs("4127130020_B_8_Grn.locs")
locsRed <- readLocs("4127130020_B_8_Red.locs")
tiffRed <- t(readTIFF("4127130020_B_8_Red.tif"))

first <- seq(1, nrow(locsGrn), nrow(locsGrn)/4)
last <- seq(nrow(locsGrn)/4, nrow(locsGrn), nrow(locsGrn)/4)


par(mar = c(2,2,2,2), mfrow = c(1,1))
#plotTIFF(tiffRed, 0, ncol(tiffRed)-1, 0, 1000, low = "white", high = "darkblue")
plotTIFF(tiffRed, low = "white", high = "darkblue", xaxt = "n", yaxt = "n")
tmp <- locsRed[first[1]:last[1], ]
true <- c(min(tmp[,1]), max(tmp[,1]))
for(i in 1:4) {
  tmp <- locsRed[first[i]:last[i], ]
  rect(min(tmp[,2]), true[1], max(tmp[,2]), true[2], lwd = 6)
  rect(min(tmp[,2]), min(tmp[,1]), max(tmp[,2]), max(tmp[,1]), col = rgb(1,0,0,0.5), lwd = 0, border = NULL)
  #rect(min(tmp[,2]), min(tmp[,1]), max(tmp[,2]), max(tmp[,1]), col = rgb(1,0,0), lwd = 3, border = NULL)
}
rect(6000, 100, 6400, 500, lwd = 4, lty = 3)


t1 <- matrix(unlist(scan(file = "4127130020_B_8.txt", what = list(0, 0, double(), double(), 0, double(), double()), skip = 1, sep = "\t", quiet = TRUE)), ncol = 7)
locs <- readLocs("4127130020_B_8_Red.locs")
comb <- BeadDataPackR:::combineFiles(t1[,c(1,5:7)], locs)
comb <- comb[order(comb[,5]),]
nSegs <- 4
nRows <- 333
beadsPerSeg <- nrow(locs) / nSegs
nCols <- beadsPerSeg / nRows
store <- list()
for(i in 1:4) {
    cat(i,"\n");
    vars <- NULL
    seg <- comb[(((i-1) * beadsPerSeg) + 1):(i * beadsPerSeg),];
    for(j in -8:8) {
        for(k in 0) {
            meanVar <- testGridShift(seg, j, k, nRows = nRows, nCols = nCols)
            vars <- c(vars, meanVar)    
        }
    }
    store[[i]] <- vars
}

par(mfrow = c(1,2), mar = c(4.3,2,2,2))
plotTIFF(tiffRed, 6000, 6400, 100, 500, low = "white", high = "darkblue", xaxt = "n", yaxt = "n")
tmp <- locsRed[first[1]:last[1], ]
true <- c(min(tmp[,1]), max(tmp[,1]))
for(i in 1:4) {
  tmp <- locsRed[first[i]:last[i], ]
  rect(min(tmp[,2]), true[1], max(tmp[,2]), true[2], lwd = 2)
  rect(min(tmp[,2]), min(tmp[,1]), max(tmp[,2]), max(tmp[,1]), col = rgb(1,0,0,0.5), lwd = 0, border = NULL)
  #rect(min(tmp[,2]), min(tmp[,1]), max(tmp[,2]), max(tmp[,1]), col = rgb(1,0,0), lwd = 3, border = NULL)
}

plot(-8:8, store[[1]], pch = 19, xlab = "Shift Size", ylab = "Mean Within Bead-type Variance", ylim = c(0,4.5), main = "")
lines(seq(-8, 8, 2), store[[1]][-which(is.na(store[[1]]))], lwd = 2)
cols <- c("magenta", "blue", "purple")
for(i in 2:4) {
    points(seq(-8, 8, 2), store[[i]][-which(is.na(store[[i]]))], col = cols[i-1], pch = 19)
    lines(seq(-8, 8, 2), store[[i]][-which(is.na(store[[i]]))], col = cols[i-1], lwd = 2)
}
legend(-8.4, 0.75, legend = c("Segment 1", "Segment 2", "Segment 3", "Segment 4"), fill = c("black", "magenta", "blue", "purple"))

setwd("../")

####################################################
## Additional File 10 - Bias across whole arrays
####################################################

setwd("Expression")

## Beads with very bright neighbours ############################

arrayNames = dir(pattern = ".txt$")[1:24]
tmp1 <- list()
NN <- neighboursMatrixForAll();
for(arrayID in 1:24) {

    arrayName = strsplit(arrayNames[arrayID], ".txt")[[1]][1];
    message(paste("Analysing", arrayName));
    
    LOCS<-readLocs(paste(arrayName, "_Grn.locs", sep = ""))
    text <- matrix(unlist(scan(file = paste(arrayName, ".txt", sep = ""), what = list(0, 0, double(), double()), skip = 1, sep = "\t", quiet = TRUE)), ncol = 4)

    combed <- BeadDataPackR:::combineFiles(text, LOCS)
    combed <- combed[order(combed[,5]),]
    
    idx <- excludeBrightBeadNeighbours(combed, paste(arrayName, "_Grn.tif", sep = ""), neighbours = NN);

    ids <- combed[idx,1]
    idx <- idx[-which((ids == 0) | (idx == 0) )]
    ids <- combed[idx,1]
    
    allIDs <- combed[-idx,1]
    intens <- log2(combed[-idx,2])
    grouped <- split(intens, allIDs)
    summary <- lapply(grouped, median, na.rm = TRUE)
    tmp1[[arrayID]] <- vector(length = length(idx))
    arrayData <- log2(combed[,2])
    for(i in seq(along = idx)) {
        tmp1[[arrayID]][i] <- arrayData[idx[i]] - summary[[paste(ids[i])]];
    }
}

## Beads whose background contains very low pixel(s) ############################

arrayNames = dir(pattern = ".txt$")[1:24]
tmp2 <- list()
for(arrayID in 1:24) {

    arrayName = strsplit(arrayNames[arrayID], ".txt")[[1]][1];
    message(paste("Analysing", arrayName));
    
    LOCS <- readLocs(paste(arrayName, "_Grn.locs", sep = ""))
    text <- matrix(unlist(scan(file = paste(arrayName, ".txt", sep = ""), what = list(0, 0, double(), double()), skip = 1, sep = "\t", quiet = TRUE)), ncol = 4)
    tiff <- readTIFF(paste(arrayName, "_Grn.tif", sep = ""))

    lowPixels <- which(tiff < 256)
    if(length(lowPixels)) {

        dimPixels <- t(sapply(lowPixels, matrixCoords, nrow(tiff))) - 1
        
        idx <- unlist(apply(dimPixels, 1, proximalBeads, text));
        ids <- text[idx,1]
        
        allIDs <- text[-idx,1]
        intens <- log2(text[-idx,2])
        grouped <- split(intens, allIDs)
        summary <- lapply(grouped, median, na.rm = TRUE)
        tmp2[[paste(arrayName)]] <- vector(length = length(idx))
        arrayData <- log2(text[,2])
        for(i in seq(along = idx)) {
            tmp2[[paste(arrayName)]][i] <- arrayData[idx[i]] - summary[[paste(ids[i])]];
        }
    }
    else {
        tmp2[[paste(arrayName)]] <- NA
    }
}

## Examine deviation from bead type median for beads within 5 ####
## invasions of a clump of nondecoded beads ######################


arrayNames = dir(pattern = ".txt$")[1:24]
tmp3 <- list()
NN <- neighboursMatrixForAll();
for(arrayID in 1:24) {

    arrayName = strsplit(arrayNames[arrayID], ".txt")[[1]][1];
    message(paste("Analysing", arrayName));
    
    LOCS<-readLocs(paste(arrayName, "_Grn.locs", sep = ""))
    text <- matrix(unlist(scan(file = paste(arrayName, ".txt", sep = ""), what = list(0, 0, double(), double()), skip = 1, sep = "\t", quiet = TRUE)), ncol = 4)

    combed <- BeadDataPackR:::combineFiles(text, LOCS)

    idx <- nearNonDecodedClusters(combed, NN)
    if(!is.na(idx)) {
        ids <- combed[idx,1]

        allIDs <- combed[-idx,1]
        intens <- log2(combed[-idx,2])
        grouped <- split(intens, allIDs)
        summary <- lapply(grouped, median, na.rm = TRUE)
        tmp3[[paste(arrayName)]] <- vector(length = length(idx))
        arrayData <- log2(combed[,2])
        for(i in seq(along = idx)) {
            tmp3[[paste(arrayName)]][i] <- arrayData[idx[i]] - summary[[paste(ids[i])]];
        }   
    }
    else {
        tmp3[[paste(arrayName)]] <- NA
    }
}

## Identify pairs of matching neighbours ##########################

arrayNames = dir(pattern = ".txt$")[1:24]
tmp4 <- list()
NN <- neighboursMatrixForAll();
for(arrayID in 1:24) {

    arrayName = strsplit(arrayNames[arrayID], ".txt")[[1]][1];
    message(paste("Analysing", arrayName));
    
    LOCS<-readLocs(paste(arrayName, "_Grn.locs", sep = ""))
    text <- matrix(unlist(scan(file = paste(arrayName, ".txt", sep = ""), what = list(0, 0, double(), double()), skip = 1, sep = "\t", quiet = TRUE)), ncol = 4)

    combed <- BeadDataPackR:::combineFiles(text, LOCS)
    revind<-cbind(1:1164798,combed[,5])[order(combed[,5]),]

    combed <- combed[order(combed[,5]),]
    newNN <- cbind(combed[1:1164798,1],NN)
    newNN[which(is.na(newNN))] <- 0


    idx <- which(apply(newNN, 1, FUN = identifyMatchingNeighbours, combed[,1]))
    ids <- combed[idx,1]
    idx <- idx[-which(ids == 0)]
    ids <- combed[idx,1]
    
    allIDs <- combed[-idx,1]
    intens <- log2(combed[-idx,2])
    grouped <- split(intens, allIDs)
    summary <- lapply(grouped, median, na.rm = TRUE)
    tmp4[[arrayID]] <- vector(length = length(idx))
    arrayData <- log2(combed[,2])
    for(i in seq(along = idx)) {
        tmp4[[arrayID]][i] <- arrayData[idx[i]] - summary[[paste(ids[i])]];
    }
}

save(tmp1, tmp2, tmp3, tmp4, file = "/home/smith68/Desktop/tmp.rda");

par(mfrow = c(2,2))
boxplot(tmp1, col = c(rep("#5D5DB5", 12), rep("purple", 12)), xaxt = "n", main = "Neighbouring large bright beads", ylab = "deviation from bead-type median log intensity")
axis(1, at = c(1:24), label = c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2", "E_1", "E_2", "F_1", "F_2", "A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2", "E_1", "E_2", "F_1", "F_2"))
abline(h = 0, lty = 3)
legend(x = -0.53, y = -5.2, legend = c("4343238066", "4343238080"), fill = c("#5D5DB5", "purple"))

boxplot(tmp2, col = c(rep("#5D5DB5", 12), rep("purple", 12)), xaxt = "n", main = "Including abnormally low pixels in background region", ylab = "deviation from bead-type median log intensity")
axis(1, at = c(1:24), label = c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2", "E_1", "E_2", "F_1", "F_2", "A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2", "E_1", "E_2", "F_1", "F_2"))
abline(h = 0, lty = 3)
legend(x = -0.53, y = -2, legend = c("4343238066", "4343238080"), fill = c("#5D5DB5", "purple"))

boxplot(tmp3, xaxt = "n", col = c(rep("#5D5DB5", 12), rep("purple", 12)), main = "Within 5 invasions of a non-decoded cluster", ylab = "deviation from bead-type median log intensity" )
axis(1, at = c(1:24), label = c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2", "E_1", "E_2", "F_1", "F_2", "A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2", "E_1", "E_2", "F_1", "F_2"))
abline(h = 0, lty = 3)
legend(x = -0.53, y = -8.2, legend = c("4343238066", "4343238080"), fill = c("#5D5DB5", "purple"))

boxplot(tmp4, xaxt = "n", col = c(rep("#5D5DB5", 12), rep("purple", 12)), main = "Neighbouring beads of the same type", ylab = "deviation from bead-type medianlog intensity")
axis(1, at = c(1:24), label = c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2", "E_1", "E_2", "F_1", "F_2", "A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2", "E_1", "E_2", "F_1", "F_2"))
abline(h = 0, lty = 3)
legend(x = -0.53, y = -5.5, legend = c("4343238066", "4343238080"), fill = c("#5D5DB5", "purple"))

setwd("../")

dev.off()
