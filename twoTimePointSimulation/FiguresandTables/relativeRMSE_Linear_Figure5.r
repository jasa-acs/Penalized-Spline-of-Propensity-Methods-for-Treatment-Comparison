rm(list=ls())
#??
DIREC="codes/twoTimePointSimulation/FiguresandTables/" ###where files are stored after combineResult_step2.R

pdf(paste(DIREC, "paperPlots/Figure5_Linear_RelativeRMSE_IPTW.pdf", sep=""))


par(mfrow = c(3,3),
    oma = c(2,2,3,2) + 0.1,
    mar = c(1,0.5,0,0) + 0.1, cex.main = 1.5)


sampleSize=500
outcome="_LinearOutcome.txt"

pencompAll=read.table(paste(DIREC, "PENCOMP_Results/sampleSize", sampleSize, outcome,sep=""), header=T, sep="\t") 
iptwAll=read.table(paste(DIREC, "IPTW_Results/sampleSize", sampleSize, outcome,sep=""), header=T, sep="\t") 
gcomputeAll=read.table(paste(DIREC, "gcompute_Results/sampleSize", sampleSize, outcome,sep=""), header=T, sep="\t") 
aipwAll=read.table(paste(DIREC, "AIPTW_Results/sampleSize", sampleSize, outcome,sep=""), header=T, sep="\t") 


###################################################################
### delta 11 - delta00

pspp1x=c(1,2,3)  ##1=both correct, 2=mispred model, 3=misweight model
aipw1x=c(1,2,3)
iptw1x=c(1,2,3)
gcompute1x=c(1,2,3)

#############low confounding
measure="RMSE"
pspp1y=pencompAll[c(1,4,7), which(names(pencompAll)==measure)]
aipw1y=aipwAll[c(1,4,7), which(names(aipwAll)==measure)]
iptw1y=iptwAll[c(1,1,4), which(names(iptwAll)==measure)]
gcompute1y=gcomputeAll[c(1,4,1), which(names(gcomputeAll)==measure)]

pspp1y=pspp1y/iptw1y[1]  ##relative to correct iptw rmse
aipw1y=aipw1y/iptw1y[1]  ##relative to correct iptw rmse
gcompute1y=gcompute1y/iptw1y[1]  ##relative to correct iptw rmse
iptw1y=iptw1y/iptw1y[1]  ##relative to correct iptw rmse


############ moderate confounding
pspp2y=pencompAll[c(2,5,8), which(names(pencompAll)==measure)]
aipw2y=aipwAll[c(2,5,8), which(names(aipwAll)==measure)]
iptw2y=iptwAll[c(2,2,5), which(names(iptwAll)==measure)]
gcompute2y=gcomputeAll[c(2,5,2), which(names(gcomputeAll)==measure)]

pspp2y=pspp2y/iptw2y[1]  ##relative to correct iptw rmse
aipw2y=aipw2y/iptw2y[1]  ##relative to correct iptw rmse
gcompute2y=gcompute2y/iptw2y[1]  ##relative to correct iptw rmse
iptw2y=iptw2y/iptw2y[1]  ##relative to correct iptw rmse



############ high confounding
pspp3y=pencompAll[c(3,6,9), which(names(pencompAll)==measure)]
aipw3y=aipwAll[c(3,6,9), which(names(aipwAll)==measure)]
iptw3y=iptwAll[c(3,3,6), which(names(iptwAll)==measure)]
gcompute3y=gcomputeAll[c(3,6,3), which(names(gcomputeAll)==measure)]


pspp3y=pspp3y/iptw3y[1]  ##relative to correct iptw rmse
aipw3y=aipw3y/iptw3y[1]  ##relative to correct iptw rmse
gcompute3y=gcompute3y/iptw3y[1]  ##relative to correct iptw rmse
iptw3y=iptw3y/iptw3y[1]  ##relative to correct iptw rmse

temp=c(0.12, 1.2)

plot(pspp1x, pspp1y,pch=17, ylim=temp, col="red", xaxt='n', yaxt='n', 
     main="", adj=0, cex=2.5, cex.axis=0.85)
points(aipw1x, aipw1y, pch=8, col="blue",  cex=2.5)
points(gcompute1x,gcompute1y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw1y, pch=22, col="green4", cex=2.5)
axis(side = 2, at = seq(0.2, 1.2, 0.2), tck = -0.03, labels = seq(0.2, 1.2, 0.2))
abline(h=1)
arrows(2, max(temp)-0.02, 2, max(temp), length=0.1)
text(2, max(temp)-0.10, paste("G-COMP=",round(gcompute1y[2],2)), col="black")
mtext(expression(bold(paste("Confounding=Low"))), line = 1.5, at=2, side=3)
mtext(expression(bold(paste(Delta, "11"))), line = 0, at=2, side=3)

plot(pspp1x, pspp2y,pch=17, ylim=c(min(temp), max(temp)),col="red", yaxt='n', xaxt='n', main="", cex=2.5)
points(aipw1x, aipw2y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute2y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw2y, pch=22, col="green4", cex=2.5)
abline(h=1)
arrows(3, max(temp)-0.02, 3, max(temp), length=0.1)
text(2.75, max(temp)-0.12, paste("IPTW=",round(iptw2y[3],2)), col="green4")
mtext(expression(bold(paste("Confounding=Moderate"))), line = 1.5, at=2, side=3)
mtext(expression(bold(paste(Delta, "11"))), line = 0, at=2, side=3)

plot(pspp1x, pspp3y,pch=17, ylim=c(min(temp), max(temp)),col="red",  yaxt='n',xaxt='n', main="",cex=2.5)
points(aipw1x, aipw3y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute3y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw3y, pch=22, col="green4", cex=2.5)
axis(side = 4, at = seq(0.2, 1.2, 0.2), tck = -0.03, labels = seq(0.2, 1.2, 0.2))
abline(h=1)
arrows(3, max(temp)-0.02, 3, max(temp), length=0.1)
text(2.75, max(temp)-0.12, paste("IPTW=",round(iptw3y[3],2)), col="green4")
mtext(expression(bold(paste("Confounding=High"))), line = 1.5, at=2, side=3)
mtext(expression(bold(paste(Delta, "11"))), line = 0, at=2, side=3)


###################################################################
### delta 10 - delta00
measure="RMSE.1"
pspp1y=pencompAll[c(1,4,7), which(names(pencompAll)==measure)]
aipw1y=aipwAll[c(1,4,7), which(names(aipwAll)==measure)]
iptw1y=iptwAll[c(1,1,4), which(names(iptwAll)==measure)]
gcompute1y=gcomputeAll[c(1,4,1), which(names(gcomputeAll)==measure)]

pspp1y=pspp1y/iptw1y[1]  ##relative to correct iptw rmse
aipw1y=aipw1y/iptw1y[1]  ##relative to correct iptw rmse
gcompute1y=gcompute1y/iptw1y[1]  ##relative to correct iptw rmse
iptw1y=iptw1y/iptw1y[1]  ##relative to correct iptw rmse


############ moderate confounding
pspp2y=pencompAll[c(2,5,8), which(names(pencompAll)==measure)]
aipw2y=aipwAll[c(2,5,8), which(names(aipwAll)==measure)]
iptw2y=iptwAll[c(2,2,5), which(names(iptwAll)==measure)]
gcompute2y=gcomputeAll[c(2,5,2), which(names(gcomputeAll)==measure)]

pspp2y=pspp2y/iptw2y[1]  ##relative to correct iptw rmse
aipw2y=aipw2y/iptw2y[1]  ##relative to correct iptw rmse
gcompute2y=gcompute2y/iptw2y[1]  ##relative to correct iptw rmse
iptw2y=iptw2y/iptw2y[1]  ##relative to correct iptw rmse



############ high confounding
pspp3y=pencompAll[c(3,6,9), which(names(pencompAll)==measure)]
aipw3y=aipwAll[c(3,6,9), which(names(aipwAll)==measure)]
iptw3y=iptwAll[c(3,3,6), which(names(iptwAll)==measure)]
gcompute3y=gcomputeAll[c(3,6,3), which(names(gcomputeAll)==measure)]


pspp3y=pspp3y/iptw3y[1]  ##relative to correct iptw rmse
aipw3y=aipw3y/iptw3y[1]  ##relative to correct iptw rmse
gcompute3y=gcompute3y/iptw3y[1]  ##relative to correct iptw rmse
iptw3y=iptw3y/iptw3y[1]  ##relative to correct iptw rmse


plot(pspp1x, pspp1y,pch=17, ylim=temp, col="red", xaxt='n', yaxt='n', 
     main="", adj=0,cex=2.5, cex.axis=0.85)
points(aipw1x, aipw1y, pch=8, col="blue",  cex=2.5)
points(gcompute1x,gcompute1y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw1y, pch=22, col="green4", cex=2.5)
axis(side = 2, at = seq(0.2, 1.2, 0.2), tck = -0.03, labels = seq(0.2, 1.2, 0.2))
abline(h=1)
arrows(2, max(temp)-0.02, 2, max(temp), length=0.1)
text(1.75, max(temp)-0.10, paste("G-COMP=",round(gcompute1y[2],2)), col="black")
arrows(3, max(temp)-0.02, 3, max(temp), length=0.1)
text(2.75, max(temp)-0.12, paste("IPTW=",round(iptw1y[3],2)), col="green4")
mtext(expression(bold(paste(Delta, "10"))), line = 0, at=2, side=3)


plot(pspp1x, pspp2y,pch=17, ylim=c(min(temp), max(temp)),col="red", yaxt='n', xaxt='n', cex=2.5)
points(aipw1x, aipw2y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute2y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw2y, pch=22, col="green4", cex=2.5)
abline(h=1)
arrows(3, max(temp)-0.02, 3, max(temp), length=0.1)
text(2.75, max(temp)-0.12, paste("IPTW=",round(iptw2y[3],2)), col="green4")
mtext(expression(bold(paste(Delta, "10"))), line = 0, at=2, side=3)

plot(pspp1x, pspp3y,pch=17, ylim=c(min(temp), max(temp)),col="red",  yaxt='n',xaxt='n',cex=2.5)
points(aipw1x, aipw3y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute3y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw3y, pch=22, col="green4", cex=2.5)
axis(side = 4, at = seq(0.2, 1.2, 0.2), tck = -0.03, labels = seq(0.2, 1.2, 0.2))
abline(h=1)
arrows(3, max(temp)-0.02, 3, max(temp), length=0.1)
text(2.75, max(temp)-0.12, paste("IPTW=",round(iptw3y[3],2)), col="green4")
mtext(expression(bold(paste(Delta, "10"))), line = 0, at=2, side=3)


###################################################################
### delta 01 - delta00

#############low confounding
measure="RMSE.2"
pspp1y=pencompAll[c(1,4,7), which(names(pencompAll)==measure)]
aipw1y=aipwAll[c(1,4,7), which(names(aipwAll)==measure)]
iptw1y=iptwAll[c(1,1,4), which(names(iptwAll)==measure)]
gcompute1y=gcomputeAll[c(1,4,1), which(names(gcomputeAll)==measure)]

pspp1y=pspp1y/iptw1y[1]  ##relative to correct iptw rmse
aipw1y=aipw1y/iptw1y[1]  ##relative to correct iptw rmse
gcompute1y=gcompute1y/iptw1y[1]  ##relative to correct iptw rmse
iptw1y=iptw1y/iptw1y[1]  ##relative to correct iptw rmse


############ moderate confounding
pspp2y=pencompAll[c(2,5,8), which(names(pencompAll)==measure)]
aipw2y=aipwAll[c(2,5,8), which(names(aipwAll)==measure)]
iptw2y=iptwAll[c(2,2,5), which(names(iptwAll)==measure)]
gcompute2y=gcomputeAll[c(2,5,2), which(names(gcomputeAll)==measure)]

pspp2y=pspp2y/iptw2y[1]  ##relative to correct iptw rmse
aipw2y=aipw2y/iptw2y[1]  ##relative to correct iptw rmse
gcompute2y=gcompute2y/iptw2y[1]  ##relative to correct iptw rmse
iptw2y=iptw2y/iptw2y[1]  ##relative to correct iptw rmse



############ high confounding
pspp3y=pencompAll[c(3,6,9), which(names(pencompAll)==measure)]
aipw3y=aipwAll[c(3,6,9), which(names(aipwAll)==measure)]
iptw3y=iptwAll[c(3,3,6), which(names(iptwAll)==measure)]
gcompute3y=gcomputeAll[c(3,6,3), which(names(gcomputeAll)==measure)]


pspp3y=pspp3y/iptw3y[1]  ##relative to correct iptw rmse
aipw3y=aipw3y/iptw3y[1]  ##relative to correct iptw rmse
gcompute3y=gcompute3y/iptw3y[1]  ##relative to correct iptw rmse
iptw3y=iptw3y/iptw3y[1]  ##relative to correct iptw rmse

plot(pspp1x, pspp1y,pch=17, ylim=temp, col="red", xaxt='n', yaxt='n', 
     main="", adj=0,cex=2.5, cex.axis=0.85)
points(aipw1x, aipw1y, pch=8, col="blue",  cex=2.5)
points(gcompute1x,gcompute1y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw1y, pch=22, col="green4", cex=2.5)
axis(side = 2, at = seq(0.2, 1.2, 0.2), tck = -0.03, labels =seq(0.2, 1.2, 0.2))
abline(h=1)
#arrows(2, max(temp)-0.02, 2, max(temp), length=0.1)
#text(1.75, max(temp)-0.10, paste("G-COMPutation=",round(gcompute1y[2],2)))
arrows(3, max(temp)-0.02, 3, max(temp), length=0.1)
text(2.75, max(temp)-0.12, paste("IPTW=",round(iptw1y[3],2)), col="green4")
axis(side = 1, at = c(1, 2, 3), labels = c("A", "B", "C"), tck = -0.01)
#mtext("A", line = 0, at=1, side=1)
#mtext("B", line = 0, at=2, side=1)
#mtext("C", line = 0, at=3, side=1)
mtext(expression(bold(paste(Delta, "01"))), line = 0, at=2, side=3)


plot(pspp1x, pspp2y,pch=17, ylim=c(min(temp), max(temp)),col="red", yaxt='n', xaxt='n', cex=2.5)
points(aipw1x, aipw2y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute2y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw2y, pch=22, col="green4", cex=2.5)
abline(h=1)
arrows(3, max(temp)-0.02, 3, max(temp), length=0.1)
text(2.75, max(temp)-0.12, paste("IPTW=",round(iptw2y[3],2)), col="green4")
axis(side = 1, at = c(1, 2, 3), labels = c("A", "B", "C"), tck = -0.01)
#mtext("A", line = 0, at=1, side=1)
#mtext("B", line = 0, at=2, side=1)
#mtext("C", line = 0, at=3, side=1)
mtext(expression(bold(paste(Delta, "01"))), line = 0, at=2, side=3)
mtext("Model Specification", line = 2, at=2, side=1)


plot(pspp1x, pspp3y,pch=17, ylim=c(min(temp), max(temp)),col="red",  yaxt='n',xaxt='n',cex=2.5)
points(aipw1x, aipw3y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute3y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw3y, pch=22, col="green4", cex=2.5)
axis(side = 4, at = seq(0.2, 1.2, 0.2), tck = -0.03, labels = seq(0.2, 1.2, 0.2))
abline(h=1)
arrows(3, max(temp)-0.02, 3, max(temp), length=0.1)
text(2.75, max(temp)-0.12, paste("IPTW=",round(iptw3y[3],2)), col="green4")
axis(side = 1, at = c(1, 2, 3), labels = c("A", "B", "C"), tck = -0.01)
#mtext("A", line = 0, at=1, side=1)
#mtext("B", line = 0, at=2, side=1)
#mtext("C", line = 0, at=3, side=1)
mtext(expression(bold(paste(Delta, "01"))), line = 0, at=2, side=3)


#title( paste("RMSE Divided by RMSE of IPTW(A)", "\n", "\n","Linear Outcome", sep=""), outer = TRUE )


#par(fig=c(0,1,0,0.35), new=TRUE)
#plot_colors <- c("red","blue", "black", "green4")
#text <- c("PENCOMP", "AIPTW", "G-COMPutation", 
#          "IPTW")
#legend("center",legend = text, text.width = max(sapply(text, strwidth)),##xpd = TRUE tells OK to plot outside the region
#       col=plot_colors, cex=1, horiz = TRUE, pch=c(17,8,1,22), xpd = TRUE, text.font=1.5, pt.cex=1.5)

##########################
#par(fig=c(0,1,0,0.45), new=TRUE)
#legend("center",legend = "Model Specification", cex=1, horiz = TRUE, bty = "n", lty=NULL, xpd=TRUE, text.font=1.5)


dev.off()










