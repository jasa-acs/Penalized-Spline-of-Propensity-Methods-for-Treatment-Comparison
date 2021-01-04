rm(list=ls())
#??
DIREC="codes/twoTimePointSimulation/FiguresandTables/" ###where files are stored after combineResult_step2.R

pdf(paste(DIREC, "paperPlots/Figure8_NonLinear_Coverage.pdf", sep=""))


par(mfrow = c(3,3),
    oma = c(2,2,3,2) + 0.1,
    mar = c(1,0.5,0,0) + 0.1, cex.main = 1.5)


sampleSize=500
outcome="_NonLinearOutcome.txt"

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
measure="X95..coverage"
pspp1y=100*(1-pencompAll[c(1,4,7), which(names(pencompAll)==measure)])
aipw1y=100*(1-aipwAll[c(1,4,7), which(names(aipwAll)==measure)])
iptw1y=100*(1-iptwAll[c(1,1,4), which(names(iptwAll)==measure)])
gcompute1y=100*(1-gcomputeAll[c(1,4,1), which(names(gcomputeAll)==measure)])

############ moderate confounding
pspp2y=100*(1-pencompAll[c(2,5,8), which(names(pencompAll)==measure)])
aipw2y=100*(1-aipwAll[c(2,5,8), which(names(aipwAll)==measure)])
iptw2y=100*(1-iptwAll[c(2,2,5), which(names(iptwAll)==measure)])
gcompute2y=100*(1-gcomputeAll[c(2,5,2), which(names(gcomputeAll)==measure)])


############ high confounding
pspp3y=100*(1-pencompAll[c(3,6,9), which(names(pencompAll)==measure)])
aipw3y=100*(1-aipwAll[c(3,6,9), which(names(aipwAll)==measure)])
iptw3y=100*(1-iptwAll[c(3,3,6), which(names(iptwAll)==measure)])
gcompute3y=100*(1-gcomputeAll[c(3,6,3), which(names(gcomputeAll)==measure)])



temp=c(0, 10)

plot(pspp1x, pspp1y,pch=17, ylim=temp, col="red", xaxt='n', yaxt='n', 
     main="", adj=0,  cex=2.5, cex.axis=0.85)
points(aipw1x, aipw1y, pch=8, col="blue",  cex=2.5)
points(gcompute1x,gcompute1y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw1y, pch=22, col="green4", cex=2.5)
axis(side = 2, at = seq(0, 10, 2), tck = -0.03, labels = seq(0, 10, 2))
abline(h=5)
arrows(2, max(temp)-0.5, 2, max(temp), length=0.1)
text(1.75, max(temp)-1, paste("G-COMP=",round(gcompute1y[2],2)), col="black")
mtext(expression(bold(paste("Confounding=Low"))), line = 1.5, at=2, side=3)
mtext(expression(bold(paste(Delta, "11"))), line = 0, at=2, side=3)


plot(pspp1x, pspp2y,pch=17, ylim=c(min(temp), max(temp)),col="red", yaxt='n', xaxt='n', main="", cex=2.5)
points(aipw1x, aipw2y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute2y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw2y, pch=22, col="green4", cex=2.5)
abline(h=5)
arrows(2, max(temp)-0.5, 2, max(temp), length=0.1)
text(1.75, max(temp)-1, paste("G-COMP=",round(gcompute2y[2],2)), col="black")

arrows(3, max(temp)-0.5, 3, max(temp), length=0.1)
text(2.75, max(temp)-1, paste("IPTW=",round(iptw2y[3],2)), col="green4")

mtext(expression(bold(paste("Confounding=Moderate"))), line = 1.5, at=2, side=3)
mtext(expression(bold(paste(Delta, "11"))), line = 0, at=2, side=3)


plot(pspp1x, pspp3y,pch=17, ylim=c(min(temp), max(temp)),col="red",  yaxt='n',xaxt='n', main="",cex=2.5)
points(aipw1x, aipw3y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute3y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw3y, pch=22, col="green4", cex=2.5)
axis(side = 4, at =  seq(0, 10, 2), tck = -0.03, labels=seq(0, 10, 2))
abline(h=5)
arrows(3, max(temp)-0.5, 3, max(temp), length=0.1)
text(2.75, max(temp)-1, paste("IPTW=",round(iptw3y[3],2)), col="green4")

arrows(2, max(temp)-0.5, 2, max(temp), length=0.1)
text(1.75, max(temp)-1, paste("G-COMP=",round(gcompute3y[2],2)), col="black")

mtext(expression(bold(paste("Confounding=High"))), line = 1.5, at=2, side=3)
mtext(expression(bold(paste(Delta, "11"))), line = 0, at=2, side=3)


###################################################################
### delta 10 - delta00

measure="X95..coverage.1"
pspp1y=100*(1-pencompAll[c(1,4,7), which(names(pencompAll)==measure)])
aipw1y=100*(1-aipwAll[c(1,4,7), which(names(aipwAll)==measure)])
iptw1y=100*(1-iptwAll[c(1,1,4), which(names(iptwAll)==measure)])
gcompute1y=100*(1-gcomputeAll[c(1,4,1), which(names(gcomputeAll)==measure)])

############ moderate confounding
pspp2y=100*(1-pencompAll[c(2,5,8), which(names(pencompAll)==measure)])
aipw2y=100*(1-aipwAll[c(2,5,8), which(names(aipwAll)==measure)])
iptw2y=100*(1-iptwAll[c(2,2,5), which(names(iptwAll)==measure)])
gcompute2y=100*(1-gcomputeAll[c(2,5,2), which(names(gcomputeAll)==measure)])


############ high confounding
pspp3y=100*(1-pencompAll[c(3,6,9), which(names(pencompAll)==measure)])
aipw3y=100*(1-aipwAll[c(3,6,9), which(names(aipwAll)==measure)])
iptw3y=100*(1-iptwAll[c(3,3,6), which(names(iptwAll)==measure)])
gcompute3y=100*(1-gcomputeAll[c(3,6,3), which(names(gcomputeAll)==measure)])


plot(pspp1x, pspp1y,pch=17, ylim=temp, col="red", xaxt='n', yaxt='n',
     main="", adj=0,cex=2.5, cex.axis=0.85)
points(aipw1x, aipw1y, pch=8, col="blue",  cex=2.5)
points(gcompute1x,gcompute1y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw1y, pch=22, col="green4", cex=2.5)
axis(side = 2, at = seq(0, 10, 2), tck = -0.03, labels = seq(0, 10, 2))
abline(h=5)
arrows(1.9, max(temp)-0.5, 1.9, max(temp), length=0.1)
text(1.85, max(temp)-1.2, paste("G-COMP=",round(gcompute1y[2],2)), col="black")
arrows(3, max(temp)-0.5, 3, max(temp), length=0.1)
text(2.75, max(temp)-1, paste("IPTW=",round(iptw1y[3],2)), col="green4")
mtext(expression(bold(paste(Delta, "10"))), line = 0, at=2, side=3)


plot(pspp1x, pspp2y,pch=17, ylim=c(min(temp), max(temp)),col="red", yaxt='n', xaxt='n', cex=2.5)
points(aipw1x, aipw2y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute2y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw2y, pch=22, col="green4", cex=2.5)
abline(h=5)

arrows(3, max(temp)-0.5, 3, max(temp), length=0.1)
text(2.75, max(temp)-1, paste("IPTW=",round(iptw2y[3],2)), col="green4")

arrows(2, max(temp)-0.5, 2, max(temp), length=0.1)
text(1.9, max(temp)-3, paste("G-COMP=",round(gcompute2y[2],2)), col="black")
text(1.9, max(temp)-2, paste("IPTW=",round(iptw2y[2],2)), col="green4")
text(1.91, max(temp)-1, paste("AIPTW=",round(aipw2y[2],2)), col="blue")

arrows(1, max(temp)-0.5, 1, max(temp), length=0.1)
text(1.25, max(temp)-1, paste("IPTW=",round(iptw2y[1],2)), col="green4")
mtext(expression(bold(paste(Delta, "10"))), line = 0, at=2, side=3)


plot(pspp1x, pspp3y,pch=17, ylim=c(min(temp), max(temp)),col="red",  yaxt='n',xaxt='n',cex=2.5)
points(aipw1x, aipw3y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute3y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw3y, pch=22, col="green4", cex=2.5)
axis(side = 4, at = seq(0, 10, 2), tck = -0.03, labels=seq(0, 10, 2))
abline(h=5)
arrows(3, max(temp)-0.5, 3, max(temp), length=0.1)
text(2.75, max(temp)-1, paste("IPTW=",round(iptw3y[3],2)), col="green4")

arrows(2, max(temp)-0.5, 2, max(temp), length=0.1)
text(1.9, max(temp)-3, paste("G-COMP=",round(gcompute3y[2],2)), col="black")
text(1.9, max(temp)-2, paste("IPTW=",round(iptw3y[2],2)), col="green4")
text(1.9, max(temp)-1, paste("AIPTW=",round(aipw3y[2],2)), col="blue")

arrows(1, max(temp)-0.5, 1, max(temp), length=0.1)
text(1.2, max(temp)-1, paste("IPTW=",round(iptw3y[1],2)), col="green4")
mtext(expression(bold(paste(Delta, "10"))), line = 0, at=2, side=3)


###################################################################
### delta 01 - delta00


#############low confounding
measure="X95..coverage.2"
pspp1y=100*(1-pencompAll[c(1,4,7), which(names(pencompAll)==measure)])
aipw1y=100*(1-aipwAll[c(1,4,7), which(names(aipwAll)==measure)])
iptw1y=100*(1-iptwAll[c(1,1,4), which(names(iptwAll)==measure)])
gcompute1y=100*(1-gcomputeAll[c(1,4,1), which(names(gcomputeAll)==measure)])

############ moderate confounding
pspp2y=100*(1-pencompAll[c(2,5,8), which(names(pencompAll)==measure)])
aipw2y=100*(1-aipwAll[c(2,5,8), which(names(aipwAll)==measure)])
iptw2y=100*(1-iptwAll[c(2,2,5), which(names(iptwAll)==measure)])
gcompute2y=100*(1-gcomputeAll[c(2,5,2), which(names(gcomputeAll)==measure)])


############ high confounding
pspp3y=100*(1-pencompAll[c(3,6,9), which(names(pencompAll)==measure)])
aipw3y=100*(1-aipwAll[c(3,6,9), which(names(aipwAll)==measure)])
iptw3y=100*(1-iptwAll[c(3,3,6), which(names(iptwAll)==measure)])
gcompute3y=100*(1-gcomputeAll[c(3,6,3), which(names(gcomputeAll)==measure)])


plot(pspp1x, pspp1y,pch=17, ylim=temp, col="red", xaxt='n', yaxt='n', 
     main="", adj=0,cex=2.5, cex.axis=0.85)
points(aipw1x, aipw1y, pch=8, col="blue",  cex=2.5)
points(gcompute1x,gcompute1y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw1y, pch=22, col="green4", cex=2.5)
axis(side = 2, at = seq(0, 10, 2), tck = -0.03, labels =seq(0, 10, 2))
abline(h=5)
arrows(3, max(temp)-0.5, 3, max(temp), length=0.1)
text(2.75, max(temp)-1, paste("IPTW=",round(iptw1y[3],2)), col="green4")
axis(side = 1, at = c(1, 2, 3), labels = c("A", "B", "C"), tck = -0.01)
mtext(expression(bold(paste(Delta, "01"))), line = 0, at=2, side=3)


plot(pspp1x, pspp2y,pch=17, ylim=c(min(temp), max(temp)),col="red", yaxt='n', xaxt='n', cex=2.5)
points(aipw1x, aipw2y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute2y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw2y, pch=22, col="green4", cex=2.5)
abline(h=5)
arrows(3, max(temp)-0.5, 3, max(temp), length=0.1)
text(2.75, max(temp)-1, paste("IPTW=",round(iptw2y[3],2)), col="green4")
axis(side = 1, at = c(1, 2, 3), labels = c("A", "B", "C"), tck = -0.01)
mtext(expression(bold(paste(Delta, "01"))), line = 0, at=2, side=3)
mtext("Model Specification", line = 2, at=2, side=1)


plot(pspp1x, pspp3y,pch=17, ylim=c(min(temp), max(temp)),col="red",  yaxt='n',xaxt='n',cex=2.5)
points(aipw1x, aipw3y, pch=8, col="blue", cex=2.5)
points(gcompute1x,gcompute3y, pch=1, col="black", cex=2.5)
points(iptw1x,iptw3y, pch=22, col="green4", cex=2.5)
axis(side = 4, at = seq(0, 10, 2), tck = -0.03, labels=seq(0, 10, 2))
abline(h=5)
arrows(3, max(temp)-0.5, 3, max(temp), length=0.1)
text(2.75, max(temp)-1, paste("IPTW=",round(iptw3y[3],2)), col="green4")
axis(side = 1, at = c(1, 2, 3), labels = c("A", "B", "C"), tck = -0.01)
mtext(expression(bold(paste(Delta, "01"))), line = 0, at=2, side=3)


#title( paste("Non-coverage Rate", "\n", "\n","Nonlinear Outcome", sep=""), outer = TRUE )


#par(fig=c(0,1,0,0.35), new=TRUE)
#plot_colors <- c("red","blue", "black", "green4")
#text <- c("PENCOMP", "AIPTW", "G-COMPutation", 
#          "IPTW")
#legend("center",legend = text, text.width = max(sapply(text, strwidth)),##xpd = TRUE tells OK to plot outside the region
#       col=plot_colors, cex=1, horiz = TRUE, pch=c(17,8,1,22), xpd = TRUE, text.font=1.5, pt.cex=1.5)

#########################
#par(fig=c(0,1,0,0.45), new=TRUE)
#legend("center",legend = "Model Specification", cex=1, horiz = TRUE, bty = "n", lty=NULL, xpd=TRUE, text.font=1.5)


dev.off()










