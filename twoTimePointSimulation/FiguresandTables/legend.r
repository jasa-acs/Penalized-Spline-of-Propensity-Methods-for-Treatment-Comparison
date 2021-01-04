rm(list=ls())
#### relative bias to RMSE of correct IPTW
#??
DIREC="codes/twoTimePointSimulation/FiguresandTables/" ###where files are stored after combineResult_step2.R

png(paste(DIREC, "paperPlots/plotsVersion5/legend.png", sep=""))
plot(1, type="n", axes=FALSE, xlab="", ylab="")

par(fig=c(0,1,0,0.35), new=TRUE)
plot_colors <- c("red","blue", "black", "green4")
text <- c("PENCOMP", "AIPTW", "G-COMP", 
          "IPTW")
legend("center",legend = text, text.width = max(sapply(text, strwidth)),##xpd = TRUE tells OK to plot outside the region
      col=plot_colors, cex=1, horiz = TRUE, pch=c(17,8,1,22), xpd = TRUE, text.font=2, pt.cex=2)

#########################
#par(fig=c(0,1,0,0.5), new=TRUE)
#legend("center",legend = "Model Specification", cex=1, horiz = TRUE, bty = "n", lty=NULL, xpd=TRUE, text.font=1.5)


dev.off()










