library(dplyr)
library(MASS)
library(gtools) 
library(matrixStats)

rm(list = ls())
cex0=2.5 ; cex1=2.2
cex_axis=2.2 ; cex_lab=2.2 ; cex_main=2.2 ; cex_leyend=1.25
cex2=2.5
dir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

df=read.table("data/probability_plot.txt",header = T)

maxis <- df %>%
  group_by(Z) %>%
  summarise(Max_prob = max(final_prob))
maxis20=maxis[which(maxis$Max_prob>-4),]

lm1=lm(maxis20$Max_prob ~ maxis20$Z)

pdf(paste("probCom.pdf",sep=""))
par(mar = c(5, 6, 4, 2) + 0.1)  # Adjusting margin
par(cex.axis = cex_axis, cex.lab = cex_lab , cex_main=cex_lab)

plot(x=df$Z,y=df$final_prob,
     cex=2.5,pch=16,col="blue",ylim = c(min(df$final_prob),0),
     xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
     

abline(lm1,lwd=2)
axis(1,cex=cex0,line=0)
axis(2,cex=cex0,line=0)
mtext(expression(paste(~tilde(K), "(x)")), side = 1, line = 4,cex=cex0) 
mtext(expression(paste(log[10],"P(x)")), side = 2, line = 3.5,cex=cex0) 


abline(lm1)
dev.off()
