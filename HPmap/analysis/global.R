library(dplyr)
rm(list = ls())

dir=dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(dir)
cex0=2.5 ; cex1=2.2
df1=read.table("data/probability_plot.txt",header = T)

maxis <- df1 %>%
  group_by(Z) %>%
  summarise(Max_prob = max(final_prob))
maxis20=maxis[which(maxis$Max_prob>-10),]

lm1=lm(maxis20$Max_prob ~ maxis20$Z)

color1="blue"#rgb(0, 0, 0, 0.5)
pdf("probCom.pdf")
par(mar = c(5, 6, 4, 2) + 0.1)  # Adjusting margin
par(cex.axis = cex1, cex.lab = cex1)
#plot(x=jitter(df1$Z),y=df1$final_prob,
ymin=min(df1$final_prob)
plot(x=df1$Z,y=df1$final_prob,
     cex=cex0,
     xaxt = 'n', yaxt = 'n', xlab = "", ylab = "",
   pch=16,col=color1,ylim = c(ymin,0))
   #xlab = expression(paste(~tilde(K), "(x)")),
   #ylab = expression(paste(log[10],"P(x)") ))
abline(lm1,lwd=2)
axis(1,cex=cex0,line=0)
axis(2,cex=cex0,line=0)
mtext(expression(paste(~tilde(K), "(x)")), side = 1, line = 4,cex=cex0)  # Increase the 'line' value to move the label further down
mtext(expression(paste(log[10],"P(x)")), side = 2, line = 3.5,cex=cex0)  # Increase the 'line' value to move the label further down


dev.off()
