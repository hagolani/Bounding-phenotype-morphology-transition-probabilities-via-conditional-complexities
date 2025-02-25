library(dplyr)


rm(list = ls())
dir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

fon = read.table("allFoN.dat",row.names = NULL)
fon$exp
for (iji in 1:nrow(fon)) {
  if(fon$exp[iji]=="circadian")fon$exp[iji]="circadian rhythm"
  if(fon$exp[iji]=="cusps")fon$exp[iji]="tooth model: cusps"
  if(fon$exp[iji]=="opc")fon$exp[iji]="tooth model: OPCr"
  if(fon$exp[iji]=="MMP")fon$exp[iji]="vector matrix map"
  if(fon$exp[iji]=="polyo")fon$exp[iji]="polyominos"
  if(fon$exp[iji]=="HP")fon$exp[iji]="HP protein map"
}

cexx=1
pdf("FoN.pdf")
par(mar = c(10, 4, 4, 2),mgp = c(3, 0.5, 0))
par(cex=cexx)
boxplot(fon$ppw ~ fon$exp,
        ylab = "",#expression(paste(~tilde(K), "(y|x)"[found], " - ", ~tilde(K), "(y|x)"[not_found])),
        xaxt = "n",xlab = "")
#abline(h=0,lwd=2,lty=2)
# Extract the levels for the x-axis labels
x_labels <- levels(as.factor(fon$exp))

# Add x-axis tick marks
axis(side = 1, at = 1:length(x_labels), labels = FALSE) # Tick marks without labels

# Add custom x-axis labels at a 45-degree angle
mtext(side = 2, line = 2, cex=cexx,
      expression(paste(~tilde(K), "(y|x)"[found], " - ", ~tilde(K), "(y|x)"[not_found])))

text(x = 1:length(x_labels), 
     y = par("usr")[3] - 0.15, # Adjust position below the x-axis
     labels = x_labels, 
     srt = 45,       # Rotate labels to 45 degrees
     adj = 1,        # Align labels properly
     xpd = TRUE,     # Allow drawing outside plot region
     cex = 0.9)      # Adjust label font size
dev.off()



