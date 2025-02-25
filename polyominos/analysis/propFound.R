library(dplyr)


rm(list = ls())
dir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

allfound=read.table("data/alldataNeigh_new_KC.txt",header = T)
numbers = unique(allfound$ParentPhenotyope)
alldata=data.frame()
for (binary in numbers) {
  
tmp.dat=allfound[which(allfound$ParentPhenotyope==binary),]
pacom=allfound[which(allfound$ParentPhenotyope==as.numeric(binary)),4][1]
#found
df1n <- tmp.dat[which(tmp.dat$proportion>0),]

#not found
df2n <- tmp.dat[which(tmp.dat$proportion==0),]


foundckc=median(df1n$Conditional_Complexity)
notfoundckc=median(df2n$Conditional_Complexity)

data.tmp=data.frame(phenotype=binary,
                    complexity=pacom,
                    foundckc,notfoundckc)
alldata=rbind(alldata,data.tmp)
}

pp=alldata$foundckc-alldata$notfoundckc
pp=pp/max(abs(pp))

ppw=pp ; ppw=as.data.frame(ppw)
ppw$exp="polyo"
write.table(ppw,"polyo_FoN.dat")

pdf("polyo_norm_foundOrNot.pdf")
par(mar = c(5, 7, 4, 2))  # c(bottom, left, top, right)
par(cex.axis = 2.2, cex.lab = 2.2)
boxplot(pp,ylab=expression(paste(~tilde(K), "(y|x)"[found], " - ",~tilde(K), "(y|x)"[not_found])))


dev.off()
length(pp)
#plot(alldata$complexity,alldata$foundckc-alldata$notfoundckc)



