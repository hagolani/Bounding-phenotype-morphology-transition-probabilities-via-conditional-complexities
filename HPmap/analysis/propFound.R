library(dplyr)
library(stringr)

rm(list = ls())
dir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)



# List all files that match the pattern "dataNeigh*_new.txt" in the "../" directory
files <- list.files(path = "data/dataAll/", pattern = "^dataToAll.*_new_KC.txt$")
# Extract the letter sequence between 'dataToAll' and '_new_KC.txt'
sequences <- gsub("dataToAll([A-Za-z]+)_new_KC.txt", "\\1", files)

allbmbs = read.table("data/wkc/allHP25.txt", header = T)
allbmbs <- allbmbs[!grepl("ParentPhenotyope", allbmbs[, 1]), ]

allbmbs$count <- as.numeric(allbmbs$count)
allbmbs$proportion <- as.numeric(allbmbs$proportion)
allbmbs$KC_col1 <- as.numeric(allbmbs$KC_col1)
allbmbs$KC_col4 <- as.numeric(allbmbs$KC_col4)
allbmbs$Conditional_Complexity <- as.numeric(allbmbs$Conditional_Complexity)


alldata=data.frame()
for (binary in sequences) {
  
allfound=read.table(paste("data/dataAll/dataToAll",binary,"_new_KC.txt",sep=""),header = T)
bmb=allbmbs[(which(allbmbs$ParentPhenotyope==binary)),]
if(nrow(bmb)==0)next

pacom=bmb$KC_col1[1]
bmb$MutantPhenotype2 <- substr(bmb$MutantPhenotype, 1, nchar(bmb$MutantPhenotype) - 1)

#found
df1n <- allfound %>%
  filter(
    sapply(Phenotype, function(pheno) any(str_detect(pheno, bmb$MutantPhenotype2)))
  )
#not found
df2n <- allfound %>%
  filter(
    sapply(Phenotype, function(pheno) !any(str_detect(pheno, bmb$MutantPhenotype2)))
  )

foundckc=median(df1n$Conditional_Complexity)
notfoundckc=median(df2n$Conditional_Complexity)

data.tmp=data.frame(phenotype=binary,
                    complexity=pacom,
                    foundckc,notfoundckc)
alldata=rbind(alldata,data.tmp)
}


write.table(alldata,"data/alldata.dat")

alldata=read.table("data/alldata.dat")


pp=alldata$foundckc-alldata$notfoundckc
pp=pp/(max(abs(pp)))

ppw=pp ; ppw=as.data.frame(ppw)
ppw$exp="HP"
write.table(ppw,"data/HP_FoN.dat")

pdf("HPmap_foundOrNot.pdf")
par(cex.axis = 2.2, cex.lab = 2.2)
par(mar = c(5, 7, 4, 2))  # c(bottom, left, top, right)
boxplot(pp,ylab=expression(paste(~tilde(K), "(y|x)"[found], " - ",~tilde(K), "(y|x)"[not_found])))#points( jitter(rep(1,length(pp)),5), pp)

dev.off()

#plot(alldata$complexity,alldata$foundckc-alldata$notfoundckc)
length(pp)


