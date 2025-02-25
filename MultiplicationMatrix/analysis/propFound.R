library(dplyr)


rm(list = ls())
dir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# List all files that match the pattern "dataNeigh*_new.txt" in the "../" directory
files <- list.files(path = "data/dataNeigh/", pattern = "^dataNeigh\\d+_new_KC.txt$")

# Extract the numeric part from each filename
numbers <- gsub("dataNeigh(\\d+)_new_KC.txt", "\\1", files)

# Convert the extracted numeric parts to a numeric vector
alldata=data.frame()
for (binary in numbers) {
if(file.exists(paste("data/inout3/in_out_count_",binary,"_KC_new2.txt",sep=""))){
allfound=read.table(paste("data/inout3/in_out_count_",binary,"_KC_new2.txt",sep=""),header = T)

bmbFound=read.table(paste("data/dataNeigh/dataNeigh",binary,"_new_KC.txt",sep=""),header = T)
pacom=bmbFound$KC_col1[1]
if(nrow(bmbFound)==0)next
#found
df1n <- allfound %>%
  filter(phenotype %in% bmbFound$MutantPhenotype)

#not found
df2n <- allfound %>%
  filter(!(phenotype %in% bmbFound$MutantPhenotype))


foundckc=median(df1n$Conditional_Complexity)
notfoundckc=median(df2n$Conditional_Complexity)

data.tmp=data.frame(phenotype=binary,
                    complexity=pacom,
                    foundckc,notfoundckc)
alldata=rbind(alldata,data.tmp)
}}



pp=alldata$foundckc-alldata$notfoundckc
pp=pp/(max(abs(pp)))

ppw=pp ; ppw=as.data.frame(ppw)
ppw$exp="MMP"
write.table(ppw,"data/MM_FoN.dat")

pdf("mmatrix_norm_foundOrNot.pdf")
par(cex.axis = 2.2, cex.lab = 2.2)
par(mar = c(5, 7, 4, 2))  # c(bottom, left, top, right)
boxplot(pp,ylab=expression(paste(~tilde(K), "(y|x)"[found], " - ",~tilde(K), "(y|x)"[not_found])))
#points( jitter(rep(1,length(pp)),5), pp)
length(pp)
dev.off()

