library(dplyr)


rm(list = ls())
dir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# List all files that match the pattern "dataNeigh*_new.txt" in the "../" directory
files <- list.files(path = "data/dataNeig_KC/", pattern = "^dataNeigh\\d+_new_KC.txt$")

# Extract the numeric part from each filename
numbers <- gsub("dataNeigh(\\d+)_new_KC.txt", "\\1", files)

# Convert the extracted numeric parts to a numeric vector
#numbers_vector <- as.numeric(numbers)
alldata=data.frame()
for (binary in numbers) {
  
allfound=read.table(paste("data/dataNeig_KC/dataToAll",binary,"_new_KC.txt",sep=""),header = T)
pacom=allfound[which(allfound$Phenotype==as.numeric(binary)),3]

bmbFound=read.table(paste("data/dataNeig_KC/dataNeigh",binary,"_new_KC.txt",sep=""),header = T)
if(nrow(bmbFound)==0)next
#found
df1n <- allfound %>%
  filter(Phenotype %in% bmbFound$MutantPhenotype)

#not found
df2n <- allfound %>%
  filter(!(Phenotype %in% bmbFound$MutantPhenotype))


foundckc=median(df1n$Conditional_Complexity)
notfoundckc=median(df2n$Conditional_Complexity)

data.tmp=data.frame(phenotype=binary,
                    complexity=pacom,
                    foundckc,notfoundckc)
alldata=rbind(alldata,data.tmp)
}

pp=alldata$foundckc-alldata$notfoundckc

pp=pp/(max(abs(pp)))

ppw=pp ; ppw=as.data.frame(ppw)
ppw$exp="circadian"
write.table(ppw,"data/circadian_FoN.dat")

pdf("circadian_foundOrNot.pdf")
par(cex.axis = 2.2, cex.lab = 2.2)
par(mar = c(5, 7, 4, 2))  # c(bottom, left, top, right)
#boxplot(pp,ylab=expression(paste(tilde(K), "(x)"[found], " - ",tilde(K), "(x)"[not_found])))
boxplot(pp,ylab=expression(paste(~tilde(K), "(y|x)"[found], " - ",~tilde(K), "(y|x)"[not_found])))


#points( jitter(rep(1,length(pp)),5), pp)
# Define color thresholds
low_threshold <- quantile(alldata$complexity, 0.33)
mid_threshold <- quantile(alldata$complexity, 0.66)

# Assign colors based on value
#point_colors <- ifelse(alldata$complexity < low_threshold, "blue",
#                       ifelse(alldata$complexity < mid_threshold, "green", "red"))

# Add jittered points to the boxplot with color-coding
#points(jitter(rep(1, length(pp)), amount = 0.09), pp, col = "black", pch = 16)

# Add a legend to the plot
#legend("topright", legend = c("Low", "Medium", "High"),title="Complexity", 
#       col = c("blue", "green", "red"), pch = 16)

dev.off()

#plot(alldata$complexity,alldata$foundckc-alldata$notfoundckc)



