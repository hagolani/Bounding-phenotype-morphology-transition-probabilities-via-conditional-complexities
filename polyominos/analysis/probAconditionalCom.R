library(dplyr)
library(MASS)
library(gtools) 
library(matrixStats)

rm(list = ls())
cex_axis=2.2 ; cex_lab=2.2 ; cex_main=2.2 ; cex_leyend=1.25
cex2=2.5
dir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

allindis <- list.files(path="data",pattern = "^phenotype")

numeric_part <- as.numeric(gsub("\\D", "", allindis))

# Order the filenames based on the numeric part
sorted_allindis <- allindis[order(numeric_part)]
allindis=sorted_allindis


cc=1
alldata=data.frame()
for (indi in allindis) {
   df.tmp=read.table(paste("data/",indi,sep=""),header = T)
   df.tmp$parent=cc ; cc=cc+2
   alldata=rbind(alldata,df.tmp)
}

df=alldata

colnames(df) = c("ConditionalComplexity","Freq","parent")


genotypes = read.table("data/Polyomino/s_28_list.txt")

genotypes_N= as.data.frame(table(genotypes))


true_indi=0
allindis=unique(df$parent)
save.stas=data.frame()
for (indi in allindis) {
    true_indi=true_indi+1
    df1=df[which(df$parent==indi),]
    df1=df1[order(df1$ConditionalComplexity),]
    
    NPheno=length(unique(df1$ConditionalComplexity))
    
    minC= min(df1$ConditionalComplexity) ; maxC= max(df1$ConditionalComplexity)
    upper_bound=df1
    upper_bound$Z2=log2(NPheno)*((df1$ConditionalComplexity-minC)/(maxC-minC))
    
    maxis <- df1 %>%
      group_by(ConditionalComplexity) %>%
      summarise(Max_prob = max(Freq))
    maxis20=maxis[which(maxis$Max_prob>-7),]
    
    # Extract x and y
    x <- maxis20$ConditionalComplexity
    y <- maxis20$Max_prob
    
    a=-1 ; b=0
    #upper_bound=df1
    upper_bound$Max_prob2= log10(2**(a*upper_bound$Z2+b))
    
    maxis21=data.frame(maxis20[which(maxis20$Max_prob==max(maxis20$Max_prob)),])
    maxis22=data.frame(maxis20[which(maxis20$Max_prob==min(maxis20$Max_prob)),])
    maxis23=rbind(maxis21,maxis22)
    
    model <- lm(Max_prob ~ ConditionalComplexity, data = maxis20)
    model1 <- lm(Max_prob2 ~ ConditionalComplexity, data = upper_bound)
    color1=rgb(0, 0, 0, 0.75)
    
    pearson.r2 = (cor(as.numeric(maxis20$ConditionalComplexity),maxis20$Max_prob))**2
    ######################    test bootstrap
    x1=maxis20$ConditionalComplexity ; y1=maxis20$Max_prob
    x2=upper_bound$ConditionalComplexity ; y2=upper_bound$Max_prob2
    
    dfmx=data.frame(x=x1,y=y1)
    dfbound=data.frame(x=x2,y=y2)
    
    fit1 <- glm(y ~ x, data=dfmx)
    ndf1 <- data.frame(x=seq(min(dfmx$x), max(dfmx$x), length.out=1e3))
    pred1 <- predict(fit1, newdata=ndf1, type='link')
    
    fit2 <- glm(y ~ x, data=dfbound)
    ndf2 <- data.frame(x=seq(min(dfbound$x), max(dfbound$x), length.out=1e3))
    pred2 <- predict(fit2, newdata=ndf2, type='link')
    
    n_boot=999
    set.seed(42)
    bf <- replicate(
      999L, {
        bdf <- dfmx[sample.int(nrow(dfmx), replace=TRUE), ]
        glm(y ~ x, data=bdf)
      },
      simplify=FALSE
    )
    
  
    bf2 <- replicate(
      999L, {
        bdf <- dfbound[sample.int(nrow(dfbound), replace=TRUE), ]
        glm(y ~ x, data=bdf)
      },
      simplify=FALSE
    )
    
    bpred <- sapply(bf, predict, newdata=ndf1, type='link')
    ci <- \(x, sd) x + as.matrix(sd*(-qt(.025, Inf))) %*% cbind(-1, 1)
    bpredci <- ci(matrixStats::rowMeans2(bpred), matrixStats::rowSds(bpred))
    
    
    
    
    
    boot_fits <- matrix(NA, nrow = n_boot, ncol = length(x))
    slope_diff <- numeric(n_boot)
    # Extract slopes and calculate the difference
    for (i in 1:n_boot) {
      boot_fits[i, ] <- predict(bf[[i]], newdata = data.frame(x = x))
      slope_diff[i] <- bf[[i]]$coefficients[2] - bf2[[i]]$coefficients[2]
    }
    
    # Calculate confidence interval
    slope_diff=na.omit(slope_diff)
    ci <- quantile(slope_diff, probs = c(0.025, 0.975))
    
    ci_lower <- apply(boot_fits, 2, quantile, probs = 0.025)
    ci_upper <- apply(boot_fits, 2, quantile, probs = 0.975)
    
    
    # Print results
    cat("Bootstrap 95% confidence interval for slope difference:", ci, "\n")
    
    # Test if 0 is within the confidence interval
    if (ci[1] > 0 || ci[2] < 0) {
      print("The slopes are significantly different.")
      a="boot: sig different"
      asv="sigdifferent"
    } else {
      print("The slopes are not significantly different.")
      a="boot: NOT sig different"
      asv="notsigdifferent"
    }
    
    
    
    
    
    #################
    
    
    
    
    pdf(paste(indi,"_c_shade_CI_v2.pdf",sep=""))
    par(mar = c(5, 6, 4, 2) + 0.1)  # Adjusting margin
    par(cex.axis = cex_axis, cex.lab = cex_lab , cex_main=cex_lab)
    
    plot(x = df1$ConditionalComplexity, y = df1$Freq,
          #main = paste("complexity", indi),
         xaxt = 'n', yaxt = 'n', xlab = "", ylab = "",
         cex = cex2, pch = 16, col = "blue", ylim = c(min(df1$Freq), 0))
          
    axis(1,cex=cex2,line=0)
    axis(2,cex=cex2,line=0)
    mtext(expression(paste(~tilde(K), "(y|x)")), side = 1, line = 4,cex=cex2)  # Increase the 'line' value to move the label further down
    mtext(expression(paste("log"[10], " P(x" %->% "y)")), side = 2, line = 3.5,cex=cex2)  # Increase the 'line' value to move the label further down
    
     # Add points for dfmx
     points(dfmx, col = "black", pch = 16, cex = 1.5)
     # Add the main fitted line
     lines(ndf1$x, fit1$family$linkinv(pred1), col = "black", lwd = 3)
     # Create a shaded area between the two prediction interval lines
     polygon(c(ndf1$x, rev(ndf1$x)),
             c(fit1$family$linkinv(bpredci[, 1]), rev(fit1$family$linkinv(bpredci[, 2]))),
             col = rgb(0, 0, 0, alpha = 0.1), border = NA) # Adjust the alpha level for shading
     # Add the regression line
     abline(model1, col = "red", lwd = 3)
     
     legend("topright",legend = c("Fitted model","Fitted 95% CI","Bound model"),
            lwd = c(2.5, 2.5,2.5),lty=c(1,1,1),cex=cex_leyend,        # Point symbols like the plot
            col = c("black","gray", "red")) # Colors for the symbols
     
     
     
        
      
      corr=cor.test(df1$Freq,df1$ConditionalComplexity,method="spearman")
      #prediction two, linear
      summary_model=summary(model)
      f_statistic <- summary_model$fstatistic[1]
      model_p_value <- pf(f_statistic, summary_model$fstatistic[2], summary_model$fstatistic[3], lower.tail = FALSE)
      #prediction three, slope
        #model  maxis20
        #model1 upper_bound
     
      group1 <- rep("Fitted", nrow(maxis20))
      group2 <- rep("Bound", nrow(upper_bound))
      
      x <- c(maxis20$ConditionalComplexity, upper_bound$ConditionalComplexity)
      y <- c(maxis20$Max_prob, upper_bound$Max_prob2)
      group <- factor(c(group1, group2))
      model3 <- lm(y ~ x * group)
      summary_model3 <- summary(model3)
      anova_model <- anova(model3)
      print(anova_model)
      # If the interaction term (x:group) is significant, the slopes are significantly different
      interaction_p_value <- summary_model3$coefficients[4, 4] # p-value of the interaction term
      
      Ncount = genotypes_N[genotypes_N$V1==true_indi,2]
      
      tmp=data.frame(spearman=corr$estimate,pvalcorr=corr$p.value,
                     pvallinear=model_p_value,pvalanova=interaction_p_value,
                     bootstrap=asv,pearson=pearson.r2,Ncount)
      save.stas=rbind(save.stas,tmp)
      
      
      
    dev.off()
    
    
    
    
    
    
    
    
}


write.table(save.stas,"data/stats_v2.dat")





save.stas=read.table("data/stats_v2.dat")
par(mfrow=c(1,1))



pdf("stats_circadian_corr_v2.pdf")
par(mar = c(5, 6, 4, 2) + 0.1)  # Adjusting margin
par(cex.axis = 2.2, cex.lab = 2.2, cex.main = 1)

sigSpear= save.stas[which(save.stas$pvalcorr<0.05 & save.stas$spearman<0 ),1]
prop1= length(sigSpear)/nrow(save.stas)
boxplot(sigSpear,
        ylim=c(-1,1),
        main=paste("Level (I) Spearman correlations. Phenotypes \n",
                   "mean=",round(mean(sigSpear),3),"median=",round(median(sigSpear),3),"var",round(var(sigSpear),3),"\n",
                   "proportion success =",round(prop1,3)),
        xlab="",ylab="Correlations")
newvals=rep(sigSpear, times = save.stas[which(save.stas$pvalcorr<0.05 & save.stas$spearman<0 ),7])
prop2= sum(save.stas[which(save.stas$pvalcorr<0.05 & save.stas$spearman<0),7])/sum(save.stas$Ncount)
boxplot(newvals,
        ylim=c(-1,1),
        main=paste("Level (I) Spearman correlations. Genotypes \n",
                   "mean=",round(mean(newvals),3),"median=",round(median(newvals),3),"var",round(var(newvals),3),"\n",
                   "proportion success =",round(prop2,3)),
        xlab="",ylab="Correlations")
abline(0,0)
dev.off()


pdf("stats_circadian_pear_r2_v2.pdf")
par(mar = c(5, 6, 4, 2) + 0.1)  # Adjusting margin
par(cex.axis = 2.2, cex.lab = 2.2, cex.main = 1)
sigPear=save.stas[which(save.stas$pvalcorr<0.05 & save.stas$pearson>0.5 ),6]
prop1= length(sigPear)/nrow(save.stas)
boxplot(sigPear,
        ylim=c(0,1),
        main=paste("Level (II) \n Pearson correlations. Phenotypes \n",
                   "mean=",round(mean(sigPear),3),"median=",round(median(sigPear),3),"var",round(var(sigPear),3),"\n",
                   "proportion success =",round(prop1,3)),
        xlab="",ylab=expression(R^2))

newvals=rep(sigPear, times = save.stas[which(save.stas$pvalcorr<0.05 & save.stas$pearson>0.5 ),7])
prop2= sum(save.stas[which(save.stas$pvalcorr<0.05 & save.stas$pearson>0.5),7])/sum(save.stas$Ncount)
boxplot(newvals,
        ylim=c(0,1),
        main=paste("Level (II) Pearson correlations. Genotypes \n",
                   "mean=",round(mean(newvals),3),"median=",round(median(newvals),3),"var",round(var(newvals),3),"\n",
                   "proportion success =",round(prop2,3)),
        xlab="",ylab="Correlations")
#abline(0,0)
#abline(0,0)
dev.off()


pdf("stats_circadian_lm_v2.pdf")
par(mar = c(5, 6, 4, 2) + 0.1)  # Adjusting margin
par(cex.axis = 2.2, cex.lab = 2.2, cex.main = 1)
boxplot(save.stas$pvallinear,
        main=paste("Level (II) \n Significant linear regressions. Phenotypes \n",
                   "mean=",round(mean(save.stas$pvallinear),3),"median=",round(median(save.stas$pvallinear),3),"var",round(var(save.stas$pvallinear),3)),
        xlab="",ylab="p-values")
abline(0.05,0)

newvals=rep(save.stas$pvallinear, times = save.stas$Ncount)
boxplot(newvals,
        ylim=c(0,1),
        main=paste("Level (II) Significant linear regressions. Genotypes \n",
                   "mean=",round(mean(newvals),3),"median=",round(median(newvals),3),"var",round(var(newvals),3)),
        xlab="",ylab="Correlations")
abline(0.05,0)

dev.off()


for (iji in 1:nrow(save.stas)) {
  if(save.stas$pearson[iji]<0.5 & save.stas$pvallinear[iji]>0.05){
    save.stas$bootstrap[iji]="sigdifferent"
  }
}



not_npredicted=nrow(save.stas[which(save.stas$bootstrap=="sigdifferent"),])
npredicted=nrow(save.stas[which(save.stas$bootstrap=="notsigdifferent"),])
total=sum(not_npredicted,npredicted)

pdf("stats_circadian_slope_v2.pdf")
par(mar = c(5, 6, 4, 2) + 0.1)  # Adjusting margin
par(cex.axis = 2.2, cex.lab = 2.2, cex.main = 1)
barplot(c(not_npredicted/total,npredicted/total),
        main=paste("Level (III) \n Predicted slopes. Phenotypes \n",
                   "predicted=",round(npredicted/total,3),"not_npredicted=",round(not_npredicted/total,3)),
        ylab="proportion of phenotypes",
        names.arg = c("Not predicted","Predicted"))

#not_npredicted=nrow(save.stas[which(save.stas$bootstrap=="sigdifferent"),])
not_npredicted=sum(save.stas[which(save.stas$bootstrap=="sigdifferent"),7])

#npredicted=nrow(save.stas[which(save.stas$bootstrap=="notsigdifferent"),])
npredicted=sum(save.stas[which(save.stas$bootstrap=="notsigdifferent"),7])

total=sum(not_npredicted,npredicted)

barplot(c(not_npredicted/total,npredicted/total),
        main=paste("Level (III) \n Predicted slopes. Genotypes \n",
                   "predicted=",round(npredicted/total,3),"not_npredicted=",round(not_npredicted/total,3)),
        ylab="proportion of phenotypes",
        names.arg = c("Not predicted","Predicted"))


dev.off()
