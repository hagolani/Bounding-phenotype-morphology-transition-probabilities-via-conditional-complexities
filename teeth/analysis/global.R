library(dplyr)
rm(list = ls())

dir=dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(dir)
cex0=2.5 ; cex1=2.2
cex_axis=2.2 ; cex_lab=2.2 ; cex_main=2.2 ; cex_leyend=1.25
cex2=2.5
df0=read.table("data/no3_varied_random_search_filtered1.txt",header = T)

df0=df0[df0$Cusps>0,]


opcORcusps=3 # 2 cusps , 3 opc

df1 <- as.data.frame(table(df0[[opcORcusps]]))
colnames(df1) = c("OPC","n")
df1$OPC <- as.numeric(as.character(df1$OPC))

df1$freq = log10(df1$n/sum(df1$n))

df1=df1[which(df1$freq>-3.5),]

maxis <- df1 %>%
  mutate(Group = ceiling(row_number() / 5)) %>%  # Create groups of 10 rows based on their row number
  group_by(Group) %>%  # Group by these sets of 10 rows
  summarise(
    Max_OPC = max(OPC),       # Get the maximum OPC value in each group
    Max_freq = max(freq)      # Get the maximum freq value in each group
  )

maxis20=maxis[which(maxis$Max_freq>-4),]



#lm1=lm(maxis20$Max_freq ~ maxis20$Max_OPC)
lm1=lm(df1$freq ~ df1$OPC)

color1="blue"#rgb(0, 0, 0, 0.5)
if(opcORcusps==2)a="cusps"
if(opcORcusps==3)a="opc"
pdf(paste("probCom_",a,".pdf",sep=""))
par(mar = c(5, 6, 4, 2) + 0.1)  # Adjusting margin
par(cex.axis = cex_axis, cex.lab = cex_lab , cex.main=cex_lab)
#plot(x=jitter(df1$Z),y=df1$final_prob,
df1=df1[df1$freq>-4,]
ymin=min(df1$freq)
plot(x=df1$OPC,y=df1$freq,
     cex=2.5,pch=16,col="blue",ylim = c(min(df1$freq),0),
     xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")


abline(lm1,lwd=2)
axis(1,cex=cex0,line=0)
axis(2,cex=cex0,line=0)
mtext(expression(paste(~tilde(K), "(x)")), side = 1, line = 4,cex=cex0)  # Increase the 'line' value to move the label further down
mtext(expression(paste(log[10],"P(x)")), side = 2, line = 3.5,cex=cex0)  # Increase the 'line' value to move the label further down

abline(lm1,lwd=2)

dev.off()
