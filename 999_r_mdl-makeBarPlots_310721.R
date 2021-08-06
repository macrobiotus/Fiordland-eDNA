rm(list=ls())
library(gridExtra)
setwd("C:/Users/Michel/Documents/fish/rscripts")
f <- readRDS(file="../data/31July2021/210703_998_r_summarize_results__data_gtestimate_accumulation_curves.Rds")
f <- data.frame(f)

dim(f)
class(f)
#f[which(f$SPECIES==""),]
#f <- f[-which(f$SPECIES==""),]

table(f$RESERVE.GROUP,f$INSIDE.RESERVE,f$SAMPLE.TYPE)
##
## 3 rows with no species. class Actinopteri (which includes salmon, trout)
##
dim(f)

level <- "SPECIES"  # One or the other
#level <- "GENUS"  

## Create barplots, which show whether we are 
##     seeing new species by adding samples?
##
b <- vector()
e <- vector()
sets <- c(1,3,4,5,7,8,9,10,11,12,17,18,19,21,22,23,24,26,27,28,29)
counter <- 0
for (i in sets) {
  counter <- counter + 1
  b[counter] <- length(unique(f[which(f$SET.ID %in% sets[1:i] & f$SAMPLE.TYPE=='BRUV'),level]))
  e[counter] <- length(unique(f[which(f$SET.ID %in% sets[1:i] & f$SAMPLE.TYPE=='eDNA'),level]))
  
}
df2 <- data.frame(bruv=b,edna=e,set=sets)
df <- df2[order(df2$set),]


library(ggplot2)

word <- 'species'
if (level == 'GENUS') { 
  word <- 'genera'
}

p1 <- ggplot(df,aes(x=seq(1:length(sets)),y=bruv)) + 
  geom_bar(stat='identity',fill = "deepskyblue4") +
  xlab('Sets') +
  ylab(paste('Total number of', word, 'discovered')) +
  ylim(0,50) +
  ggtitle('BRUV')
p2 <- ggplot(df,aes(x=seq(1:length(sets)),y=edna)) + 
  geom_bar(stat='identity',fill = "deepskyblue4") +
  xlab('Sets') +
  ylab('') +
  ylim(0,50) +
  ggtitle('e-DNA')

# library(gridExtra)
p2 <- grid.arrange(p1,p2,nrow=1)
ggsave(p2, file=paste0("../myoutput/barplots3107_",level,".pdf"), 
       device="pdf",dpi=800,
       width=200,height=100,units='mm')

