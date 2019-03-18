library(dplyr)
args<-commandArgs(TRUE)
blatnames<-c("matches","misMatches","repMatches","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd","tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStarts")

data<-read.table(args[1],header=FALSE)

colnames(data)<-blatnames
data$qName<-as.character(data$qName)
data$tName<-as.character(data$tName)

data$labelmin<-apply(as.data.frame(cbind(data$qName,data$tName)), 1, FUN=min)
data$labelmax<-apply(as.data.frame(cbind(data$qName,data$tName)), 1, FUN=max)
data<-mutate(data,combo = paste(labelmin,labelmax,sep="_"))
data<-distinct(data,combo,.keep_all=TRUE)
data$labelmin<-NULL
data$labelmax<-NULL
data$combo<-NULL
write.table(data,args[2],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
