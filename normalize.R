library(sva)
library(limma)
merge3<-as.matrix(merge3)
rownames(merge3)<-merge3[,1]
exp=merge3[,2:ncol(merge3)]
dimname<-list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimname)
batchType<-c(rep(1,32),rep(2,82),rep(3,59),rep(4,65))
modType<-c(rep("normal",32),rep('tumor',0),rep("tumor",70),rep("normal",12),rep('tumor',35),rep("normal",24),rep("tumor",27),rep("normal",38))
mod=model.matrix(~as.factor(modType))
outTab=ComBat(data,batchType,mod,par.prior = TRUE)
