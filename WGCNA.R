setwd("C:\\Users\\dell\\Downloads\\geo\\txt")
rt=read.table("normalize.txt",sep='\t',header = T,check.names = F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimname=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimname)
rt <- as.data.frame(t(rt))
#cluster
sampleTree = hclust(dist(rt), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
# 设定软阈值范围
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 获得各个阈值下的 R方 和平均连接度
sft = pickSoftThreshold(rt, powerVector = powers, verbose = 5)
# 作图：
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# 获得临近矩阵：
softPower <- sft$powerEstimate
adjacency = adjacency(rt, power = softPower);
# 将临近矩阵转为 Tom 矩阵
TOM = TOMsimilarity(adjacency);
# 计算基因之间的相异度
dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average");

# 基因多的时候使用下面的代码：
k <- softConnectivity(datE=rt,power=softPower) 
sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")


# 使用相异度来聚类为gene tree(聚类树)：
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
windows()
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
# 使用动态剪切树挖掘模块：
minModuleSize = 30;
# 动态切割树:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# 拓扑热图：
nSelect = 400 
# For reproducibility, we set the random seed 
set.seed(10); 
select = sample(13487, size = nSelect)
selectTOM = dissTOM[select, select]
dynamicColors=labels2colors(dynamicMods)
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster. 
selectTree = hclust(as.dist(selectTOM), method = "average") 
selectColors = dynamicColors[select]; 
# Open a graphical window 
sizeGrWindow(9,9) 
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot 
plotDiss = selectTOM^softPower; 
diag(plotDiss) = NA; 
TOMplot(plotDiss, 
        selectTree, 
        selectColors, 
        main = "Network heatmap plot, selected genes") 

MEList = moduleEigengenes(rt, colors = dynamicColors)
MEs = MEList$eigengenes
# 计算根据模块特征向量基因计算模块相异度：
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90) 


MEDissThres = 0.2
# 在聚类图中画出剪切线
abline(h=MEDissThres, col = "red")
# 合并模块：
merge_modules = mergeCloseModules(rt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# 合并后的颜色：
mergedColors = merge_modules$colors;
# 新模块的特征向量基因：
mergedMEs = merge_modules$newMEs;
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,cex.lab=1.2,cex.colorLabels = 1.2)


cor_ADR <- signif(WGCNA::cor(dataTrait,mergedMEs,use="p",method="pearson"),5)
p.values <- corPvalueStudent(cor_ADR,nSamples=nrow(dataTrait))
Freq_MS_max_cor <- which.max(abs(cor_ADR["condition",]))
Freq_MS_max_p <- which.min(p.values["condition",])

GS1 <- as.numeric(WGCNA::cor(dataTrait[,1],rt,use="p",method="pearson"))
# 显著性是绝对值：
GeneSignificance <- abs(GS1)
# 获得该性状在每个模块中的显著性：
ModuleSignificance <- tapply(GeneSignificance,mergedColors,mean,na.rm=T)

ADJ1=abs(cor(rt,use="p"))^softPower 
# 根据上面结果和基因所属模块信息获得连接度：
# 整体连接度 kTotal，模块内部连接度：kWithin，kOut=kTotal-kWithin， kDiff=kIn-kOut=2*kIN-kTotal 
Alldegrees1=intramodularConnectivity(ADJ1, colorh1) 

# 注意模块内基于特征向量基因连接度评估模块内其他基因： de ne a module eigengene-based connectivity measure for each gene as the correlation between a the gene expression and the module eigengene
# 如 brown 模块内：kM Ebrown(i) = cor(xi, MEbrown) ， xi is the gene expression pro le of gene i and M Ebrown is the module eigengene of the brown module
# 而 module membership 与内部连接度不同。MM 衡量了基因在全局网络中的位置。
datKME=signedKME(rt, datME, outputColumnName="MM.")
plotDendroAndColors(geneTree, dynamicColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
