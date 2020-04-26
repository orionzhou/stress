source('functions.R')
require(edgeR)
require(WGCNA)
require(pheatmap)
require(DiPALM)
dirw = file.path(dird, '14_diparm')

#{{{ read in data & filtering
data("exampleData")
cnts<-exampleData$rawCounts
colnms<-colnames(cnts)
colnms<-gsub("S1","R1_",colnms)
colnms<-gsub("S2","R2_",colnms)
colnms<-gsub("D","Drought_",colnms)
colnms<-gsub("W","Watered_",colnms)
colnames(cnts)<-colnms

safMat<-exampleData$geneAnnotations
geneLen<-setNames(abs(safMat$End-safMat$Start),nm = safMat$GeneID)
cntsDge<- DGEList(counts = cnts)
cntsDge<- calcNormFactors(cntsDge)
DroughtLog<- rpkm(cntsDge, log=T, gene.length = geneLen, prior.count=0)

geneMeans<-apply(DroughtLog,1,function(x) mean(x,na.rm = T))
hist(geneMeans,col="skyblue", breaks = seq(-12,16,0.25))
abline(v=0,col="red",lwd=3,lty=2)
dev.off()

Droughtfiltered<-DroughtLog[which(geneMeans>0),]
#
minVal<-min(Droughtfiltered[!is.na(Droughtfiltered)])-1
Droughtfiltered[is.na(Droughtfiltered)]<-minVal

tmp<-Droughtfiltered
spNms<-strsplit(x = colnames(Droughtfiltered), split = "_")
tnms<-sapply(spNms,function(x) paste(x[c(2,1)],collapse = "."))
TCsDrought<-tapply(1:length(tnms),INDEX = tnms, function(x) tmp[,x])

varFiltered<-lapply(TCsDrought,function(x) apply(x,1,function(y) var(y)>0))
varFiltered<-do.call(cbind,varFiltered)
varFiltered<-apply(varFiltered,1,function(x) all(x))
varFiltered<-names(varFiltered)[which(varFiltered)]
TCsDrought<-lapply(TCsDrought,function(x) x[varFiltered,])
TCsDrought<-lapply(TCsDrought,function(x) x[,order(as.numeric(sapply(strsplit(colnames(x),split = "_"),function(y) y[3])))])
sapply(TCsDrought,colnames)
TCsAll<-do.call(rbind,TCsDrought)
#}}}

#{{{ WGCNA
BlockModsDrought<- blockwiseModules(datExpr = t(TCsAll), power = 10, networkType = "signed", corType="bicor", TOMType="signed", minModuleSize=100, mergeCutHeight=0.2, deepSplit=1, pamRespectsDendro = F, nThreads = 4, verbose=3)
save(BlockModsDrought, file="BlockModsDrought.RData")
#}}}

load("BlockModsDrought.RData")

MEs<-BlockModsDrought[[3]]
kMEsList<-BuildModMembership(MeMat = MEs, TCsLst = TCsDrought)
Med<-sapply(TCsDrought,function(x) apply(x,1,function(y) median(y,na.rm = T)))
TCsDroughtPerm<-lapply(TCsDrought,function(x) x[sample(1:nrow(x),nrow(x),replace = T),])
kMEsPerm<-BuildModMembership(MeMat = MEs, TCsLst = TCsDroughtPerm)
MedPerm<-sapply(TCsDroughtPerm,function(x) apply(x,1,function(y) median(y,na.rm = T)))

Treat<-as.factor(c("Drought","Drought","Watered","Watered"))
design<-model.matrix(~0+Treat)
colnames(design)<-levels(Treat)
contr<-"Drought-Watered"

LimmaModskMEs<-lapply(kMEsList, function(x) BuildLimmaLM(dataMat = x, designMat = design, contrastStr = contr))
LimmaModsMed<-BuildLimmaLM(dataMat = Med, designMat = design, contrastStr = contr)
LimmaModskMEs<-do.call(cbind,lapply(LimmaModskMEs,function(x) x$t))
LimmaModsMed<-LimmaModsMed$t
gc()
LimmaModskMEsPerm<-lapply(kMEsPerm, function(x) BuildLimmaLM(dataMat = x, designMat = design, contrastStr = contr))
LimmaModsMedPerm<-BuildLimmaLM(dataMat = MedPerm, designMat = design, contrastStr = contr)
LimmaModskMEsPerm<-do.call(cbind,lapply(LimmaModskMEsPerm,function(x) x$t))
LimmaModsMedPerm<-LimmaModsMedPerm$t
gc()

TestSumskMEs<-apply(LimmaModskMEs,1, function(x) sum(abs(x),na.rm = T))
TestSumsMed<-abs(LimmaModsMed[,1])

PermSumskMEs<-apply(LimmaModskMEsPerm,1, function(x) sum(abs(x),na.rm = T))
PermSumsMed<-abs(LimmaModsMedPerm[,1])

ggPlotMultiDensities(denslist = list(Test=TestSumskMEs,Permuted=PermSumskMEs), main = "Pattern Change Scores", xlab = "Differential Pattern Score",lwidth = 1)
dev.off()
ggPlotMultiDensities(denslist = list(Test=TestSumsMed,Permuted=PermSumsMed), main = "Expression Change Scores", xlab = "Differential Expression Score",lwidth = 1)
dev.off()

AdjkMEs<-sapply(TestSumskMEs,function(x) AdjustPvalue(tVal = x, tVec = TestSumskMEs, pVec = PermSumskMEs))
AdjMed<-sapply(TestSumsMed,function(x) AdjustPvalue(tVal = x, tVec = TestSumsMed, pVec = PermSumsMed))
#
SigkMEs<-AdjkMEs[which(AdjkMEs<0.01)]
SigMed<-AdjMed[which(AdjMed<0.01)]

topgene<-"BraA07g20790R"
PlotTCs(TClst = TCsDrought,tgene = topgene, scale = T, xlab="ZT Time (hours)", xAxsLabs = c(seq(1,23,4),seq(1,23,4)), ledgeX="topleft",tcols = c("red","red","blue","blue"), tltys = c(1,2,1,2))
dev.off()

LimmaModskMEsSig<-LimmaModskMEs[names(SigkMEs),]
patternCor<-cor(t(LimmaModskMEsSig))
patternTree<-hclust(as.dist(1-patternCor),method = "complete")

expressionMat<-do.call(cbind,TCsDrought)
eMatCols<-colnames(expressionMat)
eMatCols<-gsub("^R[[:digit:]]\\_","",eMatCols)

expressionAvg<-tapply(colnames(expressionMat),INDEX = eMatCols, function(x) rowSums(as.data.frame(expressionMat[,x]),na.rm = T)/length(x))

orderVec<-strsplit(names(expressionAvg),split = "_")
orderVec<-lapply(1:2,function(x) sapply(orderVec,function(y) y[x]))
expressionAvg<-do.call(cbind,expressionAvg[order(orderVec[[1]],as.numeric(orderVec[[2]]))])

colFunc<-colorRampPalette(colors = c("darkblue","blue","lightblue","white","orange"))

pheatmap(mat = expressionAvg[names(SigkMEs),], cluster_rows = patternTree, cluster_cols = F,scale = "row", color = colFunc(25), gaps_col = 12, show_rownames = F)
dev.off()

patternClusters<-cutreeDynamic(dendro = patternTree, minClusterSize = 100, distM = 1-patternCor, deepSplit = 1)
names(patternClusters)<-patternTree$labels
table(patternClusters)
patternClusters<-tapply(X = names(patternClusters), INDEX = patternClusters, function(x) x)

clustScores<-sapply(patternClusters,function(x) mean(TestSumskMEs[x]))
topClust<-which.max(clustScores)
PlotTCsRibbon(TClst = TCsDrought, tgenes = patternClusters[[topClust]], main="Pattern-Change Cluster",xAxsLabs = c(seq(1,23,4),seq(1,23,4)), xlab="ZT Time (hours)", scale = T, tcols = c("red","red","blue","blue"), tltys = c(1,2,1,2))
dev.off()

ExpMatrixMedSig<-expressionMat[names(SigMed),]
expCor<-cor(t(ExpMatrixMedSig))
expTree<-hclust(as.dist(1-expCor),method = "complete")

colFunc<-colorRampPalette(colors = c("darkblue","white","orange"))
pheatmap(mat = expressionAvg[names(SigMed),], cluster_rows = expTree, cluster_cols = F,scale = "row", color = colFunc(25), gaps_col = 12, show_rownames = F)
expClusters<-cutreeDynamic(dendro = expTree, minClusterSize = 5, distM = 1-expCor, deepSplit = 1)
names(expClusters)<-expTree$labels
table(expClusters)
expClusters<-tapply(X = names(expClusters), INDEX = expClusters, function(x) x)
clustScores<-sapply(expClusters,function(x) mean(TestSumsMed[x]))
topClust<-which.max(clustScores)
PlotTCsRibbon(TClst = TCsDrought, tgenes = expClusters[[topClust]], main="Expression-Change Cluster",xAxsLabs = c(seq(1,23,4),seq(1,23,4)), xlab="ZT Time (hours)", scale = T, tcols = c("red","red","blue","blue"), tltys = c(1,2,1,2))
dev.off()

