##Set directory here
setwd("")
### Loads packages and my functions
source("BayesMetaFunctions.R")
dir <- getwd()

######################
#get the Tree, fix names fnalysis
Tree <- read.nexus(file.path(dir, "NewTree.nex"))
Tree <- consensus.edges(Tree)
Tree <- force.ultrametric(Tree)
Tree$tip.label <- gsub("_", " ", Tree$tip.label)
TreeAinv<-inverseA(Tree, nodes="TIPS")$Ainv

#get the Field data, flip correlates on time so that large cor = success,
#move hatch to laying b/c only one hatch sample post-processing;
#variances are not appropriately sampled for N=1
RawData <- read.csv(file.path(dir, "new data lineup.csv"), header = TRUE)
RawData[,"animal"] <- RawData$Species
flipsign <- which(RawData$MClass == "Time")
RawData$Raw.Cor[flipsign] <- RawData$Raw.Cor[flipsign]*-1
RawData$MType[which(RawData$MType == "Hatch")] <- "Laying"
RawData$MType <- droplevels(RawData$MType)

#set up the larger data sets
FullDataSet <- MakeSet(RawData)
p.varfull <- var(FullDataSet$Zs)

FunnelPlots(FullDataSet, Groups = "Rep",1)
FunnelPlots(FullDataSet, Groups = "OC",2)

AsyncTests(FullDataSet)
names(FullDataSet)

set.seed(49)
model10.1 <- BayesMeta(Data=FullDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.varfull, TreeA=TreeAinv)
ConvergePlot(model10.1, "FullPop-FullData", Type = "Pop")
print(xtable(VarianceTable2(list(model10.1), p.varfull), caption="lol"), include.rownames=FALSE,
      caption.placement = "top") 
summary(model10.1)


##For interactions X:Y (see below).  ABIGAIL If you want the intercept, leave that as is; no "-1" (linear model).
##If you are testing different groups, add a "-1"  KATE
OCDataSet <- MakeSet(RawData[which(is.na(RawData$O.C) == FALSE),])
BayesForestPlot(OCDataSet, Group="OC", title="")
model14 <- BayesMeta(Data=OCDataSet, Fix="Zs~SylRepBi:O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)






