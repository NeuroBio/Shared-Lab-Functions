########
#Coded by Cristina Robinson
#Last Modified 5-3-2018
#Written in R-Studio Version 1.1.442
#R Version 3.4.4
#JAGS v4
#metafor_v2.0-0    BEST_0.5.1   HDInterval_0.1.3
#phytools_0.6-44   dplyr_0.7.4  MCMCglmm_2.25
#maps_3.3.0        ape_5.1      coda_0.19-1
#rjags_4.3.0
########

##Set directory here
setwd("C:/Users/Karar/Documents/R/2018-Meta-Analysis")
### Loads packages and my functions
source("BayesMetaFunctions.R")
dir <- getwd()
tag <- "new"

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
OCDataSet <- MakeSet(RawData[which(is.na(RawData$O.C) == FALSE),])
p.var <- var(OCDataSet$Zs)
SylDataSet <- MakeSet(RawData[which(is.na(RawData$Syllable.Rep) == FALSE),])
p.var2 <- var(SylDataSet$Zs)


#make folder for plots
if (file.exists(paste0("Output/",tag))){
  setwd(file.path(dir, "Output",tag))
} else {
  dir.create(file.path(dir, "Output",tag))
  setwd(file.path(dir, "Output",tag))
  
}

#FunnelPlot Series
pdf("Funnel Main.pdf")
par(mfrow=c(2,2),mgp=c(1.75,.5,0), mar=c(4,4,.5,1))
FunnelPlots(FullDataSet, Groups = "Rep",1)
FunnelPlots(FullDataSet, Groups = "OC",2)
dev.off()

#Test for Asymetry
sink("Async.txt")
print("Full")
AsyncTests(FullDataSet)
print("Syl")
AsyncTests(SylDataSet)
print("OC")
AsyncTests(OCDataSet)
sink(NULL)

#Distribution plot
pdf("Distribution of Syllable Repertoire Sizes.pdf")
par(mfrow=c(1,1), mar=c(1,3,1,1), mgp=c(1.5,.5,0))
DistributionPlotV2(FullDataSet)
dev.off()


#Variance Testing, used the OC Data set, as this allows for testing all three Fixed effects models
#models run for 200K itt with burn in of 30K.  Variance estimates for prios was variances in the
#data/number of random effects tested (not incluing SE)
Threshs <- sort(unique(OCDataSet$Syllable.Rep))
Threshs <- Threshs[2:(length(Threshs)-1)]
FixedEffect <- c("Zs~1", "Zs~O.C-1", rep("Zs~SylRepBi-1", length(Threshs)))
Titles <- c("Population variance", "Song stability variance",
            paste(rep("Repertoire size variance, threhold $>$=", length(Threshs)), Threshs))
Titles <- paste(Titles, "in the song stability dataset")
Threshs <- c(38,38,Threshs)
for(i in seq_along(FixedEffect)){
  OCDataSet[,"SylRepBi"] <- cut(OCDataSet$Syllable.Rep, c(0,Threshs[i],2000), labels=c("Small", "Large"), right=FALSE)
  set.seed(49)
  model1 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~Species+idh(Se):units", Vari = p.var)
  model2 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~MType+idh(Se):units", Vari = p.var)
  model3 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~Study+idh(Se):units", Vari = p.var)
  model4 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~Species+animal+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  model5 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~MType+Study+idh(Se):units", Vari = p.var)
  model6 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~MType+Species+animal+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  model7 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~Study+Species+animal+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  model8 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  VarianceTests <- list(model1, model2, model3,
                        model4, model5, model6,
                        model7, model8)
  sink("VarTables.txt", append = TRUE)
  print(xtable(VarianceTable2(VarianceTests, p.var), caption=Titles[i]), include.rownames=FALSE,
        caption.placement = "top") 
  sink(NULL)
}



##Testing the metrics
set.seed(49)
model9 <- BayesMeta(Data=FullDataSet, Fix="Zs~MType-1", Ran="~Species+animal+Study+idh(Se):units", Vari = p.varfull, TreeA=TreeAinv)
MetricData <- ForestTable(ModelList=list(model9),
                       DataSetList=list(FullDataSet),
                       Group = sort(as.character(unique(FullDataSet$MType))))
MetricData["#Measure"] <- table(sort(FullDataSet$MType))
MetricData["#Species"] <- sapply(1:7, function(x) length(unique(FullDataSet$Species[which(FullDataSet$MType == levels(FullDataSet$MType)[x])])))
MetricData <- MetricData[order(MetricData$pMCMC),]
sink("MetricTable.txt", append = TRUE)
print(xtable(MetricData, caption="Metrics model meta-analysis in the song stability dataset"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)





#Effects in the Full Population (Dataset is what changes)
set.seed(49)
model10.1 <- BayesMeta(Data=FullDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.varfull, TreeA=TreeAinv)
model10.2 <- BayesMeta(Data=SylDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
model10.3 <- BayesMeta(Data=OCDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
ConvergePlot(model10.1, "FullPop-FullData", Type = "Pop")
ConvergePlot(model10.2, "FullPop-SylData", Type = "Pop")
ConvergePlot(model10.3, "FullPop-OCData", Type = "Pop")

PopData <- ForestTable(ModelList=list(model10.1, model10.2, model10.3),
                       DataSetList=list(FullDataSet, SylDataSet, OCDataSet),
                       Group = "")
PopData$Group <- c("Full Dataset", "Syllable Repertoire Dataset", "Song Stability Dataset")
sink("MainMetTable.txt", append = TRUE)
print(xtable(PopData, caption="Population meta-analysis"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)

#Stability testing
set.seed(49)
model11 <- BayesMeta(Data=OCDataSet, Fix="Zs~O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
ConvergePlot(model11, "SongStab", Type="Song")

sink("ForestTable.txt", append=TRUE)
print(xtable(ForestTable(ModelList=list(model11), DataSetList=list(OCDataSet), Group=c("Stable", "Plastic")),
             caption="Song stability meta-analysis in the song stability dataset"), include.rownames=FALSE,
      caption.placement = "top")
print(xtable(BESTPrinterOC(OCDataSet),
             caption="Song stability BEST Results in the song stability dataset"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)

pdf("ForestOC.pdf", width=7, height=9.5)
BayesForestPlot(OCDataSet, Group="OC", title="")
dev.off()



#All Thresholds, Seed set internally
ThresholdTest <- JackThresh(SylDataSet, Vari=p.var2, Tree=TreeAinv, Type="Thresh")
Thresholds <- sort(unique(SylDataSet$Syllable.Rep))
Thresholds <- Thresholds[2:(length(Thresholds)-1)]
for(i in seq_along(ThresholdTest)){
  ConvergePlot(ThresholdTest[[i]], "SylRep", Thresholds[i], Type="Rep")
}
sink("MainMetTable.txt", append = TRUE)
print(xtable(ForestTable(ThresholdTest, rep(list(SylDataSet), length(Thresholds)), Picks=Thresholds),
             caption="Syllable repertoire model meta-analysis in the syllable repertoire dataset"),
      include.rownames=FALSE,
      caption.placement = "top")
print(xtable(BestPrinter(SylDataSet, Thresholds[2:(length(Thresholds))]), caption="Syllable repertoire BEST analysis in the syllable repertoire dataset"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)

pdf("ForestREP.pdf", width=7, height=11)
BayesForestPlot(SylDataSet, Group="Rep", title="")
dev.off()




#Subset of thresholds in OC dataset
Picks <- c(18.5, 38, 216)
set.seed(49)
OCDataSet[,"SylRepBi"] <- cut(OCDataSet$Syllable.Rep, c(0,Picks[1],2000), labels=c("Small", "Large"), right=FALSE)
model12.1 <- BayesMeta(Data=OCDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
OCDataSet[,"SylRepBi"] <- cut(OCDataSet$Syllable.Rep, c(0,Picks[2],2000), labels=c("Small", "Large"), right=FALSE)
model12.2 <- BayesMeta(Data=OCDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
OCDataSet[,"SylRepBi"] <- cut(OCDataSet$Syllable.Rep, c(0,Picks[3],2000), labels=c("Small", "Large"), right=FALSE)
model12.3 <- BayesMeta(Data=OCDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
ConvergePlot(model12.1, "SylRepOCDataset", Picks[1], Type="Rep")
ConvergePlot(model12.2, "SylRepOCDataset", Picks[2], Type="Rep")
ConvergePlot(model12.3, "SylRepOCDataset", Picks[3], Type="Rep")

sink("SupMetaTables.txt", append = TRUE)
print(xtable(RepOC <- ForestTable(list(model12.1, model12.2, model12.2), rep(list(OCDataSet),3),Picks=Picks),
             caption="Syllable repertoire model meta-analysis in the song stability dataset"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)




#Everything with RepSize Continuous
set.seed(49)
SylDataSet["lnRep"] <- log(SylDataSet$Syllable.Rep)
model13 <- BayesMeta(Data=SylDataSet, Fix="Zs~lnRep", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)

ConvergePlot(model13, "Continuous model", Type="Cont")
sink("MainMetTable.txt", append = TRUE)
print(xtable(ForestTable(list(model13), list(SylDataSet), c("Intercept", "Larger")),
             caption="Continuous syllable repertoire model meta-analysis in the repertoire dataset"), include.rownames=FALSE,
             caption.placement = "top")
sink(NULL)






#SOng Stability and Repertoire size together
FinalTwoVarTable <- data.frame()
for(i in seq_along(Picks)){
  OCDataSet[,"SylRepBi"] <- cut(OCDataSet$Syllable.Rep, c(0,Picks[i],2000), labels=c("Small", "Large"), right=FALSE)
  set.seed(49)
  model14 <- BayesMeta(Data=OCDataSet, Fix="Zs~SylRepBi:O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  Data <- ForestTable(ModelList=list(model14), DataSetList=list(OCDataSet),
                      Group=c("SmallerStable", "SmallerPlastic","LargerStable", "LargerPlastic"))
  Data <- Tablerfor2Vars(Data=Data, Dataset=OCDataSet)
  if(i == 1){
    FinalTwoVarTable <- Data
  }else{
    FinalTwoVarTable <- rbind(FinalTwoVarTable, Data)
  }
  ConvergePlot(model14, "RepandSongStab",Picks[i], Type="Mixed")
}
sink("MainMetTable.txt", append = TRUE)
print(xtable(FinalTwoVarTable,
      caption="Combined syllable repertoire and song stability model meta-analysis in the song stability dataset"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)




###Without Offspring and minus offspring and EPP
NoOffDataSet <- SylDataSet[which(SylDataSet$MClass != "Offspring"),]
NoOfforEPPDataSet <- NoOffDataSet[which(NoOffDataSet$MType != "EPP"),]
p.var3 <- var(NoOffDataSet$Zs)
p.var4 <- var(NoOfforEPPDataSet$Zs)
FinalTwoVarTable <- data.frame()
FinalTwoVarTable2 <- data.frame()
for(i in seq_along(Picks)){
  NoOffDataSet[,"SylRepBi"] <- cut(NoOffDataSet$Syllable.Rep, c(0,Picks[i],2000), labels=c("Small", "Large"), right=FALSE)
  NoOfforEPPDataSet[,"SylRepBi"] <- cut(NoOfforEPPDataSet$Syllable.Rep, c(0,Picks[i],2000), labels=c("Small", "Large"), right=FALSE)
  set.seed(49)
  model15 <- BayesMeta(Data=NoOffDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var3, TreeA=TreeAinv)
  model16 <- BayesMeta(Data=NoOfforEPPDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var4, TreeA=TreeAinv)
  ConvergePlot(model15, "NoOff", Picks[i], Type="Rep")
  ConvergePlot(model16, "NoOfforEPP", Picks[i], Type="Rep")
  
  Data <- ForestTable(ModelList=list(model15), DataSetList=list(NoOffDataSet), Picks=Picks[i])
  Data2 <- ForestTable(ModelList=list(model16), DataSetList=list(NoOfforEPPDataSet), Picks=Picks[i])
  if(i == 1){
    FinalTwoVarTable <- Data
    FinalTwoVarTable2 <- Data2
  }else{
    FinalTwoVarTable <- rbind(FinalTwoVarTable, Data)
    FinalTwoVarTable2 <- rbind(FinalTwoVarTable2, Data2)
  }
}

sink("NoOffMetaTables.txt", append = TRUE)
print(xtable(FinalTwoVarTable,
             caption="Repertoire model meta-analysis in the dataset without Offspring Measurements"), include.rownames=FALSE,
      caption.placement = "top")
print(xtable(BestPrinter(NoOffDataSet, Picks), caption="Syllable repertoire BEST analysis in the no offspring dataset"), include.rownames=FALSE,
      caption.placement = "top")
print(xtable(FinalTwoVarTable2,
             caption="Repertoire model meta-analysis in the dataset without offspring or EPP Measurements"), include.rownames=FALSE,
      caption.placement = "top")
print(xtable(BestPrinter(NoOfforEPPDataSet, Picks), caption="Syllable repertoire BEST analysis in the no offspring or EPP dataset"), include.rownames=FALSE,
      caption.placement = "top")

sink(NULL)


###############
#Backupchecks

#Each of these has to run many meta-analysis and therefore take a long time!!!!

#Species jackknife, seed set internally
#Seed set internally
Species <- as.character(unique(SylDataSet$Species))
for(i in 1:3){
  SylDataSet[,"SylRepBi"] <- cut(SylDataSet$Syllable.Rep, c(0,Picks[i],2000), labels=c("Small", "Large"), right=FALSE)
  SPJack <- JackThresh(SylDataSet, Vari=p.var2, Tree=TreeAinv, Type="SpJack",Group = "Rep")
  SPJackTable <- ForestTable(ModelList=SPJack, DataSetList=rep(list(SylDataSet),length(SPJack)), Group=c("Smaller", "Larger"))
  SPJackTable["Removed"] <- rep(Species, each=2,1)
  SPJackTable <- SPJackTable[,c(7, 1:6)]
  sink("JackKnifeTables.txt", append = TRUE)
  print(xtable(SPJackTable),caption=paste("Syllable repertoire model meta-analysis with each species removed in turn, threshold $>$=",Picks[i]), include.rownames=FALSE,
        caption.placement = "top")
  sink(NULL)
}

#switching the plasticity categorization, seed set internally
OCJack <- JackThresh(OCDataSet, Vari=p.var, Tree=TreeAinv, Type="OCJack",Group = "OC")
Species <- as.character(unique(OCDataSet$Species))
OCJackTable <- ForestTable(ModelList=OCJack, DataSetList=rep(list(OCDataSet),length(OCJack)), Group=c("Stable", "Plastic"))
OCJackTable["Switched"] <- rep(Species, each=2,1)
OCJackTable <- OCJackTable[,c(7, 1:6)]
sink("JackKnifeTables.txt", append = TRUE)
print(xtable(OCJackTable),caption="Song stability model meta-analysis in song stability dataset with song stability categories of individual species switched", include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)






#Territory Control
#make folder for plots
if(file.exists(paste0(dir,"/Output/",tag, "/Terr"))){
  setwd(file.path(dir, "Output",tag, "Terr"))
} else {
  dir.create(file.path(dir, "Output",tag, "Terr"))
  setwd(file.path(dir, "Output",tag, "Terr"))
  
}

#Chnage out the dataset
Terr <- which(is.na(RawData$Territory.Control) == FALSE)
RawTerr <- RawData
RawTerr$Raw.Cor[Terr] <- RawTerr$Territory.Control[Terr]
Flip <- which(RawTerr$MClass[Terr] == "Time")
RawTerr$Raw.Cor[Terr][Flip] <- RawTerr$Raw.Cor[Terr][Flip]*-1

#set up the larger data set
TFullDataSet <- MakeSet(RawTerr, CUT=38)
p.varfullT <- var(TFullDataSet$Zs)
TOCDataSet <- MakeSet(RawTerr[which(is.na(RawData$O.C) == FALSE),], CUT=38)
p.var5 <- var(TOCDataSet$Zs)
TSylDataSet <- MakeSet(RawTerr[which(is.na(RawData$Syllable.Rep) == FALSE),], CUT=38)
p.var6 <- var(TSylDataSet$Zs)

#FunnelPlot Series
pdf("Funnel Territory.pdf")
par(mfrow=c(2,2),mgp=c(1.75,.5,0), mar=c(4,4,.5,1))
FunnelPlots(TFullDataSet, Groups = "Rep")
FunnelPlots(TFullDataSet, Groups = "OC",2)
dev.off()

#Test for Asymetry
sink("TAsync.txt")
print("Full")
AsyncTests(TFullDataSet)
print("Syl")
AsyncTests(TSylDataSet)
print("OC")
AsyncTests(TOCDataSet)
sink(NULL)

#POPULATION
set.seed(49)
model1.1T <- BayesMeta(Data=TFullDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.varfullT, TreeA=TreeAinv)
model1.2T <- BayesMeta(Data=TSylDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var6, TreeA=TreeAinv)
model1.3T <- BayesMeta(Data=TOCDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var5, TreeA=TreeAinv)
ConvergePlot(model1.1T, "PopTFull", Type="Pop")
ConvergePlot(model1.2T, "PopTSyl", Type="Pop")
ConvergePlot(model1.3T, "PopTOC", Type="Pop")
PopData <- ForestTable(ModelList=list(model1.1T, model1.2T, model1.3T),
                       DataSetList=list(FullDataSet, SylDataSet, OCDataSet),
                       Group = "")
PopData$Group <- c("Full Dataset", "Repertoire Dataset", "Song Stability dataset")
sink("TablesT.txt", append = TRUE)
print(xtable(PopData, caption="Population Meta-Analysis with Territory-controlled Measurements"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)

#RepSize Dataset
set.seed(49)
TSylDataSet[,"SylRepBi"] <- cut(TSylDataSet$Syllable.Rep, c(0,Picks[1],2000), labels=c("Small", "Large"), right=FALSE)
model2.1T <- BayesMeta(Data=TSylDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var6, TreeA=TreeAinv)
TSylDataSet[,"SylRepBi"] <- cut(TSylDataSet$Syllable.Rep, c(0,Picks[2],2000), labels=c("Small", "Large"), right=FALSE)
model2.2T <- BayesMeta(Data=TSylDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var6, TreeA=TreeAinv)
TSylDataSet[,"SylRepBi"] <- cut(TSylDataSet$Syllable.Rep, c(0,Picks[3],2000), labels=c("Small", "Large"), right=FALSE)
model2.3T <- BayesMeta(Data=TSylDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var6, TreeA=TreeAinv)
ConvergePlot(model2.1T,"RepSizeT", Picks[1], Type="Rep")
ConvergePlot(model2.2T,"RepSizeT", Picks[2], Type="Rep")
ConvergePlot(model2.3T,"RepSizeT", Picks[3], Type="Rep")

sink("TablesT.txt", append = TRUE)
print(xtable(ForestTable(list(model2.1T, model2.2T, model2.3T), rep(list(TSylDataSet),3),Picks=Picks),
             caption="Syllable repertoire model in the territory-controlled repertoire dataset meta-analysis"), include.rownames=FALSE,
      caption.placement = "top")
print(xtable(BestPrinter(TSylDataSet, Picks), caption="Syllable Repertoire BEST results with territory-controlled measurements"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)
pdf("ForestREPT.pdf", width=7, height=11)
BayesForestPlot(TSylDataSet, Group="Rep", title="")
dev.off()




#Stability testing
set.seed(49)
model3T <- BayesMeta(Data=TOCDataSet, Fix="Zs~O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var5, TreeA=TreeAinv)
ConvergePlot(model3T, "SongStabT", Type="Song")

sink("TablesT.txt", append=TRUE)
print(xtable(ForestTable(ModelList=list(model16), DataSetList=list(TOCDataSet), Group=c("Stable", "Plastic")),
             caption="Song stability model in the territory-controlled song stability dataset meta-analysis"), include.rownames=FALSE,
      caption.placement = "top")
print(xtable(BESTPrinterOC(TOCDataSet),
             caption="Song stability BEST results with territory-controlled measurements"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)

pdf("ForestOCT.pdf", width=7, height=9.5)
BayesForestPlot(TOCDataSet, Group="OC", title="")
dev.off()


#make folder for plots
if(file.exists(paste0(dir,"/Output/",tag, "/Misc"))){
  setwd(file.path(dir, "Output",tag, "Misc"))
} else {
  dir.create(file.path(dir, "Output",tag, "Misc"))
  setwd(file.path(dir, "Output",tag, "Misc"))
  
}

###Plots and analysis for MinMax
pdf("Distribution of MinMax Syllable Repertoire Sizes.pdf")
par(mfrow=c(1,2))
SylDataSetMAX <- SylDataSet
SylDataSetMAX$Syllable.Rep <- SylDataSet$Syl.Max
DistributionPlotV2(SylDataSetMAX, FALSE, .6)
mtext("Max Values", font=2)
mtext("A", side=2, las=2, at=7, line=1.5, font=2)

SylDataSetMIN <- SylDataSet
SylDataSetMIN$Syllable.Rep <- SylDataSet$Syl.Min
DistributionPlotV2(SylDataSetMIN, FALSE, .6)
mtext("Min Values", font=2)
mtext("B", side=2, las=2, at=7, line=1.5, font=2)
dev.off()

ThresholdTestMAX <- JackThresh(SylDataSetMAX, Vari=p.var2, Tree=TreeAinv, Type="Thresh")
ThresholdTestMIN <- JackThresh(SylDataSetMIN, Vari=p.var2, Tree=TreeAinv, Type="Thresh")

ThresholdsMAX <- sort(unique(SylDataSetMAX$Syllable.Rep))
ThresholdsMAX <- ThresholdsMAX[2:(length(ThresholdsMAX)-1)]
ThresholdsMIN <- sort(unique(SylDataSetMIN$Syllable.Rep))
ThresholdsMIN <- ThresholdsMIN[2:(length(ThresholdsMIN)-1)]

sink("MainMetTable.txt", append = TRUE)
print(xtable(ForestTable(ThresholdTestMAX, rep(list(SylDataSetMAX), length(ThresholdsMAX)), Picks=ThresholdsMAX),
             caption="Syllable repertoire model meta-analysis in the syllable repertoire dataset"),
      include.rownames=FALSE,
      caption.placement = "top")
print(xtable(ForestTable(ThresholdTestMIN, rep(list(SylDataSetMIN), length(ThresholdsMIN)), Picks=ThresholdsMIN),
             caption="Syllable repertoire model meta-analysis in the syllable repertoire dataset"),
      include.rownames=FALSE,
      caption.placement = "top")

print(xtable(BestPrinter(SylDataSetMAX, ThresholdsMAX[2:length(ThresholdsMAX)]), caption="Syllable repertoire BEST analysis in the max syllable repertoire dataset"), include.rownames=FALSE,
      caption.placement = "top")
print(xtable(BestPrinter(SylDataSetMIN, ThresholdsMIN[2:length(ThresholdsMIN)]), caption="Syllable repertoire BEST analysis in the min syllable repertoire dataset"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)


#Just the lowest versus the highest
EndsDataSet <- SylDataSet
Remove <- which((EndsDataSet$Syllable.Rep > 20 & EndsDataSet$Syllable.Rep < 100))
EndsDataSet <- EndsDataSet[-Remove,]
set.seed(49)
EndsModel <- BayesMeta(Data=EndsDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
summary(EndsModel)
EndsTable <- ForestTable(list(EndsModel), list(EndsDataSet), Picks=20)      
EndsTable[2,1] <- ">100" 
sink("JustforRev.txt")
print(xtable(EndsTable, caption="Syllable repertoire BEST analysis in the syllable repertoire dataset"), include.rownames=FALSE,
      caption.placement = "top")
print(xtable(BestPrinter(EndsDataSet, 100), caption="Syllable repertoire BEST analysis in the min syllable repertoire dataset"), include.rownames=FALSE,
      caption.placement = "top")
sink(NULL)


#ARBITRARY DATA

Nu <- 500
set.seed(49)
Arbitrary <- SylDataSet
Spec <- unique(Arbitrary$Species)
#Genertaing this data takes ~7.5 hours.
#Alternatively, load the data in about 5-10 minutes
#with the following 3 lines, and then skip the next code chunk
#load("Bestres.rDATA")
#load("ArbModels.rDATA")
#ArbClasses <- sapply(1:500, function(x) ArbModel[[x]]$X[,1])

######
#Generate the data
ArbModel <- list()
ArbClasses <- matrix(0, nrow=nrow(Arbitrary), ncol=Nu)
BestResults <- list()
for(i in 1:Nu){
  Class <- ifelse(runif(length(Spec)) < .5, "Small", "Large")
  for(j in 1:length(Spec)){
    Arbitrary$SylRepBi[which(Arbitrary$Species == Spec[j])] <- Class[j]
  }
  ArbClasses[,i] <- Arbitrary$SylRepBi 
  ArbModel[[i]] <- BayesMeta(Data=Arbitrary, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
}
#separated because of seed issues
Sm <- which(ArbClasses == 1)
ArbClasses2 <- ArbClasses
ArbClasses2[Sm] <- "Small"
ArbClasses2[-Sm] <- "Large"
for(i in 1:Nu){
  Arbitrary$SylRepBi <- ArbClasses2[,i]
  BestResults[[i]] <- BestPrinter(Arbitrary, NA)
}
######


#Data reoganization; go here post loading rDATA

ArbTable <- ForestTable(ArbModel, rep(list(Arbitrary), Nu), Picks=20)[5:7]

SmallMeasure <- colSums(ArbClasses)
SmallSpec <- sapply(1:Nu, function(x) length(unique(Arbitrary$Species[as.logical(ArbClasses[,x])])))
ArbTable["#Species"] <- c(rbind(SmallSpec, length(Spec)-SmallSpec))
ArbTable["#Measure"] <- c(rbind(SmallMeasure, nrow(Arbitrary)-SmallMeasure))
ArbTable["Threshold"] <- NULL
sig <- ceiling(which(as.numeric(ArbTable$pMCMC) < .025)/2)
percent <- as.numeric(sapply(1:Nu, function(x) BestResults[[x]]$`%<0`))
percentsig <- which(percent < 2.5 | percent > 97.5)
sigmodels <- sig[which(sig %in% percentsig)]

TabledBestResults <- matrix(unlist(BestResults), nrow=Nu, byrow = TRUE)[,2:4]

pdf("Arbitrary Runs.pdf")
par(mar=c(10,3.5,3.8,1),mfrow=c(1,2),mgp=c(1.5,.5,0))
Compare <- cbind(as.numeric(ArbTable$pMCMC[seq(1,length(ArbTable$pMCMC), by=2)]),as.numeric(ArbTable$pMCMC[seq(2,length(ArbTable$pMCMC),by=2)]))
Xaxp <- apply(Compare,1,min)
Flipper <- which(as.numeric(TabledBestResults[,3]) > 50)
TabledBestResults[Flipper,3] <- (as.numeric(TabledBestResults[Flipper,3])-100)*-1
plot(Xaxp, pch=21, bg=rgb(.4,.4,.4,.5),
     as.numeric(TabledBestResults[,3]),
     xlab="pMCMC", ylab="BEST %<0", font.lab=2,
     cex=.6, col=rgb(0,0,0,.8),
     panel.first={
       abline(v=c(.067, .013), lwd=2, lty=c(2,3), col=c("black", "grey65"))
       rect(-.2,-.2,.025,2.5, col=rgb(.5,.5,1,.8), border=NA)
     })
points(c(.067,.05,.042,.047,.029,.021,.01,.009,
         .005,.005,.003,.005,.002,.001,.002,
         .002,006,.005, .006,.029, .046),
       c(22.7,6.6,.3,rep(.1,13),1.3,.1,3.6, 7.9, 32.6),
       pch=24, bg=rgb(1,0,0,.5), cex=1)
mtext("A", side=2, las=2, at=52.5, line=1.5, font=2)
vales <- sort(as.numeric(ArbTable$pMCMC))
print(paste0(100*length(sigmodels)/Nu,"%"))
vales[floor(500*.025)]

#which species in which group
Models <- lapply(sigmodels, function(x) ArbModel[[x]])
ArbTableSig <- ForestTable(Models, rep(list(Arbitrary), length(Models)), Picks=20)[,c(2,5:7)]
Sigs <- which(ArbTableSig$pMCMC < .025)
LSig <- which(Sigs%%2==0)
Small <- sapply(1:length(Models), function(x) unique(Arbitrary$Species[which(Models[[x]]$X[,1]==1)]))
Large <- sapply(1:length(Models), function(x) unique(Arbitrary$Species[which(Models[[x]]$X[,2]==1)]))
SylOrder <- unique(Arbitrary$Species[order(Arbitrary$Syllable.Rep)])
TrueLarge <- as.vector(unlist(c(Small[ceiling(Sigs[-LSig]/2)], Large[ceiling(Sigs[LSig]/2)])))
SplitLarge <- sapply(1:length(Spec), function(x) length(which(TrueLarge==SylOrder[x])))

par(mar=c(11,3.5,0,1),mgp=c(1.5,.5,0))
hist(rep(1:length(Spec),SplitLarge), breaks=0:25, col=rgb(.5,.5,1,.8),
     ylim=c(0, 50), xaxt="n", xlab="", main="", yaxt="n",
     ylab="Percentage of Significant Groups With a Species        ",
     font.lab=2, cex.lab=.9)
rect(0,16.5,length(Spec),27.5, col=rgb(.4,.4,.4,.4), border=NA)
axis(1, seq(.5,length(Spec),by=1), paste0(SylOrder, "(",
                                          sort(Arbitrary[which(duplicated(as.character(Arbitrary$Species))==FALSE),]$Syllable.Rep),
                                          ")"), las=2, pos=0, cex.axis=.6, font.axis=3)
axis(2, seq(0,44,by=44/10), seq(0,100,by=10), las=2, pos=0)
mtext("B", side=2, las=2, at=44.5, line=1.5, font=2)
segments(0,44,25,44, lwd=1, lty=1)
segments(25,0,25,44, lwd=1, lty=1)
segments(0,22,25,22, lwd=2, lty=2)
dev.off()
