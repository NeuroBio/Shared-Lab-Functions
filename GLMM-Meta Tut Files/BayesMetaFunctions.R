####
#Coded by Cristina Robinson
#Last Modified 5-3-2018
#Written in R-Studio Version 1.1.442
#R Version 3.4.4
#JAGS v4
#metafor_v2.0-0    BEST_0.5.1   HDInterval_0.1.3
#phytools_0.6-44   dplyr_0.7.4  MCMCglmm_2.25
#maps_3.3.0        ape_5.1      coda_0.19-1
#rjags_4.3.0
####


#load packages
Packs <- c("Matrix", "coda", "ape", "maps","MCMCglmm",
           "dplyr", "metafor", "rjags",
           "BEST", "xtable", "phytools")
for(i in Packs){
  if(!eval(parse(text=paste0("require(",i,")")))){
    install.packages(i)
  }
  eval(parse(text=paste0("library(",i,")")))
}
rm(list=objects())



#convert to Zs, set of datastructure
fischerTransform <- function(Rs, Ns){
  #correction for rtoz positive bias
  Rs <- Rs-((Rs*(1-Rs^2))/(2*(Ns-3)))
  Zs <- transf.rtoz(Rs)
  return(Zs)
}

MakeSet <- function(Data, CUT=38){
  Data[,"SylRepBi"] <- cut(Data$Syllable.Rep, c(0,CUT,2000), labels=c("Small", "Large"), right=FALSE)
  Data[,"Zs"] <- fischerTransform(Data$Raw.Cor, Data$nBirds)
  Data[,"Se"] <- 1/sqrt(Data$nBirds-3)
  Data$O.C <- factor(Data$O.C, labels=c("Closed", "Open"))
  return(Data)
}

#Stats Tests
AsyncTests <- function(Data){
  print(regtest(Data$Zs, sei=Data$Se))
  print(ranktest(Data$Zs, sei=Data$Se))
}

#Main MetaAnalysis
SetPrior <- function(Ran, Var=1, Nu=.002, a=1000){
  RanVar <- unlist(strsplit(Ran,"[+]"))
  if("idh(Se):units" %in% RanVar){#vars if measurement error
    RandomN <- length(RanVar)-1
    Add <- list(V = diag(1), fix = 1)
  }else{#vars without measurement error
    RandomN <- length(RanVar)
    Add <- NULL
  }
  Gs <- list(V=matrix(Var/(RandomN+1)), nu=Nu, alpha.mu=rep(.1,1), alpha.V=diag(1)*a)
  prior<-list(R=list(V=matrix(Var/(RandomN+1)),nu=Nu),
              G=rep(list(Gs), RandomN))
  prior$G[[RandomN+1]] <- Add
  return(prior)
}
BayesMeta <- function(Data, Fix="Zs~1", Ran="~1",
                      Vari=1, NU=1, TreeA = NULL, Nitt=200000, Burnin=30000){
  #Setup basic prior
  Prior <- SetPrior(Ran,Vari, NU)
  #Prior <- list(R=list(V=diag(7)*(Vari/(1+1)),nu=1),
  #              G=list(list(V=diag(1)*(Vari/(1+1)), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))
  #G=list(list(V=diag(7)*(Vari/(1+1)), nu=7, alpha.mu=rep(0,7), alpha.V=diag(7)*1000)))
  if(!is.null(TreeA)){
    #set up tree
    drop <- which(is.na(match(rownames(TreeAinv), Data$Species)))
    TreeAinvSubSet <- TreeAinv[-drop,-drop]
    Meta <- MCMCglmm (as.formula(Fix),random=as.formula(Ran),prior=Prior,
                      ginverse=list(animal=TreeAinvSubSet), #for phylogeny
                      #mev=Data$Se^2,
                      data=Data, nitt=Nitt, burnin=Burnin, verbose=FALSE)
  }else{
    Meta <- MCMCglmm (as.formula(Fix),random=as.formula(Ran),prior=Prior,
                      #mev=Data$Se^2,
                      data=Data, nitt=Nitt, burnin=Burnin, verbose=FALSE)
  }
  return(Meta)
}

#Basic Graphs
FunnelPlots <- function(Data, Groups="Rep", pos=1){
  colors <- c(rgb(0,0,1,.75), rgb(1,0,0,.75),rgb(0,0,0,.75))
  plot(type='n', Data$Zs,1/Data$Se, xlim=c(-1.5,1.5), ylim=c(min(1/Data$Se),max(1/Data$Se)),
       xlab="Fisher's Transformed Z", ylab = "Inverse Standard Error",
       font=2, font.lab=2, cex.axis=.8, cex.lab=.9)
  abline(v=mean(Data$Zs), lty=3, lwd=2, col="grey40")
  if(Groups=="Rep"){
    log <- log(Data$Syllable.Rep)
    Black <- which(is.na(log))
    minlog <- min(log[-Black])
    maxlog <- max(log[-Black]-minlog)
    colscale <- 1-(log[-Black]-minlog)/maxlog
    points(Data$Zs[-Black], 1/Data$Se[-Black], pch=21, col="black",bg=rgb(1,colscale,colscale,1))
    points(Data$Zs[Black], 1/Data$Se[Black], pch=19)
    scale <- seq(0,1,.05)
    text(-1.35,6.75,"Smaller", adj=0, font=2)
    text(-1.35,9.5,"Larger", adj=0, font=2)
    points(rep(-1.45,length(scale)), 6.75+scale*2.75,
               pch=21, bg=rgb(1,rev(scale),rev(scale),1), col="black")
    #rect(-1.5, 7+scale*2, -1.4,7.5+scale*2, border=NA, col=rgb(1,rev(scale),rev(scale),1))
    #rect(-1.5,7,-1.4,9.5)
    points(-1.45, 6.25,pch=19)
    text(-1.35, 6.25, "Unknown", adj=0, font=2)
    mtext(LETTERS[pos], 2, 1.5, FALSE,max(1/Data$Se), las=2)
    return()
  }else if(Groups=="OC"){
    G1 <- which(Data$O.C == "Closed")
    G2 <- which(Data$O.C == "Open")
    LegTex <- c("Song-Stable", "Song-Plastic", "Unknown")
  }else if(Groups=="Metric"){
    G1 <- which(Data$MClass == "Time")
    G2 <- which(Data$MClass == "Offspring")
    LegTex <- c("Time Latency", "NumOffspring", "NumFemale")
  }
  points(Data$Zs[G1], 1/Data$Se[G1], pch=19, col=colors[1])
  points(Data$Zs[G2], 1/Data$Se[G2], pch=19, col=colors[2])
  points(Data$Zs[-c(G1,G2)], 1/Data$Se[-c(G1,G2)], pch=19, col=colors[3])
  par(font=2)
  legend("topleft", LegTex,
         cex=.9,pch=19, pt.bg="black", col=colors)
  par(font=1)
  mtext(LETTERS[pos], 2, 1.5, FALSE,max(1/Data$Se), las=2)
}
CheckPlot <- function(modelVec, title="", pos=1, SE=FALSE){
  plot.default(modelVec, type='l', ylab="", xlab= "Iterations",
               main=paste0("Trace of ", title), cex.main=1)
  if(SE != TRUE){
    mtext(LETTERS[2*pos-1], 2,1.5, las=2, FALSE, max(modelVec))
  }else{
    mtext(LETTERS[2*pos-1], 2,1.5, las=2, FALSE, 1.4)
  }
  dense <- density(modelVec)
  plot(dense, cex.main=1,
       main=paste0("Density of ", title), font.main=2)
  segments(modelVec, 0, modelVec, max(dense$y)*.03)
    mtext(LETTERS[2*pos], 2,1.5, las=2, FALSE, max(dense$y)) 
}
DistributionPlotV2 <- function(Data, lines=TRUE, specsize=.8){
  #pull data, combine, kill dups, sort, separate out closed-ended learners
  Data <- Data[which(is.na(Data$Syllable.Rep)==FALSE),]
  Data <- Data[order(Data$Syllable.Rep),]
  nodups <-Data[which(duplicated(as.character(Data$Species))==FALSE),]
  NumStud <- vector("numeric", length(nodups$Species))
  for(i in seq_along(NumStud)){
    SpecInd <- which(Data$Species == nodups$Species[i])
    NumStud[i] <- length(unique(Data$Study[SpecInd]))
  }
  Col <- rep("black", length(nodups$Species))
  Col[which(nodups$O.C==1 | nodups$O.C== "Closed")] <- "blue"
  Col[which(nodups$O.C==2| nodups$O.C== "Open")] <- "red"
  #Close <-nodups[Cind,]
  #stuff for graphing
  ticks <- length(nodups$Species)
  if(lines == TRUE){
    par(mar = c(15,5,.5,1), mgp =c(1.6,.6,0)) 
  }else{
    par(mar = c(15,3,1.5,1), mgp =c(1.6,.6,0))
  }
  
  
  #plot for ln() graph
  Syllable.Rep <- log(nodups$Syllable.Rep)
  #Close$repsize <- log(Close$repsize)
  plot(0, type="n", xlim = c(1,ticks), ylim = c(0, max(Syllable.Rep)),
       las = 2, xaxt = "n", xlab = "", ylab ="ln(Repertoire Size)",  cex = .6,
       cex.axis = .75, font.lab = 2, font.axis=2,
       #Done first so line not crossing axis
       panel.first = {
         if(lines==TRUE){
         rect(which(nodups$Syllable.Rep==17.43)-.5,-1,
              which(nodups$Syllable.Rep==216)-.5,8,
              col = "grey90", border=NA)
         rect(which(nodups$Syllable.Rep==18.15)-.5,-1,
              which(nodups$Syllable.Rep==22.5)-.5,8,
              col = "grey75", border=NA)
           abline(v=c(8.5, 15.5, 21.5), col="black", lty=2)
         }
         })
  
  segments(2:ticks-1, Syllable.Rep[2:ticks-1], 2:ticks, Syllable.Rep[2:ticks],
           lty = 1, lwd = 1.5, col = "black")
  points(1:ticks, Syllable.Rep, pch=19, col = Col, cex = 1+(NumStud-1)*.4)
  axis(1, 1:ticks, labels = paste(nodups$Species, "  (",nodups$Syllable.Rep,")", sep=""), cex.axis = specsize, las=2, font  = 4)
  par(font=2)
  legend("topleft", legend=c("Song-Stable", "Song-Plastic  ", "Unknown"),
         col=c("blue", "red", "Black"), pch = 19, cex=.8)
  par(font=1)
  mtext(side=1, "Species", line=12, font=2)
}

ThresholTestingPlot <- function(ThreshData, CUT=38){
  Thresh <- sort(unique(ThreshData$Syllable.Rep))
  Thresh <- Thresh[2:(length(Thresh)-1)]
  pvals <- data.frame(Small = numeric(), Large = numeric())
  Means <- pvals
  Cred95 <- data.frame(SmallL = numeric(), SmallU = numeric(),
                       LargeL = numeric(), LargeU = numeric())
  Cred50 <- Cred95
  for(i in seq_along(Thresh)){
    set.seed(49)
    ThreshData[,"SylRepBi"] <- cut(ThreshData$Syllable.Rep, c(0,Thresh[i],2000), labels=c("Small", "Large"), right=FALSE)
    modelThresh <- BayesMeta(Data=ThreshData, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
    summary(modelThresh)
    Means[i,] <- c(mean(modelThresh$Sol[,1]), mean(modelThresh$Sol[,2]))
    Cred95[i,] <- c(coda::HPDinterval(modelThresh$Sol[,1])[1:2], coda::HPDinterval(modelThresh$Sol[,2])[1:2])
    Cred50[i,] <- c(coda::HPDinterval(modelThresh$Sol[,1], prob=.5)[1:2],  coda::HPDinterval(modelThresh$Sol[,2], prob=.5)[1:2])
    pvals[i,] <- 2*pmax(0.5/dim(modelThresh$Sol)[1], pmin(colSums(modelThresh$Sol[,, drop = FALSE] > 0)/dim(modelThresh$Sol)[1], 
                                                          1 - colSums(modelThresh$Sol[,, drop = FALSE] > 0)/dim(modelThresh$Sol)[1]))
  }
  #cbind(Thresh,pvals)
  set1 <- 1:(length(Thresh)-1)
  set2 <- 2:length(Thresh)
  set3 <- log(Thresh[set1])
  set4 <- log(Thresh[set2])
  set5 <- 1:length(Thresh)
  set6 <- log(Thresh)
  pdf("TestIntervals.pdf")
  par(mfrow=c(1,2), mar=c(3,3,.5,.5), mgp=c(1.5,.5,0))
  plot(1,type="n", xlim=c(min(log(Thresh)), max(log(Thresh))), ylim=c(.001,1), log='y',
       ylab="pvalue", xlab="ln(Repertoire Size)")
  pMCMCFull <- pMCMC <- 2*pmax(0.5/dim(model8.1$Sol)[1], pmin(colSums(model8.1$Sol[,, drop = FALSE] > 0)/dim(model8.1$Sol)[1], 
                                                              1 - colSums(model8.1$Sol[,, drop = FALSE] > 0)/dim(model8.1$Sol)[1]))
  abline(h=pMCMCFull, lty=2, col="black")
  abline(v=log(CUT))
  segments(set3,pvals[set1,1], set4,pvals[set2,1], lwd=2, col="blue")
  segments(set3,pvals[set1,2], set4,pvals[set2,2], lwd=2, col="red")
  
  plot(1,type="n", xlim=c(min(log(Thresh)), max(log(Thresh))), ylim=c(-1,1.5),
       ylab="Posterior Mean", xlab="ln(Repertoire Size)")
  abline(h=mean(model8.1$Sol[,1]), lty=2, col="black")
  abline(h=0, lwd=2)
  abline(v=log(CUT))
  polygon(x=c(set6, rev(set6)), y=c(Cred95[set5,1],rev(Cred95[set5,2])), col=rgb(0,0,1,.1), border=rgb(0,0,1,.5))
  polygon(x=c(set6, rev(set6)), y=c(Cred95[set5,3],rev(Cred95[set5,4])), col=rgb(1,0,0,.1), border=rgb(1,0,0,.5))
  polygon(x=c(set6, rev(set6)), y=c(Cred50[set5,1],rev(Cred50[set5,2])), col=rgb(0,0,1,.5), border=rgb(0,0,1,.4), density=50)
  polygon(x=c(set6, rev(set6)), y=c(Cred50[set5,3],rev(Cred50[set5,4])), col=rgb(1,0,0,.5), border=rgb(1,0,0,.4), density=50)
  segments(set3, Means[set1,1], set4,Means[set2,1], lwd=2, col="blue")
  segments(set3, Means[set1,2], set4,Means[set2,2], lwd=2, col="red")
}
BayesForestPlot <- function(Data, title="Forest Plot", Group="Rep",
                            mainfont =.7, model1=NULL, model2=NULL){
  #if(is.null(model1) && is.null(model2)){
  #  layout(matrix(c(rep(1,6),2,3,4,5,6,7), nrow=2, ncol=6,byrow=TRUE),
  #         widths = c(.22,.06,.1,.1,.15,.36), heights = c(.05, .95)) 
    layout(matrix(c(rep(1,6),2,3,4,5,6, 7, rep(8,6)), nrow=3, ncol=6,byrow=TRUE),
           widths = c(.22,.06,.1,.1,.15,.36), heights = c(.03, 1, .07)) 

  if(Group == "OC"){
    Data <- Data[order(Data$O.C,decreasing = FALSE),]
    BoxCol <- ifelse(Data$O.C=="Closed", "blue",ifelse(Data$O.C=="Open", "red", "black"))
    BottomLab <- c("Song-Stable", "Song-Plastic")
  }else if(Group == "Rep"){
    Data <- Data[order(Data$Syllable.Rep, decreasing = FALSE),]
    BoxCol <- "red"
    BottomLab <- c("Smaller", "Larger")
  } else{
      BoxCol <- "grey50"
  }
  Start <- length(Data$Species)
  PlotSize <- Start+2
  #Title
  par(mar=rep(.01,4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.3,.3,title, adj=.3, cex=1.2, font=2)
  
  par(mar=c(3,.01,.01,.01))
  
  #species/Study labels
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(1,PlotSize))
  text(x=1,y=PlotSize, "Species (Rep) [citation]", adj=1,font=4)
  text(x=1, y=Start:1,
       paste0(Data$Species, " (", Data$Syllable.Rep, ") ", Data$Study), adj=1, font=3, cex=mainfont)
  #nBirds Labels
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(1,PlotSize))
  text(x=.5,y=PlotSize, "nBirds", adj=.5, font=2)
  text(x=.5, y=Start:1, Data$nBirds, adj=.5, cex=mainfont)
  #Weight Labels
  #plot.new()
  #plot.window(xlim=c(0,1), ylim=c(1,length(PairDate$species.name)+1))
  #text(x=.5,y=length(PairDate$species.name)+1, "Weight", adj=.5, font=2)
  #text(x=.5, y=length(PairDate$species.name):1, paste("idk", "%", sep=""), adj=.5)
  
  #Correlation Labels
  #plot.new()
  #plot.window(xlim=c(0,1), ylim=c(1,PlotSize))
  #text(x=.5,y=PlotSize, "Raw Cor", adj=.5, font=2)
  #text(x=.5, y=Start:1, Data$Raw.Cor, adj=.5, cex=mainfont)
  
  #Correlation Labels
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(1,PlotSize))
  text(x=.5,y=PlotSize, "Fisher's Z", adj=.5, font=2)
  text(x=.5, y=Start:1, round(Data$Zs, digits=3), adj=.5, cex=mainfont)
  
  #Confidence Interval
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(1,PlotSize))
  upper <- round(Data$Zs - qnorm(.05/2, lower.tail = FALSE) * Data$Se, digits=3)
  lower <- round(Data$Zs + qnorm(.05/2, lower.tail = FALSE) * Data$Se, digits=3)
  text(x=.5,y=PlotSize, "[95% CI]", adj=.5, font=2)
  text(x=.5, y=Start:1,
       paste("[",upper,";",lower,"]"), adj=.5, cex=mainfont)
  
  #Forest Plot
  par(mar=c(3,.01,2,.01), mgp=c(1.5,.5,0))
  plot(0,type='n',xlim=c(-1.5,1.5), ylim=c(1,PlotSize),
       bty='n', yaxt='n',
       xlab="Fisher's Z", font.lab=2)
  abline(v=0, lwd=2)
  abline(v=mean(Data$Zs), lty=2, lwd=2, col="grey50")
  bottom <- 1
  Fill <- .5
  rect(xright=Data$Zs-.1,xleft=Data$Zs+.1,
       ybottom=seq(PlotSize+Fill, bottom, length.out = Start),
       ytop=seq(PlotSize+Fill, bottom, length.out = Start)+.5,
       col=BoxCol, border= NA)
  segments(x0=Data$Zs,x1=Data$Zs,
           y0=seq(PlotSize+Fill, bottom, length.out = Start)+.125,
           y1=seq(PlotSize+Fill, bottom, length.out = Start)+.375)
  segments(x0=upper,x1=lower,
           y0=seq(PlotSize+Fill, bottom, length.out = Start)+.25,
           y1=seq(PlotSize+Fill, bottom, length.out = Start)+.25)
  if(Group == "Rep"){
    place <- sapply(c(18.5, 38, 216),
                    function(x) min(which(sort(Data$Syllable.Rep) == x)))
    #print((seq(PlotSize+Fill, 0, length.out = Start))[place-1])
    abline(h=c(71.7,45.9,9), lty=2, col="grey50")
  }
  plot.new()
  
  #if(is.null(model1)==FALSE){
  #  par(mar=c(.01,.01,.01,.01))
  #  plot.new()
  #  plot.window(xlim=c(0,1), ylim=c(0,1))
  #  Printout <- c(mean(model1$Sol[,1]),coda::HPDinterval(model1$Sol),
  #                effectiveSize(model1$Sol[, 1, drop = FALSE]),
  #                2*pmax(0.5/dim(model1$Sol)[1], pmin(colSums(model1$Sol[,1, drop = FALSE] > 0)/dim(model1$Sol)[1], 
  #                                                  1 - colSums(model1$Sol[, 1, drop = FALSE] > 0)/dim(model1$Sol)[1])) )
  #  Printout <- round(Printout, digits=4)
  #  if(Printout[5] == 0){
  #    Printout[5] <- "<0.0001"
  #  }
   # text(seq(.05,.6, length.out=6), .95, c("Group", "post.mean", "l-95% CI", "u-95% CI", "eff.samp",  "pMCMC"), font=2, cex=mainfont)
  #  text(seq(.05,.6, length.out=6), .70, c("Full Population",Printout), cex=mainfont)
  #  if(is.null(model2)==FALSE){
  #    i <- 2
  #    for(i in seq_along(colnames(model2$Sol))){
  #      Printout <- data.frame(mean(model2$Sol[,i]),coda::HPDinterval(model2$Sol[,i]),
  #                    effectiveSize(model2$Sol[, i, drop = FALSE]),
  #                    2*pmax(0.5/dim(model2$Sol)[1], pmin(colSums(model2$Sol[,i, drop = FALSE] > 0)/dim(model2$Sol)[1], 
  #                                                    1 - colSums(model2$Sol[, i, drop = FALSE] > 0)/dim(model2$Sol)[1])) )
  #      Printout <- round(Printout, digits=4)
  #      Printout[4] <- ceiling(Printout[4])
  #      if(Printout[5] == 0){
  #        Printout[5] <- "<0.0001*"
  #      }else if(Printout[5] < .025){
  #        Printout[5] <- paste(Printout[5], "*")
  #      }
  #      text(seq(.05,.6, length.out=6), .7-(.25*i), c(BottomLab[i],Printout), cex=mainfont)
  #    }
  #  }
  #}
}


#Misc Data pulling
VarianceTable <- function(list, Var, Size){
  #Make Table
  Table <- data.frame(Fixed=character(),Random=character(),
                      Variance=character(), DIC=character(),
                      stringsAsFactors = FALSE)
  count <- 1
  for(i in seq_along(list)){
    #get variance data, Get first row of data
    VarData <- posterior.mode(list[[i]]$VCV)
    Table[count,] <- c(deparse(list[[i]]$Fixed$formula),
                       names(VarData[1]),
                       paste0(round(VarData[1]/Var*100, digits=2), "%"),
                       round(list[[i]]$DIC, digits=2))
    #if more than one random variable (ignoring units and SE
    VarTypes <- length(VarData)-2
    if(VarTypes > 1){
      for(j in 2:VarTypes){
        Table[count+(j-1),] <- list('""',
                                    names(VarData[j]),
                                    paste0(round(VarData[j]/Var*100, digits=2), "%"),
                                    '""') 
      }
    }
    count <- count+VarTypes
  }
  
  #plotting table
  par(mar=c(.2,.2,.2,.2))
  plot.new()
  max <- length(list)+1
  leg <- length(Table$Fixed)+1
  plot.window(xlim=c(0,1), ylim=c(1,max))
  #Headers
  a <- 0
  b <- .3
  c <- .6
  d <- .9
  #text(c(a,b,c,d),max, c("Fixed", "Random", "%Variance", "DIC"), adj=0)
  #Fixed
  text(a,seq((max),1, length.out = leg),c(colnames(Table)[1],Table$Fixed), adj=0, cex=Size)
  #Random
  text(b,seq((max),1, length.out = leg),c(colnames(Table)[2],Table$Random), adj=0, cex=Size)
  #%Variance
  text(c,seq((max),1, length.out = leg),c(colnames(Table)[3],Table$Variance), adj=0, cex=Size)
  #DIC
  text(d,seq((max),1, length.out = leg),c(colnames(Table)[4],Table$DIC), adj=0, cex=Size)
}

#Extra Analyses
JackThresh <- function(dataset, Vari, Tree, Type="OCJack", Group = "OC"){
  #changes the OC grouping of a study or species
  #aternatively, tests all thresholds
  
  #Set up Data Structures
  if(Type == "Thresh"){
    Thresh <- sort(unique(dataset$Syllable.Rep))
    Thresh <- Thresh[2:(length(Thresh)-1)]
    NumItt <- length(Thresh)
    
    pvals <- data.frame(Small = numeric(NumItt), Large = numeric(NumItt), row.names=Thresh)
    Means <- pvals
    Cred95 <- data.frame(SmallL = numeric(NumItt), SmallU = numeric(NumItt),
                         LargeL = numeric(NumItt), LargeU = numeric(NumItt), row.names=Thresh)
    Cred50 <- Cred95
  }else{
    Species <- as.character(unique(dataset$Species))
    NumItt <- length(Species)
    if(Group == "Rep"){
      pvals <- data.frame(Small = numeric(NumItt), Large = numeric(NumItt), row.names=Species)
      Means <- pvals
      Cred95 <- data.frame(SmallL = numeric(NumItt), SmallU = numeric(NumItt),
                         LargeL = numeric(NumItt), LargeU = numeric(NumItt), row.names=Species)
      Cred50 <- Cred95
      FIXED <- "Zs~SylRepBi-1"
      if(Type=="OCJack"){
        warning("This thing you are doing makes no sense.  Changing the OC status will not affect the repsize results.")
      }
    }else{
      pvals <- data.frame(Closed = numeric(NumItt), Open = numeric(NumItt), row.names=Species)
      Means <- pvals
      Cred95 <- data.frame(ClosedL = numeric(NumItt), ClosedU = numeric(NumItt),
                         OpenL = numeric(NumItt), OpenU = numeric(NumItt), row.names=Species)
      Cred50 <- Cred95
      FIXED <- "Zs~O.C-1"
    }
  }
  #OC Jacknife
  if(Type == "OCJack"){
    Jack <- cbind(replicate(length(Species), dataset$O.C))
    JackKnife <- Replacer(dataset, Jack)
    models <- list()
    for(i in seq_along(Species)){
      set.seed(49)
      dataset$O.C <- JackKnife[,i]
      #This is model 8/9.3 in the main analysis!
      model <- BayesMeta(Data=dataset, Fix=FIXED,
                             Ran="~MType+Species+animal+Study+idh(Se):units",
                             Vari = Vari, TreeA=Tree)
      models[[i]] <- model
      #Means[i,] <- c(mean(model$Sol[,1]), mean(model$Sol[,2]))
      #Cred95[i,] <- c(coda::HPDinterval(model$Sol[,1])[1:2], coda::HPDinterval(model$Sol[,2])[1:2])
      #Cred50[i,] <- c(coda::HPDinterval(model$Sol[,1], prob=.5)[1:2],  coda::HPDinterval(model$Sol[,2], prob=.5)[1:2])
      #pvals[i,] <- 2*pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1], 
      #                                                      1 - colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1]))
    }
    #Species Jackknife
  }else if(Type == "SpJack"){
    HoldBack <- list()
    for(i in seq_along(Species)){
      HoldBack[[i]] <- which(dataset$Species == Species[i])
    }
    models <- list()
    for(i in seq_along(Species)){
      set.seed(49)
      Spec <- dataset[-HoldBack[[i]],]
      #This is model 8/9.2 or 8/9.3 in the main analysis!
      model <- BayesMeta(Data=Spec, Fix=FIXED,
                         Ran="~MType+Species+animal+Study+idh(Se):units",
                         Vari = Vari, TreeA=Tree)
      models[[i]] <- model
      #Means[i,] <- c(mean(model$Sol[,1]), mean(model$Sol[,2]))
     # Cred95[i,] <- c(coda::HPDinterval(model$Sol[,1])[1:2], coda::HPDinterval(model$Sol[,2])[1:2])
      #Cred50[i,] <- c(coda::HPDinterval(model$Sol[,1], prob=.5)[1:2],  coda::HPDinterval(model$Sol[,2], prob=.5)[1:2])
      #pvals[i,] <- 2*pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1], 
      #                                                1 - colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1]))
    }
    #Threshold
  }else{
    Thresh <- sort(unique(dataset$Syllable.Rep))
    Thresh <- Thresh[2:(length(Thresh)-1)]
    ThreshData <- dataset
    models <- list()
    for(i in seq_along(Thresh)){

      set.seed(49)
      #This is model 8/9.2 in the main analysis!
      ThreshData[,"SylRepBi"] <- cut(dataset$Syllable.Rep, c(0,Thresh[i],2000), labels=c("Small", "Large"), right=FALSE)
      model <- BayesMeta(Data=ThreshData, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = Vari, TreeA=TreeAinv)
      #Means[i,] <- c(mean(model$Sol[,1]), mean(model$Sol[,2]))
      #Cred95[i,] <- c(coda::HPDinterval(model$Sol[,1])[1:2], coda::HPDinterval(model$Sol[,2])[1:2])
      #Cred50[i,] <- c(coda::HPDinterval(model$Sol[,1], prob=.5)[1:2],  coda::HPDinterval(model$Sol[,2], prob=.5)[1:2])
      #pvals[i,] <- 2*pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1], 
      #                                                    1 - colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1]) )
      models[[i]] <- model
    }
  }
  #Output <- list(Means, Cred95, Cred50, pvals)
  #names(Output) <- c("Posterior.Mean", "Credibility95%", "Credibility50%", "pVals")
  #Output <- list(Output, models)
  return(models)
}
CredibilityIntPlots <- function(JakThr, model, Type){
  pvals <- JakThr$pVals
  Means <- JakThr$Posterior.Mean
  Cred95 <- JakThr$`Credibility95%`
  Cred50 <- JakThr$`Credibility50%`
  pMCMCFull<- 2*pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,1:ncol(model$Sol), drop = FALSE] > 0)/dim(model$Sol)[1], 
                     1 - colSums(model$Sol[,1:ncol(model$Sol), drop = FALSE] > 0)/dim(model$Sol)[1]))
  
  
  #breaks for segment and polygon
  set1 <- 1:(length(rownames(pvals))-1)
  set2 <- set1+1
  set5 <- 1:length(rownames(pvals))
  if(length(model$Sol[1,]) == 1){
    PostMean <- mean(model$Sol)
  }else{
    PostMean <- colMeans(model$Sol)
  }
  #X-axis points 
  if(Type == "Thresh"){
    xaxt <- as.numeric(rownames(pvals)) 
    set3 <- log(xaxt[set1])
    set4 <- log(xaxt[set2])
    set6 <- log(xaxt)
    xlabel <- "Syllable Repertoire Size"
    bmar <- 3
  }else{
    xaxt <- rownames(pvals) 
    set3 <- set1
    set4 <- set2
    set6 <- set5
    if(Type=="SPJack"){
      xlabel <- "Species Removed"
    }else{
      xlabel <- "Species Switched" 
    }
    bmar <- 10
  }
  
  
    
    #pdf("TestIntervals.pdf")
    par(mfrow=c(1,2), mar=c(bmar,4,.5,.5), mgp=c(1.5,.5,0))
    plot(1,type="n", xlim=c(min(set6), max(set6)), ylim=c(.001,1), log='y',
         xaxt="n", ylab="pMCMC", xlab="")
    axis(1, at=set6, labels = xaxt, las=2, cex.axis=.7)
    mtext(xlabel, side=1, line=bmar-1)
    abline(h=pMCMCFull, lty=2, col=c("blue","red"))
    abline(h=.025, lwd=2)
    segments(set3,pvals[set1,1], set4,pvals[set2,1], lwd=2, col="blue")
    segments(set3,pvals[set1,2], set4,pvals[set2,2], lwd=2, col="red")
    
    plot(1,type="n", xlim=c(min(set6), max(set6)), ylim=c(-1,1.5),
         xaxt="n", ylab = "Posterior Mean", xlab="")
    axis(1, at=set6, labels = xaxt, las=2, cex.axis=.7)
    mtext(xlabel, side=1, line=bmar-1)
    abline(h=PostMean, lty=2, col=c("blue","red"))
    abline(h=0, lwd=2)
    polygon(x=c(set6, rev(set6)), y=c(Cred95[set5,1],rev(Cred95[set5,2])), col=rgb(0,0,1,.1), border=NA)
    polygon(x=c(set6, rev(set6)), y=c(Cred95[set5,3],rev(Cred95[set5,4])), col=rgb(1,0,0,.1), border=NA)
    polygon(x=c(set6, rev(set6)), y=c(Cred50[set5,1],rev(Cred50[set5,2])), col=rgb(0,0,1,.5), border=NA, density=50)
    polygon(x=c(set6, rev(set6)), y=c(Cred50[set5,3],rev(Cred50[set5,4])), col=rgb(1,0,0,.5), border=NA, density=50)
    segments(set3, Means[set1,1], set4,Means[set2,1], lwd=2, col="blue")
    segments(set3, Means[set1,2], set4,Means[set2,2], lwd=2, col="red")
}
Replacer <- function(dataset, Jack){
  ind <- which(duplicated(dataset$Species) == FALSE)
  for(i in seq_along(ind)){
    ifelse(ind[i] == "Closed", Switch <- "Open", Switch <- "Closed")
    ifelse(i != length(ind), Jack[ind[i]:(ind[i+1]-1),i] <- Switch, Jack[ind[i]:length(Jack[,i]),i] <- Switch) 
  }
  return(Jack)
}


BestPrinter <- function(Dataset, Threshs){
  ColNames <- c("Threshold", "BEST Mean", "95% CredInt", "%<0")
  DataFrame <- data.frame(character(), character(),character(),character(),
                          stringsAsFactors=FALSE)
  names(DataFrame) <- ColNames
  
  for(i in seq_along(Threshs)){
    if(is.na(Threshs) == FALSE){
      CUT <- Threshs[i]
      Dataset$SylRepBi <-  cut(Dataset$Syllable.Rep, c(0,CUT,2000), labels=c("Small", "Large"), right=FALSE) 
    }
    SL <- BESTmcmc(Dataset$Zs[which(Dataset$SylRepBi == "Large")],
                   Dataset$Zs[which(Dataset$SylRepBi == "Small")],
                   numSavedSteps=16998, rnd.seed=49)
    paramSampleVec <- SL$mu1-SL$mu2
    by <- diff(hdi(paramSampleVec))/18
    breaks <- unique(c(seq(from = min(paramSampleVec), 
                           to = max(paramSampleVec), by = by), max(paramSampleVec)))
    histinfo <- hist(paramSampleVec, breaks = breaks, plot = FALSE)
    cvHt <- 0.7 * max(histinfo$density)
    pcgtCompVal <- round(100 * sum(paramSampleVec < 0)/length(paramSampleVec), 
                         1)
    if(pcgtCompVal == 0){
      pcgtCompVal <- "0.1"
    }
    HDI <- hdi(paramSampleVec, .95)
    DataFrame[i,] <- cbind(paste0(">=",Threshs[i]), round(mean(HDI), digits=3),
                           paste0("[", round(HDI[1], digits=3), ";", round(HDI[2], digits=3), "]"),
                           pcgtCompVal)
  }
  return(DataFrame)
}
VarianceTable2 <- function(list, Var, Size){
  #Make Table
  Table <- data.frame(Random=character(),
                      Variance=character(), DIC=character(),
                      stringsAsFactors = FALSE)
  count <- 1
  for(i in seq_along(list)){
    #get variance data, Get first row of data
    VarData <- posterior.mode(list[[i]]$VCV)
    Table[count,] <- c(names(VarData[1]),
                       paste0(round(VarData[1]/Var*100, digits=2), "%"),
                       round(list[[i]]$DIC, digits=2))
    #if more than one random variable (ignoring units and SE
    VarTypes <- length(VarData)-2
    if(VarTypes > 1){
      for(j in 2:VarTypes){
        VarType <- names(VarData[j])
          if(VarType == "animal"){
            VarType <- "Phylo"
          }
        Table[count+(j-1),] <- list(VarType,
                                    paste0(round(VarData[j]/Var*100, digits=2), "%"),
                                    '""') 
      }
    }
    count <- count+VarTypes
  }
  return(Table)
  #plotting table
  #par(mar=c(.2,.2,.2,.2))
  #plot.new()
  #max <- length(list)+1
  #leg <- length(Table$Fixed)+1
  #plot.window(xlim=c(0,1), ylim=c(1,max))
  #Headers
  #a <- 0
  #b <- .3
  #c <- .6
  #d <- .9
  #text(c(a,b,c,d),max, c("Fixed", "Random", "%Variance", "DIC"), adj=0)
  #Fixed
  #text(a,seq((max),1, length.out = leg),c(colnames(Table)[1],Table$Fixed), adj=0, cex=Size)
  #Random
  #text(b,seq((max),1, length.out = leg),c(colnames(Table)[2],Table$Random), adj=0, cex=Size)
  #%Variance
  #text(c,seq((max),1, length.out = leg),c(colnames(Table)[3],Table$Variance), adj=0, cex=Size)
  #DIC
  #text(d,seq((max),1, length.out = leg),c(colnames(Table)[4],Table$DIC), adj=0, cex=Size)
}
ModelDataPuller <- function(ModelList, Group){
  counter <- 1
  DataFrame <- data.frame(character(),character(),
                          character(),character(),
                          stringsAsFactors = FALSE)
  names(DataFrame) <- c("Group", "Post Mean", "95% CredInt", "pMCMC")
  for(i in seq_along(ModelList)){
    for(j in seq_along(Group)){
      credInt <- coda::HPDinterval(ModelList[[i]]$Sol[,j])
      DataFrame[counter,] <- cbind(Group[j],
                                   round(mean(ModelList[[i]]$Sol[,j]), digits=3),
                                   paste0("[", round(credInt[1], digits=3), ";", round(credInt[2], digits=3), "]"),
                                   round(2*pmax(0.5/dim(ModelList[[i]]$Sol)[1], pmin(colSums(ModelList[[i]]$Sol[,j, drop = FALSE] > 0)/dim(ModelList[[i]]$Sol)[1], 
                                                                                     1 - colSums(ModelList[[i]]$Sol[, j, drop = FALSE] > 0)/dim(ModelList[[i]]$Sol)[1])), digits=3))
      counter <- counter+1
    }
  }
  return(DataFrame)
}
ForestTable <- function(ModelList, DataSetList, Group=c("Smaller", "Larger"), Picks=0){
  ModelData <- ModelDataPuller(ModelList, Group)
  Measures <- sapply(DataSetList, nrow)
  Species <- sapply(DataSetList, function(x) length(unique(x$Species)))
  if(Group[1] %in% c("Smaller", "Stable", "Closed")){
    if(Group[1] == "Smaller"){
      G1ind <- lapply(1:length(Picks), function(x) which(DataSetList[[x]]$Syllable.Rep < Picks[x]))
    }else if(Group[1] == "Stable"){
      G1ind <- lapply(1:length(Picks), function(x) which(DataSetList[[x]]$O.C == "Closed"))
    }
    G1M <- sapply(G1ind, length)
    G1S <- sapply(1:length(G1ind), function(x) length(unique(DataSetList[[x]]$Species[G1ind[[x]]])))
    ModelData["#Species"] <-  c(rbind(G1S, Species-G1S))
    ModelData["#Measure"] <-  c(rbind(G1M, Measures-G1M))
  }else{
    ModelData["#Species"] <- Species
    ModelData["#Measure"] <- Measures
  }
  ModelData <- ModelData[,c(1,5:6,2:4)]
  if(Picks[1] != 0){
    ModelData["Threshold"] <- paste0(c("<", ">="),rep(Picks, each=2))
    ModelData <- ModelData[,c(7, 1:6)]
  }
  return(ModelData)
}
ConvergePlot <- function(model, Name, misc="", Type="Pop"){
  pdf(paste0("Convergence-",Name,"-", misc,".pdf"))
  par(mfrow=c(4,4), mgp=c(1.5,.5,0), mar=c(3,3,2,3))
  if(Type=="Pop"){
    CheckPlot(model$Sol[,1], "Population", 1)
    CheckPlot(model$VCV[,1], "Measure Type", 2)
    CheckPlot(model$VCV[,2], "Species Effect", 3)
    CheckPlot(model$VCV[,3], "Phylogenic Effects", 4)
    CheckPlot(model$VCV[,4], "Study", 5)
    CheckPlot(model$VCV[,5], "Standard Error", 6, SE=TRUE)
    CheckPlot(model$VCV[,6], "Units", 7)
  }else if(Type=="Cont"){
    CheckPlot(model$Sol[,1], "Intercept", 1)
    CheckPlot(model$Sol[,2], "Repertoire", 2)
    CheckPlot(model$VCV[,1], "Measure Type", 3)
    CheckPlot(model$VCV[,2], "Species Effect", 4)
    CheckPlot(model$VCV[,3], "Phylogenic Effects", 5)
    CheckPlot(model$VCV[,4], "Study", 6)
    CheckPlot(model$VCV[,5], "Standard Error", 7, SE=TRUE)
    CheckPlot(model$VCV[,6], "Units", 8)
  }else if(Type=="Rep"){
    CheckPlot(model$Sol[,1], "Smaller Rep", 1)
    CheckPlot(model$Sol[,2], "Larger Rep", 2)
    CheckPlot(model$VCV[,1], "Measure Type", 3)
    CheckPlot(model$VCV[,2], "Species Effect", 4)
    CheckPlot(model$VCV[,3], "Phylogenic Effects", 5)
    CheckPlot(model$VCV[,4], "Study", 6)
    CheckPlot(model$VCV[,5], "Standard Error", 7, SE=TRUE)
    CheckPlot(model$VCV[,6], "Units", 8)
  }else if(Type=="Mixed"){
    CheckPlot(model$Sol[,1], "Smaller Stable", 1)
    CheckPlot(model$Sol[,2], "Smaller Plastic", 2)
    CheckPlot(model$Sol[,3], "Larger Stable", 3)
    CheckPlot(model$Sol[,4], "Larger Plastic", 4)
    CheckPlot(model$VCV[,1], "Measure Type", 5)
    CheckPlot(model$VCV[,2], "Species Effect", 6)
    CheckPlot(model$VCV[,3], "Phylogenic Effects", 7)
    CheckPlot(model$VCV[,4], "Study", 8)
    CheckPlot(model$VCV[,5], "Standard Error", 9, SE=TRUE)
    CheckPlot(model$VCV[,6], "Units", 10)
  }else{
    CheckPlot(model$Sol[,1], "Song Stable", 1)
    CheckPlot(model$Sol[,2], "Song Plastic", 2)
    CheckPlot(model$VCV[,1], "Measure Type", 3)
    CheckPlot(model$VCV[,2], "Species Effect", 4)
    CheckPlot(model$VCV[,3], "Phylogenic Effects", 5)
    CheckPlot(model$VCV[,4], "Study", 6)
    CheckPlot(model$VCV[,5], "Standard Error", 7, SE=TRUE)
    CheckPlot(model$VCV[,6], "Units", 8)
  }
  dev.off() 
}

Tablerfor2Vars <- function(Data, Dataset){
  G1 <- which(Dataset$Syllable.Rep < Picks[i] & Dataset$O.C == "Closed")
  G2 <- which(Dataset$Syllable.Rep < Picks[i] & Dataset$O.C == "Open")
  G3 <- which(Dataset$Syllable.Rep >= Picks[i] & Dataset$O.C == "Closed")
  G4 <- which(Dataset$Syllable.Rep >= Picks[i] & Dataset$O.C == "Open")
  Data["#Measure"] <- sapply(list(G1,G2,G3,G4), length)
  Data["#Species"] <- sapply(list(G1,G2,G3,G4), function(x) length(unique(Dataset$Species[x])))
  Data["Threshold"] <- paste0(c(rep("<",2),rep(">=",2)), Picks[i])
  Data <- Data[,c(7,1:6)]
  return(Data)
}
BESTPrinterOC <- function(Dataset){
  set.seed=49
  SL <- BESTmcmc(Dataset$Zs[which(Dataset$O.C == "Open")],
                 Dataset$Zs[which(Dataset$O.C == "Closed")],
                 numSavedSteps=16998, rnd.seed=49)
  paramSampleVec <- SL$mu1-SL$mu2
  by <- diff(hdi(paramSampleVec))/18
  breaks <- unique(c(seq(from = min(paramSampleVec), 
                         to = max(paramSampleVec), by = by), max(paramSampleVec)))
  histinfo <- hist(paramSampleVec, breaks = breaks, plot = FALSE)
  cvHt <- 0.7 * max(histinfo$density)
  pcgtCompVal <- round(100 * sum(paramSampleVec < 0)/length(paramSampleVec), 
                       1)
  if(pcgtCompVal == 0){
    pcgtCompVal <- "0.1"
  }
  HDI <- hdi(paramSampleVec, .95)
  
  data <- cbind(round(mean(HDI), digits=3),
                         paste0("[", round(HDI[1], digits=3), ";", round(HDI[2], digits=3), "]"),
                         pcgtCompVal)
  colnames(data) <- c("BEST Mean", "95% CredInt", "%<0")
  return(data)
}


##Code form Liam Revel that was previously in phytools, and then got deleteted.

