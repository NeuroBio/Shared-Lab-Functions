x <- cbind(rnorm(25),rnorm(5),rnorm(25))
data <- x
#data <- x[,1]
ScatterBox <- function(data,main="",sub="",xlab="",ylab="",labels=FALSE,col1=NA,col2=rgb(0,0,0,.7)){
  #make data matrix, creat empty plot, get quantiles and get offset
  data <- as.matrix(data)
  offset <- 1:ncol(data)-1
  plot(0,type='n',xlim=c(0,ncol(data)),ylim=c(min(data),max(data)),xaxt='n',
       xlab=xlab, ylab=ylab, main=main, sub=sub)
  axis(1,at=.5+offset,labels = labels)
  quan <- apply(data, 2,quantile, c(.02,.25,.5,.75,.98))
  #box
  rect(.25+offset,quan[2,],.5+offset,quan[4,],col = col1)
  #mean
  segments(.25+offset,quan[3,],.5+offset,quan[3,], lwd=2)
  
  #get data for whiskers, scatter, and outliers
  up <- vector("numeric",ncol(data))
  down <- vector("numeric",ncol(data))
  index <- vector("list",ncol(data))
  jity <- vector("list",ncol(data))
  for(i in 1:ncol(data)){
    up[i] <- min(max(data[,i]), quan[4,i] + 1.5*(quan[4,i]-quan[2,i]))
    down[i] <- max(min(data[,i]), quan[2,i] - 1.5*(quan[4,i]-quan[2,i]))
    index[[i]] <- which(data[,i] > up[i] | data[,i] < down[i])
    jity[[i]] <- jitter(rep(.625,length(data[,i])),8)+offset[i]
  }
  
  #whisckerlength and block
  segments(.5+offset,c(quan[4,],quan[2,]),.5+offset,c(up,down),lty = 2)
  segments(.5+offset,c(up,down),.375+offset,c(up,down),lwd = 1)
  
  #Scatter
  type <- rep(16,length(x))
  if(length(unlist(index) > 0)){
    type[unlist(index)] <- 1
  }
  if(length(col2)==ncol(data)){
    collist <- vector("list",ncol(data))
    for(i in 1:ncol(data)){
      collist[[i]] <- rep(col2[i],nrow(data))
    }
    col2 <- unlist(collist)
  }
  points(unlist(jity), data, pch=type, col=col2)
}
ScatterBox(data,"Awesome new function:","ScatterBox", "thing", "random",
           c("onedata","twodata","threedata"),c("pink","cyan","limegreen"),
           c(rgb(1,0,0,.7),rgb(0,0,1,.7),rgb(0,1,0,.7)) )

#rect(x,y, x,y)

# R script by Frank Soboczenski, PhD 
# email: frank.soboczenski@kcl.ac.uk
# twitter: h21k

# Small script that checks if the required packages are installed if not it will install them !
#plotting the plot :) 
ggplot(iris) +
  theme_stata() +
  theme(line = element_blank()) + 
  stat_boxplot(aes(x = Species, y = Sepal.Length), geom='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(aes(x = Species, y = Sepal.Length, fill=Species),
               alpha = 1, size = 0.75, width = 0.25, outlier.shape = 3) + 
  annotate("rect", xmin = 1, xmax = 1.5, ymin = 0, ymax = 8, alpha = 1, fill = 'white') +
  annotate("rect", xmin = 2, xmax = 2.5, ymin = 0, ymax = 8, alpha = 1, fill = 'white') +
  annotate("rect", xmin = 3, xmax = 3.5, ymin = 0, ymax = 8, alpha = 1, fill = 'white') +
  geom_point(aes(x = as.numeric(Species) + 0.1, colour = Species, y = Sepal.Length),
             alpha = 0.5, position = position_jitter(width = 0.1))

#Notes

#outlier.shape = 3 --> display outliers as crosses
#

