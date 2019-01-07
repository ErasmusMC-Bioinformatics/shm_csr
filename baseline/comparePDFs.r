options("warn"=-1)

#from http://selection.med.yale.edu/baseline/Archive/Baseline%20Version%201.3/Baseline_Functions_Version1.3.r
# Compute p-value of two distributions
compareTwoDistsFaster <-function(sigma_S=seq(-20,20,length.out=4001), N=10000, dens1=runif(4001,0,1), dens2=runif(4001,0,1)){
#print(c(length(dens1),length(dens2)))
if(length(dens1)>1 & length(dens2)>1 ){
	dens1<-dens1/sum(dens1)
	dens2<-dens2/sum(dens2)
	cum2 <- cumsum(dens2)-dens2/2
	tmp<- sum(sapply(1:length(dens1),function(i)return(dens1[i]*cum2[i])))
	#print(tmp)
	if(tmp>0.5)tmp<-tmp-1
	return( tmp )
	}
	else {
	return(NA)
	}
	#return (sum(sapply(1:N,function(i)(sample(sigma_S,1,prob=dens1)>sample(sigma_S,1,prob=dens2))))/N)
}  


require("grid")
arg <- commandArgs(TRUE)
#arg <- c("300143","4","5")
arg[!arg=="clonal"]
input <- arg[1]
output <- arg[2]
rowIDs <- as.numeric(  sapply(arg[3:(max(3,length(arg)))],function(x){ gsub("chkbx","",x) } )  )

numbSeqs = length(rowIDs)

if ( is.na(rowIDs[1]) | numbSeqs>10 ) {
  stop( paste("Error: Please select between one and 10 seqeunces to compare.") )
}

#load( paste("output/",sessionID,".RData",sep="") )
load( input )
#input

xMarks = seq(-20,20,length.out=4001)

plot_grid_s<-function(pdf1,pdf2,Sample=100,cex=1,xlim=NULL,xMarks = seq(-20,20,length.out=4001)){
  yMax = max(c(abs(as.numeric(unlist(listPDFs[pdf1]))),abs(as.numeric(unlist(listPDFs[pdf2]))),0),na.rm=T) * 1.1

  if(length(xlim==2)){
    xMin=xlim[1]
    xMax=xlim[2]
  } else {
    xMin_CDR = xMarks[listPDFs[pdf1][[1]][["CDR"]]>0.001][1]
    xMin_FWR = xMarks[listPDFs[pdf1][[1]][["FWR"]]>0.001][1]
    xMax_CDR = xMarks[listPDFs[pdf1][[1]][["CDR"]]>0.001][length(xMarks[listPDFs[pdf1][[1]][["CDR"]]>0.001])]
    xMax_FWR = xMarks[listPDFs[pdf1][[1]][["FWR"]]>0.001][length(xMarks[listPDFs[pdf1][[1]][["FWR"]]>0.001])]
  
    xMin_CDR2 = xMarks[listPDFs[pdf2][[1]][["CDR"]]>0.001][1]
    xMin_FWR2 = xMarks[listPDFs[pdf2][[1]][["FWR"]]>0.001][1]
    xMax_CDR2 = xMarks[listPDFs[pdf2][[1]][["CDR"]]>0.001][length(xMarks[listPDFs[pdf2][[1]][["CDR"]]>0.001])]
    xMax_FWR2 = xMarks[listPDFs[pdf2][[1]][["FWR"]]>0.001][length(xMarks[listPDFs[pdf2][[1]][["FWR"]]>0.001])]
  
    xMin=min(c(xMin_CDR,xMin_FWR,xMin_CDR2,xMin_FWR2,0),na.rm=TRUE)
    xMax=max(c(xMax_CDR,xMax_FWR,xMax_CDR2,xMax_FWR2,0),na.rm=TRUE)
  }

  sigma<-approx(xMarks,xout=seq(xMin,xMax,length.out=Sample))$x
  grid.rect(gp = gpar(col=gray(0.6),fill="white",cex=cex))
  x <- sigma
  pushViewport(viewport(x=0.175,y=0.175,width=0.825,height=0.825,just=c("left","bottom"),default.units="npc"))
  #pushViewport(plotViewport(c(1.8, 1.8, 0.25, 0.25)*cex))
  pushViewport(dataViewport(x, c(yMax,-yMax),gp = gpar(cex=cex),extension=c(0.05)))
  grid.polygon(c(0,0,1,1),c(0,0.5,0.5,0),gp=gpar(col=grey(0.95),fill=grey(0.95)),default.units="npc")
  grid.polygon(c(0,0,1,1),c(1,0.5,0.5,1),gp=gpar(col=grey(0.9),fill=grey(0.9)),default.units="npc")
  grid.rect()
  grid.xaxis(gp = gpar(cex=cex/1.1))
  yticks = pretty(c(-yMax,yMax),8)
  yticks = yticks[yticks>(-yMax) & yticks<(yMax)]
  grid.yaxis(at=yticks,label=abs(yticks),gp = gpar(cex=cex/1.1))
  if(length(listPDFs[pdf1][[1]][["CDR"]])>1){
    ycdr<-approx(xMarks,listPDFs[pdf1][[1]][["CDR"]],xout=seq(xMin,xMax,length.out=Sample),yleft=0,yright=0)$y
    grid.lines(unit(x,"native"), unit(ycdr,"native"),gp=gpar(col=2,lwd=2))
  }
  if(length(listPDFs[pdf1][[1]][["FWR"]])>1){
    yfwr<-approx(xMarks,listPDFs[pdf1][[1]][["FWR"]],xout=seq(xMin,xMax,length.out=Sample),yleft=0,yright=0)$y
    grid.lines(unit(x,"native"), unit(-yfwr,"native"),gp=gpar(col=4,lwd=2))
   }

  if(length(listPDFs[pdf2][[1]][["CDR"]])>1){
    ycdr2<-approx(xMarks,listPDFs[pdf2][[1]][["CDR"]],xout=seq(xMin,xMax,length.out=Sample),yleft=0,yright=0)$y
    grid.lines(unit(x,"native"), unit(ycdr2,"native"),gp=gpar(col=2,lwd=2,lty=2))
  }
  if(length(listPDFs[pdf2][[1]][["FWR"]])>1){
    yfwr2<-approx(xMarks,listPDFs[pdf2][[1]][["FWR"]],xout=seq(xMin,xMax,length.out=Sample),yleft=0,yright=0)$y
    grid.lines(unit(x,"native"), unit(-yfwr2,"native"),gp=gpar(col=4,lwd=2,lty=2))
   }

  grid.lines(unit(c(0,1),"npc"), unit(c(0.5,0.5),"npc"),gp=gpar(col=1))
  grid.lines(unit(c(0,0),"native"), unit(c(0,1),"npc"),gp=gpar(col=1,lwd=1,lty=3))

  grid.text("All", x = unit(-2.5, "lines"), rot = 90,gp = gpar(cex=cex))
  grid.text( expression(paste("Selection Strength (", Sigma, ")", sep="")) , y = unit(-2.5, "lines"),gp = gpar(cex=cex))
  
  if(pdf1==pdf2 & length(listPDFs[pdf2][[1]][["FWR"]])>1 & length(listPDFs[pdf2][[1]][["CDR"]])>1 ){
    pCDRFWR = compareTwoDistsFaster(sigma_S=xMarks, N=10000, dens1=listPDFs[[pdf1]][["CDR"]], dens2=listPDFs[[pdf1]][["FWR"]])       
    pval = formatC(as.numeric(pCDRFWR),digits=3)
    grid.text( substitute(expression(paste(P[CDR/FWR], "=", x, sep="")),list(x=pval))[[2]] , x = unit(0.02, "npc"),y = unit(0.98, "npc"),just=c("left", "top"),gp = gpar(cex=cex*1.2))
  }
  grid.text(paste("CDR"), x = unit(0.98, "npc"),y = unit(0.98, "npc"),just=c("right", "top"),gp = gpar(cex=cex*1.5))
  grid.text(paste("FWR"), x = unit(0.98, "npc"),y = unit(0.02, "npc"),just=c("right", "bottom"),gp = gpar(cex=cex*1.5))
  popViewport(2)
}
#plot_grid_s(1)


p2col<-function(p=0.01){
  breaks=c(-.51,-0.1,-.05,-0.01,-0.005,0,0.005,0.01,0.05,0.1,0.51)
  i<-findInterval(p,breaks)
  cols = c( rgb(0.8,1,0.8), rgb(0.6,1,0.6), rgb(0.4,1,0.4), rgb(0.2,1,0.2) , rgb(0,1,0),
            rgb(1,0,0), rgb(1,.2,.2), rgb(1,.4,.4), rgb(1,.6,.6) , rgb(1,.8,.8) )
  return(cols[i])
}


plot_pvals<-function(pdf1,pdf2,cex=1,upper=TRUE){
  if(upper){
    pCDR1FWR2 = compareTwoDistsFaster(sigma_S=xMarks, N=10000, dens1=listPDFs[[pdf1]][["CDR"]], dens2=listPDFs[[pdf2]][["FWR"]])       
    pFWR1FWR2 = compareTwoDistsFaster(sigma_S=xMarks, N=10000, dens1=listPDFs[[pdf1]][["FWR"]], dens2=listPDFs[[pdf2]][["FWR"]])
    pFWR1CDR2 = compareTwoDistsFaster(sigma_S=xMarks, N=10000, dens2=listPDFs[[pdf2]][["CDR"]], dens1=listPDFs[[pdf1]][["FWR"]])       
    pCDR1CDR2 = compareTwoDistsFaster(sigma_S=xMarks, N=10000, dens2=listPDFs[[pdf2]][["CDR"]], dens1=listPDFs[[pdf1]][["CDR"]])
    grid.polygon(c(0.5,0.5,1,1),c(0,0.5,0.5,0),gp=gpar(col=p2col(pFWR1FWR2),fill=p2col(pFWR1FWR2)),default.units="npc")
    grid.polygon(c(0.5,0.5,1,1),c(1,0.5,0.5,1),gp=gpar(col=p2col(pCDR1FWR2),fill=p2col(pCDR1FWR2)),default.units="npc")
    grid.polygon(c(0.5,0.5,0,0),c(1,0.5,0.5,1),gp=gpar(col=p2col(pCDR1CDR2),fill=p2col(pCDR1CDR2)),default.units="npc")
    grid.polygon(c(0.5,0.5,0,0),c(0,0.5,0.5,0),gp=gpar(col=p2col(pFWR1CDR2),fill=p2col(pFWR1CDR2)),default.units="npc")
         
    grid.lines(c(0,1),0.5,gp=gpar(lty=2,col=gray(0.925)))
    grid.lines(0.5,c(0,1),gp=gpar(lty=2,col=gray(0.925)))

    grid.text(formatC(as.numeric(pFWR1FWR2),digits=3), x = unit(0.75, "npc"),y = unit(0.25, "npc"),just=c("center", "center"),gp = gpar(cex=cex))
    grid.text(formatC(as.numeric(pCDR1FWR2),digits=3), x = unit(0.75, "npc"),y = unit(0.75, "npc"),just=c("center", "center"),gp = gpar(cex=cex))
    grid.text(formatC(as.numeric(pCDR1CDR2),digits=3), x = unit(0.25, "npc"),y = unit(0.75, "npc"),just=c("center", "center"),gp = gpar(cex=cex))
    grid.text(formatC(as.numeric(pFWR1CDR2),digits=3), x = unit(0.25, "npc"),y = unit(0.25, "npc"),just=c("center", "center"),gp = gpar(cex=cex))
    
           
 #   grid.text(paste("P = ",formatC(pCDRFWR,digits=3)), x = unit(0.5, "npc"),y = unit(0.98, "npc"),just=c("center", "top"),gp = gpar(cex=cex))
 #   grid.text(paste("P = ",formatC(pFWRFWR,digits=3)), x = unit(0.5, "npc"),y = unit(0.02, "npc"),just=c("center", "bottom"),gp = gpar(cex=cex))
  }
  else{
  }
}


##################################################################################
################## The whole OCD's matrix ########################################
##################################################################################

#pdf(width=4*numbSeqs+1/3,height=4*numbSeqs+1/3)
pdf( output ,width=4*numbSeqs+1/3,height=4*numbSeqs+1/3) 

pushViewport(viewport(x=0.02,y=0.02,just = c("left", "bottom"),w =0.96,height=0.96,layout = grid.layout(numbSeqs+1,numbSeqs+1,widths=unit.c(unit(rep(1,numbSeqs),"null"),unit(4,"lines")),heights=unit.c(unit(4,"lines"),unit(rep(1,numbSeqs),"null")))))

for( seqOne in 1:numbSeqs+1){
  pushViewport(viewport(layout.pos.col = seqOne-1, layout.pos.row = 1))
  if(seqOne>2){ 
    grid.polygon(c(0,0,0.5,0.5),c(0,0.5,0.5,0),gp=gpar(col=grey(0.5),fill=grey(0.9)),default.units="npc")
    grid.polygon(c(1,1,0.5,0.5),c(0,0.5,0.5,0),gp=gpar(col=grey(0.5),fill=grey(0.95)),default.units="npc")
    grid.polygon(c(0,0,1,1),c(1,0.5,0.5,1),gp=gpar(col=grey(0.5)),default.units="npc")
       
    grid.text(y=.25,x=0.75,"FWR",gp = gpar(cex=1.5),just="center")
    grid.text(y=.25,x=0.25,"CDR",gp = gpar(cex=1.5),just="center")
  }
  grid.rect(gp = gpar(col=grey(0.9)))
  grid.text(y=.75,substr(paste(names(listPDFs)[rowIDs[seqOne-1]]),1,16),gp = gpar(cex=2),just="center")
  popViewport(1)
}

for( seqOne in 1:numbSeqs+1){
  pushViewport(viewport(layout.pos.row = seqOne, layout.pos.col = numbSeqs+1))
  if(seqOne<=numbSeqs){   
    grid.polygon(c(0,0.5,0.5,0),c(0,0,0.5,0.5),gp=gpar(col=grey(0.5),fill=grey(0.95)),default.units="npc")
    grid.polygon(c(0,0.5,0.5,0),c(1,1,0.5,0.5),gp=gpar(col=grey(0.5),fill=grey(0.9)),default.units="npc")
    grid.polygon(c(1,0.5,0.5,1),c(0,0,1,1),gp=gpar(col=grey(0.5)),default.units="npc")
    grid.text(x=.25,y=0.75,"CDR",gp = gpar(cex=1.5),just="center",rot=270)
    grid.text(x=.25,y=0.25,"FWR",gp = gpar(cex=1.5),just="center",rot=270)
  }
  grid.rect(gp = gpar(col=grey(0.9)))
  grid.text(x=0.75,substr(paste(names(listPDFs)[rowIDs[seqOne-1]]),1,16),gp = gpar(cex=2),rot=270,just="center")
  popViewport(1)
}

for( seqOne in 1:numbSeqs+1){
  for(seqTwo in 1:numbSeqs+1){
    pushViewport(viewport(layout.pos.col = seqTwo-1, layout.pos.row = seqOne))
    if(seqTwo>seqOne){
      plot_pvals(rowIDs[seqOne-1],rowIDs[seqTwo-1],cex=2)
      grid.rect()
    }    
    popViewport(1)
  }
}
   

xMin=0
xMax=0.01
for(pdf1 in rowIDs){
  xMin_CDR = xMarks[listPDFs[pdf1][[1]][["CDR"]]>0.001][1]
  xMin_FWR = xMarks[listPDFs[pdf1][[1]][["FWR"]]>0.001][1]
  xMax_CDR = xMarks[listPDFs[pdf1][[1]][["CDR"]]>0.001][length(xMarks[listPDFs[pdf1][[1]][["CDR"]]>0.001])]
  xMax_FWR = xMarks[listPDFs[pdf1][[1]][["FWR"]]>0.001][length(xMarks[listPDFs[pdf1][[1]][["FWR"]]>0.001])]
  xMin=min(c(xMin_CDR,xMin_FWR,xMin),na.rm=TRUE)
  xMax=max(c(xMax_CDR,xMax_FWR,xMax),na.rm=TRUE)
}



for(i in 1:numbSeqs+1){
  for(j in (i-1):numbSeqs){    
    pushViewport(viewport(layout.pos.col = i-1, layout.pos.row = j+1))
    grid.rect()
    plot_grid_s(rowIDs[i-1],rowIDs[j],cex=1)
    popViewport(1)
  }
}

dev.off() 

cat("Success", paste(rowIDs,collapse="_"),sep=":")

