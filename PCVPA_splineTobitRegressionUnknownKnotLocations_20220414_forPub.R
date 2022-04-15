##
## set-up
##

rm(list=ls())

library(rmarkdown)
library(knitr)
options(digits=2)
library(censReg) # for Tobin / censored regression models
library(rms) # for spline terms
library(boot) # to compute adjusted percentile method bootstrap CIs
library(tidyverse)
library(grid)
library(gridExtra)

st<-c("1","3","4","5","6A","6B","7F","9V","14","18C","19A","19F","23F","12F","33F") # 12F is almost exclusively censored data

wd<-paste(sep="","/Users/marc/work/MLW_LSTM/SwarthoutTodd_PCVPA/",Sys.Date(),"/") # to be changed by other users; use 'here' package instead?
if(!dir.exists(wd)){dir.create(wd)}

fileLog<-"pcvpaModelling.log"

setwd(wd)

seed<-20210121
inFile<-"../2019-06-22_data/serology_pneumo_igg_13vt_2nv2_3prot_22jun19.csv" # other users need same relative folder structure - simplify?

cat(file=fileLog,append=F,paste(sep="","This is PCVPA_splineTobitRegressionUnknownKnotLocations_20220414_forPub.R\n\n\tseed = < ",seed," >\n\tdate = < ",Sys.Date()," >\n\tinFile = < ",inFile," >\n\n"))



##
## helper functions
##

gm_mean = function(x, na.rm=TRUE){
  return(exp(sum(log(x[x > 0]), na.rm=na.rm) / sum(x>0)))
}

gm_mean_withCI = function(x, na.rm=TRUE, alpha=0.05){
  geomLog<-sum(log(x[x > 0]), na.rm=na.rm) / sum(!is.na(x) & x>0)
  geomLogSE<-sqrt( sum((log(x[x>0])-geomLog)^2, na.rm=T)/(sum(!is.na(x) & x>0)-1) ) / sqrt(sum(!is.na(x) & x>0))
  ciLog<-c(geomLog-qnorm(1-alpha/2)*geomLogSE,geomLog+qnorm(1-alpha/2)*geomLogSE)
  ci<-exp(ciLog)
  geom<-exp(geomLog)
  print(c(geomLog,geomLogSE,geomLog,ciLog))
  return(c(geom,ci))
}

knottySS<-function(pars,x=dat$ageMonths,y=log(dat[,varNameImpDL]),left=log(DL),k=3){
  # k is the number of knots, k>=1
  if(k<=0){stop("k needs to be an integer >0.")}

  # split parameters into regression coefficients and knot locations
  beta<-pars[1:(k+2)] # slope parameters
  sigma<-pars[k+3]
  gamma<-pars[(k+4):(2*k+3)] # knot locations

  res<-(-sum(censReg(y ~ lsp(x,gamma),left=left,logLikOnly=T,start=c(beta,sigma))))

  return(res)
}

bootFun<-function(data, indices, pars, xvar, yvar, left, k, initKnots){
  bootDat<-data[indices,]

  initCoefs<-coef(censReg(as.formula(paste(sep="","log(",yvar,") ~ lsp(",xvar,",initKnots)")),data=dat,left=left))

  bootMod<-optim(fn=knottySS,par=c(initCoefs,initKnots),x=bootDat[,xvar],y=log(bootDat[,yvar]),control=list(maxit=1e6),k=k,left=left)

  coefs<-bootMod$par[1:(k+2)]
  logSD<-bootMod$par[k+3]
  knots<-sort(bootMod$par[(k+4):(2*k+3)]) # occasionally the knots are in the wrong order; lsp() will order them internally, so slopes etc will not be affected by wrong order
  slopes<-cumsum(coefs[-1])

  cf<-c(coefs,logSD,knots)
  xx<-seq(0.1,60,length=250)
  logyy<-cf[1]+cf[2]*xx
  if(k>0){
    for(i in 1:k){
      logyy<-logyy+cf[2+i]*ifelse((xx-cf[2+k+1+i])>=0,xx-cf[2+k+1+i],0)
    }
  }
  yy<-exp(logyy)

  return(c(coefs,logSD,knots,slopes,yy))
}



##
## Read in the data and reformat
##

# set-up
cat(file=fileLog,append=T,paste(sep="","Reading file ",inFile,"\n\n"))
dat<-read.csv(inFile)
dat$ageMonths<-dat$age_yrs*12

levels(dat$age_cat_gmc)[levels(dat$age_cat_gmc)=="0-3m"]<-">00m"
levels(dat$age_cat_gmc)[levels(dat$age_cat_gmc)==">3-6m"]<-">03m"
levels(dat$age_cat_gmc)[levels(dat$age_cat_gmc)==">6-9m"]<-">06m"
levels(dat$age_cat_gmc)[levels(dat$age_cat_gmc)==">9m"]<-">09m"

cnames2extract<-c("labid","age_cat_gmc","ageMonths")
for(sero in c("1","3","4","5","6A","6B","7F","9V","12F","14","18C","19A","19F","23F","33F")){
  colnames(dat)[colnames(dat)==paste(sep="","st",sero,"_adj")]<-paste(sep="","res",sero,"_num_ImpHalfDL")
  idx.NRQNS<-which(is.element(el=as.character(dat[,paste(sep="","res",sero)]),set=c("NR","QNS")))
  idx.belowDL<-which(as.character(dat[,paste(sep="","res",sero)])=="<0.150") # identify left-censored observations
  idx.belowDL<-setdiff(idx.belowDL,idx.NRQNS)
  dat[,paste(sep="","res",sero,"_num")]<-as.numeric(as.character(dat[,paste(sep="","res",sero)]))
  dat[,paste(sep="","res",sero,"_num_ImpDL")]<-dat[,paste(sep="","res",sero,"_num")]
  dat[idx.belowDL,paste(sep="","res",sero,"_num_ImpDL")]<-0.15
  dat[,paste(sep="","status_",sero)]<-1
  dat[idx.belowDL,paste(sep="","status_",sero)]<-2
  dat[idx.NRQNS,paste(sep="","status_",sero)]<-NA
  cnames2extract<-c(cnames2extract,paste(sep="","res",sero),paste(sep="","res",sero,"_num"),paste(sep="","res",sero,"_num_ImpDL"),paste(sep="","res",sero,"_num_ImpHalfDL"),paste(sep="","status_",sero))
}

write.csv(dat,file="PCVPA_serology_reformatted.csv",row.names=F,quote=F)



##
## Main modelling / analysis R routine
##

doAnalysis<-function(dat,varNameOrig,varNameOrigNum,varNameImpHalfDL,varNameImpDL,DL=0.15,Nboot=1e4,seed=1234,logFile=fileLog){
  set.seed(seed)

  # filter data to only recorded values
  dat<-dat[!is.na(dat[,varNameImpHalfDL]),] # remove non-recorded values

  # compute geometric means
  cat(file=logFile,append=T,"geoM\n")
  datGeoM<-data.frame(ageCat=factor(unique(dat$age_cat_gmc)),geoM=NA,geoM_better=NA)
  for(i in 1:nrow(datGeoM)){
    datGeoM$geoM[i]<-gm_mean(dat[dat$age_cat_gmc==datGeoM$ageCat[i],varNameImpHalfDL]) # note that these use the simple DL/2 imputation for calculation
    if(sum(dat[dat$age_cat_gmc==datGeoM$ageCat[i],varNameOrig]=="<0.150")>0){
      if(sum(dat[dat$age_cat_gmc==datGeoM$ageCat[i],varNameOrig]=="<0.150")==sum(dat$age_cat_gmc==datGeoM$ageCat[i])){
        datGeoM$geoM_better[i]<-gm_mean(dat[dat$age_cat_gmc==datGeoM$ageCat[i],varNameImpDL])
      }else{
        datGeoM$geoM_better[i]<-exp(coef(censReg(as.formula(paste(sep="","log(",varNameImpDL,")~1")),data=dat[dat$age_cat_gmc==datGeoM$ageCat[i],],left=log(DL)))["(Intercept)"])
      }
    }else{
      datGeoM$geoM_better[i]<-exp(coef(lm(as.formula(paste(sep="","log(",varNameOrigNum,")~1")),data=dat[dat$age_cat_gmc==datGeoM$ageCat[i],]))["(Intercept)"])
    }
  }
  levels(datGeoM$ageCat)[levels(datGeoM$ageCat)=="0-2m"]<-"0-02m"
  levels(datGeoM$ageCat)[levels(datGeoM$ageCat)=="3-5m"]<-"03-05m"
  levels(datGeoM$ageCat)[levels(datGeoM$ageCat)=="6-8m"]<-"06-08m"
  levels(datGeoM$ageCat)[levels(datGeoM$ageCat)=="9-11m"]<-"09-11m"

  datGeoM<-datGeoM[order(as.character(datGeoM$ageCat)),]
  datGeoM$ageCatNum<-seq(1.5,58.5,by=3)

  # fit model and select best-fitting number of knots (among 0, 1, 2, 3)
  cat(file=logFile,append=T,"fit\n")
  modList<-list()
  modListKnots<-list()
  aicVect<-rep(NA,4)
  initKnots<-list(3,c(3,9),c(3,9,30))

  cat(file=logFile,append=T,"..k=0\n")
  k<-0
  modList[[k+1]]<-censReg(as.formula(paste(sep="","log(",varNameImpDL,") ~ ageMonths")),data=dat,left=log(DL))
  modListKnots[[k+1]]<-NA
  aicVect[k+1]<-AIC(modList[[k+1]])

  cat(file=logFile,append=T,"..k=1\n")
  k<-1
  modTmp<-censReg(as.formula(paste(sep="","log(",varNameImpDL,") ~ lsp(ageMonths,c(3))")),data=dat,left=log(DL)) # just to get good initial parameter values
  modListKnots[[k+1]]<-optim(fn=knottySS,par=as.numeric(c(coef(modTmp),initKnots[[k]])),control=list(maxit=1e6),k=k,left=log(DL),y=log(dat[,varNameImpDL]),x=dat$ageMonths)
  modList[[k+1]]<-censReg(as.formula(paste(sep="","log(",varNameImpDL,") ~ lsp(ageMonths,modListKnots[[k+1]]$par[k+3+1:k])")),data=dat,left=log(DL)) # we refit this, retaining the knot locations; get a slightly better fit...
  aicVect[k+1]<-2*(-summary(modList[[k+1]])$loglik)+2*length(modListKnots[[k+1]]$par)

  cat(file=logFile,append=T,"..k=2\n")
  k<-2
  modTmp<-censReg(as.formula(paste(sep="","log(",varNameImpDL,") ~ lsp(ageMonths,c(3,9))")),data=dat,left=log(DL)) # just to get good initial parameter values
  modListKnots[[k+1]]<-optim(fn=knottySS,par=as.numeric(c(coef(modTmp),initKnots[[k]])),control=list(maxit=1e6),k=k,left=log(DL),y=log(dat[,varNameImpDL]),x=dat$ageMonths)
  modList[[k+1]]<-censReg(as.formula(paste(sep="","log(",varNameImpDL,") ~ lsp(ageMonths,modListKnots[[k+1]]$par[k+3+1:k])")),data=dat,left=log(DL)) # we refit this, retaining the knot locations; get a slightly better fit...
  aicVect[k+1]<-2*(-summary(modList[[k+1]])$loglik)+2*length(modListKnots[[k+1]]$par)

  cat(file=logFile,append=T,"..k=3\n")
  k<-3
  modTmp<-censReg(as.formula(paste(sep="","log(",varNameImpDL,") ~ lsp(ageMonths,c(3,9,30))")),data=dat,left=log(DL)) # just to get good initial parameter values
  modListKnots[[k+1]]<-optim(fn=knottySS,par=as.numeric(c(coef(modTmp),initKnots[[k]])),control=list(maxit=1e6),k=k,left=log(DL),y=log(dat[,varNameImpDL]),x=dat$ageMonths)
  modList[[k+1]]<-censReg(as.formula(paste(sep="","log(",varNameImpDL,") ~ lsp(ageMonths,modListKnots[[k+1]]$par[k+3+1:k])")),data=dat,left=log(DL)) # we refit this, retaining the knot locations; get a slightly better fit...
  aicVect[k+1]<-2*(-summary(modList[[k+1]])$loglik)+2*length(modListKnots[[k+1]]$par)

  cat(file=logFile,append=T,"AIC model selection\n")
  k<-(which(aicVect==min(aicVect)))-1
  cat(file=logFile,append=T,paste(sep="",paste(collapse=" ",round(digits=2,aicVect)),"\n"))
  cat(file=logFile,append=T,paste(sep="","k = ",k,"\n"))
  if(k==0){
    cf<-coef(modK0)
  }else{
    cf<-c(as.numeric(coef(modList[[k+1]],logSigma=F)),modListKnots[[k+1]]$par[k+3+1:k])
  }

  # compute predicted values for visualisation
  cat(file=logFile,append=T,"predVals\n")
  xx<-seq(0,60,length=1000)
  logyy<-cf[1]+cf[2]*xx
  if(k>0){
    for(i in 1:k){
      logyy<-logyy+cf[2+i]*ifelse((xx-cf[2+k+1+i])>=0,xx-cf[2+k+1+i],0)
    }
  }
  yy<-exp(logyy)

  # bootstrap the whole thing to get confidence intervals and bands
  cat(file=logFile,append=T,"bootstrap\n")
  if(k==0){
    ## XXX to be done but not needed for PCVPA dataset
  }else{
    bootRuns<-boot(data=dat,statistic=bootFun,R=Nboot,xvar="ageMonths",yvar=varNameImpDL,left=log(DL),k=k,initKnots=initKnots[[k]],pars=as.numeric(c(coef(modList[[k+1]]),initKnots[[k]])))

    knotCIs<-vector("list",k)
    slopeCIs<-vector("list",k+1)

    xxCI<-seq(0,60,length=250)
    yyLow<-rep(NA,length=length(xxCI))
    yyHigh<-rep(NA,length=length(xxCI))

    for(i in 1:length(knotCIs)){
      knotCIs[[i]]<-boot.ci(bootRuns, type="perc", index=2+k+1+i) # changed this to percentile method from bca method
    }

    for(i in 1:length(slopeCIs)){
      slopeCIs[[i]]<-boot.ci(bootRuns, type="perc", index=2+2*k+1+i) # changed this to percentile method from bca method
    }

    for(i in 1:length(xxCI)){
      tmp<-boot.ci(bootRuns, type="perc", index=2+3*k+2+i) # changed this to percentile method from bca method
      yyLow[i]<-tmp$perc[4]
      yyHigh[i]<-tmp$perc[5]
    }
  }

  # output all of this
  cat(file=logFile,append=T,"ouput\n")
  res<-list(k=k,mod=modList[[k+1]],modList=modList,modListKnots=modListKnots,aicVect=aicVect,coefs=cf,xx=xx,yy=yy,xxCI=xxCI,yyLow=yyLow,yyHigh=yyHigh,knotCIs=knotCIs,slopeCIs=slopeCIs,Nboot=Nboot,dat=dat,datGeoM=datGeoM,bootRuns=bootRuns,seed=seed)
  return(res)
}

plotFun<-function(fit,strain,varNameImpDL){
par(mar=c(5, 4, 5, 2) + 0.1)
  plot(fit$datGeoM$ageCatNum,fit$datGeoM$geoM_better,log="y",xlab="age (months)",ylab=paste(sep="",strain," IgG concentration"),pch=20,ylim=c(0.15,max(c(14,max(fit$dat[[varNameImpDL]])))),cex=2)
  points(fit$dat$ageMonths,fit$dat[[varNameImpDL]],col=rgb(red=0,green=0,blue=0,alpha=50,maxColorValue=255))
  lines(fit$xx,fit$yy,col="red",lwd=2)
  abline(h=0.15,col="grey",lty=2)
  text(x=55,y=0.155,"detection limit",col="grey",adj=c(0,0))
  if(fit$k>0){
    for(i in 1:fit$k){
      abline(v=fit$coefs[2+fit$k+1+i],col="grey",lty=2)
      text(x=fit$coefs[2+fit$k+1+i]+0.1,y=14,"change point",col="grey",adj=c(0,0))
    }
  }
  polygon(x=c(fit$xxCI,fit$xxCI[length(fit$xxCI):1]),y=c(fit$yyLow,fit$yyHigh[length(fit$xxCI):1]),col=rgb(red=255,green=0,blue=0,alpha=50,maxColorValue=255),border=NA)
  if(fit$k>0){
   for(i in 1:fit$k){
     polygon(x=c(rep(fit$knotCIs[[i]]$perc[4],2),rep(fit$knotCIs[[i]]$perc[5],2)),y=c(exp(log(9.5)+0.11*(i-1)),exp(log(10.5)+0.11*(i-1)),exp(log(10.5)+0.11*(i-1)),exp(log(9.5)+0.11*(i-1))),col=rgb(red=100,green=100,blue=100,alpha=100,maxColorValue=255),border=NA)
   }
  }
  par(xpd=T)
  lgd <- legend("topleft", legend = c("data", "geometric mean per age band"), pch = c(1,20), bty="n", col=c(rgb(red=0,green=0,blue=0,alpha=50,maxColorValue=255),"black"),pt.cex=c(1,2),inset=c(0,-0.12))
  legend(lgd$rect$left+lgd$rect$w, 10^(lgd$rect$top), legend = c("model fit", "95% CI (model fit)","95% CI (change points)"), lwd=c(2,5,5), bty="n",col=c("red",rgb(red=255,green=190,blue=190,maxColorValue=255),"darkgrey")) # can't have transparent colors in the legend when in the margin; hence approximating the transparent red and gray...
  par(xpd=F)
}

ggPlotFun<-function(fit,strain,varNameImpDL,DL=0.15,colModel="red",colModelCI="red",colData="grey50",colGeoM="black",colChangepoints="grey70",alphaModelCI=50,alphaChangepoints=200,alphaData=100,cexData=2,cexGeoM=5,doLegend=TRUE,xLegPoints=0,xLegRibbons=c(30,33),yLegMin=20,plotModel=TRUE,colAnnot="black",xAnnot=60,cexAnnot=16,cexText=12,doGeoMConfInt=FALSE,axisLabels=TRUE){
  yMin<-min(c(DL,fit$datGeoM$geoM_better))
  yMax<-max(fit$dat[,varNameImpDL],na.rm=TRUE)

  if(plotModel){
    dfCI<-data.frame(x=fit$xxCI,ymin=fit$yyLow,ymax=fit$yyHigh)
    dfFit<-data.frame(x=fit$xx,y=fit$yy)
    dfKnots<-NULL
    dfKnotsCI<-NULL
    if(fit$k>0){
      dfKnots<-data.frame(x=fit$coefs[2+fit$k+1+1:fit$k])
      dfKnotsCI<-data.frame(x=numeric(),y=numeric(),k=integer())
      for(i in 1:fit$k){
        dfKnotsCI<-rbind(dfKnotsCI,data.frame(x=fit$knotCIs[[i]]$perc[4:5],ymin=rep(10^(log10((1+(i-1)*0.11)*0.75*yMax)),2),ymax=rep(10^(log10((1.1+(i-1)*0.11)*0.75*yMax)),2),k=i))
      }
    }
  }

  if(doGeoMConfInt){
    varNameOrig<-gsub(varNameImpDL,pattern="_[a-zA-Z0-9]*",replacement="")
    varNameOrigNum<-paste(sep="_",varNameOrig,"num")

    fit$dat<-fit$dat %>% mutate(
      geoMLow=NA,
      geoMUpp=NA,
      age_cat_gmc=case_when(
        as.character(age_cat_gmc)=="0-2m"~"0-02m",
        as.character(age_cat_gmc)=="3-5m"~"03-05m",
        as.character(age_cat_gmc)=="6-8m"~"06-08m",
        as.character(age_cat_gmc)=="9-11m"~"09-11m",
        TRUE~as.character(age_cat_gmc)
      )
    )

    for(i in 1:nrow(fit$datGeoM)){
      if(sum(fit$dat[fit$dat$age_cat_gmc==fit$datGeoM$ageCat[i],varNameOrig]=="<0.150")>0){
        if(sum(fit$dat[fit$dat$age_cat_gmc==fit$datGeoM$ageCat[i],varNameOrig]=="<0.150")==sum(fit$dat$age_cat_gmc==fit$datGeoM$ageCat[i])){
          tmp<-rep(DL,2)
        }else{
          tmp<-exp(confint(censReg(as.formula(paste(sep="","log(",varNameImpDL,")~1")),data=fit$dat[fit$dat$age_cat_gmc==fit$datGeoM$ageCat[i],],left=log(DL)))["(Intercept)",])
        }
      }else{
        tmp<-exp(confint(lm(as.formula(paste(sep="","log(",varNameOrigNum,")~1")),data=fit$dat[fit$dat$age_cat_gmc==fit$datGeoM$ageCat[i],]))["(Intercept)",])
      }
      fit$datGeoM$geoMLow[i]<-tmp[1]
      fit$datGeoM$geoMUpp[i]<-tmp[2]
    }
    geoM<-data.frame(age=fit$datGeoM$ageCat,geoM=paste(sep="",format(nsmall=2,round(digits=2,fit$datGeoM$geoM_better))," (",format(nsmall=2,round(digits=2,fit$datGeoM$geoMLow)),",",format(nsmall=2,round(digits=2,fit$datGeoM$geoMUpp)),")"))
  }else{
    geoM<-NULL
  }

  dfText<-data.frame(
    x=c(rep(xLegPoints+2,2),rep(xLegRibbons[2]+2,3)),
    y=c(yLegMin*1.5*1.5,yLegMin*1.5,yLegMin*1.5*1.5,yLegMin*1.5,yLegMin),
    lab=c("data","geometric mean per age band","model fit","95% CI (model fit)","95% CI (change points)")
  )
  if(!plotModel){dfText<-dfText[1:2,]}

  if(doLegend){
    mar<-c(5,0.5,0.5,0.5)
    clipOnOff<-"off"
  }else{
    mar<-rep(0.5,4)
    clipOnOff<-"on"
  }

  rgbTmp<-col2rgb(colData)
  rgbData<-rgb(red=rgbTmp[1],green=rgbTmp[2],blue=rgbTmp[3],alpha=alphaData,maxColorValue=255)
  rgbTmp<-col2rgb(colGeoM)
  rgbGeoM<-rgb(red=rgbTmp[1],green=rgbTmp[2],blue=rgbTmp[3],alpha=255,maxColorValue=255)
  rgbTmp<-col2rgb(colModel)
  rgbModel<-rgb(red=rgbTmp[1],green=rgbTmp[2],blue=rgbTmp[3],alpha=255,maxColorValue=255)
  rgbTmp<-col2rgb(colModelCI)
  rgbModelCI<-rgb(red=rgbTmp[1],green=rgbTmp[2],blue=rgbTmp[3],alpha=alphaModelCI,maxColorValue=255)
  rgbTmp<-col2rgb(colChangepoints)
  rgbChangepoints<-rgb(red=rgbTmp[1],green=rgbTmp[2],blue=rgbTmp[3],alpha=alphaChangepoints,maxColorValue=255)

  g<-ggplot() +
    geom_point(data=fit$dat,mapping=aes(x=ageMonths,y=get(varNameImpDL)),col=rgbData,pch=1,size=cexData) +
    {if(plotModel) geom_line(data=dfFit,mapping=aes(x=x,y=y),col=rgbModel,lwd=1.5)} +
    geom_point(data=fit$datGeoM,mapping=aes(x=ageCatNum,y=geoM_better),pch=20,size=cexGeoM,col=rgbGeoM) +
    {if(plotModel) geom_ribbon(data=dfCI,mapping=aes(x=x,ymin=ymin,ymax=ymax),fill=rgbModelCI)} +
    {if(plotModel) geom_vline(data=dfKnots,mapping=aes(xintercept=x),lty=2,lwd=0.75,col=rgbChangepoints)} +
    {if(plotModel) geom_ribbon(data=dfKnotsCI,mapping=aes(x=x,ymin=ymin,ymax=ymax,group=k),fill=rgbChangepoints)} +
    theme_light() +
    theme(plot.margin = unit(mar, "lines"),text=element_text(size=cexText)) +
    scale_y_log10() +
    coord_cartesian(ylim = c(yMin,yMax),clip=clipOnOff) +
    {if(axisLabels) labs(y=paste(sep="",strain," IgG concentration (ug/mL)"),x="age (months)")} +
    {if(!axisLabels) labs(y=NULL,x=NULL)} +
    {if(doLegend) geom_point(data=data.frame(x=xLegPoints,y=30),mapping=aes(x=x,y=y),pch=20,size=cexGeoM,col=rgbGeoM)} +
    {if(doLegend) geom_point(data=data.frame(x=xLegPoints,y=45),mapping=aes(x=x,y=y),pch=1,size=cexData,col=rgbData)} +
    {if(doLegend & plotModel) geom_ribbon(data=data.frame(x=xLegRibbons,ymin=rep(yLegMin*1.5*1.5,2),ymax=rep(yLegMin*1.5*1.5*1.1,2)),mapping=aes(x=x,ymin=ymin,ymax=ymax),fill=rgbModel)} +
    {if(doLegend & plotModel) geom_ribbon(data=data.frame(x=xLegRibbons,ymin=rep(yLegMin*1.5,2),ymax=rep(yLegMin*1.5*1.1,2)),mapping=aes(x=x,ymin=ymin,ymax=ymax),fill=rgbModelCI)} +
    {if(doLegend & plotModel) geom_ribbon(data=data.frame(x=xLegRibbons,ymin=rep(yLegMin,2),ymax=rep(yLegMin*1.1,2)),mapping=aes(x=x,ymin=ymin,ymax=ymax),fill=rgbChangepoints)} +
    {if(doLegend) geom_text(dat=dfText,mapping=aes(x=x,y=y,label=lab),hjust=0,vjust=0.25)} +
    {if(doGeoMConfInt) geom_errorbar(data=fit$datGeoM,mapping=aes(x=ageCatNum,ymin=geoMLow,ymax=geoMUpp),col=colGeoM,lwd=0.15,width=0.3)} +
    geom_text(data=data.frame(x=xAnnot,y=0.98*yMax,label=toupper(strain)),mapping=aes(x=x,y=y,label=label),col=colAnnot,size=cexAnnot,hjust=1,vjust=1,fontface="bold")

  return(list(g=g,fit=fit,geoM=geoM))
}




######################
## RUN THE ANALYSIS ##
######################

doAna<-FALSE

fitList<-list()

if(doAna){
  for(s in st){
    cat(file=fileLog,append=T,paste(sep="","--- STRAIN ",s," processing now ---\n"))

    fitList[[s]]<-doAnalysis(dat=dat,varNameOrig=paste(sep="","res",s),varNameOrigNum=paste(sep="","res",s,"_num"),varNameImpHalfDL=paste(sep="","res",s,"_num_ImpHalfDL"),varNameImpDL=paste(sep="","res",s,"_num_ImpDL"),DL=0.15,Nboot=1e3,seed=seed)
    pcvpaFit<-fitList[[s]]
    save(pcvpaFit,file=paste(sep="","fit",s,".RData"))

    pdf(paste(sep="","fit",s,".pdf"),width=16,height=9)
    print(ggPlotFun(fit=fitList[[s]],strain=s,varNameImpDL=paste(sep="","res",s,"_num_ImpDL"),colModel="steelblue",colModelCI="steelblue",colChangepoints="orange")$g)
    dev.off()
    cat(file=fileLog,append=T,"--- done ---\n\n")
  }

  save(list=ls(),file=paste(sep="","fit_strains",paste(collapse="-",st),".RData"))
}else{
  dateToLoad<-"2019-10-21"
  for(s in st){
    if(s != "12F"){
      load(paste(sep="","../",dateToLoad,"/fit",s,".RData"))
      fitList[[s]]<-pcvpaFit
      rm(pcvpaFit)
    }else{
      varNameImpDL<-"res12F_num_ImpDL"
      varNameImpHalfDL<-"res12F_num_ImpHalfDL"
      varNameOrig<-"res12F"
      DL<-0.15

      datTmp<-dat[!is.na(dat[,varNameImpHalfDL]),] # remove non-recorded values
      datGeoM<-data.frame(ageCat=factor(unique(datTmp$age_cat_gmc)),geoM=NA,geoM_better=NA)
      for(i in 1:nrow(datGeoM)){
        datGeoM$geoM[i]<-gm_mean(datTmp[datTmp$age_cat_gmc==datGeoM$ageCat[i],varNameImpHalfDL]) # note that these use the simple DL/2 imputation for calculation
        if(sum(datTmp[datTmp$age_cat_gmc==datGeoM$ageCat[i],varNameOrig]=="<0.150")>0){
          if(sum(datTmp[datTmp$age_cat_gmc==datGeoM$ageCat[i],varNameOrig]=="<0.150")==sum(datTmp$age_cat_gmc==datGeoM$ageCat[i])){
            datGeoM$geoM_better[i]<-gm_mean(datTmp[datTmp$age_cat_gmc==datGeoM$ageCat[i],varNameImpDL])
          }else{
            datGeoM$geoM_better[i]<-exp(coef(censReg(as.formula(paste(sep="","log(",varNameImpDL,")~1")),data=datTmp[datTmp$age_cat_gmc==datGeoM$ageCat[i],],left=log(DL)))["(Intercept)"])
          }
        }else{
          datGeoM$geoM_better[i]<-exp(coef(lm(as.formula(paste(sep="","log(",varNameOrigNum,")~1")),data=datTmp[datTmp$age_cat_gmc==datGeoM$ageCat[i],]))["(Intercept)"])
        }
      }
      levels(datGeoM$ageCat)[levels(datGeoM$ageCat)=="0-2m"]<-"0-02m"
      levels(datGeoM$ageCat)[levels(datGeoM$ageCat)=="3-5m"]<-"03-05m"
      levels(datGeoM$ageCat)[levels(datGeoM$ageCat)=="6-8m"]<-"06-08m"
      levels(datGeoM$ageCat)[levels(datGeoM$ageCat)=="9-11m"]<-"09-11m"

      datGeoM<-datGeoM[order(as.character(datGeoM$ageCat)),]
      datGeoM$ageCatNum<-seq(1.5,58.5,by=3)

      fitList[[s]]<-list(
        dat=datTmp,
        datGeoM=datGeoM
      )

      rm(datTmp,datGeoM)
    }
  }
}

# build matrix with fit parameters for each serotype
parsMat<-matrix(nrow=length(st),ncol=1+3+1+4+4)
colnames(parsMat)<-c("k","Knot1","Knot2","Knot3","Intercept","Slope1","Slope2","Slope3","Slope4","IgGTime0","IgGPeak","IgGNadir","IgGTime60")
rownames(parsMat)<-paste(sep="","serotype_",st)
parsMat<-as.data.frame(parsMat)
for(i in 1:length(st)){
  s<-st[i]
  if(!is.null(fitList[[s]]$k)){
    k<-fitList[[s]]$k
    parsMat$k[i]<-k
    parsMat$Intercept[i]<-fitList[[s]]$coefs[1]
    parsMat$Slope1[i]<-fitList[[s]]$coefs[2]
    for(j in 1:k){
      parsMat[i,paste(sep="","Knot",j)]<-fitList[[s]]$coefs[2+k+1+j]
      parsMat[i,paste(sep="","Slope",j+1)]<-parsMat[i,paste(sep="","Slope",j)]+fitList[[s]]$coefs[2+j]
    }

    xTmp<-c(0,NA,NA,60)
    yTmp<-rep(NA,4)

    if(k==1){
      if(parsMat$Slope1[i]<0 & parsMat$Slope2[i]>0){
        xTmp[3]<-parsMat$Knot1[i]
      }
    }else if(k==2){
      if(parsMat$Slope1[i]>0 & parsMat$Slope2[i]<0){
        xTmp[2]<-parsMat$Knot1[i]
      }
      if(parsMat$Slope2[i]<0 & parsMat$Slope3[i]>0){
        xTmp[3]<-parsMat$Knot2[i]
      }
    }else if(k==3){
      if(parsMat$Slope1[i]>0 & parsMat$Slope2[i]<0){
        xTmp[2]<-parsMat$Knot1[i]
      }
      if(parsMat$Slope2[i]<0 & parsMat$Slope3[i]>0){
        xTmp[3]<-parsMat$Knot2[i]
      }else if(parsMat$Slope3[i]<0 & parsMat$Slope4[i]>0){
        xTmp[3]<-parsMat$Knot3[i]
      }
    }

    for(j in 1:length(xTmp)){
      if(!is.na(xTmp[j])){
        idxTmp<-which(abs(fitList[[s]]$xx-xTmp[j])==min(abs(fitList[[s]]$xx-xTmp[j])))
        idxTmpCI<-which(abs(fitList[[s]]$xxCI-xTmp[j])==min(abs(fitList[[s]]$xxCI-xTmp[j])))
        if(j!=1){
          yTmp[j]<-paste(sep="",format(nsmall=2,round(digits=2,mean(fitList[[s]]$yy[idxTmp])))," (",paste(collapse=",",format(nsmall=2,round(digits=2,c(mean(fitList[[s]]$yyLow[idxTmpCI]),mean(fitList[[s]]$yyHigh[idxTmpCI]))))),")") # taking means just in case the value is right between 2 evaluated ones
        }else{
          yTmp[j]<-paste(sep="",format(nsmall=2,round(digits=2,exp(parsMat$Intercept[i])))," (",paste(collapse=",",format(nsmall=2,round(digits=2,c(mean(fitList[[s]]$yyLow[idxTmpCI]),mean(fitList[[s]]$yyHigh[idxTmpCI]))))),")") # taking means just in case the value is right between 2 evaluated ones
        }
      }
    }

    parsMat[i,c("IgGTime0","IgGPeak","IgGNadir","IgGTime60")]<-yTmp
    rm(xTmp,yTmp)

  }
}

parsMatTmp<-parsMat
for(i in 1:length(st)){
  s<-st[i]
  if(!is.null(fitList[[s]]$k)){
    k<-fitList[[s]]$k
    parsMatTmp$Slope1[i]<-paste(sep="",format(nsmall=2,round(digits=2,parsMat$Slope1[i]))," (",format(nsmall=2,round(digits=2,fitList[[s]]$slopeCIs[[1]]$percent[4])),",",format(nsmall=2,round(digits=2,fitList[[s]]$slopeCIs[[1]]$percent[5])),")")
    for(j in 1:k){parsMatTmp[i,paste(sep="","Slope",j+1)]<-paste(sep="",format(nsmall=2,round(digits=2,parsMat[i,paste(sep="","Slope",j+1)]))," (",format(nsmall=2,round(digits=2,fitList[[s]]$slopeCIs[[j+1]]$percent[4])),",",format(nsmall=2,round(digits=2,fitList[[s]]$slopeCIs[[j+1]]$percent[5])),")")}
  }
}

write.csv(parsMatTmp,row.names=T,quote=T,file=paste(sep="","fit_strains",paste(collapse="-",st),"_parameters.csv"))

plotList<-list()
doGeoMCI<-FALSE
for(i in 1:length(st)){
  s<-st[i]

  if(s!="12F"){
    tmp<-ggPlotFun(fit=fitList[[s]],strain=s,varNameImpDL=paste(sep="","res",s,"_num_ImpDL"),DL=0.15,colModel="orange",colModelCI="orange",colChangepoints="steelblue",alphaModelCI=75,alphaChangepoints=130,plotModel=TRUE,cexAnnot=12,cexText=16,doLegend = FALSE,doGeoMConfInt = doGeoMCI,axisLabels=FALSE)
    plotList[[s]]<-tmp$g
  }else{
    tmp<-ggPlotFun(fit=fitList[[s]],strain=s,varNameImpDL=paste(sep="","res",s,"_num_ImpDL"),DL=0.15,colModel="orange",colModelCI="orange",colChangepoints="steelblue",alphaModelCI=75,alphaChangepoints=130,plotModel=FALSE,cexAnnot=12,cexText=16,doLegend = FALSE,doGeoMConfInt = doGeoMCI,axisLabels=FALSE)
    plotList[[s]]<-tmp$g
  }

  if(doGeoMCI){
    if(i==1){
      datGeoM<-tmp$geoM
    }else{
      datGeoM<-dplyr::left_join(x=datGeoM,y=tmp$geoM,by="age")
    }

    colnames(datGeoM)[colnames(datGeoM)=="geoM"]<-paste(sep="","GMC_",s)
  }
}

if(doGeoMCI){write.csv(datGeoM,row.names=F,quote=T,file="gmcTable_PCVPA.csv")}

pdf(width=25,height=15,file="plotAllSerotypes_PCVPA_noGmcCI.pdf")
grid.arrange(grobs=plotList,nrow=3,padding = unit(0.1, "line"),bottom=grid::textGrob("age (months)", gp = grid::gpar(fontsize = 24)),left=grid::textGrob("IgG concentration (ug/mL)", rot=90, gp = grid::gpar(fontsize = 24)))
dev.off()

s<-"23F"
pdf(width=6,height=5,file="plotSerotype23F_PCVPA.pdf")
print(ggPlotFun(fit=fitList[[s]],strain=s,varNameImpDL=paste(sep="","res",s,"_num_ImpDL"),DL=0.15,colModel="orange",colModelCI="orange",colChangepoints="steelblue",alphaModelCI=75,alphaChangepoints=130,plotModel=TRUE,cexAnnot=12,doLegend = TRUE,doGeoMConfInt = TRUE)$g)
dev.off()

cat(file=fileLog,append=T,"This is the end, my friend.\n")

# extract peaks, lows and troughs
peaksLowsTroughsMat<-matrix(nrow=length(st),ncol=2+3*2)
colnames(peaksLowsTroughsMat)<-c("IgG_age0","IgG_age60","peak","peak_age","low","low_age","trough","trough_age")
rownames(peaksLowsTroughsMat)<-st
peaksLowsTroughsMat<-data.frame(peaksLowsTroughsMat)
peaksLowsTroughsMat$trough_age<-as.character(peaksLowsTroughsMat$trough_age)

evalFun<-function(x,k,cf){
    logy<-cf[1]+cf[2]*x
    if(k>0){
        for(i in 1:k){
            logy<-logy+cf[2+i]*ifelse((x-cf[2+k+1+i])>=0,x-cf[2+k+1+i],0)
        }
    }
    return(exp(logy))
}

for(s in setdiff(st,"12F")){
    peaksLowsTroughsMat[s,"IgG_age0"]<-evalFun(x=0,k=fitList[[s]]$k,cf=fitList[[s]]$coefs)
    peaksLowsTroughsMat[s,"IgG_age60"]<-evalFun(x=60,k=fitList[[s]]$k,cf=fitList[[s]]$coefs)
    if(fitList[[s]]$k==1){
        # if k==1, no peak; trough = > 1st knot, BUT only if abs(slope)<0.1
        peaksLowsTroughsMat[s,"peak"]<-NA
        peaksLowsTroughsMat[s,"peak_age"]<-NA
        peaksLowsTroughsMat[s,"low"]<-evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+1],k=fitList[[s]]$k,cf=fitList[[s]]$coefs)
        peaksLowsTroughsMat[s,"low_age"]<-round(digits=2,fitList[[s]]$coefs[2+fitList[[s]]$k+1+1])
        if(abs(parsMat[paste(sep="_","serotype",s),"Slope2"])<0.01){
            peaksLowsTroughsMat[s,"trough"]<-(evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+1],k=fitList[[s]]$k,cf=fitList[[s]]$coefs)+evalFun(x=60,k=fitList[[s]]$k,cf=fitList[[s]]$coefs))/2
            peaksLowsTroughsMat[s,"trough_age"]<-paste(collapse="-",round(digits=2,c(fitList[[s]]$coefs[2+fitList[[s]]$k+1+1],60)))
        }else{
            peaksLowsTroughsMat[s,"trough"]<-NA
            peaksLowsTroughsMat[s,"trough_age"]<-NA
        }
    }else if(fitList[[s]]$k==2){
        # if k==2, then peak = 1st knot, low = 2nd knot, trough = > 2nd knot BUT only if abs(slope)<=0.1
        peaksLowsTroughsMat[s,"peak"]<-evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+1],k=fitList[[s]]$k,cf=fitList[[s]]$coefs)
        peaksLowsTroughsMat[s,"peak_age"]<-round(digits=2,fitList[[s]]$coefs[2+fitList[[s]]$k+1+1])
        peaksLowsTroughsMat[s,"low"]<-evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+2],k=fitList[[s]]$k,cf=fitList[[s]]$coefs)
        peaksLowsTroughsMat[s,"low_age"]<-round(digits=2,fitList[[s]]$coefs[2+fitList[[s]]$k+1+2])
        if(abs(parsMat[paste(sep="_","serotype",s),"Slope3"])<0.01){
            peaksLowsTroughsMat[s,"trough"]<-(evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+2],k=fitList[[s]]$k,cf=fitList[[s]]$coefs)+evalFun(x=60,k=fitList[[s]]$k,cf=fitList[[s]]$coefs))/2
            peaksLowsTroughsMat[s,"trough_age"]<-paste(collapse="-",c(round(digits=2,fitList[[s]]$coefs[2+fitList[[s]]$k+1+2],60)))
        }else{
            peaksLowsTroughsMat[s,"trough"]<-NA
            peaksLowsTroughsMat[s,"trough_age"]<-NA
        }
    }else if(fitList[[s]]$k==3){
      # if k==3, then peak = 1st knot, low = 2nd or 3rd knot, trough = between 2nd and 3rd knot or after 3rd (check slope is flat!)
        peaksLowsTroughsMat[s,"peak"]<-evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+1],k=fitList[[s]]$k,cf=fitList[[s]]$coefs)
        peaksLowsTroughsMat[s,"peak_age"]<-round(digits=2,fitList[[s]]$coefs[2+fitList[[s]]$k+1+1])
        peaksLowsTroughsMat[s,"low"]<-min(evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+2],k=fitList[[s]]$k,cf=fitList[[s]]$coefs),
                                           evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+3],k=fitList[[s]]$k,cf=fitList[[s]]$coefs))
        idxMin<-1+which(c(evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+2],fitList[[s]]$k,cf=fitList[[s]]$coefs),evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+3],fitList[[s]]$k,cf=fitList[[s]]$coefs))==peaksLowsTroughsMat[s,"low"])
        peaksLowsTroughsMat[s,"low_age"]<-round(digits=2,fitList[[s]]$coefs[2+fitList[[s]]$k+1+idxMin])
        if(abs(parsMat[paste(sep="_","serotype",s),"Slope3"])<abs(parsMat[paste(sep="_","serotype",s),"Slope4"])){
            peaksLowsTroughsMat[s,"trough"]<-(evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+2],k=fitList[[s]]$k,cf=fitList[[s]]$coefs)+evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+3],k=fitList[[s]]$k,cf=fitList[[s]]$coefs))/2
            peaksLowsTroughsMat[s,"trough_age"]<-paste(collapse="-",round(digits=2,fitList[[s]]$coefs[2+fitList[[s]]$k+1+2:3]))
        }else if(abs(parsMat[paste(sep="_","serotype",s),"Slope4"])<0.01){
            peaksLowsTroughsMat[s,"trough"]<-(evalFun(x=fitList[[s]]$coefs[2+fitList[[s]]$k+1+3],k=fitList[[s]]$k,cf=fitList[[s]]$coefs)+evalFun(x=60,k=fitList[[s]]$k,cf=fitList[[s]]$coefs))/2
            peaksLowsTroughsMat[s,"trough_age"]<-paste(collapse="-",round(digits=2,c(fitList[[s]]$coefs[2+fitList[[s]]$k+1+3],60)))
        }

    }
}

write.csv(peaksLowsTroughsMat,row.names=T,quote=F,file=paste(sep="","fit_strains",paste(collapse="-",st),"_peaksLowsTroughs.csv"))
