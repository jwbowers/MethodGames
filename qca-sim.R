## If any of the following packages are not installed, install them.


## Now load the packages
library(QCA)
library(glmnet)
library(SIS)
library(compiler)
library(parallel)
library(data.table) ## trying this out for speed
numcores<-detectCores()

makedatamatrix<-function(ntotalvars,N){
  ## Here we allow each column of X to have different numbers of 1s versus 0s.
  ## data.table(replicate(ntotalvars,sample(c(1,0),size=N,replace=TRUE)))

  ## Here we require half of each X to be 1 versus 0.
  data.table(replicate(ntotalvars,sample(rep(c(1,0),N/2))))
}

makeoutcome<-function(X,expr){
  ## expr is a character evaluating to a logical expression
  ## for example, "V1*V2*V3 | V3*V4 | V5" meaning
  ## Y is 1 if all of X1,X2,X3 are 1, OR if X3 and X4 are 1, OR if X5 is 1.
  ## notice data.table syntax where the second element is an expression of
  ## column names
  as.numeric(X[,eval(parse(text=expr))] )
}

makemodelmatrix<-function(thedata,interaction.order){
  ## In production use Matrix(...,sparse=TRUE) from library(Matrix)
  ## No intercept
  model.matrix(as.formula(paste("~(-1+.)^",interaction.order,sep="")),data=thedata)
}


fitfn<-function(y,X,DAT){
  message(".",appendLF=FALSE)

  truthparts<-gsub("\\s","",strsplit(thetruth,"|",fixed=TRUE)[[1]])
  yF<-factor(y)

  ## Adaptive Lasso
  ### Choose lambda by minimizing misclassification
  theridge.cv <- try( cv.glmnet(X, yF, alpha=0, type.measure="class",
                                family="binomial", standardize=TRUE,
                                nfolds=min(round(nrow(X)/2), 10), grouped=FALSE))

  if(inherits(theridge.cv,'try-error')){
    lassofound<-NA
  } else {
    bhat<-as.matrix(coef(theridge.cv,s="lambda.min"))[-1,1] ## coef() is a sparseMatrix
    if(all(bhat==0)){
      ## if bhat is all zero then assign not exactly zero
      bhat<-rep(.Machine$double.eps*2,length(bhat))
    }
    adpen<-(1/pmax(abs(bhat),.Machine$double.eps))

    thelasso.cv<-try( cv.glmnet(X,yF,alpha=1,type.measure="class",
                                family="binomial",
                                exclude=which(bhat==0),
                                penalty.factor=adpen,
                                standardize=TRUE,
                                nfolds=min(round(nrow(X)/2),10),
                                grouped=FALSE) )

    if(inherits(thelasso.cv,'try-error')){
      lassofound<-NA
    } else {
      thelasso.coef<-coef(thelasso.cv,s="lambda.min")
      ### Did the adaptive lasso return non-zero coefs for the truth and only the truth?
      lassofound <- all( gsub("*",":",truthparts,fixed=TRUE) %in% row.names(thelasso.coef)[thelasso.coef[,1]!=0][-1] )
    }
  }

    ## QCA
    theqca<-try(eqmcc(DAT,outcome="Y"))
    ### Did QCA return the truth and only the truth?
    if(inherits(theqca,'try-error')){
      qcafound<-NA
    } else {
      qcafound<-all(theqca$PIs %in% truthparts)
    }

    ## SIS/ISIS
    thesis<-try( SIS(data=list(x=X,y=y),
                family=binomial(),
                folds=min(round(nrow(X)/2),10),
                maxloop=10,
                detail=FALSE) )
    if(inherits(thesis,'try-error')){
      sisfound<-NA 
    } else {
    sisfound <- all( gsub("*",":",truthparts,fixed=TRUE) %in% colnames(X)[thesis$ISISind])
    }

    return(c(qcafound=qcafound,
             lassofound=lassofound,
             sisfound=sisfound))
  } 


  myfn.maker<-function(ntotalvars,N,thetruth){
    force(ntotalvars);force(N);force(thetruth)
    function(){
      thedata<-makedatamatrix(ntotalvars,N)
      theX<-makemodelmatrix(thedata,4) ## fourway interactions
      ## theX<-unique(theX,MARGIN=2) ## delete identical columns
      theY<-makeoutcome(thedata,thetruth)
      ##while(sum(theY)==1 | any(colSums(theX) %in% c(0,N))){
      ##  ## don't allow constant Xs
      ##  thedata<-makedatamatrix(ntotalvars,N)
      ##  theX<-makemodelmatrix(thedata,4) ## fourway interactions
      ##  ## theX<-unique(theX,MARGIN=2) ## delete identical columns
      ##  theY<-makeoutcome(thedata,thetruth)
      ##}
      names(theY)<-row.names(theX)
      ## Some functions fail if we substitute data.table class objects for
      ## data.frame objects.
      thedf<-data.frame(thedata)
      thedf$Y<-theY
      fitfn(y=theY,X=theX,DAT=thedf)
    }
  }



  ## Start doing simulations
  ### Setup data. Some routines want matrices others want data.frames
  nsims<-100*numcores

  ## Easy sim: All variables in truth, p<n
  thetruth<-"V1*V2*V3 | V4*V5"
  ntotalvars<-5
  N<-30

  myfn<-myfn.maker(ntotalvars,N,thetruth)
  cmp.myfn<-cmpfun(myfn,options=list(optimize=3))

  set.seed(12345) ## for replicability
  ## mysim1<-replicate(nsims,cmp.myfn())
  mysim1<-mclapply(1:nsims,function(i){ cmp.myfn() },mc.cores=numcores)
  mysim1.arr<-simplify2array(mysim1)
  apply(mysim1.arr,1,mean,na.rm=TRUE) ## proportion of the time that qca and isis/scad and lasso identified the right terms
  apply(mysim1.arr,1,function(x){mean(is.na(x)) }) ## what proportion NAs in each method?

  save(mysim1,file="mysim1.rda")

  ## Very hard sim: most variables not in truth, p>n
  ntotalvars<-15

  myfn<-myfn.maker(ntotalvars,N,thetruth)
  cmp.myfn<-cmpfun(myfn,options=list(optimize=3))

  set.seed(12345) ## for replicability
  ## mysim2<-replicate(nsims,cmp.myfn())
  mysim2<-mclapply(1:nsims,function(i){ cmp.myfn() },mc.cores=numcores)
  mysim2.arr<-simplify2array(mysim2)
  apply(mysim2.arr,1,mean,na.rm=TRUE) ## proportion of the time that qca and isis/scad and lasso identified the right terms
  apply(mysim2.arr,1,function(x){mean(is.na(x)) }) ## what proportion NAs in each method?

  save(mysim2,file="mysim2.rda")

