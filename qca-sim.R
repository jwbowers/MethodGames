
library(brglm)
library(QCA)
##library(QCA3) ## reduce()
library(glmnet)
## library(randomForest)
library(KRLS)
library(parallel)
library(data.table) ## trying this out for speed

numcores<-detectCores()

set.seed(20130516)
makedatamatrix<-function(ntotalvars,N){
  ## Here we allow each column of X to have different numbers of 1s versus 0s.
  data.table(replicate(ntotalvars,sample(c(1,0),size=N,replace=TRUE)))
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

### Setup data. Some routines want matrices others want data.frames
nsims<-200*numcores
ntotalvars<-20
N<-15
thedata<-makedatamatrix(ntotalvars,N)
theX<-makemodelmatrix(thedata,3)
thetruth<-"V1*V2*V3 | V4*V5"
theY<-makeoutcome(thedata,thetruth)

## Some functions fail if we substitute data.table class objects for data.frame objects.
thedf<-data.frame(thedata)
thedf$Y<-theY

myfn<-function(y,X,DAT){
  message(".",appendLF=FALSE)

  ## Adaptive Lasso
  ### Choose lambda by minimizing misclassification
  theridge.cv<-cv.glmnet(theX,theY,alpha=0,type.measure="class",family="binomial",standardize=TRUE)
  bhat<-as.matrix(coef(theridge.cv,s="lambda.min"))[-1,1] ## coef() is a sparseMatrix
  if(all(bhat==0)){
	  ## if bhat is all zero then assign not exactly zero
	  bhat<-rep(.Machine$double.eps*2,length(bhat))
  }
  adpen<-(1/pmax(abs(bhat),.Machine$double.eps))

  thelasso.cv<-cv.glmnet(theX,theY,alpha=1,type.measure="class",family="binomial",
                         exclude=which(bhat==0),
                         penalty.factor=adpen,
                         standardize=TRUE)
  thelasso.coef<-coef(thelasso.cv,s="lambda.min")

  ## QCA
  theqca<-eqmcc(thedf,outcome="Y")

  ## ## Random Forest

  ## forest<-randomForest(x=thedata, y = as.factor(theY),
  ##                    ntree = 1000, importance=TRUE,
  ##                    proximity=TRUE,
  ##                    norm.votes=TRUE, keep.forest=TRUE)

  ## KRLS

  thekrls<-krls(thedata,theY,derivative=TRUE)
  
  ## Record whether the truth was found.
truthparts<-gsub("\\s","",strsplit(thetruth,"|",fixed=TRUE)[[1]])

## Did the adaptive lasso return non-zero coefs for the truth and only the truth?
  lassofound<-all(row.names(thelasso.coef)[thelasso.coef[,1]!=0][-1] %in% gsub("*",":",truthparts,fixed=TRUE) )

## Did QCA return the truth and only the truth?
qcafound<-all(theqca$PIs %in% truthparts)

## Did randomForest return the truth and only the truth?


  return(c(qcafound=qcafound,firthfound=firthfound,lassofound=lassofound))
}

cmp.myfn<-cmpfun(myfn,options=list(optimize=3))

set.seed(12345) ## for replicability
mysim<-simplify2array(mclapply(1:nsims,function(i){ cmp.myfn() },mc.cores=numcores))
apply(mysim,1,mean) ## proportion of the time that qca and brglm and lasso identified the right terms

