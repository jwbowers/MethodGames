
library(brglm)
library(QCA)
library(QCA3) ## reduce()
library(glmnet)
library(randomForest)
library(parallel)
library(data.table) ## trying this out for speed

numcores<-detectCores()
nsims<-200*numcores
ntruevars<-5
ntotalvars<-20
N<-10

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

thedata<-makedatamatrix(ntotalvars,N)
theX<-makemodelmatrix(thedata,4)
thetruth<-"V1*V2*V3 | V3*V4 | V5"
theY<-makeoutcome(thedata,thetruth)

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

  ## Random Forest

  forest<-randomForest(x=DAT[,-1],y = as.factor(Y),
                     mtry = sqrt(ncol(Dat)-1),
                     ntree = 1000,importance=TRUE,norm.votes=TRUE,keep.forest=TRUE)

  ## Record whether the truth was found.

  ## Estimate models with intercepts (for now) but exclude them
  ## from consideration [perhaps wrong].
  lassofound<-all(row.names(thelasso.coef)[thelasso.coef[,1]!=0][-1] %in% c('X1:X2','X3:X4'))
  ## Are both of the largest coefs for the two two-way interactions?
  firthfound<-all(names(sort(coef(theglm)[-1],decreasing=TRUE)[1:2]) %in% c('X1:X2','X3:X4'))
  qcafound<-all(theqca$PIs %in% c('X1*X2','X3*X4'))

  return(c(qcafound=qcafound,firthfound=firthfound,lassofound=lassofound))
}

cmp.myfn<-cmpfun(myfn,options=list(optimize=3))

set.seed(12345) ## for replicability
mysim<-simplify2array(mclapply(1:nsims,function(i){ cmp.myfn() },mc.cores=numcores))
apply(mysim,1,mean) ## proportion of the time that qca and brglm and lasso identified the right terms

