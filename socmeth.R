## File: socmeth.R Description: An R script to compare QCA, the Adaptive Lasso,
## and ISIS/SCAD Author: Jake Bowers This file lives at
## https://github.com/jwbowers/MethodGames To reproduce socmeth.Rout, run this
## file in batch mode: R CMD BATCH socmeth.R

# Setup/Install Libraries If any of the following packages are not installed,
# install them.
neededpkgs <- c("QCA", "glmnet", "SIS", "compiler", "parallel", "data.table")
localpkgs <- installed.packages()
to.install <- neededpkgs[!(neededpkgs %in% localpkgs[, "Package"])]
if (length(to.install) > 0) {
    install.packages(to.install, dependencies = TRUE)
}

## Now load the packages
library(QCA)
library(glmnet)
library(SIS)
library(data.table)  ## trying this out for speed

## These next three lines may not work on Windows.
library(compiler)
library(parallel)
numcores <- detectCores()

# Now define functions: This section might easily live in another file but is
# here for simplicity.

makedatamatrix <- function(nfeatures, N) {
    ## nfeatures: scalar integer, number of observed features of a case N: scalar,
    ## integer, number of cases Here we allow each column of X to have different
    ## numbers of 1s versus 0s but to have half 1 and half 0 on average.
    data.table(replicate(nfeatures, sample(c(1, 0), size = N, replace = TRUE)))
    ## Here we require exactly half of each X to be 1 versus 0.
    ## data.table(replicate(nfeatures,sample(rep(c(1,0),N/2))))
}

makeoutcome <- function(X, thetruth) {
    ## X: a data.table object such as that arising from makedatamatrix() thetruth: is
    ## a character evaluating to a logical thetruthession of the column names of X for
    ## example, 'V1*V2*V3 | V3*V4 | V5' meaning Y is 1 if all of X1,X2,X3 are 1, OR if
    ## X3 and X4 are 1, OR if X5 is 1.  notice data.table syntax where the second
    ## element is an thetruthession of column names
    as.numeric(X[, eval(parse(text = thetruth))])
}

makemodelmatrix <- function(thedata, interactionorder) {
    ## thedata: a data.table or data.frame object interactionorder: a scalar integer
    ## for the number of interactions to search over.  Notice: No intercept
    model.matrix(as.formula(paste("~(-1+.)^", interactionorder, sep = "")), data = thedata)
}


fitfn <- function(y, X, DAT, thetruth) {
    ## y: is the outcome X: is a matrix object arising from makemodelmatrix DAT: is a
    ## data.frame object (some functions don't like data.tables or matrices) thetruth:
    ## a character evaluting to a logical expression of the column names of X

    ### Note: All fitters are wrapped in try() environments because sometimes,
    ### especially if N is small, the fitters throw errors for some configurations of X
    ### and y. In those cases, the function returns a NA.

    message(".", appendLF = FALSE)  ## print a dot to indicate speed of script

    ## Represent thetruth in ways that can be compared to the output of the learners
    truthparts <- gsub("\\s", "", strsplit(thetruth, "|", fixed = TRUE)[[1]])
    truthpartstmp <- gsub("*", ":", truthparts, fixed = TRUE)  ## for models using ':' rather than '*'

    ## Adaptive Lasso (L1-penalized logistic regression weighted by an L2-penalized
    ## logistic regression) Choose lambda by minimizing misclassification in k-fold
    ## cross-validation
    yF <- factor(y)  ## glmnet wants y to be a factor if it is binary
    theridge.cv <- try(cv.glmnet(X, yF, alpha = 0, type.measure = "class", family = "binomial",
        standardize = TRUE, nfolds = min(round(nrow(X) / 2), 10), grouped = FALSE))

    if (inherits(theridge.cv, "try-error")) {
        lassofound <- NA
    } else {
        bhat <- as.matrix(coef(theridge.cv, s = "lambda.min"))[-1, 1]  ## coef() is a sparseMatrix
        if (all(bhat == 0)) {
            ## if bhat is all zero then assign very close to zero weight to all.  Amounts to
            ## penalizing all of the second stage to zero.
            bhat <- rep(.Machine$double.eps * 2, length(bhat))
        }
        adpen <- (1/pmax(abs(bhat), .Machine$double.eps))  ## the adaptive lasso weight

        thelasso.cv <- try(cv.glmnet(X, yF, alpha = 1, type.measure = "class", family = "binomial",
            exclude = which(bhat == 0), penalty.factor = adpen, standardize = TRUE,
            nfolds = min(round(nrow(X)/2), 10), grouped = FALSE))

        if (inherits(thelasso.cv, "try-error")) {
            lassofound <- NA
        } else {
            thelasso.coef <- coef(thelasso.cv, s = "lambda.min")
            nonzerolasso.coef <- row.names(thelasso.coef)[thelasso.coef[, 1] != 0][-1]
            ### Did the adaptive lasso return non-zero coefs for the truth and only the truth?
            lassofound <- setequal(truthpartstmp, nonzerolasso.coef)
        }
    }

    ## QCA Use both 0 and 1 outcome as to be explained.  Code 'found truth' if either
    ## finds the truth.
    theqcapos <- try(eqmcc(DAT, outcome = "Y"))
    theqcaneg <- try(eqmcc(DAT, outcome = "Y", explain = "0"))
    ### Did QCA return the truth and only the truth?
    if (inherits(theqcapos, "try-error") | inherits(theqcaneg, "try-error")) {
        qcafound <- NA
    } else {
        ## all.equal(theqcapos$solution,theqcapos$PIs)
        pospis <- gsub("[a-z][0-9]?.\\*?", "", theqcapos$PIs)
        pospis <- gsub("\\*$", "", pospis)
        negpis <- gsub("[a-z][0-9]?.\\*?", "", theqcaneg$PIs)
        negpis <- gsub("\\*$", "", negpis)
        qcafound <- all(setequal(pospis, truthparts) | setequal(negpis, truthparts))
    }

    ## ISIS/SCAD
    options(warn = -1)  ## annoying warnings from glm.fit()
    thesis <- try(SIS(x = X, y = y, family = "binomial", iter.max=100), silent = TRUE)
    options(warn = 0)  ## turn back on default warning behavior

    if (inherits(thesis, "try-error")) {
        sisfound <- NA
    } else {
        thesis.coef <- colnames(X)[thesis$ix]
        sisfound <- setequal(truthpartstmp, thesis.coef)
    }

    ## LM
    thelm <- try(lm.fit(x = X, y = y))
    if (inherits(thelm, "try-error")) {
        lmfound <- NA
    } else {
      thelm.coef <- coef(thelm)
      nonzerolm.coef <- names(thelm.coef[zapsmall(thelm.coef,digits=10) !=0 & !is.na(thelm.coef)])
      lmfound <- setequal(truthpartstmp, thelm.coef)
    }

    ## GLM
    theglm <- try(glm.fit(x = X, y = y, family=binomial(link="logit")))
    if (inherits(theglm, "try-error")) {
        glmfound <- NA
    } else {
      theglm.coef <- coef(theglm)
      nonzeroglm.coef <- names(theglm.coef[zapsmall(theglm.coef,digits=10) !=0 & !is.na(theglm.coef)])
      glmfound <- setequal(truthpartstmp, theglm.coef)
    }

    ## KRLS may be particularly promising for continuous features
    ## thekrls<-krls(thedata,theY,derivative=TRUE) derivmat<-thekrls$derivatives
    ## colnames(derivmat)<-colnames(thekrls$X) KRLS does not report interactions
    ## although they are implied. Thanks to Chad Hazlett for suggesting something like
    ## the following to assess relationships among marginals
    ## krlsfound.fn<-function(one,two,X=thekrls$X,dmat=derivmat){
    ## thelm<-lm(X[,one]~dmat[,two]) thesum<-summary(thelm)
    ## pf(thesum$fstatistic[1],thesum$fstatistic[2],thesum$fstatistic[3],
    ## lower.tail=FALSE)<=.05 } krlsfound<-all(krlsfound.fn('V1','V2'),
    ## krlsfound.fn('V4','V5'))
    return(c(qcafound = qcafound, lassofound = lassofound, sisfound = sisfound, lmfound = lmfound, glmfound = glmfound))
}


## This function makes a function that repeats the dataset and truth finding
gamefn.maker <- function(nfeatures, N, thetruth, interactionorder) {
    ## nfeatures: scalar integer, number of observed features of a case N: scalar,
    ## integer, number of cases thetruth: is a character evaluating to a logical
    ## thetruthession of the column names of X interactionorder: a scalar integer for
    ## the number of interactions to search over.
    force(nfeatures)
    force(N)
    force(thetruth)
    force(interactionorder)
    function() {
        thedata <- makedatamatrix(nfeatures, N)
        theX <- makemodelmatrix(thedata, interactionorder)  ## all 4-way interactions
        theY <- makeoutcome(thedata, thetruth)
        while (sum(theY) %in% c(0, 1, N)) {
            ## don't allow Y with only 1 positive obs or constant Y
            thedata <- makedatamatrix(nfeatures, N)
            theX <- makemodelmatrix(thedata, interactionorder)
            theY <- makeoutcome(thedata, thetruth)
        }
        names(theY) <- row.names(theX)
        thedata$Y <- theY
        fitfn(y = theY, X = theX, DAT = as.data.frame(thedata), thetruth = thetruth)
    }
}


# Specify and Run the Games
nplayers <- 100 * numcores  ## In the published article I had 8 cores, so I had 800 players

## Easy game: All variables in truth, p<n
thetruth <- "V1*V2*V3 | V4*V5"
nfeatures <- 5
N <- 40

easygamefn <- gamefn.maker(nfeatures, N, thetruth, 4)
cmp.easygamefn <- cmpfun(easygamefn, options = list(optimize = 3))  ## see if byte-compiling speeds up the runs

## Set seeds for random number generator. Since I am using a parallelize setup, I
## am going through some extra work to make these results reproducible

if (numcores > 1) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(20130501)
    mc.reset.stream()
} else {
    set.seed(20130501)
}

## easygameresults<-replicate(nplayers,cmp.easygamefn())
easygameresults <- mclapply(1:nplayers, function(i) { cmp.easygamefn() }, mc.cores = numcores, mc.set.seed = TRUE)
easygameresults.arr <- simplify2array(easygameresults)  ## turn list into array for easy summarization
apply(easygameresults.arr, 1, mean, na.rm = TRUE)  ## proportion of the time that qca and scad and lasso identified the right terms
apply(easygameresults.arr, 1, function(x) {
    mean(is.na(x))
})  ## what proportion NAs in each method?

save(easygameresults, file = "easygameresults.rda")

## Hard game: most variables not in truth, p>n
nfeatures <- 15

hardgamefn <- gamefn.maker(nfeatures, N, thetruth, 4)
cmp.hardgamefn <- cmpfun(hardgamefn, options = list(optimize = 3))

if (numcores > 1) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(20130501)
    mc.reset.stream()
} else {
    set.seed(20130501)
}

hardgameresults <- mclapply(1:nplayers, function(i) { cmp.hardgamefn() }, mc.cores = numcores, mc.set.seed = TRUE)
hardgameresults.arr <- simplify2array(hardgameresults)
apply(hardgameresults.arr, 1, mean, na.rm = TRUE)  ## proportion of the time that qca and scad and lasso identified the right terms
apply(hardgameresults.arr, 1, function(x) {
    mean(is.na(x))
})  ## what proportion NAs in each method?

save(hardgameresults, file = "hardgameresults.rda")

