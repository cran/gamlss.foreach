# this this for principle componet regression gamlss model
# ---------------------------------------------------------------------------
#--Mikis January 2019; based on a function written in 2011
# - revised August 2019  
# TO DO
#   i) prediction (OK)
#  ii) parallel fitting (probably OK)
# iii) if df are fixed or max.number is less than min(dim(X)) max.number OK
#  iv) do the svd in advance (OK)
#----------------------------------------------------------------------------
require('foreach')
require('doParallel')
#----------------------------------------------------------------------------
pc <- function (x.vars=NULL, x=NULL, x.svd=NULL,  df = NULL, center=TRUE, 
                scale=TRUE, tol = NULL, max.number=min(p,n),  k=log(n), 
                method=c("t-values", "GAIC","k-fold")) 
{
    scall <- deparse(sys.call(), width.cutoff = 500L)
    method <-  match.arg(method) 
# get where "gamlss" is in system call, it can be in gamlss() or predict.gamlss()
#----gamlss enviroment 
   rexpr <- grepl("gamlss",sys.calls()) ## 
 newData <- beta <- scaleold <- centerold <-  NULL
    for (i in length(rexpr):1) { 
      position <- i # get the position
      if (rexpr[i]==TRUE) break
    }
  gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
# get the data
    if (sys.call(position)[1]=="predict.gamlss()") 
      { # if predict is used 
      object <- get("object", envir=gamlss.env)
        beta <- getSmo(object)$beta
   centerold <- getSmo(object)$pc["center"]
    scaleold <- getSmo(object)$pc["scale"]
     newData <- get("newdata", envir=gamlss.env)
        Data <- get("data", envir=gamlss.env)
      } else if (sys.call(position)[1]=="gamlss()") { # if gamlss() is used
      if (is.null(get("gamlsscall", envir=gamlss.env)$data)) { # if no data argument but the formula can be interpreted
        Data <- data.frame(cbind(x))	
      } else {# data argument in gamlss 
        Data <- get("gamlsscall", envir=gamlss.env)$data
      }
    } else {
        Data <- get("data", envir=gamlss.env)
    }
        Data <- data.frame(eval(substitute(Data)))
# finish data 
if (is.null(x) && is.null(x.vars) && is.null(x.svd)) 
      stop("x or x.vars or x.svd has to be set in pc()")
if    (sum(c(!is.null(x), !is.null(x.vars), !is.null(x.svd)))>1L)
    stop("use only one of the arguments x, x.vars or x.svd")
if (!is.null(x.vars)) 
  {
   Xmtrx <- as.matrix(Data[, x.vars])
   Xmtrx <- scale(Xmtrx, center = center, scale = scale)
  center <- attr(Xmtrx,"scaled:center")  
   scale <- attr(Xmtrx,"scaled:scale")
       n <- dim(Xmtrx)[1]
       p <- dim(Xmtrx)[2]
       S <- La.svd(Xmtrx, nu = 0, nv= max.number)
  cNames <- colnames(Xmtrx)
 }
if (!is.null(x)) 
    {
  if ("scaled:center"%in%names(attributes(x))) 
    warning("the x matrix will be scaled again \n", "This will creates problems if prediction is used later")
   Xmtrx <- scale(x, center = center, scale = scale) 
  center <- attr(Xmtrx,"scaled:center")  
   scale <- attr(Xmtrx,"scaled:scale")
       n <- dim(Xmtrx)[1]
       p <- dim(Xmtrx)[2] 
       S <- La.svd(Xmtrx, nu = 0, nv= max.number) 
  cNames <- colnames(Xmtrx)
    } 
if (!is.null(x.svd))     
    {
     if (length(x.svd)!=6) stop("x.svd should be a list of 6, created by GetSVD()") 
       n <- dim(x.svd$u)[1]
       p <- dim(x.svd$vt)[1]
  center <- x.svd$center 
   scale <- x.svd$scale     
   Xmtrx <- x.svd[["X"]] 
  cNames <- colnames(x.svd[["X"]])
       S <- x.svd[c("d","u","vt")] 
}
if (!is.null(df)&&df>min(n,p)) stop("the df have to be less than the number of variables")
if (!is.null(df)) max.number <- df # if df is set only to few component  
## this comes from prcomp ---------------------------
if (!is.null(tol)) #-----------------------------------------
  {
              rank <- sum(S$d > (S$d[1] * tol))
        if (rank < ncol(x)) # there is no x here for some cases 
              S$vt <- S$vt[, 1:rank, drop = FALSE]
}
#-------------------------------------------------------------        
               S$d <- S$d/sqrt(max(1, n - 1))
    dimnames(S$vt) <- list(paste("PC", seq(len = nrow(S$vt)), 
                                   sep = ""), cNames)
                 r <- list(sdev = S$d, rotation = t(S$vt), center=center, scale=scale)
               r$x <- Xmtrx %*% t(S$vt) # this should be smaller if df or max.numner is set
          class(r) <- "prcomp"
                 x <- rep(0, n)
    attr(x, "PC")  <- r
   attr(x, "name") <- cNames
  attr(x, "maxNo") <- max.number
   attr(x, "call") <- substitute(gamlss.pc(data[[scall]], z, w))
     attr(x, "df") <- df
      attr(x, "k") <- k
 attr(x, "method") <- method 
attr(x, "newData") <- newData 
   attr(x, "beta") <- beta
          class(x) <- c("smooth", class(x))
    x
}
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
gamlss.pc <- function(x, y, w, xeval = NULL, ...)
{
if (is.null(xeval))
{    
     PC <-  attr(x, "PC")
 method <- attr(x, "method")
  names <- as.character(attr(x, "name"))
    edf <- as.vector(attr(x, "df"))
      k <- as.vector(attr(x, "k"))
      n <- nrow(PC$x)
      p <- ncol(PC$x)
  maxno <-  as.numeric(attr(x, "maxNo"))
#---------------------------------------------------  
    fun <- function(ind, k=k)
        {
         m <- lm(y~PC$x[,1:ind, drop = FALSE]-1, weights=w)
         AIC(m, k=k)
        }
#--------------------------------------------------  
if (is.null(edf)) # If the degrees of freedom are not fixed 
{
   MaxNo <- min(n-5, p, maxno)# 5 is added so the linear fir has at least 5 obs
if (method=="GAIC")
   
 {
   AiC <- foreach(i = 1 : MaxNo, .packages="gamlss",  
                  .export=c('fun', "MaxNo"),
                  .combine = rbind)%dopar%{
            #m <- lm(y~PC$x[,1:i, drop = FALSE]     , weights=w)
             m <- lm(y~I(PC$x[,1:i, drop = FALSE])-1, weights=w)
            AA <- AIC(m, k=k)} 
           edf <- which.min( AiC[is.finite(AiC)])
   # cat("edf", edf,"\n")
   # plot(AiC)
             T <- as.matrix(PC$x[,1:edf, drop = FALSE])
           fit <- lm(y~T-1, weights=w) 
          beta <- as.vector(PC$rotation[,1:edf]%*%matrix(coef(fit)))
  } 
if (method=="t-values")
  {
          T <- as.matrix(PC$x[,1:MaxNo, drop = FALSE])
          m <- lm(y~T-1, weights=w)
   sigterms <- abs(summary(m)$coefficients[, 3])>2
         T1 <- T[, sigterms]
        edf <- sum(sigterms)   
        fit <- if (edf==0) lm(y~1, weights=w) else lm(y~T1-1, weights=w)
        # if (edf==0) corrplot(vcov(gamlss(y~1, weights=w),"cor"))   
        # else corrplot(vcov(gamlss(y ~ T1 - 1, weights = w),"cor"))
         ii <- (1:MaxNo)*sigterms # the position on no zero beta
    expCoef <- rep(0,MaxNo)
expCoef[ii] <- coef(fit)   
       beta <- as.vector(PC$rotation[,1:MaxNo, drop = FALSE]%*%matrix(expCoef))   
        AiC <- AIC(fit, k=k)
  }
} else # df's
  {
      T <- as.matrix(PC$x[,1:edf, drop = FALSE])
    fit <- lm(y~T-1, weights=w)
    AiC <- AIC(fit, k=k)
   beta <- as.vector(PC$rotation[,1:edf]%*%matrix(coef(fit)))
  }  
   
coefSmo <- list( coef = coef(fit),
                 beta = beta,
                   pc = PC,
                  edf = edf,
                  AIC = AiC
                  )        
   class(coefSmo) <- "pc"
list(fitted.values = fitted(fit), residuals = resid(fit), var=(predict(fit, se.fit=TRUE)$se)^2, 
            nl.df = edf,  lambda=0,  coefSmo = coefSmo)
}  else # if xeval
    {
      newdata <- data.frame(attr(x, "newData"))
        names <- as.character(attr(x, "name"))
          nX  <- as.matrix(newdata[,names])
        scale <- attr(x, "scale")[[1]]
       center <- attr(x, "center")[[1]]
         beta <- attr(x, "beta")
         pred <- scale(nX, center=center,  scale=scale )%*%beta
        pred
    } 
}
#------------------------------------------------------------
#------------------------------------------------------------
plot.pc <- function(x,...)
{
  plot(x[["beta"]], type="h", xlab="knots", ylab="coefficients")
  abline(h=0)
}
#-------------------------------------------------------------
coef.pc <- function(object, ...)
{
  object[["coef"]]
}
#-------------------------------------------------------------
#-------------------------------------------------------------
print.pc  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
  cat("Principal componet fit using the gamlss function pc() \n")
  cat("Degrees of Freedom for the fit :", x$edf, "\n")
}
#-------------------------------------------------------------
#-------------------------------------------------------------
##############################################################
##############################################################
##############################################################
# This function checking a data matrix for
# variables with hight pairwise correlation between then. 
# The two arguments of the function are 
# i) a data sets and 
# ii)  a correlation value acting as a lower limit  
#      (all pair wise  correlation above this linit will be reported) 
# The function outputs an 3 columns matrix
# with the names of hight correlated pairs and their correlation 
#-------------------------------------------------------------
which.Data.Corr <- function(data, r=.90, digits=3)
{
  if (abs(r)>=1||abs(r)<=0) stop("r should be greater than  0 and lass than 1")
  
  daTa <- subset(data,  select=ifelse(sapply(data,is.factor)|sapply(data,is.character)==TRUE, FALSE, TRUE))
  Dim <- dim(daTa)
  CC <- cor(daTa)
  CC <- base::round(x = CC, digits = digits)
  CCC <- CC-diag(rep(1,Dim[2]))
  if (is.null(colnames(daTa))) colnames(daTa) <- paste0("X", seq(1:dim(data)[2]))
  if (!any(which(abs(CCC)>r))) return(cat("no correlation above", r, "\n"))
  mm <- which(abs(CCC)>r, arr.ind=T)
  nn <- mm[mm[,1]< mm[,2],]
  if (is.vector(nn))
  {
    name1 <- colnames(data)[nn[1]]
    name2 <- colnames(data)[nn[2]]
    corrs <- CCC[nn[1],nn[2]]
  } else
  { name1 <- colnames(data)[nn[,1]]
  name2 <- colnames(data)[nn[,2]]
  corrs <- CCC[nn]
  }
  cbind(name1, name2, corrs)
}
##############################################################
##############################################################
##############################################################
# This function takes the response and the x-variables
# computes their correlation and
# select variables with correlation highter that the "r"
# it is design as a way to preselect variables 
# (before a PCR analysis)
##############################################################
##############################################################
##############################################################
which.yX.Corr <- function(y, x, 
                             r =.50 , 
                          plot = TRUE, 
                  hierarchical = TRUE,
                         print = TRUE,
                        digits = 3)
{
  # get the correlations  
  CC  <- cor(y, x)
   CC <- base::round(x = CC, digits = digits)
  # if plotting is TRUE plot them 
  if (plot)
  {
    plot(as.vector(CC), pch=20, col="gray")
    abline(h=c(r, -r), col='red') 
  }
  # if   
  if(hierarchical)
  {
        DF <- as.data.frame(x)  
    nnames <- colnames(CC)[abs(CC)>r]
        ff <- as.formula(paste("~ ",paste(nnames, collapse='+')),) 
    if (any(grep(":", nnames))) ff <- eval(call("~",ff))
     # ff <- as.formula(gsub(":", "*", as.character(ff)))
        MM <- model.matrix(ff, DF)[,-1]
  } else
  {
    nnames <- colnames(CC)[abs(CC)>r]
    MM <- x[,nnames]
  } 
  if (print)
  {
    cat("cup point for correlation", r, "\n")
    cat("dimesions of new matrix", dim(MM), "\n")
  }
  invisible(MM)
}
##############################################################
##############################################################
##############################################################
#-------------------------------------------------------------
# this is a function which allows to save the SVD in 
# advance before fitting
# it saves time if all the parameters need SVD
#--------------------------------------------------------------
getSVD <- function(x=NULL,  nu=min(n, p), nv=min(n, p))
{
  n <- dim(x)[1]
  p <- dim(x)[2]
  x <- scale(x)
  center <- attr(x,"scaled:center")  
  scale <- attr(x,"scaled:scale")
  cNames <- colnames(x)
  S <- La.svd(x, nu=nu, nv=nv)
  #colnames(S$vt)  <- cNames
  S$X <- x
  S$center <- attr(x,"scaled:center")  
  S$scale <- attr(x,"scaled:scale")
  S
}
##############################################################
##############################################################
##############################################################
