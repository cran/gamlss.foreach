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
pcr <- function (x.vars=NULL, x=NULL, df = NULL, 
                 M=min(p,n),  k=log(n), r=0.2,
                 method=c("GAIC", "t-values", "SPCR")) 
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
gamlss.environment <- sys.frame(position) #gamlss or predict.gamlss
# get the data
    if (sys.call(position)[1]=="predict.gamlss()") 
      { # if predict is used 
      object <- get("object", envir=gamlss.environment)
        beta <- getSmo(object)$beta
   centerold <- getSmo(object)$pc["center"]
    scaleold <- getSmo(object)$pc["scale"]
     newData <- get("newdata", envir=gamlss.environment)
        Data <- get("data", envir=gamlss.environment)
      } else if (sys.call(position)[1]=="gamlss()") { # if gamlss() is used
      if (is.null(get("gamlsscall", envir=gamlss.environment)$data)) { # if no data argument but the formula can be interpreted
        Data <- data.frame(cbind(x))	
      } else {# data argument in gamlss 
        Data <- get("gamlsscall", envir=gamlss.environment)$data
      }
    } else {
        Data <- get("data", envir=gamlss.environment)
    }
        Data <- data.frame(eval(substitute(Data)))
# finish data 
if    (sum(c(!is.null(x), !is.null(x.vars))>1L))
    stop("use only one of the arguments x or x.vars")
if (!is.null(x.vars)) 
  {
    Xmtrx <- as.matrix(Data[, x.vars])
   } else {
    Xmtrx <- x
   }  
       n <- dim(Xmtrx)[1]
       p <- dim(Xmtrx)[2]
  cNames <- colnames(Xmtrx)
if (!is.null(df)&&df>min(n,p)) stop("the df have to be less than the number of variables")
if (!is.null(df)) M <- df # if df is set only to few component  
     sl <- sample(letters, 4)
fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
## put the starting values in the gamlss()environment
#--------
assign(startLambdaName, 0, envir=gamlss.environment)  
## this comes from prcomp ---------------------------
                      x <- rep(0, n)
        attr(x, "name") <- cNames
           attr(x, "M") <- M
           attr(x, "X") <- Xmtrx
        attr(x, "call") <- substitute(gamlss.pcr(data[[scall]], z, w))
          attr(x, "df") <- df
           attr(x, "k") <- k
           attr(x, "r") <- r
  attr(x, "gamlss.env") <- gamlss.environment
attr(x,"NameForLambda") <- startLambdaName      
      attr(x, "method") <- method 
     attr(x, "newData") <- newData 
              class(x) <- c("smooth", class(x))
    x
}
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
gamlss.pcr <- function(x, y, w, xeval = NULL, ...)
{
if (is.null(xeval))
 {    
              X <-  attr(x, "X")
              M <- attr(x, "M")    
         method <- attr(x, "method")
          names <- as.character(attr(x, "name"))
             df <- as.vector(attr(x, "df"))
              k <- as.vector(attr(x, "k"))
              r <- as.vector(attr(x, "r"))
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda"))   
#---------------------------------------------------
#
#  sup <- if (method=="SPCR") TRUE else FALSE
#  edf <- if (method=="t-values") M else df 
  #W <- if(all(w==w[1])) rep(1, length(y)) else w
  # w <- (w/ sum(w))*length(y)
   # op <- par(mfrow=c(2,1))
#cat(names(formals(get("f", envir=sys.frame(5))$valid)), "\n")
#cat("df", df,  "\n")
if (!is.null(df))
{
    fit <- fitPCR(x=X, y=y, weights=w, M=M, df=df, 
                plot=FALSE, r=r, k=k)
     fv <- fitted(fit, pc=df)
    edf <- fit$pc
    res <- fit$residuals[ ,df]
}  else
{
 if (method=="GAIC")
 {
    val <- if (method=="GAIC") get(startLambdaName, envir=gamlss.env)    
    fit <- fitPCR(x=X, y=y, weights=w, M=M, plot=FALSE, k=k)
    tdf <- fit$pc
    cc <- ceiling(M/10)
if(val!=0&&abs(val-tdf)>cc) #the difference too big
    {
      fv <- fitted(fit, pc=val)
     edf <- val
     res <- fit$residuals[ ,pc=val]   
  assign(startLambdaName, val, envir=gamlss.env) 
    }
else
    {
     fv <- fitted(fit, pc=fit$pc)
    edf <- fit$pc
    res <- fit$residuals[ ,pc=fit$pc] 
    assign(startLambdaName, tdf, envir=gamlss.env)
    }  
 } 
 if (method=="t-values")
  {
################################################################
       edf <-  M  
       fit <- fitPCR(x=X, y=y, weights=w, M=M, df=edf, 
                     plot=FALSE, r=r, k=k)       
getSigTerms <- function(obj, x=NULL, val=1.96)
  {
if (is.null(x)) stop("the function need the original x variables")
     t.val <- coef(obj)/obj$se.gamma
  whichSig <- as.numeric(abs(t.val) > val) 
     gaMma <- coef(obj, param="gamma", pc=obj$M)*whichSig 
      beta <-  obj$loadings%*%gaMma 
     scale <- obj$scale
    center <- obj$center
        fv <- as.vector(scale(x, center=center,  scale=scale)%*%beta)+obj$mean.y
       out <- list(fitted.values=fv, coefficients=as.vector(beta), gamma=gaMma) 
     out
  }
################################################################t
    tresult <- getSigTerms(fit, x=X, val=1.96)
         fv <- tresult$fitted.values
        edf <- sum( tresult$gamma!=0)
#  cat("df", df,  "\n")
        res <- y-fv
  fit$gamma <- tresult$gamma
fit$newbeta <- tresult$coefficients
   }
 if (method=="SPCR")
 {
   fit <- fitPCR(x=X, y=y, weights=w, M=M, df=df, supervised=TRUE, 
                 plot=FALSE, r=r, k=k)   
    fv <- fitted(fit)
    edf <- fit$pc
    res <- fit$residuals[ ,fit$pc]
  }
}  
   list(fitted.values = fv, residuals = res, var = fv, 
        nl.df = edf, lambda = edf,  coefSmo = fit)  
}  else # if xeval
    {
      # THIS TO BE DONE
      newdata <- data.frame(attr(x, "newData"))
        names <- as.character(attr(x, "name"))
          nX  <- as.matrix(newdata[,names])
        scale <- attr(x, "scale")[[1]]
       center <- attr(x, "center")[[1]]
         beta <- coef(fit, "beta")
         pred <- scale(nX, center=center,  scale=scale )%*%beta
        pred
    } 
}
#------------------------------------------------------------
#------------------------------------------------------------
# plot.pc <- function(x,...)
# {
#   plot(x[["beta"]], type="h", xlab="knots", ylab="coefficients")
#   abline(h=0)
# }
# #-------------------------------------------------------------
# coef.pc <- function(object, ...)
# {
#   object[["coef"]]
# }
# #-------------------------------------------------------------
# #-------------------------------------------------------------
# print.pc  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
# {   
#   cat("Principal componet fit using the gamlss function pc() \n")
#   cat("Degrees of Freedom for the fit :", x$edf, "\n")
# }
# #-------------------------------------------------------------
# #-------------------------------------------------------------
# ##############################################################
# ##############################################################
# ##############################################################
# # This function checking a data matrix for
# # variables with hight pairwise correlation between then. 
# # The two arguments of the function are 
# # i) a data sets and 
# # ii)  a correlation value acting as a lower limit  
# #      (all pair wise  correlation above this linit will be reported) 
# # The function outputs an 3 columns matrix
# # with the names of hight correlated pairs and their correlation 
# #-------------------------------------------------------------
# which.Data.Corr <- function(data, r=.90)
# {
#      if (abs(r)>=1||abs(r)<=0) stop("r should be greater than  0 and lass than 1")
#      Dim <- dim(data)
#       CC <- cor(data)
#      CCC <- CC-diag(rep(1,Dim[2]))
# if (is.null(colnames(data))) colnames(data) <- paste0("X", seq(1:dim(data)[2]))
# if (!any(which(abs(CCC)>r))) stop(cat("no correlation above", r, "\n"))
#       mm <- which(abs(CCC)>r, arr.ind=T)
#       nn <- mm[mm[,1]< mm[,2],]
#   if (is.vector(nn))
#   {
#     name1 <- colnames(data)[nn[1]]
#     name2 <- colnames(data)[nn[2]]
#     corrs <- CCC[nn[1],nn[2]]
#   } else
#   { name1 <- colnames(data)[nn[,1]]
#    name2 <- colnames(data)[nn[,2]]
#    corrs <- CCC[nn]
#   }
#   cbind(name1, name2, corrs)
# }
# ##############################################################
# ##############################################################
# ##############################################################
# # This function takes the response and the x-variables
# # computes their correlation and
# # select variables with correlation highter that the "r"
# # it is design as a way to preselect variables 
# # (before a PCR analysis)
# ##############################################################
# ##############################################################
# ##############################################################
# which.yX.Corr <- function(y, x, 
#                              r =.50 , 
#                           plot = TRUE, 
#                   hierarchical = TRUE,
#                          print = TRUE)
# {
#   # get the correlations  
#   CC  <- cor(y, x)  
#   # if plotting is TRUE plot them 
#   if (plot)
#   {
#     plot(as.vector(CC), pch=20, col="gray")
#     abline(h=c(r, -r), col='red') 
#   }
#   # if   
#   if(hierarchical)
#   {
#     DF <- as.data.frame(x)  
#     nnames <- colnames(CC)[abs(CC)>r]
#     ff <- as.formula(paste("~ ",paste(nnames, collapse='+')),) 
#     if (any(grep(":", nnames))) ff <- as.formula(gsub(":", "*", as.character(ff)))
#     MM <- model.matrix(ff, DF)[,-1]
#   } else
#   {
#     nnames <- colnames(CC)[abs(CC)>r]
#     MM <- x[,nnames]
#   } 
#   if (print)
#   {
#     cat("cup point for correlation", r, "\n")
#     cat("dimesions of new matrix", dim(MM), "\n")
#   }
#   invisible(MM)
# }
# ##############################################################
# ##############################################################
# ##############################################################
# #-------------------------------------------------------------
# # this is a function which allows to save the SVD in 
# # advance before fitting
# # it saves time if all the parameters need SVD
# #--------------------------------------------------------------
# getSVD <- function(x=NULL,  nu=min(n, p), nv=min(n, p))
# {
#   n <- dim(x)[1]
#   p <- dim(x)[2]
#   x <- scale(x)
#   center <- attr(x,"scaled:center")  
#   scale <- attr(x,"scaled:scale")
#   cNames <- colnames(x)
#   S <- La.svd(x, nu=nu, nv=nv)
#   #colnames(S$vt)  <- cNames
#   S$X <- x
#   S$center <- attr(x,"scaled:center")  
#   S$scale <- attr(x,"scaled:scale")
#   S
# }
# ##############################################################
# ##############################################################
# ##############################################################
