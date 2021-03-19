# this function is imitates the function svdpc.fit()
# of package pls
# it gets the coefficients, fitted values and residual 
# for the whole path of models 
# from the 1st PC to the M  principal componet 
############################################################
############################################################
# functions 
# 1) fitPCR(x,y,M,nu,nv) 
# 2) plot.PCR() (OK)
# 3) fitted.PCR() (OK)
# 4) coef.PCR(() (OK)
# 5) summary.PCR (OK)
# 6) print.PCR (OK)
# 7) predict (OK) does it se's?
############################################################
############################################################
# COMMENTS
#   IS THE AIC CALCULATED PROPERLY? 
#   IT SEEMS THAT IT IS SCALED RELATED 
#   FOR EXAPLE IF I CHANGED THE (CONTANT) WEIGHT 
#   THE RESULTS ARE DIFFERENT 
############################################################
############################################################
# TO DO
# i)    It needs testing on whether produce the right coefficients 
#       It looks that the coef fitted values and residuals are OK
# ii)   put stanard errors  OK 
# iii)  put local AIC       OK 
# iv)   add weights         OK needs testing
# v)    add supervised      Not yet
# vi)   add argument to do svd outside (NO maybe not appropriate) 
# vii)  add foreach        needs testing (I am not sute that spead up things) 
# viii) prediction         OK (need more checking)
#############################################################
############################################################# 
fitPCR<- function( x = NULL,      # the x matrix
                   y = NULL,      # the response
             weights = rep(1,n),  # prior weights
                   M = NULL,# maximum no of componets to fits
                  df = NULL,
      #           nu = min(n, p), # 
      #           nv = min(n, p),
          supervised = FALSE,
                   k = 2,         # For the AIC
                   r = 0.2,       # a value between 0 and 1 for supervised
                plot = TRUE)      # wheter to pplot the path
{
##############################################################
############################################################## 
weighted.sd <- function(x, w=rep(1,length(x)), na.rm = FALSE) {
    if (na.rm) {
      w <- w[i <- !is.na(x)]
      x <- x[i]
    }
    sum.w <- sum(w)
    sqrt((sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2)))
  }
##############################################################
##############################################################  
scall <- deparse(sys.call(), width.cutoff = 200L) 
  if (is.null(y)) stop("the y is needed") 
     n <- dim(x)[1] # no of observations
     p <- dim(x)[2] # no of variables
     M <- if(is.null(M)) min(n,p) else M
  if(!is.null(df)) M <- df
  if((n>=p)&&(M > p)){
       M <- p 
       warning(cat("M is reset to", p, '\n'))  
     }
  if((n <p)&&(M > n)){
       M <- n 
       warning(cat("M is reset to", n, '\n')) 
     } 
mean.y <- rep(weighted.mean(y, weights), n)
     Y <- as.matrix(sqrt(weights)*y)
 if (is.null(rownames(x))) colnames(x) <- paste0("X", seq(1:dim(x)[2]))
# IS SCALLING SHOULD BE DONE WITH WEIGHTED MEAN AND VARIANCE
# Is this  correcting SOME OF THE PROBLEMS in forecasting????
#############################################################
# (a) scale 
#   weights <- (weights/sum(weights))*n
#   always scale 
     x <- scale(x, center=apply(x,2,weighted.mean, na.rm=TRUE, w=weights),
                 scale=apply(x,2,weighted.sd, na.rm=TRUE, w=weights))
center <- attr(x,"scaled:center")   # we need this for prediction 
 scale <- attr(x,"scaled:scale")
############################################################# 
# (b) transform
     X <- sqrt(weights)*x 
############################################################# 
# (c) supervised
if (supervised) 
  {
     X <- which.yX.Corr(Y, X, r=r, plot=FALSE, print = FALSE)
cnames <- colnames(X)
     p <- dim(X)[2]
     M <- ifelse(n>=p, p, n)
  }
#############################################################       
#    to save 
      B <- matrix(0, nrow=p, ncol=M) # the coefficients of the M fits 
     Fi <- matrix(0,  nrow=n, ncol=M)# the fitted values of the M fits
    AiC <- gamMa <- seGamma <- cpT <- rep(0, M) 
 cNames <- colnames(X)               # the colunm names
############################################################
# (d) svd 
   Xsvd <- La.svd(X, nu=M, nv=M)     
      D <- Xsvd$d[1:M]
      T <- Xsvd$u[,1:M, drop = FALSE] %*% diag(D, nrow = M)# scores
      P <- t(Xsvd$vt[1:M, , drop = FALSE]) # loadings
  gamMa <- crossprod(T,Y)/D^2   # the gamma coef
for (m in 1:M)
  { 
  B[,m] <- P[, 1:m, drop = FALSE] %*% gamMa[1:m]  # the betas
 Fi[,m] <- if (supervised)  x[,cnames] %*% B[,m]+ mean.y
           else             x %*% B[,m]+ mean.y
#Res[,m] <- y-Fi[,m] 
  }
residuals <- c(y) - Fi                             # the residuals
    sig <- apply(residuals, 2, "weighted.sd", w=weights) # all estimated sigmas  
#  with foreach
if (is.null(df))
{
  AiC <- foreach(m = 1 : M, .combine=rbind)%dopar%{
    di <- -2 * dNO(y, mu=Fi[,m], sigma = sig[m], log=TRUE)
    AiC <- sum(di) + k * (m+2)}#  
  # AIC3 <- foreach(m = 1 : M, .combine=rbind)%dopar%{
  #     sigm <- sum(residuals[,m]^2) /n
  #  AIC3 <- n * log(sigm) + k * (m+2)}# this produce identical results to AiC
  }
# the se of gamma  at M  
    se.gamma <- sig[M] / D
 rownames(B) <- cNames
       gamMa <- as.vector(gamMa)
names(gamMa) <- paste0("PC", seq(1:M))
if (plot)
  {
    plot(B[1,]~seq(1,M), ylim=c(min(B),max(B)), 
         type="n", xlab="number of PC", ylab="coef")
    for (i in 1:M)
      lines(B[i,]~seq(1,M), col=i)
  if (is.null(df))  abline(v=which.min(AiC), col="gray")
  }
#  output  
out <- list(coefficients = B, 
                  scores = T, 
                loadings = P, 
                   gamma = gamMa, 
                se.gamma = se.gamma, 
                  center = center, 
                   scale = scale, 
           fitted.values = Fi,  
                    Xvar = D^2, 
                    gaic = AiC, 
                      pc = {if (is.null(df)) which.min(AiC) else df},
                       k = k,
                       M = M,
                  mean.y = mean.y[1],
                   dim.x = dim(x),
                   sigma = sig, 
               residuals = residuals,
                    call = scall)
  class(out) <- "PCR"
  out
}  
################################################################
################################################################
plot.PCR <- function(x, type=c("path", "gaic"), labels=TRUE, cex=0.8,...)
{
  type <- match.arg(type)  
if (type=="path")
{
M <- length(x$gamma)
plot( x$coefficients[1,]~seq(1,M), 
      ylim=c(min(x$coefficients),max(x$coefficients)),
      xlim=c(0, M+1), 
         type="n", xlab="number of PC", ylab="coef")
for (i in 1:M) 
  {lines(x$coefficients[i,]~seq(1,M), col=i)
  if (labels) text(x=M+1, y=x$coefficients[i,M], 
                   labels=rownames(x$coefficients)[i], cex=cex)} 
    abline(v=x$pc, col="gray")
}
else
{
plot(x$gaic, pch=20, xlab="number of componets", ylab="local AIC")
abline(v=x$pc, col="gray")  
}  
}
################################################################
################################################################
print.PCR <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\n A fitted PCR model","\n")
  #cat("\nFamily: ", deparse(x$family), "\n")
  #cat("Fitting method:", deparse(x$method), "\n")
  cat("\nCall: ", deparse(x$call, width.cutoff = 50), "\n", 
      fill = TRUE)
  cat("With gamma coefficients")
  cat(":\n")
    print.default(format(coef(x), digits = digits), 
                  print.gap = 2, quote = FALSE)
  invisible(x)
  }
################################################################
################################################################
summary.PCR <- function(object, param = c("gamma", "beta"),
                                   pc = object$pc, ...)
{
  param <- match.arg(param)  
  if (param=="gamma")
  {
    cat("The gamma coefficients of the", object$M, "Principal Componets","\n")
    return(cbind(coef=coef(object), se=object$se.gamma, t.val=coef(object)/object$se.gamma))  
  }
  else  
  {
    cat("The beta coefficients at", pc, "Principal Componets","\n")
     beta <- coef(object, param="beta", pc=pc)
  se.beta <- vcov(object, type="se", pc=pc)# the se of beta
 return(cbind(beta, se=se.beta, t.val=beta/se.beta)) 
  }  
}
 
################################################################
################################################################
# fitted.PCR
fitted.PCR <- function(object, pc=object$pc, ...)
{
object$fitted.values[, pc=pc]  
}
################################################################
################################################################  
# fitted.PCR
residuals.PCR <- function(object, pc=object$pc, ...)
{
  object$residuals[, pc=pc]  
}
################################################################
################################################################  
# coef.PCR
coef.PCR <- function(object, param = c("gamma", "beta"),
                                pc = object$pc, ...)
{
  param <- match.arg(param)  
 if (param=="gamma")
 {
   return(object$gamma)  
 }
else
 {
  return(object$coefficients[,pc] )
 }
}
################################################################
################################################################
# vcov.PCR
vcov.PCR <- function(object, pc = object$pc, 
                      type = c("vcov", "cor", "se", "all"),
                     ...) 
{
      type <- match.arg(type)
         s <- object$sigma[pc] # sqrt(sum(object$residuals[,pc]^2)/(object$dim.x[1]))
     gamma <- object$gamma[1:pc]
  se.gamma <- s / sqrt(object$Xvar[1:pc])
  # sqrt(diag(crossprod(object$scores[,1:pc], 
  #                      object$scores[,1:pc])))# the se of gamma
m.se.gamma <- if (pc==1) matrix(se.gamma, ncol=1) else diag(se.gamma)
 vcoc.beta <- object$loading[,1:pc]%*%m.se.gamma%*%t(object$loading[,1:pc])
   se.beta <- sqrt(diag(vcoc.beta))
 corr.beta <- cov2cor(vcoc.beta)
switch(type, vcov = vcoc.beta, cor = corr.beta, se = se.beta,  
       all = list(se = se.beta, vcov = vcoc.beta, cor = corr.beta))
}  

################################################################
################################################################
# predict.PCR
predict.PCR <- function(object, newdata=NULL, pc=object$pc, ...)
{
  if (is.null(newdata))
  {
    fv <- fitted(object,  pc=pc)
    return(fv)
  } else
  {
   if (dim(newdata)[2]!=object$dim.x[2])
     stop("the dimensions of the new data is not the same with the old")
    #names <- names(coef(object,"beta"))
      nX  <- as.matrix(newdata) #as.matrix(newdata[,names])
    scale <- object$scale
   center <- object$center
     beta <- coef(object, param="beta",pc=pc)
     pred <- scale(nX, center=center,  scale=scale)%*%beta
     pred <- as.vector(pred) + object$mean.y
    return(pred)
  }  
}

################################################################
################################################################
################################################################
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
################################################################
################################################################
################################################################
# test <- function(X,M=NULL) 
# {
#   
#    n <- dim(X)[1] # no of observations
#    p <- dim(X)[2] # no of variables
#    M <- if(is.null(M)) min(n,p) else M
# if((n>=p)&&(M > p)){
#   M <- p 
#   warning("M is reset to p")  
# }
# if((n <p)&&(M > n)){
#   M <- n 
#   warning("M is reset to n") 
# } 
#    cat("M=",M,"\n")
#   ll <-  Xsvd <- La.svd(X, nu=M, nv=M)   
#   list(dim(ll$u), dim(ll$vt), length(ll$d))
# }
#   
# weighted.var <- function(x, w, na.rm = FALSE) {
#   if (na.rm) {
#     w <- w[i <- !is.na(x)]
#     x <- x[i]
#   }
#   sum.w <- sum(w)
#   sum.w2 <- sum(w^2)
#   mean.w <- sum(x * w) / sum(w)
#   (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
#                                        na.rm)
# }
# 
# weighted.var2 <- function(x, w, na.rm = FALSE) {
#   if (na.rm) {
#     w <- w[i <- !is.na(x)]
#     x <- x[i]
#   }
#   sum.w <- sum(w)
#   (sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2))
# }

# weighted.var <- function(x, w = NULL, na.rm = FALSE) {
#   if (na.rm) {
#     na <- is.na(x) | is.na(w)
#     x <- x[!na]
#     w <- w[!na]
#   }
#   
#   sum(w * (x - weighted.mean(x, w)) ^ 2) / (sum(w) - 1)
# }

# CASE 0
# sig <- sqrt(sum(residuals[,m]^2)/(n))
# di <- -2 * dNO(y, mu=Fi[,m], sigma = sig, log=TRUE) 
# AiC <- sum(di) + k * (m+2)}#constant plus sigma      
# CASE 1
# sig <- weighted.sd(y-Fi[,m], w=weights) #sqrt(sum(residuals[,m]^2)/(n))
#di <- -2 * dNO(y, mu=Fi[,m], sigma = sig[m], log=TRUE)
#AiC <- sum(weights*di) + k * (m+2)}#constant plus sigma
# CASE 2    
#           AiC <- sum(weights(y-Fi[,m])^2)+ k * (m+1)}
# CASE 3
#           sig2 <- sum(weights*residuals[,m]^2)/(n)}
#           AiC  <- n*log(sig2)+  k * (m+1)
# case 4
#           sig <- weighted.sd(residuals[,m], weights)            
#           AiC  <- n*log(sig^2)+  k * (m+1)}