#########################################################################################
#########################################################################################
#########################################################################################
# A set of function for doing bayesian and not parametric boosting 
require(gamlss.foreach)
#########################################################################################
#########################################################################################
#########################################################################################
# functions
# i)    coefAll()          which probably should go in the main gamlss
# ii)   BayesianBoost()    at the moment only coef and parameters are saved 
#                          not gammas or the fitted smoothers  
#                          TO DO 
#          tidy up what to save maybe I should have 
#         options save = c("coef", "parameters", "lambda", "gammas", "smoothers")
# iii)  summary.Bayesian.boot()   summary for the coefficients
#                          TO DO
#                     there is no summary for parameters at the moments
#       plot.Bayesian.boot()      plotting the distribution of the coefficients  
#                          TO DO
#                     maybe should be on option pages=1                           
# iii)  print.Bayesian.boot()     infromation about boost
# vi)   NonParametricBoot()      This only saves coef 
#                         TO DO 
#                         saved at least parameters  
#       summary.NonParametricBoot()
#       plot.NonParametricBoot()
# v)    print.NonParametricBoot()
# vi)   ratioBoot()
# vii)  momentBaysianBoot
# viii) moment Boot
#########################################################################################
# TODO
# # save original coefficients OK DONE
# check Kriton parallelism  
# add summary on the Bayesian.boot and classical.boot classes
#########################################################################################
#########################################################################################
#########################################################################################
# FUNCTION coefAll()
# Function 1
######################################################################################### 
# this probably shoulf go in the main gamlss 
# It went but it does not get lambdas yet 
# This is with lambda
############################################################## 
# TO DO expeptions with different smoothers 
# check 
# pbm 
# pbc 
# pbz
# cs
# scs 
# ri
# pcat
# gmrf
# pvc
# ga
# ba
# nn
# tr
# lo
# ma
# pc
# pcr
# gnet
####################################################################################
####################################################################################
####################################################################################  
coefAll1 <- function(object, lambdas= TRUE, deviance=FALSE, ...)
{
  out <- list()
if ("mu" %in% object$par) #
 {
    out$mu <- coef(object, "mu")  
    if (!is.null(object$mu.s)&&lambdas)
    {
      M <- dim(object$mu.s)[2]
      outLambda <- rep(0,M)
      for (i in 1:M)
      {
        outLambda[i] <-  getSmo(object, which=i)$lambda
      }
      out$mu <- c(out$mu, lambda=outLambda) 
    }
 }
if ("sigma" %in% object$par)#
 {
  out$sigma <- coef(object, "sigma")
  if (!is.null(object$sigma.s)&&lambdas)
    {
      M <- dim(object$sigma.s)[2]
      outLambda <- rep(0,M)
      for (i in 1:M)
      {
        outLambda[i] <-  getSmo(object, parameter="sigma", which=i)$lambda
      }
      out$sigma <- c(out$sigma, lambda=outLambda) 
    }
 }  
if ("nu" %in% object$par) #
 {
    out$nu <- coef(object, "nu")  
  if (!is.null(object$nu.s)&&lambdas)
    {
      M <- dim(object$nu.s)[2]
      outLambda <- rep(0,M)
      for (i in 1:M)
      {
        outLambda[i] <-  getSmo(object, parameter="nu", which=i)$lambda
      }
      out$nu <- c(out$nu, lambda=outLambda) 
    }
 }  
if ("tau" %in% object$par) #
  {
    out$tau <- coef(object, "tau")  
  if (!is.null(object$tau.s)&&lambdas)
    {
      M <- dim(object$tau.s)[2]
      outLambda <- rep(0,M)
      for (i in 1:M)
      {
        outLambda[i] <-  getSmo(object, parameter="tau", which=i)$lambda
      }
      out$tau <- c(out$tau, lambda=outLambda) 
    }
  } 
  if ("tau" %in% object$par)
    out$tau <-  coef(object, "tau")
  if (deviance) out$deviance <- deviance(object)
  return(out)
}
#########################################################################################
#########################################################################################
#########################################################################################
# this should work with any paramatric GAMLSS models
# FUNCTION BayesianBoot
#########################################################################################
#########################################################################################
#########################################################################################
BayesianBoot <- function(obj, data=NULL, B=100, newdata=NULL)
{
# local function ------------------------------------------------------
# Bubin 1981 genreration of weights   
rFun  <- function(n)
{
     U <- runif(n-1)
    oU <- U[order(U)]
    oU <- c(0,oU,1)
    g  <- diff(oU)
    g
}
#-----------------------------------------------------------------------
#pb <- tkProgressBar(max=B)
#progress <- function(i) setTkProgressBar(pb, i)
#progress <- function(n) cat(sprintf("task %d is complete\n", n))
#---------------------------------
#opts <- list(progress=progress)
#--------------------------------- ------------------------------------- 
if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
     data <- if ("data"%in%names(obj$call)) eval(obj$call$data)
               else if (!is.null(data)) data
               else stop("data are not defined")  
if (!is.data.frame(data)) stop("data is not a data.frame")
     N <- nN <- dim(data)[1]     # length of data
if (!is.null(newdata))
{
  if(!is.data.frame(newdata))  stop("newdata should be a data.frame") 
    nN <- dim(newdata)[1]
}

  thecall <- as.call(obj$call) # the call
thecall$control <- gamlss.control(trace=FALSE) # suppressing tracing
## getting the length of the coefficients 
orig.coef <- unlist(coefAll1(obj, deviance=TRUE))
        Q <- length(orig.coef)
      res <- foreach(i = 1:B, .packages="gamlss", .errorhandling = "remove",
                     .inorder = FALSE)%dopar%{
             bootstrap_weights <- rFun(N)*N #
               thecall$weights <- bootstrap_weights
                       model_1 <- eval(thecall) # how to catch?
      list(coef=coefAll1(model_1, deviance = TRUE), par=predictAll(model_1, newdata=newdata, output="matrix"))# output
                     }
           newB <- length(res)
  if (newB!=B) warning(cat(B-newB, "simulations failed only", newB, "were finished \n"))  
      COEF <- matrix(0, ncol=Q, nrow=newB)
  if ("mu" %in% obj$par)       MU <- matrix(0, ncol=newB, nrow=nN)
  if ("sigma" %in% obj$par) SIGMA <- matrix(0, ncol=newB, nrow=nN)
  if ("nu" %in% obj$par)       NU <- matrix(0, ncol=newB, nrow=nN)
  if ("tau" %in% obj$par)     TAU <- matrix(0, ncol=newB, nrow=nN)
# maybe I should do this with foreach      
    for (i in 1:newB)
    {
      COEF[i,] <- unlist(res[[i]][["coef"]])
      if ("mu" %in% obj$par)        MU[,i] <- unlist(res[[i]][["par"]])[,"mu"]
      if ("sigma" %in% obj$par)  SIGMA[,i] <- unlist(res[[i]][["par"]])[,"sigma"]
      if ("nu" %in% obj$par)        NU[,i] <- unlist(res[[i]][["par"]])[,"nu"]
      if ("tau" %in% obj$par)      TAU[,i] <- unlist(res[[i]][["par"]])[,"tau"]
      }
parameters <- switch(length(obj$par),
                     list(mu=MU),
                     list(mu=MU, sigma=SIGMA),
                     list(mu=MU, sigma=SIGMA, nu=NU),
                     list(mu=MU, sigma=SIGMA, nu=NU, tau=TAU))
colnames(COEF) <-   names(unlist(res[[1]][["coef"]]))
# is any column of COEF has all valus NA the omit it 
if (any(is.na(COEF)))
{
  ind <- rep(0,Q)
for (i in 1:Q)  ind[i] <- ifelse(all(is.na(COEF[,i])),-i, ind[1])
 COEF <- COEF[,ind]
}
#################################################################
#                           OUTPUT                              #
#################################################################  
    out <- list(boot = COEF, 
                   B = B,
               trueB = newB,
                 par = parameters,
           orig.coef = orig.coef,   
           orig.call = obj$call) 
class(out) <- list("Bayesian.boot")
out 
}
#########################################################################################
#########################################################################################
#########################################################################################
# summary 
summary.Bayesian.boot <- function(object, ...)
{
  cat("\n Coefficients from Baysian.boot() taken from:", object$trueB, "simulations \n")  
    tcoef <-  t(apply(object$boot,2,"quantile", probs=c(0.025, 0.50, 0.975), na.rm=TRUE))
    mcoef <-  apply(object$boot,2,"mean", na.rm=TRUE )
   # out <- cbind(mean=mcoef, tcoef, mle=x$orig.coef) 
    out <- cbind(mean=mcoef, tcoef, mle=object$orig.coef[names(object$orig.coef)%in%names(mcoef)]) 
    printCoefmat(out,...)
invisible(out)
}
#########################################################################################
#########################################################################################
#########################################################################################
# plot 
plot.Bayesian.boot <- function(x, par=1, add=FALSE,  ...)
{
  if (add) {lines(density(x$boot[,par]), ... )}
  else {
    plot(density(x$boot[,par]), main=colnames(x$boot)[par], ... )
    rug(x$boot[,par])
  }
  abline(v=x$orig.coef[par], col="gray")  
  abline(v=quantile(x$boot[,par], c(0.025, .5, .975)), col=gray(.9))
}
#########################################################################################
#########################################################################################
#########################################################################################
#  iii) print.Bayesian.print
print.Bayesian.boot <- function(x,...)
{
  cat("\n Results from Baysian.boot() taken from:", x$trueB, " simulations \n")
  cat("\n The simulated values are in:", paste0(deparse(substitute(x)), "$boot"), "\n") 
  cat("\n Simulation was done from the original call: \n")
  ccc<-x$orig.call
  print(ccc)
}  
#########################################################################################
#########################################################################################
#########################################################################################
# FUNCTION classical Boot
######################################################################################### 
 NonParametricBoot <- function(obj, data=NULL, B=100, newdata=NULL)
 {
 if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
     datA <- if ("data"%in%names(obj$call)) eval(obj$call$data)
             else if (!is.null(data)) data
             else stop("data are not defined")  
 if (!is.data.frame(datA)) stop("data is not a data.frame")
               N <- nN <- dim(datA)[1]     # length of data
               R <-  dim(datA)[2] 
 if (!is.null(newdata))
  {
   if(!is.data.frame(newdata))  stop("newdata should be a data.frame") 
              nN <- dim(newdata)[1]
   }              
         thecall <- as.call(obj$call) # the call
 thecall$control <- gamlss.control(trace=FALSE) # supressing tracing
       orig.coef <- unlist(coefAll1(obj, deviance=TRUE))
      names.coef <-  names(orig.coef)
               Q <- length(orig.coef) 
          result <- foreach(i = 1:B, .packages="gamlss", .errorhandling = "remove",
                        .export = c("coefAll1"), .inorder = FALSE)%dopar%
          {
               ii <- sample(N,N,replace=T) 
        boot_data <- as.data.frame(datA[ii,])
names( boot_data) <- names(datA)
# cat(i,"\n")
             mooo <- update(obj, data= boot_data ) 
               s  <- if (is.null(newdata))   list(coef=coefAll1(mooo, deviance = TRUE)) #
                      else  list(coef=coefAll1(mooo, deviance = TRUE), par=predictAll(mooo,newdata=newdata, output="matrix")) 
            return(s)
          }
#  bootstrap_data <- data[sample(nrow(data),nrow(data),replace=T),]
#    thecall$data <- bootstrap_data
#         model_1 <- eval(thecall) # how to catch?
#         list(coef=coefAll1(model_1, deviance = TRUE), par=predictAll(model_1,newdata=newdata, output="matrix"))                   }
           newB <- length(result)
   if (newB!=B) warning(cat(B-newB, "simulation failed only", newB, "were finished \n"))        

           COEF <- matrix(0, ncol=Q, nrow=newB)
  if ("mu" %in% obj$par)       MU <- matrix(0, ncol=newB, nrow=nN)
  if ("sigma" %in% obj$par) SIGMA <- matrix(0, ncol=newB, nrow=nN)
  if ("nu" %in% obj$par)       NU <- matrix(0, ncol=newB, nrow=nN)
  if ("tau" %in% obj$par)     TAU <- matrix(0, ncol=newB, nrow=nN)
# # maybe I should do this with foreach      
 for (i in 1:newB)  COEF[i,] <- unlist(result[[i]][["coef"]])
  colnames(COEF) <- names.coef      
  if (!is.null(newdata))
  {
    for (i in 1:newB) 
    {
      if ("mu" %in% obj$par)        MU[,i] <- unlist(result[[i]][["par"]])[,"mu"]
      if ("sigma" %in% obj$par)  SIGMA[,i] <- unlist(result[[i]][["par"]])[,"sigma"]
      if ("nu" %in% obj$par)        NU[,i] <- unlist(result[[i]][["par"]])[,"nu"]
      if ("tau" %in% obj$par)      TAU[,i] <- unlist(result[[i]][["par"]])[,"tau"]
    }
    parameters <- switch(length(obj$par),
                       list(mu=MU),
                       list(mu=MU, sigma=SIGMA),
                       list(mu=MU, sigma=SIGMA, nu=NU),
                       list(mu=MU, sigma=SIGMA, nu=NU, tau=TAU))
 colnames(COEF) <-   names(unlist(result[[1]][["coef"]]))
    } else  parameters <- NULL
 # is any column of COEF has all valus NA the omit it 
 if (any(is.na(COEF)))
       {
        ind <- rep(0,Q)
        for (i in 1:Q)  ind[i] <- ifelse(all(is.na(COEF[,i])),-i, ind[1])
        COEF <- COEF[,ind]
       }
# #################################################################
# #                           OUTPUT                              #
# #################################################################  
 out <- list(boot = COEF, 
             B = B,
             trueB = newB,
             par = parameters,
             orig.coef = orig.coef,   
             orig.call = obj$call) 
 class(out) <- list("NonParametric.Boot")
 out 
 }
# #########################################################################################
# #########################################################################################
# ######################################################################################### 
 summary.NonParametric.Boot <- function(object, ...)
 {
   cat("\n Coefficients from nonparam.boot() taken from:", object$trueB, "simulations \n")  
   tcoef <-  t(apply(object$boot,2,"quantile", probs=c(0.025, 0.50, 0.975)))
   mcoef <-  apply(object$boot,2,"mean", na.rm=TRUE )
   out <- cbind(mean=mcoef, tcoef, mle=object$orig.coef[names(object$orig.coef)%in%names(mcoef)]) 
   printCoefmat(out,...)
   invisible(out)
 }
# ######################################################################################
# ######################################################################################
# ######################################################################################
# summary.NonParametric.Boot <- function(x, ...)
# {
#   cat("\n Coefficients from NonParametric.boot() taken from:", x$trueB, "simulations \n")  
#   tcoef <-  t(apply(x$boot,2,"quantile", probs=c(0.025, 0.50, 0.975)))
#   mcoef <-  apply(x$boot,2,"mean", na.rm=TRUE )
#   out <- cbind(mean=mcoef, tcoef, mle=x$orig.coef[names(x$orig.coef)%in%names(mcoef)]) 
#   printCoefmat(out)
#   invisible(out)
# }
# #######################################################################################
# #######################################################################################
# #######################################################################################
# # plot 
 plot.NonParametric.Boot <- function(x, par=1, add=FALSE, ...)
 {
   if (add) {lines(density(x$boot[,par]), ... )}
   else {
     plot(density(x$boot[,par]), main=colnames(x$boot)[par], ... )
     rug(x$boot[,par])
   }
   
   abline(v=x$orig.coef[par], col="gray")  
   abline(v=quantile(x$boot[,par], c(0.025, .5, .975)), col=gray(.9))
 }
# ######################################################################################## 
# ######################################################################################## 
# ########################################################################################
# ########################################################################################  
 print.NonParametric.Boot <- function(x,...)
 {
   cat("\n Results from classical.boot() taken from:", x$trueB, " simulations \n")
   cat("\n The simulated values are in:", paste0(deparse(substitute(x)), "$boot"), "\n") 
   cat("\n Simulation was done from the original call: \n")
   ccc<-x$orig.call
   print(ccc)
 }  
# #########################################################################################
# #########################################################################################
# #########################################################################################
# #########################################################################################
# ratioBoot <- function(y,x, B=100)
# {
#   rFun  <- function(n)
#   {
#     U <- runif(n-1)
#     oU <- U[order(U)]
#     oU <- c(0,oU,1)
#     g  <- diff(oU)
#     g
#   }
# rB <- sum(y)/sum(x)
#  N <- length(y)
# res <- foreach(i = 1:B, .packages="gamlss")%dopar%{           
#   g <- rFun(N)
#   sum(g*y)/sum(g*x)
# }
# res <- unlist(res)
# res    
# }
# #########################################################################################
# #########################################################################################
# #########################################################################################
# momentBayesianBoot <- function(y, B=100)
# {
#   rFun  <- function(n)
#   {
#     U <- runif(n-1)
#     oU <- U[order(U)]
#     oU <- c(0,oU,1)
#     g  <- diff(oU)
#     g
#   }
#       mu <- mean(y)
#    sigma <- sd(y)
#    allSK <- momentSK(y)
#     skew <- allSK$mom.skew
#     kurt <- allSK$mom.kurt-3
# original <- matrix(c(mu, sigma, skew, kurt), nrow=1, ncol=4, byrow=TRUE) 
#   colnames(original) <- c("mean", "sd", "skew", "kurt") 
#       N <- length(y)
# res <- foreach(i = 1:B, .packages="gamlss")%dopar%{           
#      g <- rFun(N)
#     m1 <- sum(g*y)
#     m2 <- sum(g*(y-m1)^2)
#     m3 <- sum(g*(y-m1)^3)
#     m4 <- sum(g*(y-m1)^4)
# sigma1 <- sqrt(m2)
#  gamma1 <- m3/(m2^1.5)
#  beta2 <- m4/(m2^2)
#     list(mean=m1, sd=sigma1, skew=gamma1, kurt=beta2-3)
#   }
#   res <- matrix(unlist(res), nrow=B, ncol=4, byrow=TRUE) 
#   colnames(res) <- c("mean", "sd", "skew", "kurt") 
#   list(original=original, boot=res)
# }
# #########################################################################################
# #########################################################################################
# #########################################################################################
# momentBoot <- function(y, B=100)
# {
#       mu <- mean(y)
#    sigma <- sd(y)
#    allSK <- momentSK(y)
#     skew <- allSK$mom.skew
#     kurt <- allSK$mom.kurt-3
# original <- matrix(c(mu, sigma, skew, kurt), nrow=1, ncol=4, byrow=TRUE) 
# colnames(original) <- c("mean", "sd", "skew", "kurt") 
#        N <- length(y)
#        g <- rep(1/N, N)
#      res <- foreach(i = 1:B, .packages="gamlss")%dopar%{   
#             bootstrap_y <- y[sample(N, replace=T)]
#                       y <-  bootstrap_y
#                      m1 <- sum(g*y)
#                      m2 <- sum(g*(y-m1)^2)
#                      m3 <- sum(g*(y-m1)^3)
#                      m4 <- sum(g*(y-m1)^4)
#                  sigma1 <- sqrt(m2)
#                  gamma1 <- m3/(m2^1.5)
#                   beta2 <- m4/(m2^2)
#     list(mean=m1, sd=sigma1, skew=gamma1, kurt=beta2-3)
#   }
#   res <- matrix(unlist(res), nrow=B, ncol=4, byrow=TRUE) 
#   colnames(res) <- c("mean", "sd", "skew", "kurt") 
#   list(original=original, boot=res)
# }
# #########################################################################################
# #########################################################################################
# #########################################################################################
# # chooseDistBayesianBoot <- function(obj, B=100)
# # {
# # #-------------------------------------  
# #   rFun  <- function(n)
# #   {
# #     U <- runif(n-1)
# #     oU <- U[order(U)]
# #     oU <- c(0,oU,1)
# #     g  <- diff(oU)
# #     g
# #   }
# # #-------------------------------------  
# #  
# #   N <- length(obj$y)
# #   res <- foreach(i = 1:B, .packages="gamlss")%dopar%{           
# #     g <- rFun(N)
# #     sum(g*y)/sum(g*x)
# #   }
# #   res <- unlist(res)
# #   res    
# # }
# #########################################################################################
# #########################################################################################
# #########################################################################################
# # Testinf whether I cam pou a progress bar 
# BayesianBootT <- function(obj, data=NULL, B=100)
# {
#   # local function ------------------------------------------------------
#   # Bubin 1981 genreration of weights   
#   rFun  <- function(n)
#   {
#     U <- runif(n-1)
#     oU <- U[order(U)]
#     oU <- c(0,oU,1)
#     g  <- diff(oU)
#     g
#   }
#   #-----------------------------------------------------------------------
#   #pb <- tkProgressBar(max=B)
#   #progress <- function(i) setTkProgressBar(pb, i)
#   #progress <- function(n) cat(sprintf("task %d is complete\n", n))
#   #---------------------------------
#   #opts <- list(progress=progress)
#   #--------------------------------- ------------------------------------- 
#   if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
#   data <- if ("data"%in%names(obj$call)) eval(obj$call$data)
#   else if (!is.null(data)) data
#   else stop("data are not defined")  
#   if (!is.data.frame(data)) stop("data is not a dataframe")
#   N <- dim(data)[1]     # length of data
#   thecall <- as.call(obj$call) # the call
#   thecall$control <- gamlss.control(trace=FALSE) # supressing tracing
#   ## getting the length of the coefficients 
#   orig.coef <- unlist(coefAll1(obj, deviance=TRUE))
#   Q <- length(orig.coef)
#   BB <- 1:B
#   p <- progressor(along = BB)
#   # my_fcn <- function(B=10) {
#   #   X<- rnorm(B)
#   #   p <- progressor(along = X)
#   #   y <- foreach(i = 1:B) %dopar% {
#   #     Sys.sleep(1.)
#   #     p(sprintf("i=%g", i))
#   #     X[i]
#   #   }
#   #   
#   res <- foreach(i = 1:B, .packages="gamlss", .errorhandling = "remove",
#                  .inorder = FALSE)%dopar%{
#                    bootstrap_weights <- rFun(N)*N #
#                    thecall$weights <- bootstrap_weights
#                    model_1 <- eval(thecall) # how to catch?
#                    p(sprintf("i=%g", i))
#                    list(coef=coefAll1(model_1, deviance = TRUE), par=predictAll(model_1,output="matrix"))# output
#                  }
#   newB <- length(res)
#   if (newB!=B) warning(cat(B-newB, "simulations failed only", newB, "were finished \n"))  
#   COEF <- matrix(0, ncol=Q, nrow=newB)
#   if ("mu" %in% obj$par)       MU <- matrix(0, ncol=newB, nrow=N)
#   if ("sigma" %in% obj$par) SIGMA <- matrix(0, ncol=newB, nrow=N)
#   if ("nu" %in% obj$par)       NU <- matrix(0, ncol=newB, nrow=N)
#   if ("tau" %in% obj$par)     TAU <- matrix(0, ncol=newB, nrow=N)
#   # maybe I should do this with foreach      
#   for (i in 1:newB)
#   {
#     COEF[i,] <- unlist(res[[i]][["coef"]])
#     if ("mu" %in% obj$par)        MU[,i] <- unlist(res[[i]][["par"]])[,"mu"]
#     if ("sigma" %in% obj$par)  SIGMA[,i] <- unlist(res[[i]][["par"]])[,"sigma"]
#     if ("nu" %in% obj$par)        NU[,i] <- unlist(res[[i]][["par"]])[,"nu"]
#     if ("tau" %in% obj$par)      TAU[,i] <- unlist(res[[i]][["par"]])[,"tau"]
#   }
#   parameters <- switch(length(obj$par),
#                        list(mu=MU),
#                        list(mu=MU, sigma=SIGMA),
#                        list(mu=MU, sigma=SIGMA, nu=NU),
#                        list(mu=MU, sigma=SIGMA, nu=NU, tau=TAU))
#   colnames(COEF) <-   names(unlist(res[[1]][["coef"]]))
#   # is any column of COEF has all valus NA the omit it 
#   if (any(is.na(COEF)))
#   {
#     ind <- rep(0,Q)
#     for (i in 1:Q)  ind[i] <- ifelse(all(is.na(COEF[,i])),-i, ind[1])
#     COEF <- COEF[,ind]
#   }
#   #################################################################
#   #                           OUTPUT                              #
#   #################################################################  
#   out <- list(boot = COEF, 
#               B = B,
#               trueB = newB,
#               par = parameters,
#               orig.coef = orig.coef,   
#               orig.call = obj$call) 
#   class(out) <- list("Bayesian.boot")
#   out 
# }
# #########################################################################################
# #########################################################################################
# #########################################################################################
# ######################################################################################### 
# # probably this will not word becaus it need random  generator
# NonParametricBootT <- function(obj, data=NULL, B=100)
# {
#   if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
#   data <- if ("data"%in%names(obj$call)) eval(obj$call$data)
#   else if (!is.null(data)) data
#   else stop("data are not defined")  
#   if (!is.data.frame(data)) stop("data is not a dataframe")
#   N <- dim(data)[1]      # length of data
#   thecall <- as.call(obj$call) # the call
#   thecall$control <- gamlss.control(trace=FALSE) # supressing tracing
#   orig.coef <- unlist(coefAll1(obj, deviance=TRUE))
#   Q <- length(orig.coef) 
#   res <- foreach(i = 1:B, .packages="gamlss", .errorhandling = "remove",
#                  .inorder = FALSE)%dopar%{           
#                    bootstrap_data <- data[sample(nrow(data),nrow(data),replace=T),]
#                    thecall$data <- bootstrap_data
#                    model_1 <- eval(thecall) # how to catch?
#                    p(sprintf("i=%g", i))
#                    list(coef=coefAll1(model_1, deviance = TRUE), par=predictAll(model_1,output="matrix"))                   }
#   newB <- length(res)
#   if (newB!=B) warning(cat(B-newB, "simulation failed only", newB, "were finished \n"))        
#   ## find which bootstrap failed
#   COEF <- matrix(0, ncol=Q, nrow=newB)
#   if ("mu" %in% obj$par)       MU <- matrix(0, ncol=newB, nrow=N)
#   if ("sigma" %in% obj$par) SIGMA <- matrix(0, ncol=newB, nrow=N)
#   if ("nu" %in% obj$par)       NU <- matrix(0, ncol=newB, nrow=N)
#   if ("tau" %in% obj$par)     TAU <- matrix(0, ncol=newB, nrow=N)
#   # maybe I should do this with foreach      
#   for (i in 1:newB)
#   {
#     COEF[i,] <- unlist(res[[i]][["coef"]])
#     if ("mu" %in% obj$par)        MU[,i] <- unlist(res[[i]][["par"]])[,"mu"]
#     if ("sigma" %in% obj$par)  SIGMA[,i] <- unlist(res[[i]][["par"]])[,"sigma"]
#     if ("nu" %in% obj$par)        NU[,i] <- unlist(res[[i]][["par"]])[,"nu"]
#     if ("tau" %in% obj$par)      TAU[,i] <- unlist(res[[i]][["par"]])[,"tau"]
#   }
#   parameters <- switch(length(obj$par),
#                        list(mu=MU),
#                        list(mu=MU, sigma=SIGMA),
#                        list(mu=MU, sigma=SIGMA, nu=NU),
#                        list(mu=MU, sigma=SIGMA, nu=NU, tau=TAU))
#   colnames(COEF) <-   names(unlist(res[[1]][["coef"]]))
#   # is any column of COEF has all valus NA the omit it 
#   if (any(is.na(COEF)))
#   {
#     ind <- rep(0,Q)
#     for (i in 1:Q)  ind[i] <- ifelse(all(is.na(COEF[,i])),-i, ind[1])
#     COEF <- COEF[,ind]
#   }
#   #################################################################
#   #                           OUTPUT                              #
#   #################################################################  
#   out <- list(boot = COEF, 
#               B = B,
#               trueB = newB,
#               par = parameters,
#               orig.coef = orig.coef,   
#               orig.call = obj$call) 
#   class(out) <- list("NonParametric.Boot")
#   out 
# }
# #########################################################################################
# 
# 
# 
# 
# #########################################################################################