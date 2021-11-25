#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# bootstap centiles
# Friday, December 8, 2006 at 17:46
# new version based on forreach on 2019 Aug 10 MS
# MS PA + BR 
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
###########################################################################
#                   'centiles.boot' arguments                             #
###########################################################################
# obj    : a fitted gamlss object from fitting a gamlss continuous distribution   
# data   : a data frame containing the variables occurring in the formula. 
#          If it is missing, then it will try to get the data frame from the GAMLSS object
#xname   : the name (as character) of the unique explanatory variable (it has to be the same as in the original fitted model) 
#xvalues : new x-variable  values in which the boootstap simulation will take place  
#cent    : a vector with elements the % centile values for which the centile curves have to be evaluated, by default is: 2.5%, 50% and 97.5%
# B     : the number of bootstraps
#...     : for extra arguments for the centiles.pred() function
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
centiles.boot <- function(
                       obj,
                      data = NULL, 
                     xname = NULL,
                   xvalues = NULL,
                     power = NULL,
                      cent = c(2.5, 50, 97.5),   
                         B = 100, 
               calibration = FALSE,
                      ...)
{
#----------------------------------------------------------------------------
##the code starts here for the centiles.boot() function
  if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
## checking the xvalues
  if (is.null(xvalues)) stop(paste("The 'xvalues' argument is not specified", "\n", ""))
##  this needs checking in case it does not cover all possibilities 
    data <- if ("data"%in%names(obj$call)) eval(obj$call$data)
            else if (!is.null(data)) data
                 else stop("data are not defined")          
    if (!is.data.frame(data)) stop("data is not a dataframe")
## checking the xname
    if (is.null(xname)) stop(paste("The xname argument is not specified", "\n", ""))
    if (!is.character(xname)) stop(paste("The xname argument is not a character", "\n", ""))
    if (!xname%in%names(data))  stop("the argument 'xname' is not included in the specified dataset.") 
############################################################################
            dDFr <- dim(data)[1]      # length of data
         thecall <- as.call(obj$call) # the call
 thecall$control <- gamlss.control(trace=FALSE) # supressing tracing
        centile  <- cent
if (calibration){
               z <- quantile(resid(obj), probs = cent/100)
               p <- pNO(z, mu = 0, sigma = 1)
         centile <- 100 * p
             }
        original <- centiles.pred(obj, xname=xname, cent=centile, xvalues=xvalues, power=power, data=data, ...)
##  foreach
  res <- foreach(i = 1:B, .packages=c("gamlss"), .errorhandling ="remove") %dopar% {           
     bootstrap_data <- data[sample(nrow(data),nrow(data),replace=T),]
       thecall$data <- bootstrap_data
                 m1 <- eval(thecall) # how to catch?
  if (calibration){
                 z <- quantile(resid(m1), probs = cent/100)
                 p <- pNO(z, mu = 0, sigma = 1)
           centile <- 100 * p
             }            
  centiles.pred(m1, xname=xname, xvalues=xvalues, cent=centile, power=power, ...  ) 
      }
## find which bootstrap failed
  trueB <- length(res)
  if (trueB<B) warning(B-trueB, " bootstraps failed", "\n" )
#  jj <- 0
#for(i in 1:length(res))  if(is.null(res[[i]])) jj=c(jj,i)
#          failed <- if (length(jj)==1) jj else jj[-1]
#         nfailed <- if (failed==0) 0 else length(failed) 
#           trueB <- B-nfailed
#############################################################################
#                                       OUTPUT                              #
#############################################################################  
  out <- list(boot0 = original, 
               boot = res, 
                  B = B, 
              trueB = trueB, 
            xvalues = xvalues,  
               cent = cent, 
      original.call = obj$call,
              yname = deparse(obj$call$formula[[2]]),
              xname = xname,
             failed = B-trueB) 
         class(out) <- list("centiles.boot")
  out 
}
#----------------------------------------------------------------------------
#   end of centiles.boot
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------      
#  centiles.boot.print
print.centiles.boot <- function(x,...)
  {
  cat("\n Results from centiles.boot() taken from:", x$trueB, " simulations \n")
  cat("\n The x-variable values are evaluated at: \n", head(x$xvalues,5),"...", tail(x$xvalues,5),"\n")
  cat("\n The centile values are evaluated at: \n", x$cent, "\n")
  cat("\n Simulation was done from the original call: \n")
  ccc<-x$original.call
  print(ccc)
  }  
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#############################################################################
#                            summary.centiles.boot                          #
#############################################################################
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
summary.centiles.boot <- function(object, fun="mean", ...)
{
         ll <- length(object$xvalues)*(length(object$cent)+1)
          B <- object$trueB
   len.cent <- object$cent
#      unres <- unlist(object$boot) 
      unres <- if(object$failed==0) unlist(object$boot) else unlist(object$boot[-object$failed])
      fact1 <- gl(ll,1,length=ll*B)
      funbb <- matrix(tapply(unres,fact1,fun,na.rm=TRUE,... ), ncol=length(object$cent)+1)
#      funbb <- matrix(tapply(unres,fact1,fun, ... ), ncol=length(object$cent)+1)

  funbb[,1] <- object$boot0[,1]
#Defines the column and rownames of the funbb matrix, by using the names of the 1st list of the obj$boot
   colnames(funbb) = colnames(object$boot[[1]]) 
   rownames(funbb) = rownames(object$boot[[1]]) 
 funbb
}
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
 plot.centiles.boot <- function(x, 
                          quantiles = c(0.025, 0.975),
                               ylab = NULL,
                               xlab = NULL,
                           location = "median", 
                           original = FALSE, 
                             scheme = c( "shaded", "lines"),
                           col.cent = "darkred",
                             col.se = "orange", 
                         col.shaded = "gray",
                         lwd.center = 1.5,   
                                ...)  
{
#----------------------------------------------------------------------------
# local function   
se.shaded <- function(x, y1, y2) 
   {
     xx <- c(x,rev(x))
     yy <- c( y1, rev(y2))
polygon(xx, yy, col = col.shaded, border = col.shaded)
}
#----------------------------------------------------------------------------
   scheme <- match.arg(scheme) 
   if (is.null(xlab)) 
      xlab <- x$xname
   if (is.null(ylab)) 
      ylab <- x$yname
# location   
 locationE <- summary(x, fun=location)
#quantiles
       c1 <- summary(x, fun="quantile", quantiles[[1]])
       c2 <- summary(x, fun="quantile", quantiles[[2]])
      miny <- min(c1[,-1])
      maxy <- max(c2[,-1])
# plotting      
plot(locationE[,1], locationE[,2], type="n" , ylim=c(miny,maxy), ylab=ylab, xlab=xlab, ... ) 
if (identical(scheme, "shaded"))
{
  for (j in 2:(length(x$cent)+1))
  {
    se.shaded(x=c1[,1],  y1=c1[,j], y2= c2[,j] )
    lines(locationE[,1], locationE[,j], col=col.cent, lwd = lwd.center)
    if (original) lines(x$boot0[,1], x$boot0[,j], col="black")
  }
}
else
{
  for (j in 2:(length(x$cent)+1))
  {  
    if (original) lines(x$boot0[,1], x$boot0[,j], col="black")
    lines(locationE[,1], locationE[,j], col=col.cent, lwd=lwd.center)
    lines(c1[,1], c1[,j], col=col.se)
    lines(c2[,1], c2[,j], col=col.se)
  }
}
invisible(list(lower=c1, upper=c2, location=locationE))
 }
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------


 