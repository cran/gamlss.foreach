#-----------------------------------------------------------#start
# to do 
#  i) what will happent if the model failed
# ii) What if we need more than one forecast?
#-----------------------------------------------------------#start
# this  is a rolling regression for time series data
fitRolling <- function(obj, data, window=365, as.time=NULL) {
  if (!is(obj,"gamlss")) stop("the obj should be a gamlss object")
  if (!grepl("$", deparse(substitute(as.time)), fixed=TRUE))
  {
    as.time <- eval(substitute(as.time), envir=as.environment(data))
  }
     pdf <- obj$family[1]
  binom  <- pdf%in%gamlss::.gamlss.bi.list # whether binomial
      fn <- eval(parse(text=pdf))()$G.dev.inc # get pdf
    lpar <- eval(parse(text=pdf))()$nopar
  if (binom) {bd <- obj$bd ; Y <- obj$y}
       N <- nrow(data)
    maxK <- N-window
       i <- 1 :  maxK
    yhat <- foreach(i = 1 :  maxK, .packages="gamlss", .combine = rbind) %dopar% {
     um1 <- try(update(obj, data = data[i:(window+(i-1)), ], 
                    mu.start = obj$mu.fv[i:(window+(i-1))],
                 sigma.start = obj$sigma.fv[i:(window+(i-1))],
                    nu.start = obj$nu.fv[i:(window+(i-1))],
                   tau.start = obj$tau.fv[i:(window+(i-1))],
                   trace=FALSE))
     if (any(class(um1)%in%"try-error")) stop("the algorithm failed in iteration", i, "\n")
     suppressMessages(pp<- predictAll(um1, data = data[i:(window+(i-1)), ],  newdata = data[window+i,], output="matrix"))
   # if) 
    if (binom)
    {
      DevIncr <- switch(lpar, 
                        fn( Y[i], mu = pp[,"mu"], bd=bd[i]),   # 1
                        fn( Y[i], mu = pp[,"mu"],              # 2
                               sigma = pp[,"sigma"], bd=bd[i]),                        
                        fn( Y[i], mu = pp[,"mu"],              # 3
                               sigma = pp[,"sigma"],
                                  nu = pp[,"nu"], bd=bd[i]),
                        fn( Y[i], mu = pp[,"mu"],              # 4
                               sigma = pp[,"sigma"],  
                                  nu = pp[,"nu"],
                                 tau = pp[,"tau"],bd=bd[i]))
    } else
    {
      DevIncr <- switch(lpar, 
                        fn( obj$y[i], mu =  pp[,"mu"]),       # 1
                        fn( obj$y[i], mu =  pp[,"mu"],        # 2
                            sigma =  pp[,"sigma"]),  # or fitted(obj, "sigma")                    
                        fn( obj$y[i], mu =  pp[,"mu"],        # 3
                            sigma =  pp[,"sigma"],
                            nu = pp[,"nu"]),
                        fn( obj$y[i], mu =  pp[,"mu"],        # 4
                                   sigma =  pp[,"sigma"],  
                                     nu = pp[,"nu"],
                                    tau = pp[,"tau"]))
    }  
    yhat <- cbind(pp, GDincr=DevIncr) 
    }
          cn <- colnames(yhat)
           M <- matrix(unlist(yhat), ncol=length(cn))
colnames(M)  <- cn
rownames(M)  <- if (!is.null(as.time)) as.time[(window+1):N]  else rownames(data[(window+1):N,])
  return(M)
}
