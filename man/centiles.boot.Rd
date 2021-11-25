\name{centiles.boot}
\alias{centiles.boot}
\alias{print.centiles.boot}
\alias{summary.centiles.boot}
\alias{plot.centiles.boot}


\title{Bootstrapping centiles curves estimated using GAMLSS}
\description{
This is a function designed for non-parametric bootstrapping centile curves (growth curves) when the fitted model is fitted using GAMLSS with a single explanatory variable (usually age). 
Non parametric  bootstrapping resample the data with replacement.  
The model is refitted for each bootstraps sample.  
Notes that if smoothing is used in the model, it is advisable (but not necessary) that the smoothing degree of freedom are fixed throughout.
}
\usage{
centiles.boot(obj, data = NULL, xname = NULL, xvalues = NULL, 
  power = NULL, cent = c(2.5, 50, 97.5), B = 100, calibration = FALSE,
  ...)
             
\method{print}{centiles.boot}(x, ...) 
\method{summary}{centiles.boot}(object, fun = "mean", ...)              
\method{plot}{centiles.boot}(x, quantiles = c(0.025, 0.975), 
       ylab = NULL, xlab = NULL, location = "median", original = FALSE, 
       scheme = c("shaded", "lines"), col.cent = "darkred", 
       col.se = "orange", col.shaded = "gray", lwd.center = 1.5, ...)  
             
}
\arguments{
  \item{obj}{a fitted gamlss object for the function \code{centiles.boot()}}
  \item{data}{a data frame containing the variables occurring in the formula. 
        If it is missing, then it will try to get the data frame from the GAMLSS object}
  \item{xname}{the name (as character) of the unique explanatory variable (it has to be the same as in the original fitted model)}
  \item{xvalues}{a vector containing the new x-variable  values for which bootstrap simulation predictions will be made}
  \item{power}{if power transformation is needed (but see example below)}
  \item{cent}{a vector of centile values for which the predicted centiles 
  have to be evaluated, by default is: 2.5, 50 and 97.5}
  \item{B}{the number of bootstraps }
  \item{calibration}{whether to calibrate the centiles, default is FALSE}
  \item{\dots}{for extra arguments, for the \code{centiles.pred()} function}
  \item{x}{an \code{a centiles.boot} object}
  \item{object}{an \code{a centiles.boot} object}
  \item{fun}{for the \code{summary()} function this is a summary statistics function. The "mean" is the  default}
  \item{quantiles}{specify  which quantiles (in the \code{plot()} function)  of the bootstrap distribution to plot} 
  \item{location}{which location parameter to plot, with default the mean} 
  \item{original}{logical if TRUE the original predicted centile values at the given \code{xvalues} are plotted (the default is FALSE)}
   \item{ylab}{y label for the plot }
  \item{xlab}{x label for the plot}
    \item{scheme}{which scheme of plotting to use  \code{"shaded"} or \code{"lines"}}
\item{col.cent}{the colour of the centile in the \code{"shaded"} scheme, with \code{"darkred"} as default}
\item{col.se}{the colour of the standard errors for the \code{"lines"} scheme with default \code{"orange"}}
\item{col.shaded}{the colour of the standard errors for the \code{"shaded"} scheme with default "gray"}
  \item{lwd.center}{the width of the central line}
}
\details{This function is designed for bootstrapping centiles curves. 
       It can be used to provide bootstrap means, standard deviations and quantiles, 
       so the variability of the centile curves 
        can be accessed (eg. by deriving confidence bands for centile curves).   
}
\value{
  The function returns an object \code{centile.boot} which has its own 
  \code{print()}, \code{summary()}, and \code{plot()} functions.  
  The object  \code{centile.boot} is a list with elements:
  \item{boot0}{Containing centile predictions from the fitted model to the original data using the \code{centile.pred()} function  on the new \code{xvalues}. This can be compared with the \code{mean} 
                of the object to assess the bias}
  \item{boot}{A list of length \code{trueB}, each containing a matrix of dimension 
               \code{length(xvalues)} by \code{(length(cent)+1)}}
  \item{B}{The number of bootstrap samples requested}
  \item{trueB}{The number of actual bootstrapping simulations performed. It is equal to 
   \code{B-number of failed simulations}}
  \item{xvalues}{The new x-variable  values for which the bootstrap simulation has taken place  }
  \item{cent}{The centile values requested}
  \item{original.call}{The call of the original gamlss fitted model}
  \item{yname}{The name of the response variable, used in the \code{plot()} function}
   \item{xname}{The name of the x-variable, used in the \code{plot()} function}
  \item{failed}{A vector containing values identifying which of the bootstrap simulations had 
               failed to converge and therefore have not included in the list \code{boot}
   }
}

\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape, (with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC, \doi{10.1201/9780429298547}. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \doi{10.18637/jss.v023.i07}.


Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  \doi{10.1201/b21973}

Stasinopoulos, M. D., Rigby, R. A., and De Bastiani F., (2018) GAMLSS: a distributional regression approach, \emph{Statistical Modelling}, Vol. \bold{18}, 
pp, 248-273, SAGE Publications Sage India: New Delhi, India. \doi{10.1177/1471082X18759144}

(see also \url{https://www.gamlss.com/}).

}

\author{Mikis Stasinopoulos, \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby 
            \email{r.rigby@londonmet.ac.uk} and Calliope Akantziliotou}
\note{See example below of how to use the function when a power transformation is used for the x-variable

Do not forget to use \code{registerDoParallel(cores = NUMBER)} or
\code{cl <- makeCluster(NUMBER)} and 
\code{registerDoParallel(cl)}
before calling the function \code{centiles.boot()}. Use   \code{closeAllConnections()} after the fits to close the connections. Where \code{NUMBER} depends on the machine used. 

}
\seealso{ \code{\link[gamlss]{gamlss}},\code{\link[gamlss]{centiles}},\code{\link[gamlss]{centiles.pred}} }
\examples{
# bring the data and fit the model
data(abdom) 
m1<-gamlss(y~poly(x,2), sigma.fo=~x, data=abdom, family=BCT)
# perform the bootstrap simulation 
# (only 10 bootstrap samples here)
registerDoParallel(cores = 2)
boC<-centiles.boot(m1,xname="x", xvalues=c(15,20,25,30,35,40,45), B=10)
stopImplicitCluster()
boC
# get summaries
summary(boC, fun="mean")
#summary(boC, "median")
#summary(boC, "quantile", 0.025)
plot(boC)
\donttest{
# with transformation in x within the formula
# unsuitable for large data set since it is slow
m2<-gamlss(y~cs(x^0.5),sigma.fo=~cs(x^0.5), data=abdom, family=BCT) 
boC<-centiles.boot(m2,xname="x", xvalues=c(15,20,25,30,35,40,45), B=10)
summary(boC)
plot(boC)
#
# now with x-variable previously transformed
# better for large data set as it is faster 
nx<-abdom$x^0.5
newd<-data.frame(abdom, nx=abdom$x^0.5)
m3<-gamlss(y~cs(nx),sigma.fo=~cs(nx), data=newd, family=BCT)
boC <- centiles.boot(m3, xname="nx", xvalues=c(15,20,25,30,35,40,45), data=newd, power=0.5, B=10)
summary(boC)
#plot(boC)
# the original variables can be added in
#points(y~x,data=abdom)
}}
\keyword{regression}
