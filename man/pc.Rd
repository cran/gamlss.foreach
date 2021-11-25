\name{pc}
\alias{pc}
\alias{pcr}
\alias{gamlss.pc}
\alias{gamlss.pcr}
\alias{getSVD}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to Fit Principal Component Regression in GAMLSS}
\description{
The functions \code{pcr()} and \code{pc()} can be use to fit  principal component regression (PCR) within a GAMLSS model. They  can be used as an extra additive term (together with other additive terms for example \code{pb()}) but the idea is mainly to be used on their own as a way of reducing the dimensionality of the the (scaled) x-variables. The functions can be used even when the number of the explanatory variables say \code{p} is greater than the number of observations \code{n}. 

The two functions differ on the way PCR is implemented within the GAMLSS algorithm see for example Stasinopoulos et.al (2021).
In the function \code{pc()} the singular value decomposition of the scaled x's is done in the beginning and different re-weighted linear models are fitted on the PC scores see  algorithm  1 in Stasinopoulos et al. (2021). In the function \code{pcr()} at each iteration a new weighted PCR is performed using the function \code{fitPCR()} see algorithm  2 in Stasinopoulos et al. (2021). 


The functions \code{gamlss.pcr()} and \code{gamlss.pc()} are supporting functions. The  are not intended to be called directly by users. The function \code{gamlss.pc()} is using the linear model function \code{lm()} to fit the first principal components while the function code{gamlss.pcr()} uses \code{fitPCR()}.

The function \code{getSVD()} creates a singular value decomposition of a design matrix \code{X} using the \R{} function \code{La.svd()}. 
}
\usage{
pc(x.vars = NULL, x = NULL, x.svd = NULL, df = NULL, 
   center = TRUE, scale = TRUE, tol = NULL, 
   max.number = min(p, n), k = log(n), 
   method = c( "t-values","GAIC","k-fold"))

pcr(x.vars = NULL, x = NULL, df = NULL, 
    M = min(p, n), k = log(n), 
    r = 0.2, method = c("GAIC", "t-values", "SPCR"))

gamlss.pc(x, y, w, xeval = NULL, ...)

gamlss.pcr(x, y, w, xeval = NULL, ...)

getSVD(x = NULL, nu = min(n, p), nv = min(n, p))

}

\arguments{
  \item{x.vars}{A character vector showing the names of the x-variables. The variables should exist in the original \code{data} argument declared in the \code{gamlss()} function}
  \item{x}{For the function \code{pc()} and \code{getSVD()} \code{x} is a design matrix of dimensions \code{n x p} contain all the explanatory variables terms. 
  
For the function \code{gamlss.pc()}, \code{x} is a vector of zeros which curries all in information needed for the principal components fits in its attributes
}
  \item{x.svd}{A list created by the function \code{getSVD()}. This will speed up the time of fitting, (especial for large data sets), since all  the principal components calculation are done in advance. Also if all the parameters of the distribution are modelled by principal components the calculation needed to be done only once.}
  \item{df}{(if is not \code{NULL}) the number of principal components to be fitted. If it is \code{NULL} the number of principal components is automatically calculated  using a GAIC criterion.}
  \item{center}{whether to center the explanatory variables with default \code{TRUE}}
  \item{scale}{whether to scale the explanatory variables with default \code{TRUE}}
  \item{r}{the cut point for correlation coefficient to be use \code{SPCR}} 
  \item{tol}{CHECK THIS?????}
  \item{max.number, M}{The maximum number of principal component to be used in the fit.}
   \item{method}{method used for choosing the number of components}
  \item{k}{the penalty for GAIC}
  \item{y}{the iterative response variable}
  \item{w}{the iterative weights }
  \item{xeval}{used in prediction}
  \item{\dots}{for extra arguments}  
   \item{nu}{the number of left singular vectors to be computed. This must between 0 and n = nrow(x).}
  \item{nv}{	
the number of right singular vectors to be computed. This must be between 0 and p = ncol(x).}
}
\details{
There are three different ways of declaring the list of x-variables (two for the function \code{pcr()}):

\emph{x.vars}: this should be a character vector having the names of the explanatory variables. The names  should be  contained in the names  of variables of the \code{data} argument of the function \code{gamlss()}, see example below.

\emph{x}: This should be a design matrix (preferable unscaled since this could create problems when try to predict), see examples.

\emph{x.svd}: This should be a list created by the function \code{getSVD()} which is used as an argument a design matrix, see examples.
}

\value{
For the function \code{pc()} returns an object \code{pc} with elements \code{"coef"}, \code{"beta"}, \code{"pc"},   \code{"edf"},  \code{"AIC"}. The object \code{pc} has methods \code{plot()},  \code{coef()} and \code{print()}. 

For the function \code{pcr()} returns an object \code{PCR} see for the help for function \code{fitPCR}.  
}

\references{
Bjorn-Helge Mevik, Ron Wehrens and Kristian Hovde Liland
  (2019). pls: Partial Least Squares and Principal Component
  Regression. R package version 2.7-2.
  https://CRAN.R-project.org/package=pls

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

Stasinopoulos, M. D., Rigby, R. A., Georgikopoulos N., and De Bastiani F., (2021) Principal  component regression in GAMLSS  applied to Greek-German government bond yield spreads, \emph{Statistical Modelling} 
\doi{10.1177/1471082X211022980}.  


(see also \url{https://www.gamlss.com/}).  

}

\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby 
}
\note{ Do not forget to use \code{registerDoParallel(cores = NUMBER)} or
\code{cl <- makeCluster(NUMBER)} and 
\code{registerDoParallel(cl)}
before calling the function \code{pc()} without specifying the degrees of freedom.  Use \code{closeAllConnections()} after the fits to close the connections. The \code{NUMBER} depends on the machine used. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{centiles.boot}}, \code{\link{fitRolling}}
}
\examples{
# the pc() function
# fitting the same model using different arguments
# using x.vars
p1 <- gamlss(y~pc(x.vars=c("x1","x2","x3","x4","x5","x6")), data=usair)
registerDoParallel(cores = 2)
t1 <- gamlss(y~pcr(x.vars=c("x1","x2","x3","x4","x5","x6")), data=usair)
# using x
X <- model.matrix(~x1+x2+x3+x4+x5+x6, data=usair)[,-1]
p2 <- gamlss(y~pc(x=X), data=usair)
t2 <- gamlss(y~pcr(x=X), data=usair)
# using x.svd
svdX <- getSVD(X)
p3 <- gamlss(y~pc(x.svd=svdX), data=usair)
# selecting the componets 
p3 <- gamlss(y~pc(x.svd=svdX, df=3), data=usair)
stopImplicitCluster()
plot(getSmo(t2))
plot(getSmo(t2), "gaic")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
