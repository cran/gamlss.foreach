\name{fitted.PCR}
\alias{fitted.PCR}
\alias{coef.PCR}
\alias{plot.PCR}
\alias{summary.PCR}
\alias{predict.PCR}
\alias{vcov.PCR}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods for PCR objects}
\description{The functions below are methods for PCR objects}
\usage{
\method{fitted}{PCR}(object, pc = object$pc, ...)
\method{plot}{PCR}(x, type = c("path", "gaic"),
                   labels = TRUE, cex = 0.8, ...)
\method{coef}{PCR}(object, param = c("gamma", "beta"), 
              pc = object$pc, ...)
\method{predict}{PCR}(object, newdata = NULL, 
                pc = object$pc, ...)
\method{vcov}{PCR}(object, pc = object$pc, 
        type = c("vcov", "cor", "se", "all"),
                      ...)
}

\arguments{
  \item{object, x}{an PCR object}
  \item{pc}{the number of PC  (by default the one minimising the local GAIC)}
  \item{type}{for \code{plot()} whether to plot the path of coefficients or the path of GAIC and for \code{vcov()} whether variance covariance, correlation or se's}
  \item{param}{getting the gamma or the beta coefficients}
  \item{newdata}{new data for prediction}
  \item{labels}{whether to plot the labels of the variables}
  \item{cex}{the size of the text when plotting the labels of the variables}
   \item{\dots}{for extra arguments}
}

\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape, (with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).

}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby and Fernada De Bastiani  }



%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{  \code{\link{fitPCR}} }


\keyword{regression}% use one of  RShowDoc("KEYWORDS")

