\name{gamlss.foreach-package}
\alias{gamlss.foreach-package}
\alias{gamlss.foreach}
\docType{package}
\title{Computational Intensive Functions within GAMLSS} 
\description{This package is intended for functions needed parallel computations  provided by the package \pkg{foreach}. 

At the moment the following functions exist:

\code{centiles.boot()}, which is designed get bootstrap confidence intervals for centile curves

\code{fitRolling()}, rolling regression  which is common in time series analysis when one step ahead forecasts is required.

\code{fitPCR()}, for univariate  principal component regression. I


} 

\details{
The DESCRIPTION file:
\packageDESCRIPTION{gamlss.foreach}
\packageIndices{gamlss.foreach}
}
\author{
Mikis Stasinopoulos, \email{d.stasinopoulos@londonmet.ac.uk},and Bob Rigby 
            \email{r.rigby@londonmet.ac.uk}  
            
Maintainer: Mikis Stasinopoulos, \email{d.stasinopoulos@londonmet.ac.uk} 
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

(see also \url{https://www.gamlss.com/}).

}

\keyword{ package }
\seealso{
\code{\link[gamlss]{gamlss}},\code{\link[gamlss]{centiles}},\code{\link[gamlss]{centiles.pred}}
}
\examples{
\donttest{
library(gamlss.foreach)
# fixed degrees of freedom
cl <- makePSOCKcluster(2)
registerDoParallel(2)
data(db)
nage <- with(db, age^0.33)
ndb <- data.frame(db, nage)
m1 <- gamlss(head~cs(nage, 12), sigma.fo=~cs(nage,4), nu.fo=~nage, 
             tau.fo=~nage, family=BCT, data=ndb)
test1 <- centiles.boot(m1, xname="nage", xvalues=seq(0.01,20,0.2),B=10, power=0.33)
test1
plot(test1)
# degrees of freedom varying
m2 <- gamlss(head~pb(nage), sigma.fo=~pb(nage), nu.fo=~pb(nage), 
             tau.fo=~pb(nage), family=BCT, data=ndb)
test2 <- centiles.boot(m2, xname="nage", xvalues=seq(0.01,20,0.2),B=10, power=0.33)
stopImplicitCluster()
test2
plot(test2)
}
}
