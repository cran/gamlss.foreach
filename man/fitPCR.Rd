\name{fitPCR}
\alias{fitPCR}

\title{
Function to fit simple Principal Component Regression.
}
\description{
This function is a univariate (i.e. one response) version of a principal component regression. It is based on the function \code{svdpc.fit()}
 of package \pkg{pls} but it has been generalised to take prior weights. It gets a (single) response variable \code{y} (n x 1) and a matrix of  explanatory variables of dimensions n x p  and fits different  principal component regressions up to principal component M. Note that M can be less or equal to p (if \eqn{n > p}{n > p}) or less or equal to n if \eqn{n <p}{n < p}, that is, when there they are less observations than variables.
 
The function is used by the GAMLSS additive term function \code{pcr()} to fit a principal component regression model within \code{gamlss()}.
}
\usage{

fitPCR(x = NULL, y = NULL, weights = rep(1, n), 
       M = NULL, df = NULL, supervised = FALSE, 
       k = 2, r = 0.2, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix of explanatory variables}
  \item{y}{the response variable}
  \item{weights}{prior weights}
  \item{M}{if set  specifies the maximum number of components to be considered}
  \item{df}{if set specifies the number of components}
  \item{supervised}{whether supervised PCR should be used or not, \code{default=FALSE}}
  \item{k}{the penalty of GAIC}
  \item{r}{a correlation value (between zero and one) used smoothing parameter when \code{supervised=TRUE}}
  \item{plot}{Whether to plot the coefficient path}
}

\details{ 
More details here
}
\value{
It returns a object \code{PCR} which can be used with methods \code{print()}, 
\code{summary()}, \code{plot()}, \code{fitted()} and \code{coef()}. The object has elements:

  \item{coefficients}{The beta coefficients for 1 to M principal components}
  \item{scores}{the n x M dimensional matrix T o=f scores} 
  \item{loadings}{the p x M dimensional matrix P of loadings}  
  \item{gamma}{the first M principal component coefficients}                 \item{se.gamma}{the standard errors of the  M principal component coefficients}
  \item{center}{the location parameters used to scale the x's} 
  \item{scale}{the scale parameters used to scale the x's} 
  \item{fitted.values}{matrix of n x M dimensions}  
  \item{Xvar}{sum of squares  of the scores i.e. diag(T'T)} 
  \item{gaic}{The GAIC values} 
  \item{pc}{number  of PC i.e. which value of  GAIC has the minimum}
  \item{k}{which penalty for GAIC}
  \item{M}{the maximum of PC tried}
  \item{sigma}{The estimated sigma from the M fitted components}
  \item{residuals}{The n x M matrix of the residuals}
  \item{call}{the function call}
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
\author{Mikis Stasinopoulos, Robert Rigby and Fernanda De Bastiani.}


\seealso{
\code{\link{pc}}
}
\examples{
library(glmnet)
data(QuickStartExample)
attach(QuickStartExample)
hist(y, main="(a)")
if (is.null(rownames(x))) colnames(x) <- paste0("X", 
    seq(1:dim(x)[2]))
############################################################
# fitPCR
############################################################
# fitting
registerDoParallel(cores = 2)
MM<- fitPCR(x,y, k=log(100))
stopImplicitCluster()
points(MM$coef[,16]~rep(16,20))
names(MM)
MM
#----------------------------------------------------------
# plotting
plot(MM)
plot(MM, "gaic")
#----------------------------------------------------------
print(MM)
#----------------------------------------------------------
coef(MM)                        # the gammas
coef(MM, param="beta")          # the betas
coef(MM, param="beta", pc=1)  # at position 1
#----------------------------------------------------------
# plotting y and and fitted balues at different points
plot(y)
points(fitted(MM, pc=3), col=2)
points(fitted(MM, pc=20), col=3)
#----------------------------------------------------------
# variance covariance 
vcov(MM, type="se", pc=1) 
vcov(MM, type="se", pc=2)
vcov(MM, type="se", pc=20)
# library(corrplot)
# corrplot(vcov(MM, type="cor", pc=10))
# corrplot(vcov(MM, type="cor", pc=20))
#----------------------------------------------------------
summary(MM)
summary(MM, param="beta", pc=15)
summary(MM, param ="beta", pc=3) 
summary(MM, param ="beta") # at default
#----------------------------------------------------------
predict(MM, newdata= x[1:5,])
fitted(MM)[1:5]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{rgression}% use one of  RShowDoc("KEYWORDS")
