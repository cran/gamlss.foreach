\name{fitRolling}
\alias{fitRolling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to Fit Rolling Regression in gamlss}
\description{
Rolling regression is common in time series analysis when one step ahead forecasts is required. The function \code{fitRolling()} works as follows: A model is fitted first to the whole data set using \code{gamlss()}. Then the function \code{fitRolling()} can be used. The function uses a fixed size rolling widow  i.e.  365 days. The model is refitted repeatedly for the different windows over time (like a local regression in smoothing).  Each time one step ahead forecast of distribution parameters are saved together with the prediction global deviance. The result is presented as a matrix  with time as rows and parameters and the prediction deviance as columns.
}
\usage{
fitRolling(obj, data, window = 365, as.time = NULL)
}

\arguments{
  \item{obj}{a gamlss fitted model}
  \item{data}{the original data of the fitted model}
  \item{window}{the number of observation to include in the window (typically this will be a year)}
  \item{as.time}{if a column indicating time exist in the data set this can be specified here}
}
\details{
If the total observations are \code{N} and the window size \code{n} then we will need \code{N-n} different fits. The parallelization of the fits is achieved using the function \code{foreach()} from the package \pkg{foreach}.  
}
\value{Returns a matrix containing as columns the one ahead prediction parameters of the distribution as well as the  prediction global deviance.}
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
\author{
Mikis Stasinopoulos, \email{d.stasinopoulos@londonmet.ac.uk}
}
\note{
Do not forget to use \code{registerDoParallel(cores = NUMBER)} or
\code{cl <- makeCluster(NUMBER)} and 
\code{registerDoParallel(cl)}
before calling the function \code{fitRolling()} and \code{closeAllConnections()} after the fits. Where \code{NUMBER} depends on the machine used. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{ \code{\link[gamlss]{gamlss}}}
\examples{
\donttest{
# fitting the aids data 45 observations
m1 <- gamlss(formula = y ~ pb(x) + qrt, family = NBI, data = aids) 
# get rolling regression with a window of 30
# there are 45-40=15 fits to do
# declaring cores (not needed for small data like this)
registerDoParallel(cores = 2)
FF <- fitRolling(m1, data=aids, window=30)
FF
stopImplicitCluster()
# check the first prediction
m30_1 <-update(m1, data=aids[1:30,])
predictAll(m30_1, newdata=aids[31,],output="matrix")
FF[1,]
# plot all the data 
plot(y~x, data=aids, xlim=c(0,45), ylim=c(0, 700), col=gray(.8))
# the first 30 observations
points(y~x, data=aids[1:30,], xlim=c(0,45))
# One step ahead forecasts
lines(FF[,"mu"]~as.numeric(rownames(FF)), col="red")
lines(fitted(m1)~aids$x, col="blue")
}
}

\keyword{regression}
