\name{BayesianBoot}
\alias{BayesianBoot}
\alias{NonParametricBoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non parametric and Bayesian Bootstrapping for  GAMLSS models
}
\description{
The function takes a GAMLSS fitted model and bootstrap it to create \code{B} bootstrap samples.
}
\usage{
NonParametricBoot(obj, data = NULL, B = 100, newdata = NULL)

BayesianBoot(obj, data = NULL, B = 100, newdata = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a \code{gamlss} fitted model}
  \item{data}{a data frame}
  \item{B}{the number of boostrap samples}
  \item{newdata}{new data for \code{predictAll()}}
}
\details{
The function  \code{NonParametric()} perform non-parametric bootstraping, Efron and Tibshirani (1993) while the function  \code{BayesianBoot()} perform Bayesian bootstrap
Rubin (1981)}
\value{
An \code{Bayesian.boot} object with elements
 \item{boot}{the bootstrap samples}
  \item{B}{the required number of boostraps}
  \item{trueB}{the actual number of boostraps}
  \item{par}{the distribution parameters}
  \item{orig.coef}{the fitted coeficients from the GAMLSS model}
  \item{orig.call}{the call from the GAMLSS model}
}
\references{
Efron, B. and Tibshirani, R, (1993), \emph{An introduction to the bootstrap},
  Chapman and Hall New York, Monographs on statistics and applied probability,
  vulume 57.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape, (with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC, \doi{10.1201/9780429298547}. An older version can be found in \url{https://www.gamlss.com/}.

Rubin, D. B. (1981) the bayesian bootstrap. \emph{The annals of statistics},
pp. 130-134.

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
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{ \code{\link[gamlss]{gamlss}}}

\examples{
m1 <-gamlss(y~x+qrt, data=aids, family=NBI())
registerDoParallel(cores = 2)
B1 <- BayesianBoot(m1)
summary(B1)
plot(B1)
B2 <- NonParametricBoot(m1)
stopImplicitCluster()
summary(B2)
plot(B2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
 \keyword{ regression}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
