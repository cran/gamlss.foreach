\name{which.Data.Corr}
\alias{which.Data.Corr}
\alias{which.yX.Corr}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Detecting Hight Pair-Wise Correlations in Data}
\description{There are two function here. 

The function \code{which.Data.Corr()} is taking as an argument a \code{data.frame} or a data matrix and it reports the pairs of variables which have higher correlation  than \code{r}.

The function \code{which.yX.Corr()} it takes as arguments a continuous response variable, \code{y}, and a set of continuous explanatory variables, \code{x}, (which may include first order interactions), 
 and it creates a \code{data.frame} containing all variables with a pair-wise correlation above \code{r}. If the set of the continuous explanatory variables contains first order interactions, then by default, (\code{hierarchical = TRUE}), the main effects of the first order interactions will be also included so hierarchy will be preserved.  
}
\usage{
which.Data.Corr(data, r = 0.9, digits=3)
which.yX.Corr(y, x, r = 0.5, plot = TRUE, 
              hierarchical = TRUE, print = TRUE, digits=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} or a matrix}
  \item{r}{a correlation values (acting as a lower limit)}
  \item{y}{the response variable (continuous)}
  \item{x}{the (continuous) explanatory variables}
  \item{plot}{whether to plot the results or not}
  \item{print}{whether to print the dim of the new matrix or not}
  \item{hierarchical}{This is designed for make sure that if first order interactions are included in the list the main effects will be also included}
  \item{digits}{the number of digits to print.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
The function \code{which.Data.Corr()} creates a matrix with three columns. The first two columns contain the names of the variables having pair-wise correlation higher than \code{r} and the third column show their correlation.

The function \code{which.yX.Corr()} creats a design matrix which containts variables which have 
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
. 
}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby and Fernada De Bastiani.
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pc}}
}
\examples{
data(oil, package="gamlss.data")
dim(oil)
# which variables are highly correlated?
CC<- which.Data.Corr(oil, r=0.999)
head(CC)
# 6 of them
# get the explanatory variables
form1 <- as.formula(paste("OILPRICE ~ ",
          paste(names(oil)[-1],collapse='+'))) 
# no interactions
X <- model.matrix(form1, data=oil)[,-1]
dim(X)
sX <- which.yX.Corr(oil$OILPRICE,x=X, r=0.4)
dim(sX)
# first order interactions
form2 <- as.formula(paste("OILPRICE ~ ",
        paste0(paste0("(",paste(names(oil)[-1], 
        collapse='+')), ")^2"))) 
form2
XX <- model.matrix(form2, data=oil)[,-1]
dim(XX)
which.yX.Corr(oil$OILPRICE,x=XX, r=0.4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% use one of  RShowDoc("KEYWORDS")

