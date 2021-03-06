% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pace.R
\name{pace}
\alias{pace}
\title{Multiple robust estimators of average causal effects within principal strata (pace)}
\usage{
pace(X, Z, S, Y, family.Y, nboot = 0)
}
\arguments{
\item{X}{is a matrix of pre-treatment covariates without intercept  (n x p).}

\item{Z}{is a vector of treatment  (n x 1).}

\item{S}{is a vector of binary intermediate outcome  (n x 1).}

\item{Y}{is a vector of outcome  ( n x 1).}

\item{family.Y}{specifies the family for the outcome model.
\code{"gaussian"}: a linear regression model for the continuous outcome.

\code{"binomial"}: a logistic regression model for the binary outcome.}

\item{nboot}{is the number of bootstrap samples.}
}
\value{
\itemize{

\item \code{tau10w}: a principal score weighting estimator of tau10
\item \code{tau10sw}: a stabilized weighting estimator of tau10
\item \code{tau10reg}: a regression estimator of tau10 using outcome mean regression and inverse probability of treatment weighting
\item \code{tau10reg2}: a regression estimator of tau10 using outcome mean and principal score regression
\item \code{tau10aw}: a triply robust estimator of tau10
\item bootstrap variance estimator
\itemize{
\item \code{ve.tau10w}
\item \code{ve.tau10sw}
\item \code{ve.tau10reg}
\item \code{ve.tau10reg2}
\item \code{ve.tau10aw}
}

\item \code{tau00w}: a principal score weighting estimator of tau00
\item \code{tau00sw}: a stabilized weighting estimator of tau00
\item \code{tau00reg}: a regression estimator of tau00 using outcome mean regression and inverse probability of treatment weighting
\item \code{tau00reg2}: a regression estimator of tau00 using outcome mean and principal score regression
\item \code{tau00aw}: a triply robust estimator of tau00
\item bootstrap variance estimator
\itemize{
\item \code{ve.tau00w}
\item \code{ve.tau00sw}
\item \code{ve.tau00reg}
\item \code{ve.tau00reg2}
\item \code{ve.tau00aw}
}

\item \code{tau11w}: a principal score weighting estimator of tau11
\item \code{tau11sw}: a stabilized weighting estimator of tau11
\item \code{tau11reg}: a regression estimator of tau11 using outcome mean regression and inverse probability of treatment weighting
\item \code{tau11reg2}: a regression estimator of tau11 using outcome mean and principal score regression
\item \code{tau11aw}: a triply robust estimator of tau11
\item bootstrap variance estimator
\itemize{
\item \code{ve.tau11w}
\item \code{ve.tau11sw}
\item \code{ve.tau11reg}
\item \code{ve.tau11reg2}
\item \code{ve.tau11aw}
}
}
}
\description{
Implements various estimators of PACE from randomized experiments and observational studies
}
\details{
Details will be provided in the reference paper.
}
\examples{

library(stats)

set.seed(1)
#X consists of 5 covariates
n <- 500
X <- rnorm(n,0.25,1)
for(k in 2:4){
 X <- cbind(X,rnorm(n,0.25,1))
}
X <- cbind(X,rbinom(n,1,0.5))

# treatment assignment
Xtilde1 <- (X-0.25)/1
theta <- 0
eta <- c(0,1,1,1,1,theta)/2.5
px0 <- exp(cbind(1,Xtilde1)\%*\%eta)
pix  <- px0/(1+px0)
Z   <- rbinom(n,1,pix)

# S-model
eta1 <- c( 2,-1,1,-1,1,theta)/2.5
eta0 <- c(-2,1,-1,1,-1,theta)/2.5
p1xtemp <- exp(cbind(1,Xtilde1)\%*\%eta1)
p0xtemp <- exp(cbind(1,Xtilde1)\%*\%eta0)
p1x <- p1xtemp/(1+p1xtemp)
p0x <- p0xtemp/(1+p0xtemp)
s1<-rbinom(n,1,p1x)
s0<-rbinom(n,1,p0x)
S<-s1*Z+s0*(1-Z)

# Y-model (continuous outcome)
Y1<- Xtilde1\%*\%rep(1,5)*(1+S+Z) +rnorm(n)
# Y-model (binary outcome)
leYtemp <- Xtilde1\%*\%rep(1,5)*(1+S+Z)/4
eY<-exp(leYtemp)/(1+exp(leYtemp) )
Y2<-rbinom(n,1,eY)

out1<-pace(X,Z,S,Y1,family.Y="gaussian",nboot=50)
out1$tau10aw
out1$ve.tau10aw

out2<-pace(X,Z,S,Y2,family.Y="binomial",nboot=50)
out2$tau10aw
out2$ve.tau10aw
}
\references{
A reference paper will come up soon.
}
