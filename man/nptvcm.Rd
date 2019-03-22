\name{nptvcm}
\alias{nptvcm}
\title{NonParametric Time-Varying Coefficient Mmodel for Panel Data}
\description{
  Nonparametric time-varying coefﬁcients model based penalty for longitudinal data
    with pre-speciﬁed ﬁnite time points, also known as panel data.
}
\usage{
  nptvcm (Xit, Yit, lambda=NULL, tuning.select='CV', nfolds=10, estimateSE=F, penalty='SCAD',nlambda=50, lambda.max=5)
}
\arguments{
  \item{Xit}{a 3-dim n x T x p array, the observational value of covariates on panel data, the 1st index is sample index, the 2nd index is time index, the 3rd variable.}
  \item{Yit}{a n x T matrix, the observational value of response variable on panel data.}
  \item{lambda}{a positive value, user-specified penalty parameter, (optional) default as NULL.}
  \item{tuning.select}{a string, the method to select the penalty parameter, "CV" and "BIC" are supported in current version, (optional) default as "CV".}
  \item{nfolds}{a positive integer, specify the folds of cross validation when tuning.select='CV', (optional) default as 10.}
  \item{estimateSE}{a logical value, specify whether estimate the SE of estimated curve.}
  \item{penalty}{a string, specify the penalty function, supporting "SCAD", "MCP" and "laso", (optional) default as "SCAD".}
  \item{nlambda}{a positive integer, specify the grid of candidates on lambda, (optional) default as 50.}
  \item{lambda.max}{a positive value, specify the maximum value of lambda grids. (optional) default as 5.}
}
\details{
  For specific model and methodology, see reference.
}
\value{
  return a correlation matrix with a type of specified structure.
}
\references{
  Huazhen Lin, Hyokyoung G. Hong, Baoying Yang, Wei Liu, Yong Zhang, Gang-Zhi
Fan, Yi Li (2019). Nonparametric time-varying coeﬃcient models for panel data: Study of collection rate of public pension contributions. Statistics in Biosciences.
}
\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  plot.nptvcm
}
\examples{
  N=100;TT=5;
  datalist <- gendata(N, TT, 1)
  Yit <- datalist$Yit
  Xit <- datalist$Xit
  res <- nptvcm(Xit, Yit, estimateSE = T, tuning.select='CV')
}
\keyword{penalized regression on panel data}
