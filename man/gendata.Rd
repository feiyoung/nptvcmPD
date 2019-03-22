\name{gendata}
\alias{gendata}
\title{Generate the panel data for two simulation settings}
\description{
  Generate the panel data with 2 covariates for two simulation settings
}
\usage{
  gendata(N, TT, seed, set = 1, r= 0.5)
}
\arguments{
  \item{N}{a positive integer, the sample size of panel data.}
  \item{TT}{a positive integer, the length of observational time.}
  \item{seed}{an integer, specify the random seed for data reproducibility.}
  \item{set}{an integer with 1 or 2, specify which setting is to be used, (optional ) default as 1.}
  \item{r}{a value between 0 and 1, specify the autoregression coefficient of X_{jt}, (optional) default as 0.5.}
}
\details{
  The specific simulation settings are referred in the reference.
}
\value{
  return a list including the following components,
  \item{Yit }{a N x TT matrix, the observational value of response variable on panel data.}
  \item{Xit }{a 3-dim N x TT x 2 array, the observational value of covariates on panel data, the 1st index is sample index, the 2nd index is time index, the 3rd variable index.}
  \item{gt0}{an integer, specify the random seed for data reproducibility.}
  \item{gamma1t0}{an integer with 1 or 2, specify which setting is to be used, (optional ) default as 1.}
  \item{gamma2t0}{a value between 0 and 1, specify the autoregression coefficient of X_{jt}, (optional) default as 0.5.}
  \item{ft0}{}
  \item{beta1t0}{}
  \item{beta2t0}
}
\references{
  Huazhen Lin, Hyokyoung G. Hong, Baoying Yang, Wei Liu, Yong Zhang, Gang-Zhi
Fan, Yi Li (2019). Nonparametric time-varying coeﬃcient models for panel data: Study of collection rate of public pension contributions. Statistics in Biosciences.
}
\author{
  Wei Liu
}
\note{
  Nothing
}
\seealso{
  Nothing
}
\examples{
  N <- 100; TT<-5;
  datalist <- gendata(N, TT, seed=1)
  ##
  N <- 100; TT<-10;
  datalist <- gendata(N, TT, seed=1)
  Yit <- datalist$Yit
  Xit <- datalist$Xit
  ###
  ##
  N <- 100; TT<-10;
  datalist <- gendata(N, TT, seed=1, set=2)
  Yit <- datalist$Yit
  Xit <- datalist$Xit
}
\keyword{Panel data generating}
