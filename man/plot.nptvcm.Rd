\name{plot.nptvcm}
\alias{plot.nptvcm}
\title{Visualize the result for object nptvcm.}
\description{
  Plot the figures about BIC/CV vs bandwidth, time-varing regression coefficient estimate and the confidence band.
}
\usage{
## S3 method
  plot(obj)
}
\arguments{
  \item{obj}{a nptvcm object obtained from function \link{nptvcm}.}
}
\details{
  Nothing.
}
\value{
  No value.
}
\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  \link{nptvcm},\link{plot}
}
\examples{
  N <- 100; TT<-5;
  datalist <- gendata(N, TT, seed=1)
  Yit <- datalist$Yit
  Xit <- datalist$Xit
  res <- nptvcm(Xit, Yit, estimateSE = T, tuning.select='BIC')
  plot(res)
  ##
  N <- 100; TT<-10;
  datalist <- gendata(N, TT, seed=1)
  Yit <- datalist$Yit
  Xit <- datalist$Xit
  res <- nptvcm(Xit, Yit, estimateSE = T, tuning.select='CV')
  plot(res)
  ###
  ##
  N <- 100; TT<-10;
  datalist <- gendata(N, TT, seed=1, set=2)
  Yit <- datalist$Yit
  Xit <- datalist$Xit
  res <- nptvcm(Xit, Yit, estimateSE = T, tuning.select='CV')
  plot(res)
}
\keyword{Panel data generating}
