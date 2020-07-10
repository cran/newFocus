\name{ctbab}
\alias{ctbab}
\title{
Closed testing with branch and bound 
}
\description{
Closed testing with branch and bound algorithm specifically for globaltest
}
\usage{
ctbab(y, Cm, Tm, upnode, level, lownode, tmin, ctrue, lf, ls, alpha, count = 0, maxIt = 0)
}

\arguments{
  \item{y}{
  The response variable
}
  \item{Cm}{
  The matrix for calculating critical values of globaltest
}
  \item{Tm}{
  The matrix for calculating test statistics of globaltest
}
  \item{upnode}{
  The upper node that is used to bound critical values
}
  \item{level}{
  The level that the GO term of interest 
}
  \item{lownode}{
  The lower node that is used to bound critical values
}
  \item{tmin}{
  The minimum test statistic 
}
  \item{ctrue}{
  The true critical value corresponding to the minimum test statistic
}
  \item{lf}{
  The lambda vector corresponding to the upper node
}
  \item{ls}{
  The lambda vector corresponding to the lower node
}
  \item{alpha}{
  The significance level
}
  \item{count}{
  An integer stores the repetitions of the branch and bound, i.e. how many time branch and bound is implemented
}
  \item{maxIt}{
  The maximal number of repetitions prespecified by user
}
}

\value{
It will retrun the rejection indicator by closed testing with branch and bound algorithm.
}

\author{
Ningning Xu
}

\references{
Xu, N., & Goeman, J. (2020). Closed testing with Globaltest with applications on metabolomics data. arXiv preprint arXiv:2001.01541.
}

