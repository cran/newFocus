\name{discov}
\alias{discov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
True discoveries
}
\description{
True discoveries calculated by the partial closed testing
}
\usage{
discov(response, alternative, null, data, maxit = 0, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{
  The response variable
}
  \item{alternative}{
  The alternative hypothesis, which is a character vector, i.e. a set of genes
}
  \item{null}{
  The null hypothesis
}
  \item{data}{
  A data frame with response and all covariates included
}
  \item{maxit}{
The maximal number of repetitions prespecified by user
}
  \item{alpha}{
  The significance level
}
}

\value{
It will return a non-negative integer: the lower bound for the number of true discovereis of the alternative gene set.
}
\author{
Ningning Xu
}
 
