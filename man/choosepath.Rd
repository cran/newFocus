\name{choosepath}
\alias{choosepath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A set of focus set index
}
\description{
The function aims to find out the focus set index for which the true discoveries is the most and all other focus sets that are disjoint with it .
}
\usage{
choosepath(startingindex = 1, fsets, lowdv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{startingindex}{
  The index of focus set that has the first largest number of true discovereis 
}
  \item{fsets}{
  A list of focus level gene sets,or GO (Gene Ontology) terms  
}
  \item{lowdv}{
  A non-negative integer vector, which are the number of true discovereis, the length of the vector is the same as the list of focus level sets
}
}

\value{
The function will return an integer or a numeric vector.
}

\author{
Ningning Xu
}

