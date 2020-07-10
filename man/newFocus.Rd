\name{newFocus}
\alias{newFocus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The new focus level procedure
}
\description{
The new focus level procedure for calculating true discoveries for focus level nodes
}
\usage{
newFocus(response, fsets, null, data, maxit = 0, alpha = 0.05, adj = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{
  The response variable
}
  \item{fsets}{
  A list of focus level sets
}
  \item{null}{
  The null hypothesis
}
  \item{data}{
  The data frame with response and all covariates included
}
  \item{maxit}{
  The maximal number of repetitions prespecified by user
}
  \item{alpha}{
  The significance level
}
  \item{adj}{
  The number of focus sets that are fully rejected by partial closed testing, which is used to adjust the number of focus sets, The dafault value is 0. 
}
}

\value{
The function will return a focus subject with the lower bound for each focus level node and the focus level nodes
}

\author{
Ningning Xu
}

\references{
Goeman, J. J., & Mansmann, U. (2008). Multiple testing on the directed acyclic graph of gene ontology. Bioinformatics, 24(4), 537-544.
}

\examples{
## example data set
n= 100
m = 5
X = matrix(0, n, m,byrow = TRUE )
for ( i in 1:n){
  set.seed(1234+i)
  X[i,] =  as.vector(arima.sim(model = list(order = c(1, 0, 0), ar = 0.2), n = m) )
}
y = rbinom(n,1,0.6)
X[which(y==1),1:3] = X[which(y==1),1:3] + 0.8
xs = paste("x",seq(1,m,1),sep="") 
colnames(X) = xs

mydata = as.data.frame(cbind(X,y))

## focus level sets
fl = list(c("x1", "x2"), c("x3", "x4"), "x5")
names(fl) = c("12", "34", "5")

## get td for focus level sets
focus_subject = newFocus(response = y, fsets = fl, data = mydata)

## get td for any set of interest given the focus subject
setofinterest = c("x1", "x2","x3", "x4")
pick(focus_subject, setofinterest)

}
