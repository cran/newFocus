\name{pick}
\alias{pick}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
True discoveries for non-focus level node
}
\description{
The number of true discoveries for the non-focus level GO terms is calculated given the focus subject.
}
\usage{
pick(focus_obj, setofinterest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{focus_obj}{
  The focus subject from function \link{newFocus}
}
  \item{setofinterest}{
  A gene set or GO term of interest
}
}

\value{
It will return an integer: the lower bound for the number of true discoveries in the set of interest
}

\author{
Ningning Xu
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