\name{cholInvArray}
\alias{cholInvArray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Precision matrices to variance matrices for Winbugs output}
\description{
  Given an array containing simulations from the posterior of a precision matrix, each individual precision matrix is converted to variances, covariances, and correlations.
}
\usage{
cholInvArray(x, prefix = "T", chol=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An array of winbugs output, with precision matrix entries of the form "T[1,3]"}
  \item{prefix}{The name of the precision matrix in winbugs, the "T" in "T[1,2  ]"}
  \item{chol}{If TRUE, the cholesky decomposition is returned instead of the inverse}
}
\details{
  Inverts the matrices with the cholesky decomposition, but operating on all matrices simultaneously using array arithmetic.
}
\value{
An array with the third dimension's precision matrix entries changed to 
  \item{"sdT[i,i]"}{for the standard deviation of component i}
  \item{"covT[i,j]"}{for the covariance between i and j}
    \item{"corrT[i,j]"}{for the correlations between i and j}
}
 
\examples{
# create a random positive definite matrix by 
# generating a lower triangle
  N=4
  lmat = diag(runif(N, 1, 10))
  thetri = lower.tri(lmat)
  lmat[thetri] = rnorm(sum(thetri), 0, 2)
#  precmat = solve(lmat \%*\% t(lmat))
  precmat = solve(lmat \%*\% t(lmat))
 
# put this matrix into an array 
  precarray = array(c(precmat), dim=c(1,1,length(precmat)))
  dimnames(precarray) = list(NULL, NULL, 
    paste("T[", rep(1:N, N), ",", rep(1:N, rep(N,N)), "]",sep="") )

# invert it with cholInvArray and the solve function
  cholInvArray(precarray)[1,1,]
  # the off diagonals of solve(precmat) should be 
  # the covT elements of cholInvArray(precarray)
  solve(precmat)
  # the standard deviations in cholInvArray(precarray) should be the 
  # root of the diagonals of solve(precmat)
  sqrt(diag(solve(precmat)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
 

