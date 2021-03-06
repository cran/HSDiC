% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HSDiC_ADMM.R
\name{HSDiC_ADMM}
\alias{HSDiC_ADMM}
\title{Homogeneity Detection Incorporating Prior Constraint Information by ADMM}
\usage{
HSDiC_ADMM(X, Y, A.eq, A.ge, A.lbs, A.ubs, b.eq, b.ge, b.lbs, b.ubs,
  penalty = c("MCP", "SCAD", "adlasso", "lasso"), lambda2,
  admmScale1 = 1/nrow(X), admmScale2 = 1, admmAbsTol = 1e-04,
  admmRelTol = 1e-04, nADMM = 2000, admmVaryScale = FALSE)
}
\arguments{
\item{X}{n-by-p design matrix.}

\item{Y}{n-by-1 response matrix.}

\item{A.eq}{equality constraint matrix.}

\item{A.ge}{inequality constraint matrix.}

\item{A.lbs}{low-bounds matrix on variables, see examples.}

\item{A.ubs}{upper-bounds matrix on variables, see examples.}

\item{b.eq}{equality constraint vector.}

\item{b.ge}{inequality constraint vector.}

\item{b.lbs}{low-bounds on variables, see details.}

\item{b.ubs}{upper-bounds on variables, see details.}

\item{penalty}{The penalty to be applied to the model. Either "lasso" (the default), "SCAD",
or "MCP".}

\item{lambda2}{penalty tuning parameter for thresholding function.}

\item{admmScale1}{first ADMM scale parameter, 1/nrow(X) is default.}

\item{admmScale2}{second ADMM scale parameter, 1 is default.}

\item{admmAbsTol}{absolute tolerance for ADMM, 1e-04 is default.}

\item{admmRelTol}{relative tolerance for ADMM, 1e-04 is default.}

\item{nADMM}{maximum number of iterations for ADMM, 2000 is default.}

\item{admmVaryScale}{dynamically chance the ADMM scale parameter, FALSE is default}
}
\value{
\item{betahat}{solution vector.}

\item{stats.ADMM_inters}{number of iterations.}
}
\description{
simultaneous homogeneity detection and variable selection incorporating prior constraint by ADMM
algorithm. The problem turn to solving quadratic programming problems of the form
\emph{min(-d^T b + 1/2 b^T D b)} with the constraints \emph{A^T b >= b_0}. The penalty is the pairwise
fusion with p(p-1)/2 number of penalties.
}
\references{
{'Pairwise Fusion Approach Incorporating Prior Constraint Information' by Yaguang Li}
}
\seealso{
\code{\link{solve.QP}}
}
\examples{

## data generation
set.seed(111)
n=100
p=50
r <- 1 #0.5, 0.8, 1

beta <- r*c(sample(rep(1:2, each = 10)), rep(0,10), -sample(rep(1:2, each = 10)) )
X <- matrix(rnorm(n*p),nrow = n)
sigma = 1
Y <- X \%*\% beta + sigma * rnorm(n, 0, 1)


# equalities
A.eq <- rbind(rep(1,p))
b.eq <- c(0)

# inequalities
A.ge <- diag( c(rep(1,30), rep(-1,20)) )
b.ge <- rep(0,p)

# lower-bounds
A.lbs <- diag(1, p)
b.lbs <- rep(-2, p)

# upper-bounds on variables
A.ubs <- diag(-1, p)
b.ubs <- rep(-2, p)

ptm <- proc.time()
fit <- HSDiC_ADMM(X, Y, A.eq, A.ge, A.lbs, A.ubs, b.eq, b.ge, b.lbs, b.ubs,
                 penalty = "adlasso", lambda2 = 0.8, admmScale2 = 1)
proc.time() - ptm

## table(round(fit$beta,1))

plot(beta, type="p", pch = 20, cex = 1)
points(fit$beta, col = 3)


}
