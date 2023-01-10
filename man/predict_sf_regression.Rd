\name{predict_sf_regression}
\alias{predict_sf_regression}
\title{Prediction for a scalar-on-function regression model}
\description{This function is used to make prediction for a new set of functional and scalar (if any) predictors based upon a fitted scalar-on-function regression model in the output of \code{\link{rob.sf.reg}}.
}
\usage{
predict_sf_regression(object, Xnew, Xnew.scl = NULL)
}
\arguments{
  \item{object}{An output object obtained from \code{\link{rob.sf.reg}}.}
  \item{Xnew}{A list of matrices consisting of the new observations of functional predictors. The argument \code{Xnew} must have the same length and the same structure as the input \code{X} of \code{\link{rob.sf.reg}}.}
  \item{Xnew.scl}{A matrix consisting of the new observations of scalar predictors. The argument \code{Xnew.scl} must have the same length and the same structure as the input \code{X.scl} of \code{\link{rob.sf.reg}}.}
}
\value{An \eqn{n_{test} \times 1}-dimensional matrix of predicted values of the scalar response variable for the given set of new functional and scalar (if any) predictors \code{Xnew} and \code{Xnew.scl}, respectively. Here, \eqn{n_{test}}, the number of rows of the matrix of predicted values, equals to the number of rows of \code{Xnew} and and \code{Xnew.scl} (if any).
}
\author{Ufuk Beyaztas and Han Lin Shang}
\examples{
set.seed(2022)
sim.data <- generate.sf.data(n = 400, n.pred = 5, n.gp = 101, out.p = 0.1)
out.indx <- sim.data$out.indx
indx.test <- sample(c(1:400)[-out.indx], 120)
indx.train <- c(1:400)[-indx.test]
Y <- sim.data$Y
X <- sim.data$X
Y.train <- Y[indx.train,]
Y.test <- Y[indx.test,]
X.train <- X.test <- list()
for(i in 1:5){
  X.train[[i]] <- X[[i]][indx.train,]
  X.test[[i]] <- X[[i]][indx.test,]
}
gp <- rep(list(seq(0, 1, length.out = 101)), 5) # grid points of Xs

model.tau <- rob.sf.reg(Y.train, X.train, emodel = "robust", fmodel = "tau", gp = gp)
pred.tau <- predict_sf_regression(object = model.tau, Xnew = X.test)
round(mean((Y.test - pred.tau)^2), 4)        # 1.868 (tau method)
}
