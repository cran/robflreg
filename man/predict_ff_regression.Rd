\name{predict_ff_regression}
\alias{predict_ff_regression}
\title{Prediction for a function-on-function regression model}
\description{This function is used to make prediction for a new set of functional predictors based upon a fitted function-on-function regression model in the output of \code{\link{rob.ff.reg}}.
}
\usage{
predict_ff_regression(object, Xnew)
}
\arguments{
  \item{object}{An output object obtained from \code{\link{rob.ff.reg}}.}
  \item{Xnew}{A list of matrices consisting of the new observations of functional predictors. The argument \code{Xnew} must have the same length and the same structure as the input \code{X} of \code{\link{rob.ff.reg}}.}
}
\value{
An \eqn{n_{test} \times p}-dimensional matrix of predicted functions of the response variable for the given set of new functional predictors \code{Xnew}. Here, \eqn{n_{test}}, the number of rows of the matrix of predicted values, equals to the number of rows of \code{Xnew}, and \eqn{p} equals to the number of columns of \code{Y}, the input in the \code{\link{rob.ff.reg}}.
}
\author{
Ufuk Beyaztas and Han Lin Shang
}
\examples{
set.seed(2022)
sim.data <- generate.ff.data(n.pred = 5, n.curve = 200, n.gp = 101, out.p = 0.1)
out.indx <- sim.data$out.indx
Y <- sim.data$Y
X <- sim.data$X
indx.test <- sample(c(1:200)[-out.indx], 60)
indx.train <- c(1:200)[-indx.test]
Y.train <- Y[indx.train,]
Y.test <- Y[indx.test,]
X.train <- X.test <- list()
for(i in 1:5){
  X.train[[i]] <- X[[i]][indx.train,]
  X.test[[i]] <- X[[i]][indx.test,]
}
gpY = seq(0, 1, length.out = 101) # grid points of Y
gpX <- rep(list(seq(0, 1, length.out = 101)), 5) # grid points of Xs

model.MM <- rob.ff.reg(Y = Y.train, X = X.train, model = "full", emodel = "robust",
                       fmodel = "MM", gpY = gpY, gpX = gpX)
pred.MM <- predict_ff_regression(object = model.MM, Xnew = X.test)
round(mean((Y.test - pred.MM)^2), 4)        # 0.5925 (MM method)
}

