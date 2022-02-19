
fit.ordinal <- function(x, y, iter=200, lambda = 0, step = 1){
    p <- ncol(x)
    n <- nrow(x)
    k <- ncol(y) - 1

    beta <- matrix(0,p,k)
    coef0 <- rep(0,k+1)
    weight <- rep(1,n)
    step <- step
    coef <- fit_ordinal(x,y,weight,beta,coef0,lambda,iter,step)

    return (coef)
}

one.hot <- function(y){
    y <- model.matrix(~ factor(as.numeric(y) - 1) + 0)
    colnames(y) <- NULL
    y
}

#' @export
loss <- function(x,y,beta,intercept){
    y <- one.hot(y)
    p <- ncol(x)
    n <- nrow(x)
    k <- ncol(y) - 1

    beta0 <- matrix(0,p,k)
    beta0[,1] <- beta
    coef0 <- c(intercept,0)

    weight <- rep(1,n)
    lambda <- 0

    loss_function(x,y,weight,beta0,coef0,lambda)
}

#' @export
grad <- function(x,y,beta,intercept){
    y <- one.hot(y)
    p <- ncol(x)
    n <- nrow(x)
    k <- ncol(y) - 1

    beta0 <- matrix(0,p,k)
    beta0[,1] <- beta
    coef0 <- c(intercept,0)

    grad_cpp(x,y,beta0,coef0)
}

#' @export
predict <- function(x,beta,intercept){
    p <- ncol(x)
    n <- nrow(x)
    k <- length(intercept)

    beta0 <- matrix(0,p,k)
    beta0[,1] <- beta
    coef0 <- c(intercept,0)

    pro <- pro(x,beta0,coef0)

    return(apply(pro,1,which.max) - 1)
}

#' @export
Brier.score  <- function(x,y,beta,intercept){
    y <- one.hot(y)
    p <- ncol(x)
    n <- nrow(x)
    k <- ncol(y) - 1

    beta0 <- matrix(0,p,k)
    beta0[,1] <- beta
    coef0 <- c(intercept,0)

    pro <- pro(x,beta0,coef0)

    sum((y-pro)^2)/n
}


#' @export
TPR <- function(true_beta_idx,fit_beta_idx,p=1){
    length(intersect(true_beta_idx,fit_beta_idx)) / length(true_beta_idx)
}

#' @export
TNR <- function(true_beta_idx,fit_beta_idx,p){
    (p - length(union(true_beta_idx,fit_beta_idx))) / (p - length(true_beta_idx))
}

#' @export
ReErr <- function(true_beta,fit_beta){
    sum((true_beta - fit_beta)^2) / sum(true_beta^2)
}

#' @export
SLE <- function(true_beta_idx,fit_beta_idx){
    length(fit_beta_idx) - length(true_beta_idx)
}

#' @export
MR <- function(x,y,beta,intercept){
    sum(predict(x,beta,intercept) != y) / length(y)
}

#' @export
MAE <- function(x,y,beta,intercept){
    sum(abs(predict(x,beta,intercept) - y)) / length(y)
}
