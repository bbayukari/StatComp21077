#' @title Fit the gamma model with approximate newton method.
#' @description Use approximate newton method fit the gamma model which is an example of generalized linear model.
#' @param x Input matrix, of dimension \eqn{n \times p}; each row is an observation
#' vector and each column is a predictor/feature/variable.
#' @param y The response variable, of \code{n} observations.
#' @param weight Observation weights. When \code{weight = NULL},
#' we set \code{weight = 1} for each observation as default.
#' @param start the start point of fitness.
#' @param lambda A single lambda value for L2 regularized. Default is 0.
#'
#' @return Intercept and linear coefficient.
#' @export
#'
#' @examples
#' \dontrun{
#' n <- 1000
#' p <- 10
#' dataset <- generate.data(n, p, p,family = "gamma", seed = 1)
#' coef <- fit.gamma.appro(dataset[["x"]],dataset[["y"]])
#' }
fit.gamma.appro <- function(x, y, weight = NULL, start = NULL, lambda = 0){
   # check x
    if (is.data.frame(x)) {
        x <- as.matrix(x)
    }
    if (!is.numeric(x)) {
        stop("x must be a *numeric* matrix/data.frame!")
    }
    if (anyNA(x) || any(is.infinite(x))) {
        stop("x has missing value or infinite value!")
    }
    nvars <- ncol(x)
    nobs <- nrow(x)

    # check weight
    if (is.null(weight)) {
        weight <- rep(1, nobs)
    }
    stopifnot(is.vector(weight))
    if (length(weight) != nobs) {
        stop("Rows of x must be the same as length of weight!")
    }
    stopifnot(all(is.numeric(weight)), all(weight >= 0))

    # check y:
    if (anyNA(y)) {
        stop("y has missing value!")
    }
    if (any(is.infinite(y))) {
        stop("y has infinite value!")
    }
    if (any(y < 0)) {
      stop("y must be positive value when family = 'gamma'.")
    }

    # check lambda:
    if(lambda < 0){
        stop("lambda must be positive.")
    }

    # check start:
    if(is.null(start)) {
        start <- rep(0,nvars + 1)
    }
    stopifnot(is.vector(weight))
    if(length(start) != nvars + 1) {
        stop("start must be match with the number of variable data.")
    }

    # compute
    coef <- gamma_fit_approximate_newton_method(x,y,weight,start[-1],start[1],lambda)

    return (coef)
}

#' @title Fit the gamma model with IWLS method.
#' @description Use IWLS method fit the gamma model which is an example of generalized linear model.
#' @param x Input matrix, of dimension \eqn{n \times p}; each row is an observation
#' vector and each column is a predictor/feature/variable.
#' @param y The response variable, of \code{n} observations.
#' @param weight Observation weights. When \code{weight = NULL},
#' we set \code{weight = 1} for each observation as default.
#' @param start the start point of fitness.
#' @param lambda A single lambda value for L2 regularized. Default is 0.
#'
#' @return Intercept and linear coefficient.
#' @export
#'
#' @examples
#' \dontrun{
#' n <- 1000
#' p <- 10
#' dataset <- generate.data(n, p, p,family = "gamma", seed = 1)
#' coef <- fit.gamma.IWLS(dataset[["x"]],dataset[["y"]])
#' }
fit.gamma.IWLS <- function(x, y, weight = NULL, start = NULL, lambda = 0){
  # check x
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x)) {
    stop("x must be a *numeric* matrix/data.frame!")
  }
  if (anyNA(x) || any(is.infinite(x))) {
    stop("x has missing value or infinite value!")
  }
  nvars <- ncol(x)
  nobs <- nrow(x)

  # check weight
  if (is.null(weight)) {
    weight <- rep(1, nobs)
  }
  stopifnot(is.vector(weight))
  if (length(weight) != nobs) {
    stop("Rows of x must be the same as length of weight!")
  }
  stopifnot(all(is.numeric(weight)), all(weight >= 0))

  # check y:
  if (anyNA(y)) {
    stop("y has missing value!")
  }
  if (any(is.infinite(y))) {
    stop("y has infinite value!")
  }
  if (any(y < 0)) {
    stop("y must be positive value when family = 'gamma'.")
  }

  # check lambda:
  if(lambda < 0){
    stop("lambda must be positive.")
  }

  # check start:
  if(is.null(start)) {
    start <- rep(0,nvars + 1)
  }
  stopifnot(is.vector(weight))
  if(length(start) != nvars + 1) {
    stop("start must be match with the number of variable data.")
  }

  # compute
  coef <- gamma_fit_IWLS_method(x,y,weight,start[-1],start[1],lambda)

  return (coef)
}
