#' @export
demo <- function(x,y,beta,string,n,m){
    Rcpp::sourceCpp(code = string)
    loss <- get_loss_function()
    gradient <- get_gradient_function()
    list(
        loss_API(n,x,y,beta,loss),
        gradient_API(m,x,y,beta,gradient)
    )
}
