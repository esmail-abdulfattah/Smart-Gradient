#' makeSmart Function
#'
#' This function allows you to convert a non-smart numerical gradient to a SMART Gradient
#' @param fn A function to be minimized (or maximized),
#' @param gr is an optional generic gradient function of type: gr(fun, x, ...) and it returns the gradient for the "BFGS"
#' @param ... Further arguments to be passed to fn and gr.
#' @param gr.args Specific arguments to 'gr' needs to be passed.
#' @keywords Model
#' @export
#' @examples
#'
#' myfun <- function(x) {
#' res <- 0.0
#' for(i in 1:(length(x)-1))
#' res <- res + 100*(x[i+1] - x[i]^2)^2 + (1-x[i])^2
#' return(res)
#' }
#' mygrad <- function(wrapper.myfun,x){
#'    h = 1e-3
#'    grad <- numeric(length(x))
#'    for(i in 1:length(x)){
#'      e = numeric(length(x))
#'      e[i] = 1
#'      grad[i] <- (wrapper.myfun(x+h*e) - wrapper.myfun(x-h*e))/(2*h)
#'    }
#'    return(grad)
#' }
#' 
#'
#' library("stats")
#' library("smartGrad")
#' x_dimension = 5
#' x_initial = rnorm(x_dimension)
#' result <- optim(par = x_initial,
#'                 fn = myfun,
#'                 gr = makeSmart(fn = myfun,gr = mygrad),
#'                 method = c("BFGS"))

#--------------------------------------------------------------------->


makeSmart <- function(fn, gr, ..., gr.args = list())
{
  stopifnot(!is.null(fn))
  stopifnot(!is.null(gr))

  g1 <- .gr.wrapper(fn, gr = gr, .enable = TRUE, .verbose = FALSE, ..., gr.args = gr.args)
  return(g1)
}

