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
#' library(lbfgs)
#' library(stats)
#'
#' myfun <- function(x) { 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2}
#'
#' mygrad  <- function(fun,x) {
#'   h = 0.001
#'   grad <- numeric(2)
#'   grad[1] <- (fun(x + c(h,0)) - fun(x - c(h,0))) / (2 * h)
#'   grad[2] <- (fun(x + c(0,h)) - fun(x - c(0,h))) / (2 * h)
#'   return (grad)
#' }
#'
#' mySmartgrad = smartGrad::makeSmart(fn = myfun,
#'                                    gr = mygrad)
#'
#' x_initial = rnorm(2)
#' result1 <- stats::optim(par = x_initial,
#'                         fn = myfun,
#'                         gr = mySmartgrad,
#'                         method = c("BFGS"),
#'                         hessian = FALSE)
#'
#' result2 <- lbfgs::lbfgs(myfun, mySmartgrad, vars = x_initial)

#--------------------------------------------------------------------->


makeSmart <- function(fn, gr, ..., gr.args = list())
{
  stopifnot(!is.null(fn))
  stopifnot(!is.null(gr))

  g1 <- .gr.wrapper(fn, gr = gr, .enable = TRUE, .verbose = FALSE, ..., gr.args = gr.args)
  return(g1)
}

