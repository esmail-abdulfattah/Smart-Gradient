#' optim Function
#'
#' @description The function to compute gradient for, is a wrapper function containing 'fn' hence do not need these extra arguments. For every call to 'gr', the function to compute the gradient of will change, and the point to compute the gradient is always rep(0, length(x)).
#'
#' @param fn A function to be minimized (or maximized),
#' @param gr is an optional generic gradient function of type: gr(fun, x, ...) and it returns the gradient for the "BFGS"
#'
#' The default gradient function use a central differences with fixed step.size (argument step.size=0.001)
#' @param gr.args Specific arguments to 'gr' needs to be passed.
#' @param smart Logical. Do you want to applay Smart Gradient Technique on the gradient. Default is TRUE.
#' @param '...' are optional arguments to 'fn': fn(x, ...)
#' @param control a list of control parameters. See stats::optim function for details.
#' @param hessian Logical. Do you want to return the hessian?
#' @keywords Model
#' @export
#' @examples
#' myfun <- function(x) { 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2}
#'
#' mygrad  <- function(fun,x) {
#'  h = 0.001
#'  grad <- numeric(2)
#'  grad[1] <- (fun(x + c(h,0)) - fun(x - c(h,0))) / (2 * h)
#'  grad[2] <- (fun(x + c(0,h)) - fun(x - c(0,h))) / (2 * h)
#'  return (grad)
#'}
#'
#' x_initial = rnorm(2)
#' result  <- smartGrad::optim(par = x_initial,
#'                            fn=myfun,
#'                            gr = mygrad,
#'                            method = c("BFGS"),
#'                            smart = TRUE)
#' @details More details about this function can be found in stats::optim.

#--------------------------------------------------------------------->

optim <- function(par, fn, gr = NULL, ...,
      control = list(), hessian = FALSE, smart = TRUE, gr.args = list())
{

  library("stats")
  method = "BFGS"
  lower = -Inf; upper = Inf
  g1 <- .gr.wrapper(fn, gr = gr, .enable = smart, .verbose = FALSE, ..., gr.args = gr.args)
  r.opt <- stats::optim(par=par, fn=fn, gr = g1,
                 method = method,
                 lower = lower, upper = upper, ...,
                 control = control, hessian = hessian)
  .reset.gr.wrapper(fn, gr = gr, .enable = smart, .verbose = FALSE, ..., gr.args = gr.args)
  return(r.opt)
}




