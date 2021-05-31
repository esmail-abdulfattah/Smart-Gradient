#' optim Function
#'
#' This function allows you to
#' @param
#' @keywords Model
#' @export
#' @examples
#' optim()

#--------------------------------------------------------------------->

optim <- function(par, fn, gr = NULL, ...,
      method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                 "Brent"),
      lower = -Inf, upper = Inf,
      control = list(), hessian = FALSE, smart = TRUE, gr.args = list())
{

  library("stats")
  g1 <- .gr.wrapper(fn, gr = gr, .enable = smart, .verbose = FALSE, ..., gr.args = gr.args)
  r.opt <- stats::optim(par=par, fn=fn, gr = g1,
                 method = method,
                 lower = lower, upper = upper, ...,
                 control = control, hessian = hessian)
}

