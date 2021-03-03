#' optim Function
#'
#' This function allows you to
#' @param #love Do you love cats? Defaults to TRUE.
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
  #lapply(sys.call()[-1], as.character)$fn
  #if (!is.null(gr)) formals(gr)$fun = fn

  library("stats")
  g1 <- .gr.wrapper(fn, gr = gr, .enable = smart, .verbose = FALSE, ..., gr.args = gr.args)
  r.opt <- stats::optim(par=par, fn=fn, gr = g1,
                 method = method,
                 lower = lower, upper = upper, ...,
                 control = control, hessian = hessian)
}

