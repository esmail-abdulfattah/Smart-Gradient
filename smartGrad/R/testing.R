#' testing Function
#'
#' This function allows you to
#' @param
#' @keywords Model
#' @export
#' @examples
#' testing()

#--------------------------------------------------------------------->


testing <- function(fn ="",  verbose = list(), xdim = NULL)
{
  setwd("~/Documents/smartGrad/R")
  source("verbose.R")
  source("wrapper.R")
  source("funcs.R")

  if(is.null(xdim)) xdim = 2
  if(is.null(verbose$MSE)) verbose$MSE = FALSE
  if(is.null(verbose$iters)) verbose$iters = FALSE
  if(is.null(verbose$wheel$show)) verbose$wheel$show = FALSE
  if(is.null(verbose$wheel$wheelparts)) verbose$wheel$wheelparts = 0.5
  if(is.null(verbose$wheel$position)) verbose$wheel$position = 6

  res.fns <- .myfuncs(fn)
  myfn <- res.fns$myfn
  myegr <- res.fns$myegr


  ## I add an argument here, just to make sure it passes through correctly...
  f1 <- myfn

  ## use simple estimates
  g1.new <- .gr.wrapper(f1, .enable = TRUE, .verbose = verbose$iters, gr.args = list(step.size = 0.00099))
  g1.plain <- .gr.wrapper(f1, .enable = FALSE)

  g1 <- function(x, ...) {
    n <- length(x)
    g <- myegr(x)

    err.new <- mean(abs(g - g1.new(x)))
    err.default <- mean(abs(g - g1.plain(x)))
    ##print(round(dig = 6, c(err.new = err.new, err.default = err.default, ratio = err.new/err.default)))

    G <- Global  #get("Global", envir = .GlobalEnv)
    G$err.trace <- c(G$err.trace, err.new - err.default)
    G$default.trace <- c(G$default.trace, err.default)
    G$new.trace <- c(G$new.trace, err.new)
    G$x_trace_x1 = c(G$x_trace_x1,x[1])
    G$x_trace_x2 = c(G$x_trace_x2,x[2])
    G$y_trace = c(G$y_trace,myfn(x))

    #assign("Global", G, envir = .GlobalEnv)
    Global <<- G
    return (g)
  }

  Global <<- list(err.trace = c(), default.trace = c(), new.trace = c(),x_trace_x1 = c(),x_trace_x2 = c(), y_trace = c())
  x_initial = rnorm(xdim, mean = 1, sd = 2)
  r.opt <- stats::optim(x_initial, f1, g1, method = "BFGS", control = list(maxit = 100000))

  cat("-> x* is : ", r.opt$par, "\n")
  cat("-> value of fn(x*) is: ", r.opt$value)


  if(verbose$MSE)
  {
    addTitle = fn
    print(.plot_MSE(c(1:length(Global$default.trace)),Global$default.trace,Global$new.trace, addTitle,FALSE))
    #plot(Global$new.trace, pch = 19, log = "y", type = "l", lwd = 3, col = "blue")
    #lines(Global$default.trace, lwd = 3, lty = 2, col = "red")

  }

  if(verbose$wheel$show)
  {
    #0.5: from 0 to pi/2 and 2 is for full wheel
    #It: getting the value of x at Iteration "It" using Global
    It = verbose$wheel$position
    resrot <- .rotate_the_Wheel(It, Global, verbose$wheel$wheelparts,f1,g1)
    print(.plot_rotated_track(resrot$errors,"",Global$new.trace[It],Global$default.trace[It]))
    print(.get_the_wheel(resrot$vects,resrot$errors))
  }


}
