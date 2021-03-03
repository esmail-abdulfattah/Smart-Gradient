.gr.wrapper <- function(fn = NULL, gr = NULL, gr.args = list(), ...,
                       .enable = TRUE, .verbose = FALSE)
{
  ## '...' are optional arguments to 'fn': fn(x, ...)
  ## 'gr' is an optional generic gradient function of type: gr(fun, x, ...) where spesific
  ## arguments to 'gr' needs to be passed in 'gr.args'. The '...' arguments to 'fn' in
  ## 'gr.wrapper' are also passed to 'gr' since this is how its done in 'optim'. The function
  ## to compute gradient for, is a wrapper function containing 'fn' hence do not need these
  ## extra arguments. For every call to 'gr', the function to compute the gradient of will
  ## change, and the point to compute the gradient is always rep(0, length(x)).

  ## the default gradient function use a central differences with fixed step.size (argument
  ## step.size=0.001)

  stopifnot(!is.null(fn))
  is.empty.list <- function(a) {
    if (!(is.null(a) || is.list(a))) {
      ## wrong format, ignore
      return (TRUE)
    }
    return (is.null(a) || (is.list(a) && length(a) == 0))
  }

  fun <- local({

    gr.grad.default <- function(fun, x, step.size = 0.001, ...) {
      not.used.args <- list(...)
      n <- length(x)
      grad <- numeric(n)
      e <- rep(0, n)
      for(i in 1:n) {
        e[] <- 0
        e[i] <- 1
        grad[i] <- (fun(x + step.size * e) - fun(x - step.size * e)) / (2 * step.size)
      }
      return (grad)
    }

    grw <- list()
    grw$fn <- fn
    grw$fn.arg.name <- names(formals(fn))[1]
    grw$fn.args <- list(...)
    grw$gr <- if (is.null(gr)) gr.grad.default else gr
    grw$gr.arg.name <- names(formals(grw$gr))[2]
    grw$gr.args <- if (is.empty.list(gr.args)) list() else gr.args

    grw$par.prev <- c()
    grw$n <- 0
    grw$A <- matrix()
    grw$AA <- matrix()
    grw$enable <- .enable
    grw$verbose <- .verbose
    grw$step.len <- 0.001
    grw$eps.sd <- .Machine$double.eps^(1/5)
    grw$iter <- 0
    grw.par <- grw

    MGS <- function(G) {
      n <- dim(G)[1]
      q <- numeric(n)
      for (i in 1:n) {
        r <- sqrt(sum(G[,i]*G[,i]))
        q <- G[,i]/r
        G[,i] <- q

        if ((i + 1) <= n) {
          for (j in (i + 1):n) {
            r <- sum(q * G[,j])
            G[,j] <- G[,j] - r*q
          }
        }
      }
      return(G)
    }

    gradient <- function(par, ...) {
      gr.args <- list(...)
      grw <<- grw.par
      grw$iter <- grw$iter + 1

      first.time <- FALSE
      if (grw$n == 0) {
        first.time <- TRUE
        grw$par.prev <- par
        grw$n <- length(par)
        if (grw$enable) {
          grw$A <- matrix(rnorm(grw$n^2, sd = grw$eps.sd), grw$n, grw$n)
          diag(grw$A) <- diag(grw$A) + 1
        } else {
          grw$A <- diag(grw$n)
        }
        grw$AA <- grw$A
      }

      if (!first.time && grw$enable) {
        grw$A[, 2:grw$n] <- grw$A[, 1:(grw$n-1)]
        dpar <- scale(par - grw$par.prev)
        grw$A[, 1] <- dpar + rnorm(grw$n, sd = grw$eps.sd)
        grw$par.prev <- par
        grw$AA <- MGS(grw$A)
        if (grw$verbose) {
          print(paste0("Iteration: ", grw$iter))
          rownames(grw$AA) <- paste0("par", 1:grw$n)
          colnames(grw$AA) <- paste0("dir", 1:grw$n)
          print(round(dig = 2, grw$AA))
        }
      }

      grw$par <- par
      tmp.fn <- function(x) {
        dpar <- drop(grw$AA %*% x)
        x <- grw$par + dpar
        args <- c(list(x), grw$fn.args)
        nm <- names(args)
        nm[1] <- grw$fn.arg.name
        names(args) <- nm
        return (do.call(grw$fn, args = args))
      }
      x <- rep(0, grw$n)
      args <- c(tmp.fn, list(x), grw$gr.args, ...)
      nm <- names(args)
      nm[1:2] <- c("", grw$gr.arg.name)
      names(args) <- nm
      gg <- do.call(grw$gr, args = args)
      grad <- solve(t(grw$AA), gg)
      grw.par <<- grw
      return(grad)
    }
    return(gradient)
  })
  return (fun)
}
