gr.wrapper <- function (fn = NULL, gr = NULL, enable = TRUE, verbose = FALSE, ...) {
    stopifnot(!is.null(fn))

    fun <- local({

        gr.grad.default <- function(fun, x, h = 0.001) {
            n <- length(x)
            grad <- numeric(n)
            e <- rep(0, n)
            for(i in 1:n) {
                e[] <- 0
                e[i] <- 1
                grad[i] <- (fun(x + h * e) - fun(x - h * e)) / (2 * h)
            }
            return (grad)
        }

        grw <- list()
        grw$fn <- fn
        grw$gr <- if (is.null(gr)) gr.grad.default else gr
        grw$par.prev <- c()
        grw$n <- 0
        grw$A <- matrix()
        grw$AA <- matrix()
        grw$enable <- enable
        grw$verbose <- verbose
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
                
                if((i + 1) <= n) {
                    for (j in (i + 1):n) {
                        r <- sum(q * G[,j]) 
                        G[,j] <- G[,j] - r*q
                    }
                }
            }
            return(G)
        }
        
        gradient <- function(par, ...) {
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


            if (TRUE) {
                grw$par <- par
                tmp.fn <- function(x) {
                    n <- length(x)
                    dpar <- rep(0, n)
                    for(i in 1:n) {
                        dpar <- dpar + x[i] * grw$AA[, i]
                    }
                    return (grw$fn(grw$par + dpar))
                }
                ## yes, always evaluate the gradient in zero(-vector)
                gg <- grw$gr(tmp.fn, rep(0, grw$n))
            } else {
                ## old code
                gg <- numeric(grw$n)
                for(i in 1:grw$n) {
                gg[i] <- (grw$fn(par + grw$step.len * grw$AA[, i]) -
                          grw$fn(par - grw$step.len * grw$AA[, i])) / (2 * grw$step.len)
                }
            }
            grad <- solve(t(grw$AA), gg)
            grw.par <<- grw
            return(grad)
        }
        return(gradient)
    })
    return (fun)
}
f1 <- function(x) {   ## Rosenbrock Banana function with higher dimension 
  res = 0.0
  for(i in 1:(length(x)-1))
      res = res + 100*(x[i+1] - x[i]^2)^2 + (1-x[i])^2
  return(res)
}

if (TRUE) {
    ## use simple estimates
    g1.new <- gr.wrapper(f1, enable = TRUE, verbose = FALSE)
    g1.plain <- gr.wrapper(f1, enable = FALSE)
} else {
    ## use good estimates
    library(numDeriv)
    g1.new <- gr.wrapper(f1, gr = grad, enable = TRUE, verbose = FALSE)
    g1.plain <- gr.wrapper(f1, gr = grad, enable = FALSE)
}

g1 <- function(x) {
    n <- length(x)
    g <- numeric(n)
    for(i in 1:(n-1)) {
        g[i] <- g[i] -400 * x[i] * (x[i+1] - x[i]^2) - 2 * (1 - x[i])
        g[i+1] <- g[i+1] + 200 * (x[i+1] - x[i]^2)
    }

    err.new <- mean(abs(g - g1.new(x)))
    err.default <- mean(abs(g - g1.plain(x)))
    ##print(round(dig = 6, c(err.new = err.new, err.default = err.default, ratio = err.new/err.default)))

    G <- get("Global", envir = .GlobalEnv)
    G$err.trace <- c(G$err.trace, err.new - err.default)
    G$default.trace <- c(G$default.trace, err.default)
    G$new.trace <- c(G$new.trace, err.new)
    assign("Global", G, envir = .GlobalEnv)

    return (g)
}

Global <- list(err.trace = c(), default.trace = c(), new.trace = c())
dim <- 5
x_initial = rnorm(dim, mean = 1, sd = 2)
r.opt <- optim(x_initial, f1, g1, method = "BFGS", control = list(maxit = 100000))

print(r.opt$value)
print(r.opt$par)

plot(Global$new.trace, pch = 19, log = "y", type = "l", lwd = 3, col = "blue")
lines(Global$default.trace, lwd = 3, lty = 2, col = "red")
