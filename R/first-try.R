gr.wrapper <- function (fn = NULL, enable = TRUE, verbose = FALSE, ...) {
    stopifnot(!is.null(fn))

    fun <- local({
        grw <- list()
        grw$fn <- fn
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
                dpar <- scale(par - grw$par.prev, center = FALSE)
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

            gg <- numeric(grw$n)
            for(i in 1:grw$n) {
                gg[i] <- (grw$fn(par + grw$step.len * grw$AA[, i]) -
                          grw$fn(par - grw$step.len * grw$AA[, i])) / (2 * grw$step.len)
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

g1.new <- gr.wrapper(f1, enable = TRUE, verbose = TRUE)
g1.plain <- gr.wrapper(f1, enable = FALSE)

g1 <- function(x) {
    n <- length(x)
    g <- numeric(n)
    for(i in 1:(n-1)) {
        g[i] <- g[i] -400 * x[i] * (x[i+1] - x[i]^2) - 2 * (1 - x[i])
        g[i+1] <- g[i+1] + 200 * (x[i+1] - x[i]^2)
    }

    err.new <- mean(abs(g - g1.new(x)))
    err.default <- mean(abs(g - g1.plain(x)))
    print(round(dig = 6, c(err.new = err.new, err.default = err.default, ratio = err.new/err.default)))

    G <- get("Global", envir = .GlobalEnv)
    G$r.trace <- c(G$r.trace, err.new - err.default)
    assign("Global", G, envir = .GlobalEnv)

    return (g)
}

Global <- list(r.trace = c())
dim = 5
x_initial = rnorm(dim, sd = 2)
r.opt <- optim(x_initial, f1, g1, method = "BFGS", control = list(maxit = 100000))

print(r.opt$value)
print(r.opt$par)

plot(Global$r.trace, pch = 19)
abline(h = 0)
