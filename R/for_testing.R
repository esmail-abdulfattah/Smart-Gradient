

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
      
      # gg <- numeric(grw$n)
      # for(i in 1:grw$n) {
      #   gg[i] <- (grw$fn(par + grw$step.len * grw$AA[, i]) -
      #               grw$fn(par - grw$step.len * grw$AA[, i])) / (2 * grw$step.len)}
      
      trans_fn <- function(x) {Gx = grw$AA%*%x; grw$fn(Gx)}
      gg = mygrad(trans_fn, solve(grw$AA,par))

      grad <- solve(t(grw$AA), gg)
      grw.par <<- grw
      return(grad)
    }
    return(gradient)
  })
  return (fun)
}

#### ---> Functions Names
# Rosenbrock_Banana
# Extended_Trigonometric_dim5
# Extended_Freudenstein_Roth
# Perturbed_Quadratic #D
# Raydan_1 #D
# Raydan_2
# FLETCHCR #D
# COSINE #D
# Generalized_Quartic #D

myFun <- COSINE ### "Function Name"
getExactGrad <- gr_COSINE ### gr_"Function Name"

myGrad.new <- gr.wrapper(myFun, enable = TRUE, verbose = FALSE)
myGrad.plain <- gr.wrapper(myFun, enable = FALSE)

myGrad <- function(x) {
  n <- length(x)
  g <- getExactGrad(x)
  
  err.new <- mean(abs(g - myGrad.new(x)))
  err.default <- mean(abs(g - myGrad.plain(x)))
  ##print(round(dig = 6, c(err.new = err.new, err.default = err.default, ratio = err.new/err.default)))
  
  G <- get("Global", envir = .GlobalEnv)
  G$err.trace <- c(G$err.trace, err.new - err.default)
  G$default.trace <- c(G$default.trace, err.default)
  G$new.trace <- c(G$new.trace, err.new)
  
  G$x_trace_x1 = c(G$x_trace_x1,x[1]) 
  G$x_trace_x2 = c(G$x_trace_x2,x[2]) 
  G$y_trace = c(G$y_trace,Rosenbrock_Banana(x))
  assign("Global", G, envir = .GlobalEnv)
  
  return (g)
}

iterations = 1
new_avg_MSE <- numeric(iterations)
default_avg_MSE <- numeric(iterations)

mygrad <- function(f,x)
{
  h = 1e-3
  res = numeric(length(x))
  for(i in 1:length(x))
  {
    e = numeric(length(x))
    e[i] = 1
    res[i] = (f(x+h*e) - f(x-h*e))/(2*h)
  }
  return(res)
}

for(iter in 1:iterations)
{
  Global <- list(err.trace = c(), default.trace = c(), new.trace = c(), x_trace_x1 = c(),x_trace_x2 = c(), y_trace = c())
  dim = 2
  x_initial = rnorm(dim, mean = 1, sd = 2)
  r.opt <- optim(x_initial, fn = myFun, gr = myGrad, method = "BFGS", control = list(maxit = 100000))
  
  print(r.opt$value)
  print(r.opt$par)
  
  addTitle = "Rosenbrock Banana Function dimension 5"
  plot_MSE(c(1:length(Global$default.trace)),Global$default.trace,Global$new.trace, addTitle,FALSE)
  
  default_avg_MSE[iter] <- mean(Global$default.trace)
  new_avg_MSE[iter] <- mean(Global$new.trace)
}

print(round(dig = 6,mean(default_avg_MSE)))
print(round(dig = 6,mean(new_avg_MSE)))


