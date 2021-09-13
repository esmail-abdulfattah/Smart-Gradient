

smartGrad::testing(fn="Rosenbrock_Banana",verbose = list(MSE = TRUE), xdim = 50)
smartGrad::testing(fn="Rosenbrock_Banana",verbose = list(iters = TRUE))
smartGrad::testing(fn="Rosenbrock_Banana",verbose = list(wheel= list(show = TRUE, wheelparts = 2,position = 6)))


myfun <- function(x) {
  res = 0.0
  for(i in 1:(length(x)-1))
    res = res + 100*(x[i+1] - x[i]^2)^2 + (1-x[i])^2
  return(res)
}

mygrad <- function(x) {
  step.size = 0.001
  n <- length(x)
  grad <- numeric(n)
  e <- rep(0, n)
  for(i in 1:n) {
    e[] <- 0
    e[i] <- 1
    grad[i] <- (myfun(x + step.size * e) - myfun(x - step.size * e)) / (2 * step.size)
  }
  return (grad)
}

par = rnorm(2,0,2)
df = stats::optim(par, fn=myfun,
                 gr = mygrad,
                 method = c("BFGS"),
                 hessian = FALSE)


gr.grad.default <- function(fun,x) {
  step.size = 0.001
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

st = smartGrad::optim(par, fn=myfun,
                        gr = gr.grad.default,
                        method = c("BFGS"),
                        hessian = FALSE,
                        smart = TRUE)
df$par
st$par





