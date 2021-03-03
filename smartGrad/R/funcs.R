.myfuncs <- function(fn = "")
{
  myfn <- NULL
  myegr <- NULL

  setwd("~/Documents/smartGrad/R")
  source("funcs.R")

  if(fn=="Rosenbrock_Banana"){

    myfn <- .Rosenbrock_Banana
    myegr <- .gr_Rosenbrock_Banana

  }else if(fn=="Extended_Trigonometric_dim5"){

    myfn <- .Extended_Trigonometric_dim5
    myegr <- .gr_Extended_Trigonometric_dim5

  }else if(fn=="Extended_Freudenstein_Roth"){

    myfn <- .Extended_Freudenstein_Roth
    myegr <- .gr_Extended_Freudenstein_Roth

  }else if(fn=="Perturbed_Quadratic"){

    myfn <- .Perturbed_Quadratic
    myegr <- .gr_Perturbed_Quadratic

  }else if(fn=="Raydan_1"){

    myfn <- .Raydan_1
    myegr <- .gr_Raydan_1

  }else if(fn=="Raydan_2"){

    myfn <- .Raydan_2
    myegr <- .gr_Raydan_2

  }else if(fn=="FLETCHCR"){

    myfn <- .FLETCHCR
    myegr <- .gr_FLETCHCR

  }else if(fn=="COSINE"){

    myfn <- .COSINE
    myegr <- .gr_COSINE

  }else if(fn=="Generalized_Quartic"){

    myfn <- .Generalized_Quartic
    myegr <- .gr_Generalized_Quartic

  }else if(fn=="Schumer_Steiglitz"){

    myfn <- .Schumer_Steiglitz
    myegr <- .gr_Schumer_Steiglitz

  }else if(fn=="ackley2"){

    myfn <- .ackley2
    myegr <- .gr_ackley2

  }else if(fn=="ackley3"){

    myfn <- .ackley3
    myegr <- .gr_ackley3

  }else{
    print(paste(fn, "is not added yet"))
  }

  return(list(myfn = myfn,myegr=myegr))
}

.gr.grad.default <- function(fun, x, step.size = 0.001, ...) {
  gr.args <- list(...)
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
.Rosenbrock_Banana <- function(x) {
  res <- 0.0
  for(i in 1:(length(x)-1))
    res <- res + 100*(x[i+1] - x[i]^2)^2 + (1-x[i])^2
  return(res)
}

.gr_Rosenbrock_Banana <- function(x){
  n <- length(x)
  g <- numeric(n)
  for(i in 1:(n-1)) {
    g[i] <- g[i] -400 * x[i] * (x[i+1] - x[i]^2) - 2 * (1 - x[i])
    g[i+1] <- g[i+1] + 200 * (x[i+1] - x[i]^2)
  }
  return(g)
}

.Extended_Trigonometric_dim5 <- function(x){
  res <- 0.0
  n <- length(x)
  for(i in 1:n) {
    c <- (n - sum(cos(x))) + i*(1 - cos(x[i])) - sin(x[i])
    res = res + (c*c)
  }
  return(res)
}

.gr_Extended_Trigonometric_dim5 <- function(x){
  gr <- numeric(5)
  gr[1] <- 4*cos(x[1])^2 + 2*cos(x[1])*(-6 + cos(x[2]) + cos(x[3]) +
                                          cos(x[4]) + cos(x[5]) - 7*sin(x[1])) - 2*sin(x[1])*(-46 + 8*cos(x[2]) +
                                                                                                9*cos(x[3]) + 10*cos(x[4]) + 11*cos(x[5]) + 2*sin(x[1]) + sin(x[2]) +
                                                                                                sin(x[3]) + sin(x[4])+ sin(x[5]))
  gr[2] <- 6*cos(x[2])^2 + 2*cos(x[2])*(-7 + cos(x[1]) + cos(x[3]) + cos(x[4]) +
                                          cos(x[5]) ) - 2*sin(x[2])*(-54 + 8*cos(x[1]) + 12*cos(x[2]) +
                                                                       10*cos(x[3]) + 11*cos(x[4]) + 12*cos(x[5]) + sin(x[1]) +
                                                                       3*sin(x[2]) + sin(x[3]) + sin(x[4]) + sin(x[5]))
  gr[3] <- 8*cos(x[3])^2 + 2*cos(x[3])*(-8 + cos(x[1]) + cos(x[2]) + cos(x[4]) +
                                          cos(x[5]) - 19*sin(x[3])) - 2*sin(x[3])*( -64 + 9*cos(x[1]) + 10*cos(x[2]) +
                                                                                      12*cos(x[4]) + 13*cos(x[5]) + sin(x[2]) + 4*sin(x[3]) + sin(x[4]) + sin(x[5])  )
  gr[4] <- 2*( 5*cos(x[4])^2 + cos(x[4])*(-9 + cos(x[1]) + cos(x[2]) + cos(x[3]) +
                                            cos(x[5])) - sin(x[4])*(-76 + 10*cos(x[1]) + 11*cos(x[2]) + 12*cos(x[3]) +
                                                                      28*cos(x[4]) + 14*cos(x[5]) + sin(x[1]) + sin(x[2]) + sin(x[3]) + 5*sin(x[4]) + sin(x[5])))
  gr[5] <- 2*(6*cos(x[5])^2 + cos(x[5])* (-10 + cos(x[1]) + cos(x[2]) + cos(x[3]) + cos(x[4])
                                          - 39*sin(x[5])) -sin(x[5])*(-90 + 11*cos(x[1]) + 12*cos(x[2]) + 13*cos(x[3]) +
                                                                        14*cos(x[4]) + sum(sin(x[1:4])) + 6*sin(x[5])))

  return(gr)
}

.Extended_Freudenstein_Roth <- function(x){
  res <- 0.0
  n <-length(x)
  for(i in 1:(n/2))
  {
    c1 <- (-13 + x[(2*i -1)] + ((5 - x[2*i])*x[2*i] - 2)*x[2*i])^2
    c2 <- (-29 + x[(2*i-1)] + ((x[2*i] + 1)*x[2*i] -14)*x[2*i])^2
    res <- res + (c1 + c2)
  }
  return(res)
}

.gr_Extended_Freudenstein_Roth <- function(x){
  n <-length(x)
  gr <- numeric(n)
  gr[1] <- 4*(-21 +x[1] - 8*x[2] + 3*x[2]^2)
  for(i in 2:n)
  {
    if((i %% 2) == 0){
      gr[i] <- 4*(216 -8*x[i-1] + 6*x[i] + 6*x[i-1]*x[i] - 60*x[i]^2 + 2*x[i]^3 -10*x[i]^4 + 3*x[i]^5)
    } else {
      if(n>=(i+1)) gr[i] <- 4*(-21 +x[i] - 8*x[i+1] + 3*x[i+1]^2)
    }
  }

  return(gr)
}

.Perturbed_Quadratic <- function(x){
  res <- 0.0
  n <-length(x)
  c1 <- 1:n
  c2 <- x^2
  return(sum(c1*c2) + (1/100)*sum(x)^2)
}

.gr_Perturbed_Quadratic <- function(x){
  gr <- numeric(length(x))
  n <-length(x)
  c <- sum(x)
  for(i in 1:n) gr[i] <- 2*i*x[i] + (1/50)*c
  return(gr)
}

.Raydan_1 <- function(x){
  n <-length(x)
  return(sum(((1:n)/10)*(exp(x) -x)))
}

.gr_Raydan_1 <- function(x){
  n <-length(x)
  gr <- numeric(n)
  for(i in 1:n) gr[i] <- (i/10)*(-1 + exp(x[i]))
  return(gr)
}

.Raydan_2 <- function(x){
  n <-length(x)
  return(sum((exp(x) -x)))
}

.gr_Raydan_2 <- function(x){
  n <-length(x)
  gr <- numeric(n)
  for(i in 1:n) gr[i] <- (-1 + exp(x[i]))
  return(gr)
}

.FLETCHCR <- function(x){
  res <- 0.0
  n <-length(x)
  for(i in 1:(n-1)) res <- res + 100*(x[i+1] - x[i] + 1 - x[i]^2)^2
  return(res)
}

.gr_FLETCHCR <- function(x){
  n <-length(x)
  gr <- numeric(n)
  gr[1] <- 200*(-1 - x[1] + 3*x[1]^2 + 2*x[1]^3 - x[2]  - 2*x[1]*x[2])

  for(i in 2:n)
  {
    if(i==2){
      if(n < 3){
        gr[i] <- -200*(1 + x[1] -x[2] + x[1]^2 )
      }else {
        gr[i] <- -200*(x[1] + x[1]^2 - 3*x[2]^2 - 2*x[2]^3 + x[3] + 2*x[2]*x[3])
      }
    } else {
      if(n>=(i+1)) {
        gr[i] <- -200*(x[i-1] + x[i-1]^2 - 3*x[i]^2 - 2*x[i]^3 + x[i+1] + 2*x[i]*x[i+1])
      } else {
        gr[i] <- -200*(-1 + x[i-1] + x[i-1]^2 -x[i]  )
      }
    }
  }
  return(gr)
}

.COSINE <- function(x){
  res <- 0.0
  n <-length(x)
  for(i in 1:(n-1)) res <- res + cos(-0.5*x[i+1] + x[i]^2)
  return(res)
}

.gr_COSINE <- function(x){
  n <-length(x)
  gr <- numeric(n)
  gr[1] <- -2*sin(x[1]^2 - 0.5*x[2])*x[1]

  for(i in 2:n)
  {
    if(n>=(i+1)) {
      gr[i] <- 0.5*(sin(x[i-1]^2 - 0.5*x[i]) - 4*sin(x[i]^2 - 0.5*x[i+1])*x[i])
    } else {
      gr[i] <- 0.5*sin(x[i-1]^2 - 0.5*x[i])
    }
  }
  return(gr)
}

.Generalized_Quartic <- function(x){
  res <- 0.0
  n <-length(x)
  for(i in 1:(n-1)) res <- res + (x[i]^2 + (x[i+1] + x[i]^2)^2)
  return(res)
}

.gr_Generalized_Quartic <- function(x){
  n <-length(x)
  gr <- numeric(n)
  gr[1] <- 2*(x[1] + 2*x[1]^3 + 2*x[1]*x[2])
  for(i in 2:n)
  {
    if(n>=(i+1)) {
      gr[i] <- 2*(x[i-1]^2 + 2*x[i] + 2*x[i]^3 + 2*x[i]*x[i+1])
    } else {
      gr[i] <- 2*(x[i-1]^2 + x[i])
    }
  }
  return(gr)
}

.Schumer_Steiglitz <- function(x){
  res <- 0.0
  n <-length(x)
  return(sum(x[1:n]^4))
}

.gr_Schumer_Steiglitz <- function(x){
  n <-length(x)
  gr <- numeric(n)
  gr[1:n] = 4*x[1:n]^3
  return(gr)

}

.ackley2 <- function(x){
  return(-200*exp(-0.02*sqrt(x[1]^2+ x[2]^2)))
}

.gr_ackley2 <- function(x){
  gr <- numeric(2)
  c1 <- sqrt(x[1]^2+ x[2]^2)
  c2 <- exp(-0.02*c1)
  gr[1] = 4*x[1]*c2/c1
  gr[2] = 4*x[2]*c2/c1

  return(gr)
}

.ackley3 <- function(x){
  return(200*exp(-0.02*sqrt(x[1]^2 + x[2]^2)) + 5*exp(cos(3*x[1])+sin(3*x[2])))
}

.gr_ackley3 <- function(x){
  gr <- numeric(2)
  c1 <- sqrt(x[1]^2+ x[2]^2)
  c2 <- -4*exp(-0.02*c1)/c1
  c3 <- 15*exp(cos(3*x[1]) + sin(3*x[2]))
  gr[1] = -c3*sin(3*x[1]) + x[1]*c2
  gr[2] = c3*cos(3*x[2]) + x[2]*c2

  return(gr)
}

