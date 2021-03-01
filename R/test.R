fn <- function(x, ...) sum(x^2)
gr <- function(x, ...) {
    print(list(...))
    return (2*x)
}

x <- c(1, 2, 3)
a <- optim(x, fn, gr, method = "BFGS")
a <- optim(x, fn, gr, method = "BFGS", AA = 3)
