plot_MSE <- function(x1,y1,y2,s,addLineType=FALSE){
  library(ggplot2)
  library(viridis)
  
  df.merged <- rbind(data.frame(x = x1, y = y1), data.frame(x = x1, y = y2))
  df.merged$clas = c(rep(c("Plain","New Approach"), each = length(x1)))
  p = ggplot(NULL)
  if(addLineType) {p = p +geom_line(data = df.merged,aes(x = x, y = y, color = clas,linetype = clas))
  } else {p = p +geom_line(data = df.merged,aes(x = x, y = y, color = clas))}
  
  p = p + scale_size_manual(values = c(1,2)) +
    #scale_color_manual(values = plasma(2)) +
    theme(legend.position="bottom", legend.direction="horizontal", 
          legend.title = element_blank(),panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white")) +
    labs(y="Mean Absolute Error", x = "Iterations") + 
    ggtitle(s)
  
    p = p + theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "grey"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "grey"))
  
  return(p)
}
get_the_wheel <- function(vects,errors){
  play_with_colors <- function(c)
  {
    c = 10000*c
    max_c = max(c)
    min_c = min(c)
    mid_c = 0.5*(max_c + min_c)
    mid_c_2 =  0.5*(min_c + mid_c)
    mid_c_5 =  0.5*(mid_c + max_c)
    mid_c_1 =  0.5*(min_c + mid_c_2)
    mid_c_3 =  0.5*(mid_c_2 + mid_c)
    mid_c_4 =  0.5*(mid_c + mid_c_5)
    mid_c_6 =  0.5*(mid_c_5 + max_c)
    
    # val1 = paste("<=",round(mid_c_1,8))
    # val2 = paste(">",round(mid_c_1,8), "& <", round(mid_c_2,8))
    # val3 = paste(">",round(mid_c_2,8), "& <", round(mid_c_3,8))
    # val4 = paste(">",round(mid_c_3,8), "& <", round(mid_c,8))
    # val5 = paste(">",round(mid_c,8), "& <", round(mid_c_4,8))
    # val6 = paste(">",round(mid_c_4,8), "& <", round(mid_c_5,8))
    # val7 = paste(">",round(mid_c_5,8), "& <", round(mid_c_6,8))
    # val8 = paste(">=",round(mid_c_6,8))
    
    res = c()
    for(i in 1:length(c))
    {
      if(c[i] <= mid_c_1) {
        res = c(res, "1")
      } else if(c[i] > mid_c_1 && c[i] <= mid_c_2){
        res = c(res, "2")
      } else if(c[i] > mid_c_2 && c[i] <= mid_c_3){
        res = c(res, "3")
      } else if(c[i] > mid_c_3 && c[i] <= mid_c){
        res = c(res, "4")
      } else if(c[i] > mid_c && c[i] <= mid_c_4){
        res = c(res, "5")
      } else if(c[i] > mid_c_4 && c[i] <= mid_c_5){
        res = c(res, "6")
      } else if(c[i] > mid_c_5 && c[i] <= mid_c_6){
        res = c(res, "7")
      } else{ 
        res = c(res,"8")}
    }
    
    return(res)
  }
  
  d=data.frame(x=rep(0,dim(vects)[1]), y=rep(0,dim(vects)[1]), vx=vects[,1], vy=vects[,2],Error_Increasing = play_with_colors(errors))
  p = ggplot() + 
    geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy,color = Error_Increasing), arrow=arrow(), size=1) + 
    geom_point(data=d, mapping=aes(x=x, y=y), size=4, shape=21, fill="white")   +
    theme(legend.position="bottom", legend.direction="horizontal",panel.grid.major = element_line(size = 0.5,
      linetype = 'solid',colour = "white"),panel.background = element_rect(fill = "white", colour = "black")) +
    labs(y="", x = "Vectors Directions")  + 
    scale_color_manual(values = inferno(8)) + labs(color='Higher number represents higher error')
  
  return(p)
  
}
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
rotate_the_Wheel <- function(It = 4, Global, wheel_parts = 2){
  trans_fn <- function(x) {Gx = G%*%x; myFun(Gx)}
  xnode = c(Global$x_trace_x1[It],Global$x_trace_x2[It])
  gexact = getExactGrad(xnode)
  errors = c()
  a = 0; i =1
  vects = c(); step_angle = 0.01
  if(wheel_parts<2) wheel_parts = wheel_parts + step_angle
  while(a<wheel_parts)
  {
    angle = pi*a
    ROT = matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),2,2)
    rot_line = ROT%*%c(1,0)
    a = a + step_angle
    G = MGS(matrix(c(rot_line[1],rot_line[2],0,1),2,2))
    gest = solve(t(G), mygrad(trans_fn, solve(G,xnode)))
    errors = c(errors,mean(abs(gexact - gest)))
    vects = rbind(vects,c(rot_line))
    i = i +1
  }
  
  return(list(errors=errors,vects=vects))
}
plot_rotated_track <- function(y,s,smart,default){
  library(ggplot2)
  library(viridis)
  
  df.merged <- rbind(data.frame(x = c(1:length(y)), y = y))
  p = ggplot(NULL) +geom_line(data = df.merged,aes(x = x, y = y)) + scale_size_manual(values = c(1,2)) +
    theme(legend.position="bottom", legend.direction="horizontal",
          legend.title = element_blank(),panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white")) +
    labs(y="Mean Absolute Error", x = "Directions") + 
    ggtitle(s)
  p = p + geom_hline(yintercept=smart, linetype="dashed", color = "red") +
      geom_hline(yintercept=default, linetype="dashed", color = "green")
  p = p + theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "grey"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "grey"))
  
  return(p)
}

Rosenbrock_Banana <- function(x) {  
  res <- 0.0
  for(i in 1:(length(x)-1))
    res <- res + 100*(x[i+1] - x[i]^2)^2 + (1-x[i])^2
  return(res)
}

gr_Rosenbrock_Banana <- function(x){
  n <- length(x)
  g <- numeric(n)
  for(i in 1:(n-1)) {
    g[i] <- g[i] -400 * x[i] * (x[i+1] - x[i]^2) - 2 * (1 - x[i])
    g[i+1] <- g[i+1] + 200 * (x[i+1] - x[i]^2)
  }
  return(g)
}

Extended_Trigonometric_dim5 <- function(x){
  res <- 0.0
  n <- length(x)
  for(i in 1:n) {
    c <- (n - sum(cos(x))) + i*(1 - cos(x[i])) - sin(x[i])
    res = res + (c*c)
  }
  return(res)
}

gr_Extended_Trigonometric_dim5 <- function(x){
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

Extended_Freudenstein_Roth <- function(x){
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

gr_Extended_Freudenstein_Roth <- function(x){
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

Perturbed_Quadratic <- function(x){
  res <- 0.0
  n <-length(x)
  c1 <- 1:n
  c2 <- x^2
  return(sum(c1*c2) + (1/100)*sum(x)^2)
}

gr_Perturbed_Quadratic <- function(x){
  gr <- numeric(length(x))
  n <-length(x)
  c <- sum(x)
  for(i in 1:n) gr[i] <- 2*i*x[i] + (1/50)*c
  return(gr)
}

Raydan_1 <- function(x){
  n <-length(x)
  return(sum(((1:n)/10)*(exp(x) -x)))
}

gr_Raydan_1 <- function(x){
  n <-length(x)
  gr <- numeric(n)
  for(i in 1:n) gr[i] <- (i/10)*(-1 + exp(x[i]))
  return(gr)
}

Raydan_2 <- function(x){
  n <-length(x)
  return(sum((exp(x) -x)))
}

gr_Raydan_2 <- function(x){
  n <-length(x)
  gr <- numeric(n)
  for(i in 1:n) gr[i] <- (-1 + exp(x[i]))
  return(gr)
}

FLETCHCR <- function(x){
  res <- 0.0
  n <-length(x)
  for(i in 1:(n-1)) res <- res + 100*(x[i+1] - x[i] + 1 - x[i]^2)^2
  return(res)
}

gr_FLETCHCR <- function(x){
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

COSINE <- function(x){
  res <- 0.0
  n <-length(x)
  for(i in 1:(n-1)) res <- res + cos(-0.5*x[i+1] + x[i]^2)
  return(res)
}

gr_COSINE <- function(x){
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

Generalized_Quartic <- function(x){
  res <- 0.0
  n <-length(x)
  for(i in 1:(n-1)) res <- res + (x[i]^2 + (x[i+1] + x[i]^2)^2)
  return(res)
}

gr_Generalized_Quartic <- function(x){
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

default.gr <- function(fn,x){
  s = 1e-3
  res = numeric(length(x))
  for(i in 1:length(x))
  {
    e = numeric(length(x))
    e[i] = 1
    res[i] = (fn(x+s*e) - fn(x-s*e))/(2*s)
  }
  return(res)
}

Schumer_Steiglitz <- function(x){
  res <- 0.0
  n <-length(x)
  return(sum(x[1:n]^4))
}

gr_Schumer_Steiglitz <- function(x){
  n <-length(x)
  gr <- numeric(n)
  gr[1:n] = 4*x[1:n]^3
  return(gr)
  
}


ackley2 <- function(x){
  return(-200*exp(-0.02*sqrt(x[1]^2+ x[2]^2)))
}

gr_ackley2 <- function(x){
  gr <- numeric(2)
  c1 <- sqrt(x[1]^2+ x[2]^2)
  c2 <- exp(-0.02*c1)
  gr[1] = 4*x[1]*c2/c1
  gr[2] = 4*x[2]*c2/c1
  
  return(gr)
}

ackley3 <- function(x){
  return(200*exp(-0.02*sqrt(x[1]^2 + x[2]^2)) + 5*exp(cos(3*x[1])+sin(3*x[2])))
}

gr_ackley3 <- function(x){
  gr <- numeric(2)
  c1 <- sqrt(x[1]^2+ x[2]^2)
  c2 <- -4*exp(-0.02*c1)/c1
  c3 <- 15*exp(cos(3*x[1]) + sin(3*x[2]))
  gr[1] = -c3*sin(3*x[1]) + x[1]*c2
  gr[2] = c3*cos(3*x[2]) + x[2]*c2
  
  return(gr)
}

#### ------> List of Functions
# Rosenbrock_Banana
# Extended_Trigonometric_dim5
# Extended_Freudenstein_Roth
# Perturbed_Quadratic
# Raydan_1
# Raydan_2
# FLETCHCR
# COSINE
# Generalized_Quartic
# Schumer_Steiglitz 
# ackley2
# ackley3

###Comparing exact gradient and numerical gradients

fn <- ackley3  
gr <- gr_ackley3  #gr_"functionName"
x = 1:2
default.gr(fn,x)              
numDeriv::grad(fn,x)
gr(x)

###Testing Plot Function
plot_MSE(c(1:100),rnorm(100),rnorm(100),FALSE)
