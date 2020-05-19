#### Simulação e ajuste do modelo dinamico para um sarima ####

### pacotes ###
library(sn)

set.seed(999)
x <- rep(0, 1200)
x[1] <- rnorm(1)
for(i in 2:1200){
  x[i] <- 0.7*x[i-1] + 0.3 * ifelse(i>12, x[i-12], 0) + rnorm(1,0,2)
}

ts.plot(x)

N <- 1000#length(x)

y <- x[1:N]

# modelo
p <- 12
w <- 2*pi/p

Gt <- matrix(c(1,1,0,0,
               0,1,0,0,
               0,0,cos(w),sin(w),
               0,0,-sin(w),cos(w)),4,4, byrow = T)

Ft <- c(1,0,1,0)

n <- length(Ft)

m0 <- rep(0, n)
C0ast <- diag(100, n, n)
d0 <- 0.01
n0 <- 0.01
S0 <- d0/n0
a <- matrix(NA, length(y), n) 
m <- matrix(NA, length(y), n) 

f <- c() 
e <- c() 
Cast <- list() -> C
Rast <- list() -> R
Qast <- list() -> Q
A <- list()
d <- c()
n <- c()
S <- c()
V <- c()

delta = 0.98

a[1, ]    <- Gt %*% m0
Rast[[1]] <- Gt %*% C0ast %*% t(Gt) / delta
R[[1]]    <- S0 * Rast[[1]]
f[1]      <- t(Ft) %*% a[1, ]
Qast[[1]] <- t(Ft) %*% Rast[[1]] %*% Ft + 1
Q[[1]]    <- S0*Qast[[1]]
e[1]      <- y[1] - f[1]
n[1]      <- n0 + 1
d[1]      <- d0 + e[1]^2*solve(Qast[[1]])
S[1]      <- d[1]/n[1]
A[[1]]    <- Rast[[1]] %*% Ft %*% solve(Qast[[1]])
m[1, ]    <- a[1, ] + A[[1]] %*% e[1]
Cast[[1]] <- Rast[[1]] - as.vector(Qast[[1]]) * (A[[1]] %*% t(A[[1]]))
C[[1]]    <- S[1]*Cast[[1]]

cred_f <- matrix(NA, length(y), 2)
cred_nivel <- matrix(NA, length(y), 2)
cred_beta <- matrix(NA, length(y), 2)
cred_saz1 <- matrix(NA, length(y), 2)
cred_saz2 <- matrix(NA, length(y), 2)

for(t in 2:length(y)){
  a[t, ]      <- Gt %*% m[t-1, ]
  Rast[[t]]   <- Gt %*% Cast[[t-1]] %*% t(Gt) / delta
  R[[t]]      <- S[t-1] * Rast[[t]]
  f[t]        <- t(Ft) %*% a[t, ]
  Qast[[t]]   <- t(Ft) %*% Rast[[t]] %*% Ft + 1
  Q[[t]]      <- S[t-1]*Qast[[t]]
  e[t]        <- y[t] - f[t]
  n[t]        <- n[t-1] + 1
  d[t]        <- d[t-1] + e[t]^2*solve(Qast[[t]])
  S[t]        <- d[t]/n[t]
  A[[t]]      <- Rast[[t]] %*% Ft %*% solve(Qast[[t]])
  m[t, ]      <- a[t, ] + A[[t]] %*% e[t]
  Cast[[t]]   <- Rast[[t]] - as.vector(Qast[[t]]) * (A[[t]] %*% t(A[[t]]))
  C[[t]]      <- S[t]*Cast[[t]]
  cred_f[t, ] <- qst(c(0.025, 0.975), xi = f[t], omega = as.numeric(sqrt(Qast[[t]])),alpha=0,nu=n[t-1])
}

N <- length(y)
B = list()
m_s <- matrix(NA, N, length(Ft)) 
C_s <- list()
m_s[N, ] <- m[N, ]
C_s[[N]] <- C[[N]]

for(t in (N-1):1){
  B[[t]] <- C[[t]] %*% t(Gt) %*% solve(R[[t+1]])
  m_s[t, ] <- m[t, ] + B[[t]] %*% (m_s[t+1, ] - a[t+1, ])
  C_s[[t]] <- C[[t]] - B[[t]] %*% (R[[t+1]] - C_s[[t+1]]) %*% t(B[[t]]) 
  cred_nivel[t, ] <- qst(c(0.025, 0.975), 
                         xi = m_s[t,1], 
                         omega = as.numeric(sqrt(C_s[[t]][1,1])),
                         alpha=0,nu=n[t])
  cred_beta[t, ] <- qst(c(0.025, 0.975), 
                        xi = m_s[t,2], 
                        omega = as.numeric(sqrt(C_s[[t]][2,2])),
                        alpha=0,nu=n[t])  
  cred_saz1[t, ] <- qst(c(0.025, 0.975), 
                        xi = m_s[t,3], 
                        omega = as.numeric(sqrt(C_s[[t]][3,3])),
                        alpha=0,nu=n[t])
  cred_saz2[t, ] <- qst(c(0.025, 0.975), 
                        xi = m_s[t,4], 
                        omega = as.numeric(sqrt(C_s[[t]][4,4])),
                        alpha=0,nu=n[t])
}


plot(y, type = "l")
lines(f, col = 2, lty = 2)


eqm_ajuste = mean((y-f)^2)
eam_ajuste = mean(abs(y-f))
mape_ajuste = mean(abs((y - f)/y))*100 



h=150

a_prev = matrix(NA, h, length(Ft))
f_prev = vector()
R_prev = list()
Q_prev = vector()
theta_prev = matrix(NA, h, length(Ft))

a_prev[1,] = m[N,]
R_prev[[1]] = C[[N]]
f_prev[1] <- t(Ft) %*% a_prev[1, ]
Q_prev[1] <- t(Ft) %*% R_prev[[1]] %*% Ft + 1

cred_f_prev <- matrix(NA,length(y),2)
cred_f_prev[1, ] <- qst(c(0.025, 0.975),
                        xi = f_prev[1],
                        omega = as.numeric(sqrt(Q_prev[1])),
                        alpha = 0, nu = n[t])

for (t in 2:h){
  a_prev[t, ] <- Gt %*% a_prev[t-1, ]
  R_prev[[t]] <- Gt %*% R_prev[[t-1]] %*% t(Gt) / delta  
  f_prev[t] <- t(Ft) %*% a_prev[t, ]
  Q_prev[t] <- t(Ft) %*% R_prev[[t]] %*% Ft + 1
  
  cred_f_prev[t, ] <- qst(c(0.025, 0.975),
                          xi = f_prev[t],
                          omega = as.numeric(sqrt(Q_prev[t])),
                          alpha = 0, nu = n[t])
}


y_prev = x[1001:1150]

mse = mean((y_prev-f_prev)^2)
mae = mean(abs(y_prev-f_prev))
mape = mean(abs((y_prev - f_prev)/y_prev))*100



# library(keras)
# model = 



