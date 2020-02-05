library(rstan)

# generate random matrix
n <- sample(1:5, size = 1)
p <- sample(1:5, size = 1)
A <-  matrix(runif(n), nrow=n, ncol=p) 
A

# using R function
R_svd <- svd(A, nu = n, nv = p)
#R_svd

# using Stan function
stan_args <- list( 
    p = nrow(A),
    q = ncol(A),
    M = A
    )

SVD_stan <- stan( 
    file = "Stan/SVD.stan",
    data = stan_args,
    chains = 1,
    iter = 1,
    warmup = 0,
    algorithm = "Fixed_param",
    cores = 1
)

Stan_A <- extract(SVD_stan, pars = "M_svd")
Stan_A <- as.matrix(Stan_A$M_svd[, , 1:p], nrow = n, ncol = p)
if(n == 1){
    Stan_A <- t(Stan_A)
}
Stan_A
A

num_digits <- 10
matrices_equal <- rep(FALSE, 10)
for(i in 1:num_digits){
    matrices_equal[i] <- all(round(Stan_A, digits = i) == round(A, digits = i))
    if(!matrices_equal[i] | i == 10){
        message("Matrices equal up to ", i, " digits")
        break
    }
}


U <- extract(SVD_stan, pars = "U")
V <- extract(SVD_stan, pars = "V")
D <- extract(SVD_stan, pars = "D")

U <- as.matrix(U$U[, , 1:n], nrow = n, ncol = p)
V <- as.matrix(V$V[, , 1:p], nrow = p, ncol = min(n, p))
D <- diag(as.vector(D$D[, , 1:p]), nrow = n, ncol = p)
D
U
V
