
# generate random matrix of random size
n <- sample(2:10, size = 1)
p <- sample(2:10, size = 1)
A <- matrix(runif(n), nrow=n, ncol=p) 

# using R function
R_svd <- svd(A)
R_svd

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

U <- extract(SVD_stan, pars = "U")
V <- extract(SVD_stan, pars = "V")
D <- extract(SVD_stan, pars = "D")

U <- U$U[, , 1:nrow(A)]
V <- V$V[, , 1:ncol(A)]
D <- D$D[, , 1:ncol(A)]

Stan_A <- U %*% D %*% t(V)

Stan_A
A
