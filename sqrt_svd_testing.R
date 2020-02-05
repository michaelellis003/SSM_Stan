library(rstan)

# generate random matrix
n <- sample(1:5, size = 1)
A <- matrix(runif(n^2)*2-1, ncol=n) 
A <- t(A) %*% A

# using Stan function
stan_args <- list( 
    p = nrow(A),
    q = ncol(A),
    M = A,
    calc_sqrt = 1
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

sqrt_A <- extract(SVD_stan, pars = "sqrt_M")
sqrt_A <- as.matrix(sqrt_A$sqrt_M[, , 1:n], nrow = n, ncol = n)
Stan_A <- t(sqrt_A) %*% sqrt_A
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
