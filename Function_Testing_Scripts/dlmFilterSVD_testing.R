library(data.table)
library(dlm)
library(rstan)

dat <- data.table( 
    y = as.vector(Nile),
    Year = 1871:1970
)

ggplot(dat, aes(x = Year, y = y)) + 
    geom_point() + 
    geom_line() + 
    theme_minimal() + 
    ylab("Nile")

# Model 1 filter using functions in dlm package
mod1 <- dlmModPoly( 
    order = 1,
    dV = 15100,
    dW = 755
)

NileFilt1 <- dlmFilter( 
    Nile,
    mod1
)

dat[, mod1_filtered := NileFilt1$m[-1]]

# Model 1 in Stan
stan_args <- list( 
    N = nrow(dat),    # Number of observations
    M = 1,            # Dimension of observation
    P = 1,            # Dimension of state vector
    
    # Observation Stuff
    y = matrix(dat$y, nrow = nrow(dat), ncol = 1),         # Observations
    FF = mod1$FF,      # Observation equation matrix
    V = mod1$V,        # Known constant covariance matix for observation equation
    
    # State Stuff
    GG = mod1$GG,       # State/transition equation matrix
    W = mod1$W,         # Known constant covariance matix for state transition equation
    
    m0 = as.array(mod1$m0),       # Initial state vector prior mean
    C0 = mod1$C0        # Initial state vector prior variance-covariance matix
)

mod1_Stan <- stan( 
    file = "Stan/dlmFilter.stan",
    data = stan_args,
    chains = 1,
    iter = 1,
    warmup = 0,
    algorithm = "Fixed_param",
    cores = 1
)

dat[, mod1_stan_filtered := as.vector(unlist(extract(mod1_Stan, pars = "m")))[-1]]

## compare parameter estimates
as.vector(unlist(extract(mod1_Stan, pars = "m")))
as.vector(unlist(NileFilt1$m))

as.vector(unlist(extract(mod1_Stan, pars = "U_C")))
as.vector(unlist(NileFilt1$U.C))

sqrt(as.vector(unlist(extract(mod1_Stan, pars = "D_C"))))
as.vector(unlist(NileFilt1$D.C))

sqrt(as.vector(unlist(extract(mod1_Stan, pars = "U_R"))))
as.vector(unlist(NileFilt1$U.R))

sqrt(as.vector(unlist(extract(mod1_Stan, pars = "D_R"))))
as.vector(unlist(NileFilt1$D.R))

as.vector(unlist(extract(mod1_Stan, pars = "a")))
as.vector(unlist(NileFilt1$a))


# Model 2 in dlm
mod2 <- dlmModPoly( 
    order = 1,
    dV = 15100,
    dW = 7550
)

NileFilt2 <- dlmFilter( 
    Nile,
    mod2
)

dat[, mod2_filtered := NileFilt2$m[-1]]

# Model 2 in Stan
stan_args <- list( 
    N = nrow(dat),    # Number of observations
    M = 1,            # Dimension of observation
    P = 1,            # Dimension of state vector
    
    # Observation Stuff
    y = matrix(dat$y, nrow = nrow(dat), ncol = 1),         # Observations
    FF = mod2$FF,      # Observation equation matrix
    V = mod2$V,        # Known constant covariance matix for observation equation
    
    # State Stuff
    GG = mod2$GG,       # State/transition equation matrix
    W = mod2$W,         # Known constant covariance matix for state transition equation
    
    m0 = as.array(mod2$m0),       # Initial state vector prior mean
    C0 = mod2$C0        # Initial state vector prior variance-covariance matix
)


mod2_Stan <- stan( 
    file = "Stan/dlmFilter.stan",
    data = stan_args,
    chains = 1,
    iter = 1,
    warmup = 0,
    algorithm = "Fixed_param",
    cores = 1
)

dat[, mod2_stan_filtered := as.vector(unlist(extract(mod2_Stan, pars = "m")))[-1]]

ggplot(dat) + 
    geom_point(aes(x = Year, y = y)) + 
    geom_line(aes(x = Year, y = mod1_filtered, color = "Model 1 - DLM")) + 
    geom_line(aes(x = Year, y = mod2_filtered, color = "Model 2 - DLM")) + 
    geom_line(aes(x = Year, y = mod1_stan_filtered, color = "Model 1 - Stan"), linetype="dashed") + 
    geom_line(aes(x = Year, y = mod2_stan_filtered, color = "Model 2 - Stan"), linetype="dashed") + 
    theme_minimal() + 
    ylab("Nile")
