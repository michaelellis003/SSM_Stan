library(data.table)
library(dlm)
library(rstan)

y <- as.vector(Nile)

dat <- data.table( 
    Nile = y,
    Year = 1871:1970
)

ggplot(dat, aes(x = Year, y = Nile)) + 
    geom_point() + 
    geom_line() + 
    theme_minimal() + 
    ylab("Nile")

# Filter using functions in dlm package
mod <- dlmModPoly(order = 1)
dlm_fit <- dlmGibbsDIG(y, mod, 
                       shape.y = 0.01, rate.y = 0.01, 
                       shape.theta = 0.01, rate.theta = 0.01, 
                       n.sample = 5000)

# Model 1 in Stan
stan_args <- list( 
    N = length(y),             # Number of observations
    M = 1,             # Dimension of observations
    P = nrow(mod$GG),        # Dimension of state vector
    
    # Observation Stuff
    y = matrix(y, nrow = length(y), ncol = 1),                   # Observations
    y_missing = matrix(0, nrow = length(y), ncol = 1),    
    num_missing = rep(0, length(y)),
    FF = mod$FF,             # Observation equation matrix
    GG = mod$GG,             # State/transition equation matrix
    m0 = as.array(mod$m0),   # Initial state vector prior mean
    C0 = mod$C0              # Initial state vector prior variance-covariance matix
)

mod1_Stan <- stan( 
    file = "Stan/stan_dlm.stan",
    data = stan_args,
    chains = 2,
    iter = 5000,
    warmup = 2500,
    cores = 2,
    control = list( 
        adapt_delta = 0.8,
        max_treedepth = 10
        )
)

stan_filtered <- extract(mod1_Stan, pars = "m")
dat[,
    ':='(
        mod1_stan_filtered = as.vector(stan_filtered$m[, -1, 1]),
        mod2_stan_filtered = as.vector(stan_filtered$m[, -1, 2])
    )
    ]
