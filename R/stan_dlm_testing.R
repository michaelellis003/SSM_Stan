library(data.table)
library(dlm)
library(rstan)

y <- matrix(0, nrow = length(Nile), ncol = 2)
y[, 1] <- y[, 2] <- as.vector(Nile)

dat <- data.table( 
    y1 = y[, 1],
    y2 = y[, 2],
    Year = 1871:1970
)

plot_dat <- melt(
    dat,
    id.vars = "Year",
    variable.name = "y"
)

ggplot(plot_dat, aes(x = Year, y = value, col = y)) + 
    geom_point() + 
    geom_line() + 
    facet_wrap(~y, ncol = 1) +
    theme_minimal() + 
    ylab("Nile")

# Filter using functions in dlm package
mod <- dlmModPoly(order = 1)
mod$FF <- mod$FF %x% diag(2)
mod$GG <- mod$GG %x% diag(2)
mod$V <- bdiag(15100, 15100)
mod$W <-bdiag(755, 7550)
mod$m0 <- rep(0, nrow(mod$GG))
mod$C0 <- diag(1e7, nrow(mod$GG))

# Model 1 in Stan
stan_args <- list( 
    N = nrow(y),             # Number of observations
    M = ncol(y),             # Dimension of observations
    P = nrow(mod$GG),        # Dimension of state vector
    
    # Observation Stuff
    y = y,                   # Observations
    y_missing = matrix(0, nrow = nrow(dat), ncol = ncol(y)),    
    num_missing = rep(0, nrow(dat)),
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
        adapt_delta = 0.9,
        max_treedepth = 15
        )
)

stan_filtered <- extract(mod1_Stan, pars = "m")
dat[,
    ':='(
        mod1_stan_filtered = as.vector(stan_filtered$m[, -1, 1]),
        mod2_stan_filtered = as.vector(stan_filtered$m[, -1, 2])
    )
    ]
