library(data.table)
library(dlm)
library(rstan)

y <- matrix(0, nrow = length(Nile), ncol = 2)
y[, 1] <- y[, 2] <- as.vector(Nile)

rand_rows_case1 <- sample(nrow(y), size = 3)
y[rand_rows_case1, 1:ncol(y)] <- NA

rand_rows_case2 <- sample(nrow(y), size = 2)
rand_cols_case2 <- sample(ncol(y), size = 2)
y[rand_rows_case2[1], rand_cols_case2[1]] <- NA
y[rand_rows_case2[2], rand_cols_case2[2]] <- NA

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

NileFilt <- dlmFilter( 
    y,
    mod
)

dat[,
    ':='(
        mod1_filtered = NileFilt$m[-1, 1],
        mod2_filtered = NileFilt$m[-1, 2]
    )
    ]

# Model 1 in Stan
y_missing <- matrix(0, nrow = nrow(y), ncol = ncol(y))

y[rand_rows_case1, 1:ncol(y)] <- Inf
y_missing[rand_rows_case1, 1:ncol(y)] <- 1

y[rand_rows_case2[1], rand_cols_case2[1]] <- Inf
y[rand_rows_case2[2], rand_cols_case2[2]] <- Inf
y_missing[rand_rows_case2[1], rand_cols_case2[1]] <- 1
y_missing[rand_rows_case2[2], rand_cols_case2[2]] <- 1
num_missing <- rowSums(y_missing)

stan_args <- list( 
    N = nrow(y),    # Number of observations
    M = ncol(y),            # Dimension of observation
    P = nrow(mod$GG),            # Dimension of state vector
    
    # Observation Stuff
    y = y,         # Observations
    y_missing = y_missing,
    num_missing = num_missing,
    FF = mod$FF,      # Observation equation matrix
    V = mod$V,        # Known constant covariance matix for observation equation
    
    # State Stuff
    GG = mod$GG,       # State/transition equation matrix
    W = mod$W,         # Known constant covariance matix for state transition equation
    
    m0 = as.array(mod$m0),       # Initial state vector prior mean
    C0 = mod$C0        # Initial state vector prior variance-covariance matix
)

mod1_Stan <- stan( 
    file = "Stan/dlmFilter_missingvalues.stan",
    data = stan_args,
    chains = 1,
    iter = 1,
    warmup = 0,
    algorithm = "Fixed_param",
    cores = 1
)

stan_filtered <- extract(mod1_Stan, pars = "m")
dat[,
    ':='(
        mod1_stan_filtered = as.vector(stan_filtered$m[, -1, 1]),
        mod2_stan_filtered = as.vector(stan_filtered$m[, -1, 2])
    )
    ]
