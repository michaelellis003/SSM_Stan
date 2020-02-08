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

# Smooth using functions in dlm package
mod1 <- dlmModPoly( 
    order = 1,
    dV = 15100,
    dW = 755
)

mod2 <- dlmModPoly( 
    order = 1,
    dV = 15100,
    dW = 7550
)

NileSmooth <- dlmSmooth( 
    Nile,
    mod1
)

hwid_mod1 <- qnorm(0.025, lower = FALSE) * sqrt(unlist(dlmSvd2var(NileSmooth$U.S, NileSmooth$D.S)))

dat[, ':='( 
    mod1_smooth = NileSmooth$s[-1],
    mod1_smooth_upper = NileSmooth$s[-1] + hwid_mod1[-1],
    mod1_smooth_lower = NileSmooth$s[-1] - hwid_mod1[-1]
    )]

NileSmooth <- dlmSmooth( 
    Nile,
    mod2
)

hwid_mod2 <- qnorm(0.025, lower = FALSE) * sqrt(unlist(dlmSvd2var(NileSmooth$U.S, NileSmooth$D.S)))

dat[, ':='( 
    mod2_smooth = NileSmooth$s[-1],
    mod2_smooth_upper = NileSmooth$s[-1] + hwid_mod2[-1],
    mod2_smooth_lower = NileSmooth$s[-1] - hwid_mod2[-1]
)]

ggplot(dat) + 
    geom_point(aes(x = Year, y = y)) + 
    geom_line(aes(x = Year, y = mod1_smooth, color = "Model 1")) + 
    geom_ribbon(aes(x = Year, ymin = mod1_smooth_lower, ymax = mod1_smooth_upper, fill = "Model 1"), alpha = 0.25) + 
    geom_line(aes(x = Year, y = mod2_smooth, color = "Model 2")) + 
    geom_ribbon(aes(x = Year, ymin = mod2_smooth_lower, ymax = mod2_smooth_upper, fill = "Model 2"), alpha = 0.25) +
    theme_minimal() + 
    ylab("Nile")

# Smooth using functions in Stan
# Model 1 in Stan
stan_args_mod1 <- list( 
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
    file = "Stan/dlmSmooth.stan",
    data = stan_args_mod1,
    chains = 1,
    iter = 1,
    warmup = 0,
    algorithm = "Fixed_param",
    cores = 1
)

dat[, mod1_stan_smoothed := as.vector(unlist(extract(mod1_Stan, pars = "s")))[-1]]

# Model 2 in Stan
stan_args_mod2 <- list( 
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
    file = "Stan/dlmSmooth.Stan",
    data = stan_args_mod2,
    chains = 1,
    iter = 1,
    warmup = 0,
    algorithm = "Fixed_param",
    cores = 1
)

dat[, mod2_stan_smoothed := as.vector(unlist(extract(mod2_Stan, pars = "s")))[-1]]

ggplot(dat) + 
    geom_point(aes(x = Year, y = y)) + 
    geom_line(aes(x = Year, y = mod1_smooth, color = "Model 1 - DLM")) + 
    geom_ribbon(aes(x = Year, ymin = mod1_smooth_lower, ymax = mod1_smooth_upper, fill = "Model 1 - DLM"), alpha = 0.25) + 
    geom_line(aes(x = Year, y = mod2_smooth, color = "Model 2 - DLM")) + 
    geom_ribbon(aes(x = Year, ymin = mod2_smooth_lower, ymax = mod2_smooth_upper, fill = "Model 2 - DLM"), alpha = 0.25) +
    geom_line(aes(x = Year, y = mod1_stan_smoothed, color = "Model 1 - Stan"), linetype="dashed") + 
    geom_line(aes(x = Year, y = mod2_stan_smoothed, color = "Model 2 - Stan"), linetype="dashed") + 
    theme_minimal() + 
    ylab("Nile")




















