args <- commandArgs(trailingOnly = TRUE)
grid <- args[1] # e.g. "grid.213.7.Rdata"
env.name <- args[2] # e.g. "env.scale.2.sd.PC.Rdata"
columns <- args[3] # Arguments have type character, so coerce to numeric
# see https://wiki.csiro.au/display/ASC/Run+R+script+as+a+batch+job+on+a+Linux+cluster
len.ind <- as.numeric(args[4]) # length index: 0 if env. data has all scales; 1 if the smallest scale is missing 
iters <- as.numeric(args[5]) # markov iterations, e.g. 25000
outpath <- args[6]

### for an example:
## DO NOT RUN:
# 1) grid <- "grid.213.7.Rdata"
# 2) env.name <- "env.scale.2.sd.PC.Rdata"
# 3) columns <- "-c(26,29,30)"
# 4) len.ind <- as.numeric("1")
# 5) iters <- as.numeric("10")

# s <- 5
####

library(vegan)
library(rosalia)
library(mistnet)
library(doParallel)
# Find out how many cores are available (if you don't already know)
detectCores()
# Create cluster with desired number of cores
cl <- makeForkCluster() # !!! set to 2, since the first process (20,000 sites) takes ages on one core, while the other can finish on the second core 
# Register cluster
registerDoParallel(cl)
# Find out how many cores are being used
getDoParWorkers()

########
# # outside the loop
`%plus%` = mistnet:::`%plus%`
logistic = binomial()$linkinv # logistic inverse link
rbern = function(p){rbinom(length(p), size = 1, prob = p)}

# load("env.scale.2.sd.PC.Rdata") ### !!! change to "mean"
env2 <- get(load(env.name))
# load("grid.213.7.Rdata") # !!! change to "size"
grid2 <- get(load(grid))

######################
out <- foreach (s = 1:length(env2), .packages=c("rosalia", "mistnet","vegan")) %dopar% { 
# !!! change to "1:7" if means are used; 1:6 if sd is used
y <- as.matrix(grid2[[s+len.ind]]) # !!! for "sd" change to [[s+1]]...
# y <- as.matrix(grid.213.7[[s]])                                     # !!! change to just "s", if means instead of sd are used
# transform into presence-absence matrix:
y <- decostand(y, "pa")

###################################################################
####### to be changed depending on which env. data should be used
env = true_env <- x <- soil.13 <- scale(env2[[s]][,eval(parse(text=columns))]) #!!!! change to the right variable set
    # 1 pH
    # 2 elevation
    # 3 shade
    # c(1,2) pH & elevation 
    # 4 PC1
    # c(4,5) PC1 & PC2
###################################################################

set.seed(1)
n_spp = ncol(y)   # number of species
n_loc = nrow(y)    # number of locations
n_env = ncol(x)     # number of environmental predictors

n_gibbs = 5000  # number of Gibbs sampling iterations
######################

y_stats = crossprod(y)
y_env_stats = t(true_env) %*% y
y_sim = matrix(0.5, nrow = nrow(y), ncol = ncol(y))
env = true_env
alpha_env = delta_alpha_env = matrix(0, nrow = n_env, ncol = n_spp)
alpha_species = qlogis((colSums(y) + 1) / (nrow(y) + 2))
delta_alpha_species = rep(0, n_spp)
beta = delta_beta = matrix(0, nrow = n_spp, ncol = n_spp)
alpha = matrix(0, nrow = n_spp, ncol = n_spp) 

alpha_env_prior = rosalia::make_logistic_prior(scale = 2)$log_grad 
alpha_species_prior = rosalia::make_logistic_prior(scale = 2)$log_grad
beta_prior = rosalia::make_logistic_prior(scale = 0.5)$log_grad

initial_learning_rate = 1 # step size at start of optimization

##############################################
### needs to be changed
maxit = iters        
                                        # !!! 25000, Number of rounds of optimization # try 50000 and 100000 as well
########################################
start_time = as.integer(Sys.time())
times = integer(maxit)
        for(i in 1:maxit){
  ##############################
  # Gibbs sampling for predicted species composition
  ##############################

  # Update alpha
  alpha = env %*% alpha_env %plus% alpha_species

  # Sample entries in y_sim from their conditional 
  # distribution (Gibbs sampling)
    for(j in sample.int(n_spp)){
    y_sim[,j] = rbern(logistic(alpha[ , j] + y_sim %*% beta[ , j]))
  }

  ##############################
  # Stochastic approximation for updating alpha and beta
  ##############################

  # Update learning rate and momentum
 
    learning_rate = initial_learning_rate * 1000 / (998 + 1 + i)
 momentum = .9 * (1 - 1/(.1 * i + 2))

  # Calculate sufficient statistics
  y_sim_stats = crossprod(y_sim)
  y_sim_env_stats = t(env) %*% y_sim


  stats_difference = y_stats - y_sim_stats
  beta_grad = (stats_difference + beta_prior(beta)) / n_loc

  alpha_species_grad = (
    diag(stats_difference) + 
      alpha_species_prior(alpha_species)
  ) / n_loc
  diag(beta_grad) = 0 # beta_ii is 0 by convention
  y_env_difference = y_env_stats - y_sim_env_stats
  alpha_env_grad = (y_env_difference + 
                      alpha_env_prior(alpha_env))  / n_loc


  delta_beta = beta_grad * learning_rate + momentum * delta_beta
  delta_alpha_species = alpha_species_grad * learning_rate + 
    momentum  * delta_alpha_species
  delta_alpha_env = alpha_env_grad * learning_rate +
    momentum  * delta_alpha_env

  beta = beta + delta_beta
  alpha_species = alpha_species + delta_alpha_species
  alpha_env = alpha_env + delta_alpha_env
  times[i] = as.integer(Sys.time()) - start_time
}
        list("alpha_species"=alpha_species, "alpha_env"=alpha_env, "beta"=beta)

}

##  give name here
#save(out, file = paste(grid,env.name,columns,iters,"Rdata", sep = "."))
save(out, file = outpath)
