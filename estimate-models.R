source("functions.R")
install.packages(c("plyr", "geoR", "mnormt", "MHadaptive"))
daten.model <- filter.data(data)

X <- 40*h**daten.model$x
S <- 40*h**daten.model$s
n <- daten.model$n

f <- 1 -daten.model[,"stopprob"]

modes.EU <- Mode(c(1,5,0.9), n, f, X, S, prior = "yes", hyperpars = c(0.7,0.3))

modes.Regret <- Mode(c(1,0.1,5,0.9), n, f, X, S, prior = "yes", hyperpars = c(0.7,0.3))

# As standard: Take inverse Hessian as covariance matrix of proposal density
# Caution: this neglects restrictions on sigma and kappa. This does not cause
# trouble in the Markov Chain algorithm, but makes the inverse hessian
# difficult to interpret as covariance matrix.
# Hence, simulated quantiles from the Markov chain of the parameters
# should be taken for inference.
HESS.EU <- solve(-modes.EU$hessian)
HESS.Regret <- solve(-modes.Regret$hessian)

# Draw dispersed starting values to run multiple chains form different
# points to assess robustness of convergence.
set.seed(196)
starters.EU <- rbind(modes.EU$par, rmnorm(n = 3, mean = modes.EU$par, varcov = 5 * HESS.EU))
starters.Regret <- rbind(modes.Regret$par, rmnorm(n = 3, mean = modes.Regret$par, varcov = 5 * HESS.Regret))

metros.EU <- lapply(1:nrow(starters.EU), function(x) MH(starters.EU[x, ], HESS.EU, 5500, 0, n, f, X, S, prior = "yes", mode= modes.EU))

# Run the chains for Regret
metros.Regret <-  lapply(1:nrow(starters.Regret), function(x) MH(starters.Regret[x, ], HESS.Regret, 5500, 0, n, f, X, S, prior = "yes", mode= modes.Regret))

######## End of pooled models - per-subject models follow ##########

library(parallel)
cores <- detectCores()-1
clust <- makeCluster(cores)
# Inherit the workspace to each node
work.env <- c(setdiff(ls(), c("metros.EU.sub", "metros.Regret.sub")))
#work.env <- ls()
clusterExport(clust, work.env)
clusterEvalQ(clust, {library(plyr)
                     library(geoR)
                     library(MHadaptive)})
# Broadcast random seeds for child processes of R. Makes results reproducible.
clusterSetRNGStream(cl = clust, iseed = 190684)

# Find the modes for per-subject models 
EU.modes.wider.prior<- clusterApplyLB(cl=clust, x=1:14, fun = estimate.model, start.val = c(1,5,0.9), round=1, prior = "yes", hyperpars=c(0.7,1))
Regret.modes.wide.prior <- clusterApplyLB(cl=clust, x=1:44, fun = estimate.model, start.val = c(1,0.1,5,0.9), round=1, prior = "yes", hyperpars=c(0.7,1))

stopCluster(clust)
rm(clust)

#EU.modes <- lapply(1:44, estimate.model, start.val = c(1,5,0.9), prior = "yes", hyperpars=c(0.7,0.3))
#Regret.modes <- lapply(1:44, estimate.model, start.val = c(1,0.1,5,0.9), round = 1, prior = "yes", hyperpars=c(0.7,0.3))

posterior.modes.EU.wider <- unlist(sapply(1:14, function(x) EU.modes.wider.prior[[x]]$par))
posterior.modes.Regret.wide <- unlist(sapply(1:44, function(x) Regret.modes.wide.prior[[x]]$par))

posterior.hess.Regret <- lapply(1:44, function(x) solve(-Regret.modes[[x]]$hessian))
posterior.hess.EU <- lapply(1:44, function(x) solve(-EU.modes[[x]]$hessian))

# Here we can run the MCMC algorithm in parallel
clust <- makeCluster(cores)
# Inherit the workspace to each node
work.env <- c(setdiff(ls(), c("metros.EU", "metros.Regret", "metros.EU.sub", "metros.Regret.sub", 
                              "simulated.decisions.EU", 
                              "simulated.decisions.Regret",
                              "stopping.frequencies.EU", 
                              "stopping.frequencies.Regret",
                              "EUs",
                              "Regs")))
clusterExport(clust, work.env)
clusterEvalQ(clust, {library(plyr)
                     library(geoR)
                     library(mnormt)})
# Broadcast random seeds for child processes of R. Makes results reproducible.
clusterSetRNGStream(cl = clust, iseed = 190684)

metros.Regret.sub <- clusterApplyLB(cl=clust, x=1:44, fun = MH.ind, iterations = 10000, round=1, warmup=2500, model = "Regret", start.val='dispersed')

#test <- lapply(20:30, MH.ind, iterations = 600, warmup=0, model = "Regret")

#posterior(c(0.2557, 0.0196, 0, 0.9507), n, f, X, S, prior="yes", hyperpars=c(0.7,0.3))

metros.EU.sub <- clusterApplyLB(cl=clust, x=1:44, fun = MH.ind, iterations = 10000, round=1, warmup=2500, model = "EU", start.val='dispersed')

warmup <- 2500
posterior.means.regret <- sapply(1:44, function(x) unlist(apply(metros.Regret.sub[[x]]$trace[(warmup+1):5500,], 2, mean)))
posterior.means.EU <- sapply(1:44, function(x) unlist(apply(metros.EU.sub[[x]]$trace[(warmup+1):5500,], 2, mean)))


DIC.EU <- sapply(1:44, function(x) metros.EU.sub[[x]]$DIC)
DIC.Regret <- sapply(1:44, function(x) metros.Regret.sub[[x]]$DIC)

plot(sort(DIC.Regret-DIC.EU)); abline(h=c(5,0,-5))
