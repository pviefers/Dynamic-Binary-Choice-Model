source("functions.R")
install.packages(c("plyr", "geoR", "mnormt", "MHadaptive"))

library(parallel)
cores <- detectCores()
clust <- makeCluster(cores)
# Inherit the workspace to each node
work.env <- c(setdiff(ls(), c("metros.EU", "metros.Regret", "metros.EU.sub", "metros.Regret.sub", 
                              "simulated.decisions.EU", 
                              "simulated.decisions.Regret",
                              "stopping.frequencies.EU", 
                              "stopping.frequencies.Regret",
                              "EUs",
                              "Regs")))#work.env <- ls()
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
MLEs.Regret <- unlist(sapply(1:44, function(x) Regret.MLEs[[x]]$par))


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

metros.Regret.sub <- clusterApplyLB(cl=clust, x=1:44, fun = MH.ind, iterations = 10000, round=1, warmup=2500, model = "Regret")

#test <- lapply(20:30, MH.ind, iterations = 600, warmup=0, model = "Regret")

#posterior(c(0.2557, 0.0196, 0, 0.9507), n, f, X, S, prior="yes", hyperpars=c(0.7,0.3))

metros.EU.sub <- clusterApplyLB(cl=clust, x=1:44, fun = MH.ind, iterations = 10000, round=1, warmup=2500, model = "EU", start.val='dispersed')

warmup <- 0
posterior.means.regret <- sapply(1:44, function(x) unlist(apply(metros.Regret.sub[[x]]$trace[(warmup+1):7500,], 2, mean)))
posterior.means.EU <- sapply(1:44, function(x) unlist(apply(metros.EU.sub[[x]]$trace[(warmup+1):7500,], 2, mean)))


DIC.EU <- sapply(1:44, function(x) metros.EU.sub[[x]]$DIC)
DIC.Regret <- sapply(1:44, function(x) metros.Regret.sub[[x]]$DIC)

plot(sort(DIC.Regret-DIC.EU)); abline(h=c(5,0,-5))

length(which(DIC.Regret-DIC.EU < 0))
