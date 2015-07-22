# Workhorse function
fct <- function(x, s, j, theta, k){
  # This function computes part of the regret value function
  (z**j) * (1 - z) * max(utils(s, theta), utils(x*h**j, theta))
}

transform <- function(x){
  y <- x$trace
  if(ncol(y)==3){
    y[,1] <- 1 - y[,1]
    y[,3] <- 0.5/(1-y[,3]) - 0.5
  } else {
    y[,1] <- 1 - y[,1]
    y[,4] <- 0.5/(1-y[,4]) - 0.5
  }
  return(y)
}

HPD <- function(x, prob){
  y <- mcmc(x)
  out <- HPDinterval(y)
  return(out)
}

# Workhorse function
V <- function(x, s, m, theta, k){
  # This function returns the expected (regret) utility at any point (x,s) for the 
  # cut-off strategy that stops when m upticks occur.
  a <- (z**m)*utils(x*h**m, theta)
  # sapply(x, FUN, ...) applies the function 'FUN' to every element of x and takes additional
  # arguments supplied after FUN. Below mapply() works the same way, but that
  # it applies the function to all elements of a matrix instead of a vector.
  # Look at function 'fct' below; it simply outsources part of the computations for
  # better readability. I try to follow your suggestion of functional programming.
  b <- k*sum(sapply(0:(m-1), fct, x=x, s=s, theta=theta, k=k))
  c <- k * (z**m) * max(utils(s, theta), utils(x*h**m, theta))
  return(a-b-c)
}


# Workhorse function
qstar <- function(x, s, theta, k){
  # Function returns the (minimal) optimal number of upticks before stopping
  res <- which.max(sapply(0:30, V, x=x, s=s, theta = theta, k=k))-1
  # This is a workhorse function. You may rather want to evaluate bstar().
  return(res)
}

# Returns the (minimal) optimal cut-off
bstar <- function(x,s,theta,k){
  # Function returns the (minimal) optimal cut-off
  res <- x*h**qstar(x,s,theta,k)
  return(res)
}

# Utility function
utils <- function(x, theta){
  if(type=="crra"){
    # utility function for CRRA
    # u(x - K) equals:
    u <- ((x/K)**theta - 1)*(K/theta)
    if(theta==0){
      u <- 0
    }
  }
  
  if(type == "cara"){
    # utility function for CARA
    # u(x - K) equals:
    u <- (1 - exp(-theta*(x-K)))
  }
  
  return(u)
}


# Computes the Stopping value
sv.crra <- function(x, s, theta, k){
  # computes the stopping value at a given point (x,s) for a given parameter (theta,k).
  stop.val <- utils(x, theta) - k * utils(s, theta)
  return(stop.val)
}

con.val <- function(y, s, theta, kappa, steps){
  theta <- theta 
  k <- kappa
  xhi <- outer(array(y, length(y)),h**(0:steps), FUN = "*")
  u.xhi <- apply(xhi, 2, utils, theta)
  zs <- z**(0:steps)
  immediate <- t(apply(u.xhi, 1, function(x) x*zs))
  tmp <- immediate
  
  if(k > 0){
    benchmark <- utils(s, theta)
    
    # Each row below is a timepoint and each column is an uptick in m
    # Matrix with max(u(s), u(x*h**i)): each row is for fixed t and varying m
    uvss2 <- array(0, c(length(y), ncol(u.xhi)))
    uvss <- t(sapply(1:length(y), function(x) pmax(benchmark[x], u.xhi[x,])))
    uvss2[, 2:ncol(uvss2)] <- t(sapply(1:length(y), function(x) pmax(benchmark[x], u.xhi[x,1:(ncol(u.xhi)-1)])))
    # interim is now the interim for date t for all m=0,..,100 in a vector, i.e. the l-th row of b
    interim <- k * (1-z) * (t(apply(t(apply(uvss2, 1, function(x) x*c(0,zs[1:(length(zs)-1)]))), 1, cumsum)))
    #interim <- k * (1-z) * (t(apply(t(apply(uvss, 1, function(x) x*zs)), 1, cumsum)))
    # a is now all a's in a vector, i.e. the l-th row of a
    reach <- k * t(apply(uvss, 1, function(x) zs*x))
    tmp <- immediate - interim - reach
  }
  
  res <- apply(tmp, 1, max)
  return(res)
}

cv.crra <- function(x, s, theta, k){
  # Expected utility from waiting one more period and continuing with the optimal cut-off strategy from there
  tmp <- cbind(x, s, seq_along(x))
  df1 <- unique(tmp[,1:2])
  conval <- delta*((p*(df1[,1]<df1[,2]))*con.val(h*df1[,1], df1[,2], theta, k, 120) + p*(df1[,1]==df1[,2])*con.val(h*df1[,1], h*df1[,2], theta, k, 120) + (1-p)*con.val(df1[,1]/h, df1[,2], theta, k, 120)) -  (1-delta)*(k*utils(df1[,2], theta))
  df2 <- cbind(df1, conval)
  ContinuationValue <- merge(tmp, df2, by=c('x', 's'))
  return(ContinuationValue[order( ContinuationValue[,3]),4])  
}

# Computes the stopping probability at a given point and for given parameters.
stop.prob <- function(x,s,theta,kappa,sigma,w){
  out <- 1 - ((1-w)*pnorm(cv.crra(x, s, theta, kappa) - sv.crra(x, s, theta, kappa), mean = 0, sd = sigma) +  w)
  return(out)
}

filter.data.old <- function(dat, i=NULL){
  require(plyr)
  
  # Do restrict attention to subject i?
  if(!is.null(i)){
    # A problem arises for subjects who did not stop in the first round.
    # Therefore sort the dat such that the rounds invested come first.
    daten <- dat[dat$individuals==i, ]
    daten <- daten[order(daten$invested, decreasing=TRUE), ]
    # Now iterate
    for(i in 1:nrow(daten)){
      if(daten[i,]$invested==TRUE){
        # The path until a subject stopped
        tmp <- process[(daten[i,]$rowIndex+1),1:(daten[i,]$investmentIdx+1),1:2]
        # The point where it stopped
        tmp2 <- process[(daten[i,]$rowIndex+1),(daten[i,]$investmentIdx+1),1:2]
      } else {
        # The whole path, where we need to delete NA's at the end.
        tmp <- process[(daten[i,]$rowIndex+1),,1:2]
        tmp <- tmp[!is.na(tmp)[,1], ]
      }
      
      # Dynamically build the objects
      if(i==1){
        master <- tmp
        stops <- tmp2
      } else {
        master <- rbind(master, tmp)
        if(daten$invested[i]==TRUE){
          stops <- rbind(stops, tmp2)
        }
      }
    }
    
    # Do we take the pooled dat else?
  } else {
    daten <- dat
    daten <- daten[order(daten$invested, decreasing=TRUE), ]
    for(i in 1:nrow(daten)){
      if(daten[i,]$invested==TRUE){
        tmp <- process[(daten[i,]$rowIndex+1),1:(daten[i,]$investmentIdx+1),1:2]
        tmp2 <- process[(daten[i,]$rowIndex+1),(daten[i,]$investmentIdx+1),1:2]
      } else {
        tmp <- process[(daten[i,]$rowIndex+1),,1:2]
        tmp <- tmp[!is.na(tmp)[,1], ]
      }
      if(i==1){
        master <- tmp
        # The following line makes sure this function also works with simulated data 
        if(class(tmp2)=='numeric'){ tmp2 <- t(as.matrix(tmp2))}
        stops <- tmp2
      } else {
        master <- rbind(master, tmp)
        if(daten$invested[i]==TRUE){
          if(class(tmp2)=='numeric'){ tmp2 <- t(as.matrix(tmp2))}
          stops <- rbind(stops, tmp2)
        }
      }
    }
  }
  
  #x <-master[,1]
  
  # Delete rownames - otherwise an ugly warning pops up.
  row.names(stops) <- NULL
  row.names(master) <- NULL
  
  
  # Give meaningful names and make columns integer for matching below.
  colnames(master) <- colnames(stops) <- c('x', 's')
  master[,1] <- round(master[,1])
  master[,2] <- round(master[,2])
  stops[,1] <- round(stops[,1])
  stops[,2] <- round(stops[,2])
  
  # Function count counts the number of occurrences of the tuples (x, s)
  master.count <- as.matrix(count(data.frame(master), c("x", "s")))
  stop.count <- as.matrix(count(data.frame(stops), c("x", "s")))
  
  # Since at many points (x, s) we do not observe a stopping decision
  # find the rows in the set of points subjects saw that have at least
  # one stopping decision and insert values there
  indx <- as.matrix(match(data.frame(t(stop.count[,1:2])), data.frame(t(master.count[,1:2]))))
  
  master.count <- data.frame(master.count)
  
  # Prob is the column with the relative frequency of stopping decisions
  # to visits: initialize as zero
  master.count$st <- 0
  
  # Insert the number of stopping decisions we can compute
  master.count$st[indx] <- stop.count[,3]
  
  # Compute the distance in ticks between x and s 
  # (rounded to the nearest integer to make sure it does not become a float)
  master.count$dist <- round(master.count$s - master.count$x)
  
  # Stopping probabilities: variable freq is n(x,s) in the notation of
  # the paper.
  master.count$probs <- master.count$st / master.count$freq
  
  colnames(master.count) <- c("x", "s", "n", "st", "dist", "stopprob")
  return(master.count)
}


# Filter the sufficient statistics from the data. For i = NULL (the default) do this across all
# subjects
filter.data <- function(dat, i=NULL){
  require(plyr)
  
  # Do restrict attention to subject i?
  if(!is.null(i)){
    # A problem arises for subjects who did not stop in the first round.
    # Therefore sort the dat such that the rounds invested come first.
    daten <- dat[dat$individuals==i, ]
    daten <- daten[order(daten$invested, decreasing=TRUE), ]
    rowindices <- as.vector(daten$rowIndex+1)
    invIdx <- as.vector(daten$investmentIdx+1)
    investeds <- as.vector(daten$invested)
    l <- nrow(daten)
    all.points <- vector(mode = "list", length = l)
    all.stops <- vector(mode = "list", length = l)
    for(i in 1:l){
      k <- rowindices[i]
      a <- invIdx[i]
      if(investeds[i]==TRUE){
        # The path until a subject stopped
        tmp <- process[k, 1:a, 1:2]
        # The point where it stopped
        tmp2 <- process[k, a, 1:2]
        if(class(tmp2)=='numeric'){ tmp2 <- t(as.matrix(tmp2))}
        all.stops[[i]] <- tmp2
      } else {
        # The whole path, where we need to delete NA's at the end.
        tmp <- process[k,,1:2]
        tmp <- tmp[!is.na(tmp)[,1], ]
      }
      
      all.points[[i]] <- tmp
    }
    
    # Do we take the pooled dat else?
  } else {
    daten <- dat
    daten <- daten[order(daten$invested, decreasing=TRUE), ]
    rowindices <- as.vector(daten$rowIndex+1)
    invIdx <- as.vector(daten$investmentIdx+1)
    investeds <- as.vector(daten$invested)
    l <- nrow(daten)
    all.points <- vector(mode = "list", length = l)
    all.stops <- vector(mode = "list", length = l)
    for(i in 1:l){
      k <- rowindices[i]
      a <- invIdx[i]
      if(investeds[i]==TRUE){
        tmp <- process[k,1:a,1:2]
        tmp2 <- process[k,a,1:2]
        if(class(tmp2)=='numeric'){ tmp2 <- t(as.matrix(tmp2))}
        all.stops[[i]] <- tmp2
      } else {
        tmp <- process[k,,1:2]
        tmp <- tmp[!is.na(tmp)[,1], ]
      }
      
      all.points[[i]] <- tmp
    }
  }
  
  master <- quick.matrix.2.df(do.call("rbind", all.points), c('x', 's'))
  stops <- quick.matrix.2.df(do.call("rbind", all.stops), c('x', 's'))
  
  # Delete rownames - otherwise an ugly warning pops up.
  row.names(stops) <- NULL
  row.names(master) <- NULL
  
  
  # Give meaningful names and make columns integer for matching below.
  colnames(master) <- colnames(stops) <- c('x', 's')
  master[,1] <- round(master[,1])
  master[,2] <- round(master[,2])
  stops[,1] <- round(stops[,1])
  stops[,2] <- round(stops[,2])
  
  # Function count counts the number of occurrences of the tuples (x, s)
  master.count <- as.matrix(count(data.frame(master), c("x", "s")))
  stop.count <- as.matrix(count(data.frame(stops), c("x", "s")))
  
  # Since at many points (x, s) we do not observe a stopping decision
  # find the rows in the set of points subjects saw that have at least
  # one stopping decision and insert values there
  indx <- as.matrix(match(data.frame(t(stop.count[,1:2])), data.frame(t(master.count[,1:2]))))
  
  master.count <- quick.matrix.2.df(master.count, c('x', 's', 'freq'))
  
  # Prob is the column with the relative frequency of stopping decisions
  # to visits: initialize as zero
  master.count$st <- 0
  
  # Insert the number of stopping decisions we can compute
  master.count$st[indx] <- stop.count[,3]
  
  # Compute the distance in ticks between x and s 
  # (rounded to the nearest integer to make sure it does not become a float)
  master.count$dist <- round(master.count$s - master.count$x)
  
  # Stopping probabilities: variable freq is n(x,s) in the notation of
  # the paper.
  master.count$probs <- master.count$st / master.count$freq
  
  colnames(master.count) <- c("x", "s", "n", "st", "dist", "stopprob")
  return(master.count)
}

simulate.stops <- function(M, theta, kappa, sigma, w,i=NULL){
  # Carrier for the return object
  master <- array(, c(M*dim(process)[1], 3))
  # Set of indices to write into master below
  a <- seq(1, nrow(master), M)
  b <- seq(M, nrow(master), M)
  for(i in 1:dim(process)[1]){
    tmp <- process[i,,1:2]
    tmp <- tmp[!is.na(tmp)[,1], ]
    tmp2 <- stop.prob(40*h**tmp[,1], 40*h**tmp[,2], theta, kappa, sigma, w)
    daten <-  matrix(c(rep(-2, M),rep(i-1, M), rep(0, M)), nrow=M)
    #daten <- data.frame(investmentIdx = rep(-2, M),
    #                    rowIndex = rep(i-1, M),
    #                    invested = rep(FALSE, M))
    for(j in 1:M){
      tmp3 <- vapply(tmp2, function(x) sample(c(0,1), 1, prob = c(1-x, x)), 1)
      if(max(tmp3)!=0){
        #daten$investmentIdx[j] <- min(which(tmp3==1))-1
        #daten$invested[j] <- TRUE
        daten[j, 1] <- min(which(tmp3==1))-1
        daten[j, 3] <- TRUE
      }
    }
    master[a[i]:b[i], ] <- daten
    
    #if(i==1){
    #  master <- daten
    #} else {
    #  master <- rbind(master, daten)
    #}
  }
  
  # Hat tip: Hadley Wickham
  out <- quick.matrix.2.df(master, c("investmentIdx", "rowIndex", "invested"))
  out[,3] <- as.logical(out[,3])
  #return(master)
  return(out)
}

simulate.stops.old <- function(M, theta, kappa, sigma, w,i=NULL){
  # Carrier for the return object
  #master <- array(, c(M*dim(process)[1], 3))
  # Set of indices to write into master below
  #a <- seq(1, nrow(master), M)
  #b <- seq(M, nrow(master), M)
  for(i in 1:dim(process)[1]){
    tmp <- process[i,,1:2]
    tmp <- tmp[!is.na(tmp)[,1], ]
    tmp2 <- stop.prob(40*h**tmp[,1], 40*h**tmp[,2], theta, kappa, sigma, w)
    #daten <-  matrix(c(rep(-2, M),rep(i-1, M), rep(0, M)), nrow=M)
    daten <- data.frame(investmentIdx = rep(-2, M),
                        rowIndex = rep(i-1, M),
                        invested = rep(FALSE, M))
    for(j in 1:M){
      tmp3 <- vapply(tmp2, function(x) sample(c(0,1), 1, prob = c(1-x, x)), 1)
      if(max(tmp3)!=0){
        daten$investmentIdx[j] <- min(which(tmp3==1))-1
        daten$invested[j] <- TRUE
        #daten[j, 1] <- min(which(tmp3==1))-1
        #daten[j, 3] <- TRUE
      }
    }
    #master[a[i]:b[i], ] <- daten
    
    if(i==1){
      master <- daten
    } else {
      master <- rbind(master, daten)
    }
  }
  
  # Hat tip: Hadley Wickham
  #out <- quickdf(master, c("InvestmentIdx", "rowIndex", "invested"))
  
  return(master)
  #return(out)
}

quick.matrix.2.df <- function(x, y) {
  d <- dim(x)
  nrows <- d[1L]
  ir <- seq_len(nrows)
  ncols <- d[2L]
  ic <- seq_len(ncols)
  #dn <- dimnames(x)
  #if (is.null(row.names)) 
  #  row.names <- dn[[1L]]
  #collabs <- dn[[2L]]
  #if (any(empty <- !nzchar(collabs))) 
  #  collabs[empty] <- paste0("V", ic)[empty]
  value <- vector("list", ncols)
  #if (mode(x) == "character" && stringsAsFactors) {
  #  for (i in ic) value[[i]] <- as.factor(x[, i])
  #}
  #else {
  for (i in ic) value[[i]] <- as.vector(x[, i])
  #}
  #if (is.null(row.names) || length(row.names) != nrows) 
  #  row.names <- .set_row_names(nrows)
  #if (length(collabs) == ncols) 
  names(value) <- y
  #else if (!optional) 
  #  names(value) <- paste0("V", ic)
  attr(value, "row.names") <- .set_row_names(nrows)
  class(value) <- "data.frame"
  value
}

quick.numeric.2.df <- function(x, y){
  nrows <- length(x)
  row.names <- .set_row_names(nrows)
  value <- list(x)
  names(value) <- y
  attr(value, "row.names") <- row.names
  class(value) <- "data.frame"
  value
}

quick.list.2.df <- function(l) {
  class(l) <- "data.frame"
  attr(l, "row.names") <- .set_row_names(length(l[[1]]))
  l
}

classifier <- function(x, by, max.ticks){
  if(max.ticks > max(x)){
    warning('max.ticks larger than largest element of x -> I will use max(x) instead of max.ticks')
    max.ticks <- max(x)
  }
  if(max.ticks<=0){stop('max.ticks must be strictly greater zero!')}
  labels <- rep(NA, length(x))
  a <- seq(0, max.ticks, by)
  b <- seq(by, max.ticks, by)
  l.a <- length(a)
  l.b <- length(b)
  indx <- rep(NA, length(x))
  aa <- a[seq(1, l.a, 2)]
  bb <- b[seq(1, l.b, 2)]
  l.aa <- length(aa)
  l.bb <- length(bb)
  if(l.aa > l.bb){
    labs <- c(paste(aa[1:l.bb], 'to', bb, sep=' '), paste(aa[l.aa], "or more", sep=' '))
    c <- cbind(aa[1:l.bb], bb)
    tmp <- sapply(x, function(x) apply(c, 1, function(y) y[1] <= x & x <= y[2]))
    indx <- which(apply(tmp, 2, sum)==0)
    labels[indx] <- labs[length(labs)]
    indx2 <- apply(tmp[, -indx], 2, function(x) which(x==TRUE))
    labels[-indx] <- labs[indx2]
    indx.zero <- which(x==0)
    labels[indx.zero] <- labs[1]
  } else {
    labs <- paste(a, 'to', b, sep=' ')
    c <- cbind(a[1:l.b], b)
    tmp <- sapply(x, function(x) apply(c, 1, function(y) y[1] < x & x <= y[2]))
    indx <- which(apply(tmp, 2, sum)==0)
    labels[indx] <- labs[length(labs)]
    indx2 <- apply(tmp[, -indx], 2, function(x) which(x==TRUE))
    labels[-indx] <- labs[indx2]
  }
  return(labels)
}

classifier.2.step <- function(y){
  if(y<=1){
    out <- '0 to 1'
  } 
  if(y>1 & y<=3){
    out <- '2 to 3'
  } 
  if(y>3 & y <=5){
    out <- '4 to 5'
  }
  if(y>5){
    out <- '6 or more'
  } 
  #if(y > 7){
  #  out <- '8 or more'
  #}
  return(out)
}

classifier.3.step <- function(y){
  if(y<=2){
    out <- '0 to 2'
  } 
  if(y>2 & y<=5){
    out <- '3 to 5'
  } 
  if(y>5){
    out <- '6 or more'
  }
  #if(y > 8){
  #  out <- 'more than 8'
  #}
  return(out)
}


# Return the probability to continue at (x, s) given the parameters.
Psi.rho <- function(x, s, theta, k, w, sigma){
  rho <- cv.crra(x, s, theta, k) - sv.crra(x, s, theta, k)
  Psi.rho <- (1 - w) * pnorm(rho, mean=0, sd=sigma) + w
  return(Psi.rho)
}

# Calculates the log-posterior given a prior. If prior is NULL, returns log-likelihood.
posterior <- function(Theta, n, f, x, s, prior=NULL, hyperpars=NULL){
  Theta <- Theta * parFac
  
  if(length(Theta)==4){
    Lower <- c(1e-5,0,0,0) ; Upper <- c(Inf,1,1e6,1)
    penFac <- 1 + sum(pmax(0, Lower-Theta)^1.1) + sum (pmax(0, Theta-Upper)^1.1)
    
    Theta <- pmax(Lower, pmin(Upper, Theta))
    
    theta <- Theta[1]
    kappa <- Theta[2]
    sigma <- Theta[3]
    w <- Theta[4]
    
    cat ("theta =", round(c(theta), digits=4), ", kappa = ", c(round(kappa, digits=4)), ", sigma = ", round(sigma,4), ", w = ", c(round(w, digits=4)))
  }
  
  if(length(Theta)==3){
    Lower <- c(1e-5,0,0) ; Upper <- c(Inf,1e6,1)
    penFac <- 1 + sum(pmax(0, Lower-Theta)^1.1) + sum (pmax(0, Theta-Upper)^1.1)
    
    Theta <- pmax(Lower, pmin(Upper, Theta))
    
    theta <- Theta[1]
    kappa <- 0
    sigma <- Theta[2]
    w <- Theta[3]
    
    cat ("theta =", round(c(theta), digits=4), ", kappa = ", c(round(kappa, digits=4)), ", sigma = ", round(sigma,4), ", w = ", c(round(w, digits=4)))
  }
  cont.prob <- Psi.rho(x, s, theta, kappa, w, sigma)
  # Avoid singularities here:
  cont.prob[cont.prob==1] <- 0.9999
  log.stop.prob <- log(1 - cont.prob)
  tmp <- f * log(cont.prob/(1-cont.prob)) + log.stop.prob
  out <- sum(n * (tmp))
  if(!is.null(prior)){
    out <- out + prior.theta(Theta, hyperpars)
  }
  if(sigma==0){ out = -Inf}
  if(is.nan(out)){ out = -Inf }
  cat (" => ln(pd) = ", out, "\n")
  return(out*penFac)
}

# Returns the log-prior probability of the parameters given hyperparameters
prior.theta <- function(Theta, hyperpars){
  require(geoR)
  if(length(Theta)==3){
    theta <- Theta[1]
    sigma <- Theta[2]
    tremble <- Theta[3]
    m.theta <- hyperpars[1] 
    s.theta <- hyperpars[2]
    a <- dnorm(theta, mean=m.theta, sd=s.theta, log=TRUE)
    b <- dbeta(tremble, shape1=5, shape2=1, log=TRUE)#dunif(tremble, min=0, max=1, log=TRUE)
    c <- dinvchisq(sigma, 1, 1, log=TRUE)
    if(sigma==0){c <- -Inf}
    out <- a + b + c
  }
  if(length(Theta)==4){
    theta <- Theta[1]
    kappa <- Theta[2]
    sigma <- Theta[3]
    tremble <- Theta[4]
    m.theta <- hyperpars[1] 
    s.theta <- hyperpars[2]
    a <- dnorm(theta, mean=m.theta, sd=s.theta, log=TRUE)
    b <- dbeta(tremble, shape1=5, shape2=1, log=TRUE)#dunif(tremble, min=0, max=1, log=TRUE)
    c <- dunif(kappa, min=0, max=1, log=TRUE)
    d <- dinvchisq(sigma, 1, 1, log=TRUE)
    if(sigma==0){c <- -Inf}
    out <- a + b + c + d
  }
  return(out)
}

# A wrapper around optim() that finds the posterior mode.
Mode <- function(sp,n,f,X,S,prior,hyperpars){
  parFac <<- pmax(1e-2, abs(sp))
  mode <- optim(sp/parFac,
                fn=posterior, method="BFGS", 
                n=n, f=f, x=X, s=S, 
                hessian = TRUE,
                prior = prior,
                hyperpars = hyperpars,
                control=list(fnscale=-1))
  if (0) {
    cat ("Restart at first optimum -- try to refine solution\n")
    sp <- mode$par * parFac
    parFac <<- pmax(1e-2, abs(sp))
    mode <- optim(sp/parFac, 
                  fn=posterior, method="Nelder-Mead", 
                  n=n, f=f, x=X, s=S, 
                  hessian = TRUE,
                  prior = prior,
                  hyperpars = hyperpars,
                  control=list(fnscale=-1))
  }
  mode$par <- mode$par * parFac
  parFac <<- 1
  if(!is.null(prior)){
    mode$prior <- prior
    mode$hyperpars <- hyperpars
  }
  return(mode)
}

# A wrapper around the Metro_Hastings function from package 'MHadaptive'.
# IMPORTANT: The multivariate normal proposal density of the sampler uses
# the function mvrnorm() to draw proposals. The fct. mvrnorm() has the 
# covariance --- not its 'square root' -- as its second argument.
MH <- function(start, HESS, iterations, warmup, n, f, X, S, prior, mode){
  require(MHadaptive)
  
  parFac <<- 1
  
  prior = mode$prior
  hyperpars = mode$hyperpars
  
  mha <- Metro_Hastings(posterior, pars=start, 
                        n = n, 
                        f = f,
                        x = X,
                        s = S,
                        prior = prior,
                        hyperpars = hyperpars,
                        prop_sigma = HESS,
                        iterations=iterations,
                        burn_in=warmup)
  return(mha)
}


MH.ind <- function(subject, iterations, warmup, round, model, start.val=NULL){
  require(MHadaptive)
  require(mnormt)
  
  
  
  cat ("Subject", subject)
  
  if(is.null(start.val)){
    if(model == "EU"){
      start = EU.modes[[subject]]$par
      HESS = makePositiveDefinite(solve(-EU.modes[[subject]]$hessian))
      prior = EU.modes[[subject]]$prior
      hyperpars = EU.modes[[subject]]$hyperpars
    }
    
    if(model == "Regret"){
      start = Regret.modes[[subject]]$par
      HESS = makePositiveDefinite(solve(-Regret.modes[[subject]]$hessian))
      prior = Regret.modes[[subject]]$prior
      hyperpars = Regret.modes[[subject]]$hyperpars
    }
  } else {
    if(model == "EU"){
      Lower <- c(1e-5,0.5,0.09) ; Upper <- c(Inf,1e3,0.95)
      HESS = makePositiveDefinite(solve(-EU.modes[[subject]]$hessian))
      start = rmnorm(n=1, mean=EU.modes[[subject]]$par, varcov=makePositiveDefinite(2*HESS))
      start <- pmax(Lower, pmin(Upper, start))
      prior = EU.modes[[subject]]$prior
      hyperpars = EU.modes[[subject]]$hyperpars
    }
    
    if(model == "Regret"){
      Lower <- c(1e-5,0,0.5,0.09) ; Upper <- c(Inf,1,1e3,0.95)
      HESS = makePositiveDefinite(solve(-Regret.modes[[subject]]$hessian))
      start = rmnorm(n=1, mean=Regret.modes[[subject]]$par, varcov=makePositiveDefinite(2*HESS))
      start <- pmax(Lower, pmin(Upper, start))
      prior = Regret.modes[[subject]]$prior
      hyperpars = Regret.modes[[subject]]$hyperpars
    }
  }  
  
  parFac <<- 1
  
  daten.per.subject <- filter.data(data[data$serieIndex>round, ], subject)
  
  X <- 40*h**daten.per.subject$x
  S <- 40*h**daten.per.subject$s
  n <- daten.per.subject$n
  
  #f in the notation below is the _continuation probability_
  # whereas it is the stopping probability below - unfortunate!
  f <- 1 - daten.per.subject$stopprob
  
  
  
  mha <- Metro_Hastings(posterior, pars=start, 
                        n = n, 
                        f = f,
                        x = X,
                        s = S,
                        prior = prior,
                        hyperpars = hyperpars,
                        prop_sigma = HESS,
                        iterations=iterations,
                        burn_in=warmup)
  return(mha)
}

estimate.model <- function(start.val, subject, round, prior=NULL, hyperpars=NULL){
  
  daten.per.subject <- filter.data(data[data$serieIndex>round, ], subject)
  
  X <- 40*h**daten.per.subject$x
  S <- 40*h**daten.per.subject$s
  n <- daten.per.subject$n
  
  #f in the notation below is the _continuation probability_
  # whereas it is the stopping probability below - unfortunate!
  f <- 1 -daten.per.subject$stopprob
  
  out <- Mode(start.val, n, f, X, S, prior, hyperpars)
  return(out)  
}

xstar <- function(theta){
  
  # The optimal number of upticks from the starting value x0 
  # We take a grid of value for the number of upticks q, 
  # ranging from zero upticks (the starting value)
  # up to 30 upticks (which is 229.7396 for our parameters).
  q=which.min(sapply(x0*h**seq(0:30),   # this is the set of possible values of the process we defined as chi in the paper
                     function(x){
                       sign(utils(x*h,theta)-h**alpha*utils(x,theta))   # find the first tick q* for which stopping becomes optimal
                     }
  )
  )
  
  # Return the value of the process after q* upticks
  return(x0*h**q)
}

Gelman.diag.ind <- function(chain1, chain2, bin.width, Model){
  # No of iterations
  n <- nrow(chain1[[1]]$trace)
  # No. of parameters
  k <- ncol(chain1[[1]]$trace)
  # Intervals to consider
  ints <- seq(bin.width, n, bin.width)
  # No. of intervals
  m <- length(ints)
  
  # Carrier
  psrf <- array(, c(m, k))
  
  # Take fct. off the shelve
  for(x in 1:length(chain1)){
    for(i in 1:m){
      chains <- list(mcmc(chain1[[x]]$trace[1:ints[i], ]), mcmc(chain2[[x]]$trace[1:ints[i], ]))
      tmp <- gelman.diag(chains, multivariate=FALSE)
      psrf<- tmp$psrf[, 1]
    }
    if(x == 1){
      PSRF <- psrf                         
    } else {
      PSRF <- c(PSRF, psrf)
    }
  }
  
  if(Model == "EU"){
    PSRF <- data.frame(PSRF, 
                       Parameter = rep(c("theta", "sigma", "w"), times = m*44),
                       Subject = rep(1:44, each=3*m),
                       Iterations = rep(ints, each=3, times=44))
  } else {
    PSRF <- data.frame(PSRF, 
                       Parameter = rep(c("theta", "kappa", "sigma", "w"), times = m*44),
                       Subject = rep(1:44, each=4*m),
                       Iterations = rep(ints, each=4, times=44))
  }
  
  
  return(PSRF)
}


Gelman.diag.pool <- function(x, bin.width){
  # No of iterations
  n <- nrow(x[[1]]$trace)
  # No. of parameters
  k <- ncol(x[[1]]$trace)
  # Intervals to consider
  ints <- seq(bin.width, n, bin.width)
  # No. of intervals
  m <- length(ints)
  
  # Carrier
  psrf <- array(, c(m, k))
  ci <- array(, c(m, k))
  
  # Take fct. off the shelve
  for(i in 1:m){
    chains <- lapply(1:length(x), function(y) mcmc(x[[y]]$trace[1:ints[i], ]))
    tmp <- gelman.diag(chains, multivariate=FALSE)
    psrf[i, ] <- tmp$psrf[, 1]
    ci[i, ] <- tmp$psrf[, 2]
  }
  
  convergence.idx <- ints[max(apply(psrf,2, function(x) min(which(x < 1.1))))]
  
  return(list(psrf = psrf, CI = ci, conv = convergence.idx))
}

# Creates the data matrix necessary to plot our barplots.
# x: output from filter.data
# EU, Regret: output from filter.data based on simulated choices.
# If items EU and Regret are provided this also binds the model predictions to the final output.
plot.data <- function(actual, block.size, cap.dist, dist.block.size, from=NULL, to = NULL, max.dist=NULL, EU=NULL, Regret=NULL){
  # Be careful: indx contains NA's and deleting them breaks the exact matching!
  
  # 1. Find rows in plot.dat.sim that have a corresponding row in actual
  if(!is.null(EU)){
    indx.EU <- as.matrix(match(data.frame(t(EU[,1:2])), data.frame(t(actual[,1:2]))))
    kills.EU <- which(is.na(indx.EU[,1]))
    index.EU <- indx.EU[-kills.EU,]
    tmp <- array(0, c(nrow(actual), 3))
    tmp[index.EU,] <- as.matrix(EU[-kills.EU,c(3,4,6)])
    
    
    actual$predicted.n.EU <- tmp[,1]
    actual$predicted.st.EU <- tmp[,2]
    actual$predicted.stopprob.EU <- tmp[,3]
  }
  if(!is.null(Regret)){
    indx.Regret <- as.matrix(match(data.frame(t(Regret[,1:2])), data.frame(t(actual[,1:2]))))
    # 2. Which do NOT have one
    kills.Regret <- which(is.na(indx.Regret[,1]))
    # 3. Take these out, as we cannot merge them with the original data (we could by expanding the original, data, but we omit them)
    index.Regret <- indx.Regret[-kills.Regret,]
    # 4. Write simulated ct and st into an array that we then bind to the actual item below    
    tmp <- array(0, c(nrow(actual), 3))
    tmp[index.Regret,] <- as.matrix(Regret[-kills.Regret,c(3,4,6)])
    
    actual$predicted.n.Regret <- tmp[,1]
    actual$predicted.st.Regret <- tmp[,2]
    actual$predicted.stopprob.Regret <- tmp[,3]
  }
  
  # Determine blocks
  ref.seq<-unique(actual$x)
  counter <- 0
  idx <- array(, c(length(ref.seq)))
  indx <- array(, c(nrow(actual)))
  for(i in 1:length(ref.seq)){
    if(ref.seq[i]%%block.size==0){counter <- counter + 1}
    idx[i] <- counter 
    tmp <- which(actual$x==ref.seq[i])
    indx[tmp] <- counter
  }
  actual$xblocks <- indx
  
  for(i in unique(indx)){
    name <- paste(round(40*h**actual$x[min(which(indx==i))],1), " to ", round(40*h**actual$x[max(which(indx==i))], 1))
    actual$xblocks[which(indx==i)] <- name
  }
  
  actual$x <- 40*h**actual$x
  actual$s <- 40*h**actual$s
  
  if(!is.null(max.dist)){
    plot.dat2 <- actual[actual$dist <= max.dist, ]
  } else {
    plot.dat2 <- actual
  }
  check.dat <- plot.dat2
  if(!is.null(from)){
    plot.dat2 <- plot.dat2[plot.dat2$x >= from & plot.dat2$x <= to, ]
  }
  
  # Construct blocks of ticks below the maximum:
  #plot.dat2$distblock <- classifier(plot.dat2$dist, dist.block.size, cap.dist)
  
  if(dist.block.size == 2){
    plot.dat2$distblock <- sapply(plot.dat2$dist, classifier.2.step)
  }
  if(dist.block.size == 3){
    plot.dat2$distblock <- sapply(plot.dat2$dist, classifier.3.step)
  }
  
  # Bar plot as we had it before: you can change the variables here to get different
  # graphs with different grouping etc.
  # Now aggregate the stopping probabilities
  tmp <- unique(plot.dat2$xblocks)
  n <- length(tmp)
  k <- length(unique(plot.dat2$distblock))
  bars <- array(0, c(n, k))
  if(!is.null(EU)){
    preds.EU <- array(0, c(n, k))
  }
  if(!is.null(Regret)){
    preds.Regret <- array(0, c(n, k))
  }
  counter <- 1
  for(i in tmp){
    tmp1 <- tapply(plot.dat2$st[plot.dat2$xblocks==i], plot.dat2$distblock[plot.dat2$xblocks==i], sum)
    tmp2 <- tapply(plot.dat2$n[plot.dat2$xblocks==i], plot.dat2$distblock[plot.dat2$xblocks==i], sum)
    bars[counter, 1:length(tmp1)] <- tmp1/tmp2
    if(length(tmp1)==k){  colnames(bars) <- names(tmp1)}#unique(plot.dat2$distblock) }
    if(!is.null(EU)){
      tmp3 <- tapply(plot.dat2$predicted.st.EU[plot.dat2$xblocks==i], plot.dat2$distblock[plot.dat2$xblocks==i], sum)
      tmp4 <- tapply(plot.dat2$predicted.n.EU[plot.dat2$xblocks==i], plot.dat2$distblock[plot.dat2$xblocks==i], sum)
      preds.EU[counter, 1:length(tmp1)] <- tmp3/tmp4
    }
    if(!is.null(Regret)){
      tmp5 <- tapply(plot.dat2$predicted.st.Regret[plot.dat2$xblocks==i], plot.dat2$distblock[plot.dat2$xblocks==i], sum)
      tmp6 <- tapply(plot.dat2$predicted.n.Regret[plot.dat2$xblocks==i], plot.dat2$distblock[plot.dat2$xblocks==i], sum)
      preds.Regret[counter, 1:length(tmp1)] <- tmp5/tmp6
    }
    counter <- counter + 1
  }
  
  require(reshape2)
  #colnames(bars) <- names(tmp1)#unique(plot.dat2$distblock)
  rownames(bars) <- tmp#unique(plot.dat2$xblocks)
  plot.dat4 <- melt(bars)
  if(!is.null(EU)){
    colnames(preds.EU) <- unique(plot.dat2$distblock)
    rownames(preds.EU) <- unique(plot.dat2$xblocks)
    tmp <- melt(preds.EU)
    plot.dat4$predicted.prob.EU <- tmp[,3]
  }
  if(!is.null(Regret)){
    colnames(preds.Regret) <- unique(plot.dat2$distblock)
    rownames(preds.Regret) <- unique(plot.dat2$xblocks)
    tmp2 <- melt(preds.Regret)
    plot.dat4$predicted.prob.Regret <- tmp2[,3]
  }
  
  if(!is.null(EU) & !is.null(Regret)){
    colnames(plot.dat4) <- c("x", "dist", "stopprob", "predicted.prob.EU", "predicted.prob.Regret")
  } else {
    if(!is.null(EU)){
      colnames(plot.dat4) <- c("x", "dist", "stopprob", "predicted.prob.EU")
    }
    if(!is.null(Regret)){
      colnames(plot.dat4) <- c("x", "dist", "stopprob", "predicted.prob.Regret")
    }
    if(is.null(EU) & is.null(Regret)){
      colnames(plot.dat4) <- c("x", "dist", "stopprob")
    }
  }
  #return(list(plot.dat = plot.dat4, check.dat = check.dat))
  return(plot.dat4)
}

posterior.sample <- function(x, N, M){
  rows <- sample(M:nrow(x), size = N, replace = TRUE, prob = rep(1/(nrow(x)-M+1), (nrow(x)-M+1)))
  out <- x[rows, ]
  return(out)
}


HLprior <- function(theta){
  r <- 1 - theta
  p <- ifelse(r<(-0.95), exp(log(0.01)/0.95)**(-r),
              ifelse(r>=(-0.95) & r<(-0.49), 0.01, 
                     ifelse(r>=(-0.49) & r < (-0.15), 0.06,
                            ifelse(r >=(-0.15) & r < 0.15, 0.26,
                                   ifelse(r>=0.15 & r < 0.41, 0.26,
                                          ifelse(r >= 0.41 & r < 0.68, 0.23, 
                                                 ifelse(r>=0.68 & r < 0.97, 0.13,
                                                        ifelse(r>=0.97 & r < 1.37, 0.03, 0.03**r
                                                        ))))))))
  return(p)
}