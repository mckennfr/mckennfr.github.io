# GSIM 4x4 model 

# elasticities
Em <- c(US = -1.25, JP = -1.25, EU = -1.25, RW = -1.25) # import demand

Ex <- c(US = 1.5, JP = 1.5, EU = 1.5, RW = 1.5) # export supply

Es <- c(US = 5, JP = 5, EU = 5, RW = 5) # substitution


# trade flows
M <- matrix(c(  0,  50, 200, 300,
              500,   0, 150, 200,
              300, 100, 200, 200,
               50, 100, 110,  20),
            nrow = 4, ncol=4, byrow = TRUE)

row.names(M) <- c("US", "JP", "EU", "RW")
colnames(M)  <- c("US", "JP", "EU", "RW")

# initial import tariffs 
T0 <- matrix(c(  1.00, 1.21, 1.41, 1.22,
                 1.37, 1.00, 1.31, 1.23,
                 1.32, 1.36, 1.00, 1.18,
                 1.57, 1.41, 1.25, 1.15),
              nrow = 4, ncol=4, byrow = TRUE)

row.names(T0) <- c("US", "JP", "EU", "RW")
colnames(T0)  <- c("US", "JP", "EU", "RW")

# final import tariffs
T1 <- T0
T1["EU","US"] <- T1["US","EU"] <- 1 # US-EU trade agreement

# a vector of ones
i_agg <- matrix(1, nrow=nrow(M), ncol=1)

# import shares at internal prices
mult <- diag( c(1 / (t(i_agg) %*% (M * T0))))
row.names(mult) <- row.names(M)
colnames(mult) <- colnames(M)

theta <- (M * T0) %*% mult

# export shares at world prices
mult <- diag(1/c(M %*% i_agg))
row.names(mult) <- row.names(M)
colnames(mult) <- colnames(M)

phi <- mult %*% M

# own price elasticities
eqn3 <- theta %*% diag( Em ) - (1 - theta) %*% diag( Es )
colnames(eqn3) <- colnames(M)

# cross price elasticities
eqn2 <- theta %*% diag( Em + Es )
colnames(eqn3) <- colnames(M)

# cross price effects
cross_price_effects <- function(x = rep(0,4), eq2, dT) {
  n <- ncol(eq2)
  rslt <- matrix(NA_real_, nrow=n, ncol = n)
  for(i in 1:n) {
    tmp <- matrix(1, nrow=n, ncol=1) %*% eq2[,i]
    tmp[n*1:n + 1:n - n] <- 0
    rslt[,i] <- tmp %*% (dT[,i] + x)
  }
  rslt
}


# change in demand
demand_effects <- function(x = rep(0,4), phi, eq3, chi, dT) {
  n <- length(x)
  rslt <- rep(NA_real_, n)
  for(i in 1:n) {
    rslt[i] <- (phi[i,] * eqn3[i,]) %*% (dT[i,] + x[i])
  }
  rslt + rowSums(phi * chi)
}

# reported prices
# x <- c(0.0792, -0.0316, 0.0480, -0.0184)
#
# try from zero...
x <- rep(0, 4)

repeat {
  chi <- cross_price_effects(x, eqn2, T1/T0 - 1)
  
  dSupply <- x * Ex
  
  dDemand <- demand_effects(x, phi, eqn3, chi, T1/T0 - 1)
  
  xDemand <- dDemand - dSupply

  Jacob <- matrix(NA, nrow=4, ncol=4)
  h <- 1/32768
  
  for(i in 1:4) {
    x1 <- x
    x1[i] <- x[i] + h
    
    chi1 <- cross_price_effects(x1, eqn2, T1/T0 - 1)
    
    dSupply1 <- x1 * Ex
    
    dDemand1 <- demand_effects(x1, phi, eqn3, chi1, T1/T0 - 1)
    
    xDemand1 <- dDemand1 - dSupply1
    Jacob[i,] <- (xDemand1 - xDemand ) / h
  }
  
  x1 <- x - solve(Jacob) %*% xDemand / 1024 # 2^10 
  
  chi1 <- cross_price_effects(x1, eqn2, T1/T0 - 1)
  
  dSupply1 <- x1 * Ex
  
  dDemand1 <- demand_effects(x1, phi, eqn3, chi1, T1/T0 - 1)
  
  xDemand1 <- dDemand1 - dSupply1
  
  SSE <- sum(xDemand * xDemand)
  SSE1 <- sum(xDemand1 * xDemand1)
  if (SSE1 < SSE) x <- x1
  else break;
}

