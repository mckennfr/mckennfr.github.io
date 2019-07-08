# GSIM 4x4 model 

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

# calc numerical Jacobian
Jacobian <- function(x, eqn2, eqn3, phi, dT, xDemand) {
  J <- matrix(NA, nrow=4, ncol=4)
  h <- 1/32768
  
  for(i in 1:4) {
    x1 <- x
    x1[i] <- x[i] + h
    
    chi1 <- cross_price_effects(x1, eqn2, dT)
    
    dSupply1 <- x1 * Ex
    
    dDemand1 <- demand_effects(x1, phi, eqn3, chi1, dT)
    
    xDemand1 <- dDemand1 - dSupply1
    J[i,] <- (xDemand1 - xDemand ) / h
  }
  J
}

calcEqm <- function(x0, phi, eqn2, eqn3, dT) {
  x <- x0
  repeat {
    chi <- cross_price_effects(x, eqn2, dT)
    
    dSupply <- x * Ex
    
    dDemand <- demand_effects(x, phi, eqn3, chi, dT)
    
    xDemand <- dDemand - dSupply
    
    Jacob <- Jacobian(x, eqn2, eqn3, phi, dT, xDemand)
    
    x1 <- x - solve(Jacob) %*% xDemand / 1024 # 2^10 
    
    chi1 <- cross_price_effects(x1, eqn2, dT)
    
    dSupply1 <- x1 * Ex
    
    dDemand1 <- demand_effects(x1, phi, eqn3, chi1, dT)
    
    xDemand1 <- dDemand1 - dSupply1
    
    SSE <- sum(xDemand * xDemand)
    SSE1 <- sum(xDemand1 * xDemand1)
    if (SSE1 < SSE) x <- x1
    else break;
  }
  x
}

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
#T1["EU","US"] <- T1["US","EU"] <- 1 # US-EU trade agreement
T1["JP","US"] <- T1["US","JP"] <- 1 # US-JP agreement

# US alone...
#T1["JP","US"] <- T1["EU","US"] <- T1["RW","US"] <- 1

# JP alone...
#T1["US","JP"] <- T1["EU","JP"] <- T1["RW","JP"] <- 1

# a vector of ones
i_agg <- matrix(1, nrow=nrow(M), ncol=1)

# Total imports and exports
Tot_imp <- t(i_agg) %*% M
row.names(Tot_imp) <- "Imports"

Tot_exp <- M %*% i_agg
colnames(Tot_exp) <- "Exports"

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

# reported prices
# x <- c(0.0792, -0.0316, 0.0480, -0.0184)
#
# try from zero...
x <- calcEqm(rep(0, 4), phi, eqn2, eqn3, T1/T0 - 1)

chi <- cross_price_effects(x, eqn2, T1/T0 - 1)

dPCT <- (x %*% t(i_agg) + T1/T0 - 1) * eqn3 + chi

M1 <- M * (1 + dPCT) * (1 + x %*% t(i_agg))

dM <- M1 - M

dP <- (1 + x %*% t(i_agg)) * T1/T0 -1

dCompP <- matrix(NA, nrow=1, ncol = ncol(dP))
row.names(dCompP) <- "Composite Price"
colnames(dCompP) <- c("US", "JP", "EU", "RW")
for(i in 1:ncol(dCompP)) {
  dCompP[1,i] <- sum(c(dP[,i] * theta[,i]))
}

dTariff <- M1 * (T1 - 1) - M * (T0 - 1)

dRev <- colSums(dTariff)

dCS <- c((0.5 * dCompP ^ 2) * colSums(M * T0) * sign(dCompP) * Em - colSums(M * T0) * dCompP)
names(dCS) <- c("US", "JP", "EU", "RW")

dPS <- c(x * Tot_exp *(1 + 0.5 * Ex * x))
names(dPS) <- c("US", "JP", "EU", "RW")

dW <- dPS + dCS + dRev


# make report

df_report <- data.frame(Region = c("United States", "Japan", "European Union", "Rest of World"),
                        ProducerSurplus = format(round(dPS,1),nsmall=1),
                        ConsumerSurplus = format(round(dCS,1),nsmall=1),
                        TariffRevenue = format(round(dRev,1),nsmall=1),
                        NetWelfare = format(round(dW,1),nsmall=1))

df_plot <- data.frame(Region = c("United States", "Japan", "European Union", "Rest of World"), 
                      Component = "Producer Surplus", Value = dPS)

df_plot <- rbind(df_plot, 
                 data.frame(Region = c("United States", "Japan", "European Union", "Rest of World"), 
                            Component = "Consumer Surplus", Value = dCS))

df_plot <- rbind(df_plot, 
                 data.frame(Region = c("United States", "Japan", "European Union", "Rest of World"), 
                            Component = "Tariff Revenue", Value = dRev))

df_netW <- data.frame(Region=c("United States", "Japan", "European Union", "Rest of World"),
                      Component=NA,
                      Value=dW)


library(ggplot2)
g <- ggplot(df_plot, aes(x=Region, y=Value, fill = Component)) + 
     geom_bar(stat = "identity") +
     scale_x_discrete(limits=c("United States", "Japan", "European Union", "Rest of World")) +
     geom_hline(yintercept = 0) + 
     geom_point(data=df_netW, size=2, show.legend = FALSE) +
     labs(x="", y="Welfare", title="Total Welfare Effects", subtitle = "By Component")

plot(g)
