library(SDEtools)
library(Matrix)
library(fields)

setwd('C:/Users/utilisateur/Documents/DTU/Première Année 2023-2024/Cours/Stochastic Differential Equations/KolmogorovEquations')
source('fvade.R')

#### Global variables
## Spatial grid
xmax <<- 10
xv <<- seq(0,xmax,length=101)
dx <<- diff(xv)
xc <<- xv[-1] - 0.5*dx

forwardkolmogorov <- function(u,D){
  
  ## Discretize the generator
  G <- fvade(u,D,xv,'r')
  
  pi <- StationaryDistribution(G)
  phi <- pi / dx
  
  return(list(
    "G"=G,
    "phi"=phi))
}


solve_forwardEquation <- function (tv,phi0,G) {
  # phi0 is the initial condition
  return(sapply(tv,function(t) as.numeric(phi0 %*% expm(G*t)))/dx)
}

solve_backwardEquation <- function(tv,h,G) {
  # h is the initial terminal condition
  return(sapply(tv,function(t) as.numeric(expm(G*(T-t)) %*% h)))
}

computeMoment <- function(PHI,k) apply(PHI*dx*xc^k,2,sum)
computeVariance <- function(PHI) computeMoment(PHI,2) - computeMoment(PHI,1)^2
computeStandardDeviation <- function(PHI) sqrt(computeVariance(PHI))



 