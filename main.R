setwd('C:/Users/utilisateur/Documents/DTU/Première Année 2023-2024/Cours/Stochastic Differential Equations/KolmogorovEquations')

source('fvade.R')
source('kolmogorov-equations.R')


# Parameters
xi <- 2
gamma <- 1/2
lambda <- 1/2

## Drift and intensity
f <- function(x) lambda * (xi - x)
g <- function(x) gamma*sqrt(abs(x))

## Diffusivity and its spatial derivative
D <- function(x) 0.5*gamma^2*x
Dp <- function(x) 0.5*gamma^2

## Advective flow field
u <- function(x) f(x) - Dp(x)


result <- forwardkolmogorov(u,D)
G <- result$G
phi <- result$phi

x11()
plot(xc, phi, main="Forward Kolmogorov Equation Solution", xlab="xc", ylab="phi")
curve(dgamma(x, rate=2*lambda/gamma^2, shape=2*lambda*xi/gamma^2), 
      from=0, to=xmax, add=TRUE, col="red", lty=2, lwd=2)

# Ajout d'une légende
legend("topright", legend=c("Numerical Solution", "Theoretical Solution"),
       col=c("black", "red"), lty=c(1, 2), lwd=c(1, 2))

################################################################################
#### Solve the forward Kolmogorov equation
## Initial condition for the SDE
x0 <- xi/4

## Initial condition for the FKE is a dirac
phi0 <- numeric(length(xc))
phi0[sum(xc<x0)] <- 1 

## Time grid
tv <- seq(0,10,0.1)

## Solve the Forward Kolmogorov Equation
PHI <- solve_forwardEquation(tv,phi0,G)

x11()
image.plot(tv,xc,t(PHI))

x11()
CDF <- apply(PHI*dx,2,cumsum)
image.plot(tv,xc,t(CDF))

x11()
plot(xc,PHI[,length(tv)])
lines(xc,phi)

################################################################################
#### Mean, variance and standard variation

EX <- computeMoment(PHI,1)
VX <- computeVariance(PHI)
sX <- computeStandardDeviation(PHI)

x11()
plot(tv,EX,lwd=2,ylim=c(0,max(EX+sX)))
lines(tv,EX+sX,lty="dashed")
lines(tv,EX-sX,lty="dashed")

# True stationnary 
EXa <- xi+(x0-xi)*exp(-lambda*tv)
lines(tv,EXa,col=2)
sXa <- sqrt(gamma^2*xi/2/lambda)
abline(h=xi,col=2,lty=2)
abline(h=xi+sXa,col=2,lty=2)
abline(h=xi-sXa,col=2,lty=2)

################################################################################
### Backward Equation
h <- (xc>=2)
T <- tail(tv,1)
psi <- solve_backwardEquation(tv,h,G)

x11()
image.plot(tv,xc,t(psi))