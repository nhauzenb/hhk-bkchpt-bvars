###############################################################################
###-------------------------- Auxiliary functions --------------------------### 
###-------------------  for non-conjugate BVAR function --------------------###
###############################################################################
#Function creating the lags of y (alternative to embed)
mlag <- function(X,lag)
{
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)  
}

# Get posteriors for the horseshoe prior (see Makalic & Schmidt, 2015)
get.hs.min <- function(bdraw,lam.hs,nu.hs,tau.hs,zeta.hs,update.ls = TRUE){
  k <- length(bdraw)
  # Local shrinkage scalings
  if(update.ls){
    lam.hs <- 1/rgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
    nu.hs <-  1/rgamma(k,shape=1,rate=1+1/lam.hs)
  }else{
    lam.hs <- lam.hs
    nu.hs <- nu.hs
  }
  # Global shrinkage parameter
  tau.hs  <- 1/rgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lam.hs)/2) 
  zeta.hs <- 1/rgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lam.hs*tau.hs),"lam"=lam.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

# Get posteriors for the NG prior
get.ng <- function(bdraw,psi.ng,a_tau){
  k <- length(bdraw)
  psi.ng <- matrix(0,k,1)
  # Global shrinkage parameter
  tau.ng <- rgamma(1,0.01 + a_tau*k,0.01 + (a_tau * sum(psi.ng))/2)
  # Local shrinkage scalings
  for(kk in 1:k){
    psi.ng[kk] <- GIGrvg::rgig(n=1, lambda=a_tau-0.5, chi= bdraw[kk]^2+1e-15, psi=tau.ng*a_tau)
  }
  return(list("psi"=psi.ng, "tau" = tau.ng))
}

###--------------------------------------------------------------------------###
###------------------- Non-conjugate BVAR function --------------------------###
###--------------------CCCM (2022, JoE) algorithm----------------------------###
###--------------------------------------------------------------------------###
nconjVAR.func <- function(data.setup = data.setup, 
                          mdl.setup  = mdl.setup, 
                          mcmc.setup = mcmc.setup){

list2env(data.setup, .GlobalEnv)
list2env(mdl.setup , .GlobalEnv)
list2env(mcmc.setup, .GlobalEnv)

ntot     <- nsave*nthin + nburn
save.set <- seq(nthin,nsave*nthin,nthin) + nburn
save.ind <- 0
  
###--------------------------------------------------------------------------###
###------------------------ Data set-up -------------------------------------###
###------------- conditional on the first p observations --------------------###
###--------------------------------------------------------------------------###
Yinit <- Y[1:p,,drop=F] # Initial observations 
X <- mlag(Y,p)
Y <- Y[(p+1):nrow(Y),]  
X <- X[(p+1):nrow(X),]
Ylbl <- colnames(Y)
Xlbl <-  paste0(rep(Ylbl,p),rep(paste0("_t-",1:p),each=M))
if(cons){
  X <- cbind(X,1) 
  Xlbl <- c(Xlbl, "cons")
}
colnames(X) <- Xlbl

M <- ncol(Y)   # No. of endogenous variables
N <- nrow(Y)   # No. of observations
K <- ncol(X)   # No. of coefficients per equation
Mp <- M*p      # No. of coefficients related to lags
v <- (M+1)*M/2 # No. of free elements in the variance-covariance matrix
k <- M*K       # k = (M*K) = (M*(M*p+1)) # No. of all coefficients of the system 
  
###--------------------------------------------------------------------------###
###----------------- Estimate a VAR(p) model with OLS -----------------------###
###--------------------------------------------------------------------------###
if(K > N){
  tXX <- crossprod(X) + 1e-3*diag(1,K)
  df <- N
}else{
  tXX <- crossprod(X) + 1e-6*diag(1,K)
  df <- (N-K+10)
}
  
A.ols     <- solve(tXX)%*%crossprod(X,Y)        # Reduced-form OLS coefficients
eps.ols   <- eps_draw <- Y - X %*%A.ols         # Reduced-form shocks
SSR.ols   <- crossprod(eps.ols)                 # SSR based on OLS
SIG.ols   <- SSR.ols/df                         # OLS variance-covariance matrix
L.inv.ols <- t(chol(SIG.ols))                   # Lower Cholesky factor of SIG.ols
L.inv.ols <- L.inv.ols*t(matrix(1/diag(L.inv.ols), M, M))
L.ols     <- solve(L.inv.ols)
B.ols     <- A.ols%*%t(L.ols)                   # Stuctural-form OLS coefficients
eta.ols   <- eta_draw <- eps.ols%*%t(L.ols)     # Stuctural-form shocks
  
OLS.est <- list("A" = A.ols,"SIG" = SIG.ols)

# Run a set of AR(p) models
s2.ARp <- matrix(NA,M,1)
for(mm in 1:M)
{
  ar.var <- try(arima(Y[,mm],order = c(p,0,0))$sigma2, silent = TRUE)
  if(is(ar.var, "try-error")){
    ar.cor  <- cor(Y[2:N,mm], Y[1:(N-1),mm])
    ar.var  <- var(Y[2:N,mm] - ar.cor*Y[1:(N-1),mm])
  } 
  s2.ARp[mm,1] <- ar.var
}
  
###--------------------------------------------------------------------------###
###------------------ Prior and initialization set-up -----------------------###
###--------------------------------------------------------------------------###
  
###----------------- Prior mean of VAR coefficients -------------------------###
prmean <- rep(pr.mean,M) 
A_pr <- matrix(0,K,M) 
diag(A_pr[1:M,1:M]) <- prmean
A_draw <- A.ols
  
###--------- Initialization of shrinkage prior on VAR coefficients ----------###
###----------------------- with Minnesota moments ---------------------------###
own.slct <- matrix(FALSE,K, M)
other.slct <- matrix(FALSE,K, M)
for(mm in 1:M){
  own.slct[seq(mm, M*p, by = M),mm]   <- TRUE       # identify own lags
  other.slct[setdiff(1:(M*p), seq(mm, M*p, by = M)),mm] <- TRUE # identify other lags
}  
  
# Global scalings (differentiate between own and other lags for Minnesota)
tau_A   <- rep(1, M)
tau.1_A <- rep(1, M)
tau.2_A <- rep(0.5, M)
zeta_A  <- zeta.1_A <- zeta.2_A <- rep(1, M)
  
# Local scalings accoriding to Minnesota moments 
lam_A <- matrix(1,K, M)
nu_A  <- matrix(1,K, M)
# Elements of prior variance-covariance matrix
psi_A <- theta_A  <- matrix(1,K,M)
  
for(mm in 1:M){
  own.slct.mm <- own.slct[,mm]
  other.slct.mm <- other.slct[,mm]
  lam_ai <- rep(1, K)
  lam_ai[1:(M*p)] <- 1/rep((1:p)^2, each = M) 
    
  if(cons){
    lam_ai <- lam_ai*s2.ARp[mm]/c(rep(s2.ARp, p),0.01^2) 
  }else{
    lam_ai <- lam_ai*s2.ARp[mm]/rep(s2.ARp, p)
  }
  lam_A[,mm] <- nu_A[,mm]   <- lam_ai
  theta_A[own.slct.mm,mm]   <- psi_A[own.slct.mm,mm]   <- tau.1_A[mm]*lam_A[own.slct.mm,mm]
  theta_A[other.slct.mm,mm] <- psi_A[other.slct.mm,mm] <- tau.2_A[mm]*lam_A[other.slct.mm,mm]
}
psi_A[K,] <- theta_A[K,] <- lam_A[K,]
  
###-------- Prior mean and variance of contemporaneous relationships ---------###
L_pr       <- matrix(0,M,M)
L_draw     <- L.ols
L.inv_draw <- L.inv.ols
  
psi_L   <- theta_L <- matrix(1,M,M)
lam_L   <- matrix(1, M, M)
nu_L    <- matrix(1, M, M)
tau_L   <- 1
zeta_L  <- 1
  
id.ind_L <- matrix(1:M^2, M, M)
id.slct_L  <- lower.tri(id.ind_L, diag = F)
  
###------------- Prior and initialization of SV processes -------------------###
# Prior of SV processes
sv_pr <- specify_priors(
  mu  = sv_normal(0,10),     # prior on unconditional mean in the state equation
  phi = sv_beta(shape1 = 25, shape2 = 1.5), #informative prior to push the model towards a random walk in the state equation (for comparability)
  sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*0.1)), # Gamma prior on the state innovation variance
  nu = sv_infinity(),
  rho = sv_constant(0)
)
  
# Initialization of SV processes
Ht_draw     <- matrix(0,N,M)
svpara_draw <- matrix(NA,3,M)
rownames(svpara_draw) <- c("mu", "phi", "sigma")
sv_draw<- list()
for (mm in 1:M) sv_draw[[mm]] <- list(mu = 0, phi = 0.95, sigma = 0.1, nu = Inf, rho = 0, beta = NA, latent0 = 0)
  
###--------------------------------------------------------------------------###
###---------------------------- Store objects -------------------------------###
###--------------------------------------------------------------------------###
A_store                <- array(NA, c(nsave,K,M))
dimnames(A_store)      <- list(1:nsave,Xlbl, Ylbl)
L_store                <- array(NA, c(nsave,M,M))
dimnames(L_store)      <- list(1:nsave,Ylbl, Ylbl)
L.inv_store            <- array(NA, c(nsave,M,M))
dimnames(L.inv_store)  <- list(1:nsave,Ylbl, Ylbl)
Ht_store               <- array(NA, c(nsave,N,M))
dimnames(Ht_store)     <- list(1:nsave,insmp.dates, Ylbl)
SIGt_store             <- array(NA, c(nsave,N,M,M))
dimnames(SIGt_store)   <- list(1:nsave,insmp.dates, Ylbl, Ylbl)
svpara_store           <- array(NA, c(nsave,3, M))
dimnames(svpara_store) <- list(1:nsave,c("mu", "phi", "sigma"), Ylbl)
tau_pr_store              <- array(NA, c(nsave,M, 2))
dimnames(tau_pr_store)    <- list(1:nsave,Ylbl,c("tau.1", "tau.2"))
lam_pr_store              <- array(NA, c(nsave,K, M))
dimnames(lam_pr_store)    <- list(1:nsave,Xlbl, Ylbl)
logDet_pr_store           <- array(NA, c(nsave, 1))
dimnames(logDet_pr_store) <- list(1:nsave, "logDet")

pb <- txtProgressBar(min = 0, max = ntot, style = 3)
irep <- 1
start <- Sys.time()

X.mm <- XX.mm  <- kronecker(X,rep(1,M))

for(irep in 1:ntot){
###------------------------------ Step 1: -----------------------------------###
###------------------- Sample contemp. relationships ------------------------###
###--------------------- and prior hyperparameters --------------------------###
###--------------------------------------------------------------------------###
eps_draw <- (Y-X%*%A_draw) # Reduced-form shocks
    
# Sample equation-specific contemp. relationships
for(mm in 2:M){
  norm.mm <- exp(-Ht_draw[,mm]/2)   
  eps.norm <- eps_draw*norm.mm
  Y.eps.mm <- eps.norm[,mm,drop = F]
  X.eps.mm <- eps.norm[,1:(mm-1),drop = F]
      
  l_pr.mm     <- as.numeric(L_pr[mm,1:(mm-1)])
  if(mm == 2){
    Vlinv_pr.mm <- as.numeric(1/theta_L[mm,1:(mm-1)])
    Vl_po.mm    <- 1/(sum((X.eps.mm)^2) + Vlinv_pr.mm)
    l_po.mm     <- Vl_po.mm*(sum(X.eps.mm*Y.eps.mm) + Vlinv_pr.mm*l_pr.mm)
    l_draw.mm   <- rnorm(1,l_po.mm, sqrt(Vl_po.mm))
  }else{
    Vlinv_pr.mm <- diag(1/theta_L[mm,1:(mm-1)])
    Vl_po.mm    <- solve(crossprod(X.eps.mm) + Vlinv_pr.mm)
    l_po.mm     <- Vl_po.mm%*%(crossprod(X.eps.mm,Y.eps.mm) + Vlinv_pr.mm%*%l_pr.mm)
    l_draw.mm   <- l_po.mm + t(chol(Vl_po.mm))%*%rnorm(mm-1)
  }
  L_draw[mm,1:(mm-1)] <- -l_draw.mm
}
    
# Update overall shrinkage on contemp. relationships
if(prior == "MIN"){
  hs_L <- get.hs.min(bdraw=as.numeric(L_draw[id.slct_L]-L_pr[id.slct_L]),lam.hs = lam_L[id.slct_L], nu.hs = nu_L[id.slct_L], tau.hs = tau_L, zeta.hs = zeta_L, update.ls = FALSE)
  psi_L[id.slct_L]          <- hs_L$psi
  tau_L <- hs_L$tau; zeta_L <- hs_L$zeta
}else if(prior == "HS"){
  hs_L <- get.hs.min(bdraw=as.numeric(L_draw[id.slct_L]-L_pr[id.slct_L]),lam.hs = lam_L[id.slct_L], nu.hs = nu_L[id.slct_L], tau.hs = tau_L, zeta.hs = zeta_L, update.ls = TRUE)
  psi_L[id.slct_L]          <- hs_L$psi
  lam_L[id.slct_L]          <- hs_L$lam 
  nu_L[id.slct_L]           <- hs_L$nu
  tau_L <- hs_L$tau; zeta_L <- hs_L$zeta
}else if(prior == "NG"){
  ng_L <- get.ng(bdraw=as.numeric(L_draw[id.slct_L]-L_pr[id.slct_L]), psi.ng = psi_L[id.slct_L], a_tau = 0.6)
  psi_L[id.slct_L]  <- ng_L$psi
}else if(prior == "lasso"){
  ng_L <- get.ng(bdraw=as.numeric(L_draw[id.slct_L]-L_pr[id.slct_L]), psi.ng = psi_L, a_tau = 1)
  psi_L[id.slct_L]  <- ng_L$psi
}
psi_L[psi_L < 1e-10] <- 1e-10
psi_L[psi_L > 1e2]   <- 1e2
    
theta_L <- psi_L
###------------------------------ Step 2: -----------------------------------###
###----------------------- Sample VAR coefficients --------------------------###
###--------------------- and prior hyperparameters --------------------------###
###--------------------------------------------------------------------------###
norm.vec <- exp(c(t(-Ht_draw))/2)
# Draw VAR coefficients equation by equation 
for(mm in 1:M){
  A_draw[,mm] <- 0              # Zero out coefficients that we want to update
  eps.mm      <- Y - X%*%A_draw # Equation-specific reduced-form shocks
  Y.mm        <- c(L_draw %*% t(eps.mm))*norm.vec
  #X.mm       <- kronecker(X,L_draw[,mm])*norm.vec
  X.mm[]      <- XX.mm*rep(L_draw[,mm], N)*norm.vec
      
  a_pr.mm    <- A_pr[,mm]
  Vinv_pr.mm <- diag(1/theta_A[,mm])
  V_po.mm    <- solve(crossprod(X.mm) + Vinv_pr.mm)
  a_po.mm    <- V_po.mm%*%(crossprod(X.mm,Y.mm) + Vinv_pr.mm%*%a_pr.mm)
  a_draw.mm  <- a_po.mm + t(chol(V_po.mm))%*%rnorm(K)
      
  A_draw[,mm] <- a_draw.mm
  # Update equation-specific shrinkage on VAR coefficients 
  if(prior == "MIN"){
    own.slct.mm <- own.slct[,mm]
    other.slct.mm <- other.slct[,mm]
        
    min_A.own   <- get.hs.min(bdraw     = (a_draw.mm[own.slct.mm] - a_pr.mm[own.slct.mm]),
                              lam.hs    = lam_A[own.slct.mm,mm], 
                              nu.hs     = nu_A[own.slct.mm ,mm], 
                              tau.hs    = tau.1_A[mm], 
                              zeta.hs   = zeta.1_A[mm], 
                              update.ls = FALSE)
    psi_A[own.slct.mm,mm] <- min_A.own$psi
    tau.1_A[mm]           <- min_A.own$tau; zeta.1_A[mm] <- min_A.own$zeta
        
    min_A.other   <- get.hs.min(bdraw     = (a_draw.mm[other.slct.mm] - a_pr.mm[other.slct.mm]),
                                lam.hs    = lam_A[other.slct.mm,mm], 
                                nu.hs     = nu_A[other.slct.mm ,mm], 
                                tau.hs    = tau.2_A[mm], 
                                zeta.hs   = zeta.2_A[mm], 
                                update.ls = FALSE)
    psi_A[other.slct.mm,mm] <- min_A.other$psi
    tau.2_A[mm]             <- min_A.other$tau; zeta.2_A[mm] <- min_A.other$zeta
  }else if(prior == "HS"){
    hs_A <- get.hs.min(bdraw     = (a_draw.mm - a_pr.mm),
                       lam.hs    = lam_A[,mm], 
                       nu.hs     = nu_A[,mm], 
                       tau.hs    = tau_A[mm], 
                       zeta.hs   = zeta_A[mm], 
                       update.ls = TRUE)
    psi_A[,mm] <- hs_A$psi
    lam_A[,mm] <- hs_A$lam; nu_A[,mm]  <- hs_A$nu
    tau_A[mm]  <- hs_A$tau; zeta_A[mm] <- hs_A$zeta
  }else if(prior == "NG"){
    ng_A <- get.ng(bdraw  = (a_draw.mm - a_pr.mm), 
                   psi.ng = psi_A[,mm], 
                   a_tau  = 0.1)
    psi_A[,mm] <- ng_A$psi
    lam_A[,mm] <- ng_A$psi*(ng_A$tau/2)
    tau_A[mm]  <- 2/ng_A$tau
  }else if(prior == "lasso"){
    ng_A       <- get.ng(bdraw  = (a_draw.mm - a_pr.mm), 
                         psi.ng = psi_A[,mm], 
                         a_tau  = 1)
    psi_A[,mm] <- ng_A$psi
    lam_A[,mm] <- ng_A$psi*(ng_A$tau/2)
    tau_A[mm]  <- 2/ng_A$tau
  }
}
    
psi_A[psi_A < 1e-10] <- 1e-10
psi_A[psi_A > 1e2]   <- 1e2

theta_A <- psi_A
    
###------------------------------ Step 3: -----------------------------------###
###------------- Sample structural error variance processes -----------------###
###--------------------------------------------------------------------------###
eta_draw <- (Y - X%*%A_draw)%*%t(L_draw) # Structural-form shocks
    
if(sv == "homo"){
  for(mm in 1:M){
    eht.mm <- 1/rgamma(1, 0.01 + N/2, 0.01+ sum(eta_draw[,mm]^2)/2)
    Ht_draw[,mm]     <- log(eht.mm)  
    svpara_draw[,mm] <- c(log(eht.mm), 0,0)
  }
}else if(sv == "SV"){
  for(mm in 1:M){
    sv_mm <- svsample_general_cpp(eta_draw[,mm], startpara = sv_draw[[mm]], startlatent = Ht_draw[,mm], priorspec = sv_pr)
    svpara_mm <- sv_mm$para[, c("mu", "phi","sigma")]
    sv_draw[[mm]][c("mu", "phi","sigma")] <- as.list(svpara_mm)
        
    svpara_draw[,mm] <- svpara_mm
    Ht_draw[,mm]     <- as.numeric(sv_mm$latent)
  }
}
    
###------------------------------ Step 4: -----------------------------------###
###----------------------------- Storage  -----------------------------------###
###--------------------------------------------------------------------------###
if(irep %in% save.set)
{
  save.ind <- save.ind + 1
  A_store[save.ind,,]    <- A_draw
  L_store[save.ind,,]    <- L_draw
  Ht_store[save.ind,,]   <- Ht_draw
  svpara_store[save.ind,,] <- svpara_draw

  L.inv_draw <- solve(L_draw)
  L.inv_store[save.ind,,] <- L.inv_draw
      
  SIGt_draw <- array(NA,c(N,M,M))
  for(tt in 1:N){
    SIGt_draw[tt,,] <- L.inv_draw %*% diag(exp(Ht_draw[tt,])) %*% t(L.inv_draw)
  }
  SIGt_store[save.ind,,,] <- SIGt_draw
  
  if(prior %in% c("MIN", "MINh")){
    tau_pr_store[save.ind,,"tau.1"]  <-  tau.1_A
    tau_pr_store[save.ind,,"tau.2"]  <-  tau.2_A
  }else{
    tau_pr_store[save.ind,,]  <-  tau_A
  }
  lam_pr_store[save.ind,,]  <- lam_A 
  logDet_pr_store[save.ind,] <- sum(log(theta_A))
}
setTxtProgressBar(pb, irep)
}
end <- Sys.time()
# Time in minutes
est.time        <- c(as.numeric((ts(end)-ts(start))/60), ntot)
names(est.time) <- c("Minutes", "Draws")
  
est.mod <- list("A"       = A_store,
                "L"       = L_store, 
                "L.inv"   = L.inv_store, 
                "Ht"      = Ht_store, 
                "svpara"  = svpara_store, 
                "SIGt"    = SIGt_store, 
                
                "tau_pr"    = tau_pr_store,
                "lam_pr"    = lam_pr_store,
                "logDet_pr" = logDet_pr_store, 
                
                "OLS.est"  = OLS.est, 
                "est.time" = est.time, 
                
                "Ylbl"     = Ylbl,
                "Xlbl"     = Xlbl, 
                "X"        = X, 
                "Y"        = Y,
                "N"        = N,
                "Mp"       = Mp,
                "K"        = K)
}
