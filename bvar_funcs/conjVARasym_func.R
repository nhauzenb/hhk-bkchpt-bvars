###############################################################################
###-------------------------- Auxiliary functions --------------------------### 
###-------------------- for asymmetric conjugate BVAR ----------------------### 
###############################################################################
###-------- Function creating the lags of y (alternative to embed) ---------###
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

###----------------------- Moments of gamma distribution --------------------###
gammacoef <- function(mode, sd){
  k.shape <- (2+mode^2/sd^2+sqrt((4+mode^2/sd^2)*mode^2/sd^2))/2
  tau.scale <- sqrt(sd^2/k.shape)
  return(data.frame(shape = k.shape, scale = tau.scale))
}

###---- Get standard posterior moments of a regression for m-th equation ----###
get.post.mom.mm <- function(y.mm = y.mm, x.mm = x.mm, N = N, pr_mm.list){
  bl_pr.mm     <- pr_mm.list[["bl_pr.mm"]] 
  v_pr.mm      <- pr_mm.list[["v_pr.mm"]] 
  v_pr.inv.mm  <- pr_mm.list[["v_pr.inv.mm"]] 
  s_pr.mm      <- pr_mm.list[["s_pr.mm"]] 
  S_pr.mm      <- pr_mm.list[["S_pr.mm"]] 
  
  # Posterior moments
  V_po.inv.mm     <- (crossprod(x.mm) + diag(v_pr.inv.mm)) # Posterior precision
  Vchol_po.inv.mm <- t(chol(V_po.inv.mm))                  # Cholesky of posterior precision
  V_po.mm         <- solve(V_po.inv.mm)                    # Posterior variance-covariance
  Vchol_po.mm     <- t(chol(V_po.mm))                      # Cholesky of posterior variance-covariance
  bl_po.mm        <- V_po.mm%*%(crossprod(x.mm, y.mm) + diag(v_pr.inv.mm)%*%bl_pr.mm) # Posterior mean
  
  S_po.mm         <- (S_pr.mm + crossprod(y.mm) + t(bl_pr.mm)%*%diag(v_pr.inv.mm)%*%bl_pr.mm -  t(bl_po.mm)%*%V_po.inv.mm%*%bl_po.mm)/2
  s_po.mm         <- (N + s_pr.mm)/2
  
  return(list("bl_po.mm"        = bl_po.mm, 
              "V_po.mm"         = V_po.mm, 
              "Vchol_po.mm"     = Vchol_po.mm,
              "V_po.inv.mm"     = V_po.inv.mm,
              "Vchol_po.inv.mm" = Vchol_po.inv.mm,
              "s_po.mm"         = s_po.mm, 
              "S_po.mm"         = S_po.mm
  ))
}

###------- Marginal likelihood (simplyfied: without constant terms) ---------###
get.logML.mm <- function(po_mm.list = po_mm.list, pr_mm.list = pr_mm.list){
  v_pr.mm         <- pr_mm.list[["v_pr.mm"]]
  s_pr.mm         <- pr_mm.list[["s_pr.mm"]]
  S_pr.mm         <- pr_mm.list[["S_pr.mm"]]
  
  Vchol_po.inv.mm <- po_mm.list[["Vchol_po.inv.mm"]] 
  S_po.mm         <- po_mm.list[["S_po.mm"]]
  s_po.mm         <- po_mm.list[["s_po.mm"]]
  
  # Evaluate marginal likelihood
  cons.logML.mm <- 0#-N/2*log(2*pi) + log(gamma(s_po.mm[jj])) - log(gamma(s_pr.mm)) + s_pr.mm*log(S_pr.mm)
  logML.mm <- cons.logML.mm - s_po.mm*log(S_po.mm) - 1/2*(sum(log(v_pr.mm)) +  2*sum(log(diag(Vchol_po.inv.mm))))
  
  if (is.infinite(logML.mm)) logML.mm <- -10^10
  return(logML.mm)
}  
  
###---------------- Objective function for ML maximization ------------------###
obj_func.asym <- function(Y = Y, X = X, own.slct = own.slct, other.slct = other.slct, N = N, M = M, K = K,
                        lam_B = lam_B, psi_L = psi_L, B_pr = B_pr, L_pr = L_pr, s_pr = s_pr, S_pr = S_pr, 
                        tau = tau, tau.1_pr = tau.1_pr, tau.2_pr = tau.2_pr){
  
  tau.1 <- tau[1]
  tau.2 <- tau[2]
  logML <- c()
  for (mm in 1:M){
    # Set-up for m-th equation
    own.slct.mm   <- own.slct[,mm]
    other.slct.mm <- other.slct[,mm]
    K.mm          <- K + (mm - 1) # Get dimension of m-th equation 
    y.mm          <- Y[,mm]
    
    # Prior variance-covariance  elements of m-th equation that change with respect to tau
    psi_B.mm <- lam_B.mm    <- lam_B[,mm] 
    psi_B.mm[own.slct.mm]   <- tau.1*lam_B.mm[own.slct.mm]
    psi_B.mm[other.slct.mm] <- tau.2*lam_B.mm[other.slct.mm] 
    
    if(mm == 1){ 
      x.mm     <- X            # Design matrices for the 1-st equation
      bl_pr.mm <- c(B_pr[,mm]) # Prior mean for 1-st equation
      psi_L.mm <- NULL         # For 1-st equation we do not have contemp. relations 
    }else{
      x.mm     <- cbind(X, -Y[,1:(mm-1)])         # Design matrices for the m-th equation
      bl_pr.mm <- c(B_pr[,mm], L_pr[mm,1:(mm-1)]) # Prior mean for m-th equation
      psi_L.mm <- psi_L[mm,1:(mm-1)]              # Prior variances for (m-1)-th contemp. relations 
    }
    
    v_pr.mm <- c(psi_B.mm, psi_L.mm) # Prior variances for m-th equation
    v_pr.inv.mm  <- 1/v_pr.mm
    # Prior DoF and scaling for m-th equation 
    s_pr.mm <-  s_pr + mm - M
    S_pr.mm  <- S_pr[mm,mm]
    
    pr_mm.list <- list(
      "bl_pr.mm"     = bl_pr.mm, 
      "v_pr.mm"     = v_pr.mm, 
      "v_pr.inv.mm" = v_pr.inv.mm,
      "s_pr.mm"     = s_pr.mm, 
      "S_pr.mm"     = S_pr.mm)
    
    po_mm.list <- get.post.mom.mm(y.mm = y.mm, x.mm = x.mm, N = N, pr_mm.list = pr_mm.list)

    logML[mm] <- get.logML.mm(po_mm.list = po_mm.list, pr_mm.list = pr_mm.list)
    
  }
  
  logML <- sum(logML)
  logPR  <- dgamma(tau.1, shape = tau.1_pr[[1]], scale = tau.1_pr[[2]], log = TRUE) + dgamma(tau.2, shape = tau.2_pr[[1]], scale = tau.2_pr[[2]], log = TRUE)
  
  obj_value <- -(logML + logPR)
  return(obj_value)

}

###--------------------------------------------------------------------------###
###---------------- Asymmetric conjugate BVAR function ----------------------###
###--------------------------------------------------------------------------###
conjVARasym.func <- function(data.setup = data.setup, 
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
Yinit <- Y[1:p,,drop=F] # Initial observations (e.g., for SOC and SUR prior)
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
Mp <- M*p       # No. of coefficients related to lags
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
  
A.ols     <- solve(tXX)%*%crossprod(X,Y)          # Reduced-form OLS coefficients
eps.ols   <- eps_draw <- Y - X %*%A.ols           # Reduced-form shocks
SSR.ols   <- crossprod(eps.ols)                   # SSR based on OLS
SIG.ols   <- SSR.ols/df                           # OLS variance-covariance matrix
L.inv.ols <- t(chol(SIG.ols))                     # Lower Cholesky factor of SIG.ols
L.inv.ols <- L.inv.ols*t(matrix(1/diag(L.inv.ols), M, M))
L.ols     <- solve(L.inv.ols)
B.ols     <- A.ols%*%t(L.ols)                     # Stuctural-form OLS coefficients
eta.ols   <- eta_draw <- eps.ols%*%t(L.ols)       # Stuctural-form shocks
  
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
###--------------------------- Prior set-up ---------------------------------###
###------------ consistent with asymmetric Minnesota prior ------------------###
###--------------------------------------------------------------------------###

###---------------- Identifiers for own and other lags ----------------------###
own.slct <- matrix(FALSE,K, M)
other.slct <- matrix(FALSE,K, M)
for(mm in 1:M){
  own.slct[seq(mm, M*p, by = M),mm]   <- TRUE       # identify own lags
  other.slct[setdiff(1:(M*p), seq(mm, M*p, by = M)),mm] <- TRUE # identify other lags
}

###------------- Prior mean of structural VAR coefficients ------------------###
prmean <- rep(pr.mean,M) 
B_pr <- matrix(0,K,M) 
diag(B_pr[1:M,1:M]) <- prmean

###------------- Local scalings according to Minnesota moments --------------### 
###------------- Elements of prior variance-covariance matrix ---------------### 
lam_B             <- matrix(1,K, M)
lam_B.mm          <- rep(1, K)
lam_B.mm[1:(M*p)] <- 1/rep((1:p)^2, each = M)
if(cons){
  lam_B.mm        <- lam_B.mm/c(rep(s2.ARp, p),0.01^2) 
}else{
  lam_B.mm        <- lam_B.mm/rep(s2.ARp, p)
}
lam_B[] <-  lam_B.mm

###-------- Prior mean and variance of contemporaneous relationships ---------###
L_pr       <- matrix(0,M,M)
L_draw     <- L.ols

tau_L   <- 1 # Fixed global shrinkage parameter for contemporaneous relations
psi_L   <- matrix(0,M,M) 
lam_L   <- matrix(0, M, M)
for(mm in 2:M){
  lam_L[mm,1:(mm-1)] <- c(1/s2.ARp[1:(mm-1)]) 
  psi_L[mm,1:(mm-1)] <-  tau_L*lam_L[mm,1:(mm-1)] 
}  

###----------- Prior moments of structural error variances ------------------###
s_pr <-  M + 2
S_pr <- diag(as.numeric(s2.ARp))

###------- Global scalings (differentiate between own and other lags) -------###
tau.1_pr <- gammacoef(mode = 1, sd = 0.4) 
tau.2_pr <- gammacoef(mode = 1, sd = 0.4) 
tau.1_lb <- tau.2_lb <- 1/k 
tau.1_ub <- tau.2_ub <- 10/M
###--------------------------------------------------------------------------###
###-------------------------- Posterior set-up ------------------------------###
###--- for changing global scalings of Minnesota based on a Griddy-Gibbs ----###
###--------------------------------------------------------------------------###
if(tau_gs > 1){
tau.1_gs <- tau.2_gs <- sqrt(tau_gs) 
tau.1_grid <- seq(tau.1_lb, tau.1_ub, length.out = tau.1_gs)
tau.2_grid <- seq(tau.2_lb, tau.2_ub, length.out = tau.2_gs)
tau.2D_grid <- expand.grid(tau.1_grid, tau.2_grid)
colnames(tau.2D_grid) <- c("tau.1", "tau.2")
  
tau.2D_grid$logPR <- dgamma(tau.2D_grid[,"tau.1"], shape = tau.1_pr[[1]], scale = tau.1_pr[[2]], log = TRUE) + dgamma(tau.2D_grid[,"tau.2"], shape = tau.2_pr[[1]], scale = tau.2_pr[[2]], log = TRUE)

# Grid-specific posterior store objects for m-th equation
po_grid.list <- vector("list", tau_gs)
logML        <- rep(0, tau_gs)
logDet_V.pr  <- rep(0, tau_gs)
jj <- mm <- 1
for(jj in 1:tau_gs){
# Grid-specific posterior store objects for j-th combination
po.jj_grid.list <- vector("list", M); names(po_grid.list) <- var.slct
logML.mm <- logDet_V.pr.mm <- c()
# Variance-covariance elements of j-th combination
for (mm in 1:M){
  # Set-up for m-th equation
  own.slct.mm       <- own.slct[,mm]
  other.slct.mm     <- other.slct[,mm]
  K.mm              <- K + (mm - 1) # Get dimension of m-th equation 
  y.mm <- Y[,mm]
  
  # Prior variance-covariance  elements of m-th equation that change with respect to tau
  psi_B.jj.mm <- lam_B.mm    <- lam_B[,mm] 
  psi_B.jj.mm[own.slct.mm]   <- tau.2D_grid[jj,"tau.1"]*lam_B.mm[own.slct.mm]
  psi_B.jj.mm[other.slct.mm] <- tau.2D_grid[jj,"tau.2"]*lam_B.mm[other.slct.mm] 
  
  if(mm == 1){ 
    x.mm <- X    # Design matrices for the 1-st equation
    bl_pr.mm <- c(B_pr[,mm]) # Prior mean for 1-st equation
    psi_L.mm   <- NULL      # For 1-st equation we do not have contemp. relations 
  }else{
    x.mm <- cbind(X, -Y[,1:(mm-1)])  # Design matrices for the m-th equation
    bl_pr.mm <- c(B_pr[,mm], L_pr[mm,1:(mm-1)]) # Prior mean for m-th equation
    psi_L.mm   <- psi_L[mm,1:(mm-1)] # Prior variances for (m-1)-th contemp. relations 
  }
  
  v_pr.jj.mm <- c(psi_B.jj.mm, psi_L.mm) # Prior variances for m-th equation
  v_pr.inv.jj.mm  <- 1/v_pr.jj.mm
  # Prior DoF and scaling for m-th equation 
  s_pr.mm <-  s_pr + mm - M
  S_pr.mm  <- S_pr[mm,mm]
  
  pr_jj.mm.list <- list(
    "bl_pr.mm"     = bl_pr.mm, 
    "v_pr.mm"     = v_pr.jj.mm, 
    "v_pr.inv.mm" = v_pr.inv.jj.mm,
    "s_pr.mm"     = s_pr.mm, 
    "S_pr.mm"     = S_pr.mm)
  
  po_jj.mm.list <- get.post.mom.mm(y.mm = y.mm, x.mm = x.mm, N = N, pr_jj.mm.list)
  po.jj_grid.list[[mm]] <- po_jj.mm.list
  
  logML.mm[mm]     <- get.logML.mm(po_mm.list = po_jj.mm.list, pr_mm.list = pr_jj.mm.list)
  logDet_V.pr.mm[mm] <- sum(log(v_pr.jj.mm))
  
}   
 
po_grid.list[[jj]] <- po.jj_grid.list
logML[[jj]]        <- sum(logML.mm)
logDet_V.pr[[jj]]    <- mean(logDet_V.pr.mm) 
}
 
tau.2D_grid$logML     <- logML
tau.2D_grid$logPO     <- tau.2D_grid$logPR + tau.2D_grid$logML 
tau.2D_grid$logDet_V.pr <- logDet_V.pr

# Define weights of hyperparameter combinations according to ML 
tau.2D_grid[,"logPOresc"] <- tau.2D_grid[,"logPO"] - max(tau.2D_grid[,"logPO"])
tau.2D_grid[,"weight"] <- exp(tau.2D_grid[,"logPOresc"])/sum(exp(tau.2D_grid[,"logPOresc"])) 
tau.1_map <- tau.2D_grid[which(tau.2D_grid[,"weight"] == max(tau.2D_grid[,"weight"])),"tau.1"]
tau.2_map <- tau.2D_grid[which(tau.2D_grid[,"weight"] == max(tau.2D_grid[,"weight"])),"tau.2"]
  
}else{
  if(tau_gs == 0){ # Flat prior
      tau_optim <- list("par" = c(100,100), "value" = NA)
  }else{ # Optimize ML 
      tau_optim <- optim(par=c(0.5, 0.5), fn = obj_func.asym, Y = Y, X = X, own.slct = own.slct, other.slct = other.slct, N = N, M = M, K = K, 
                         lam_B = lam_B, psi_L = psi_L, B_pr = B_pr, L_pr = L_pr, s_pr  = s_pr, S_pr  = S_pr, 
                         tau.1_pr = tau.1_pr, tau.2_pr = tau.2_pr,
                         lower= c(tau.1_lb,tau.2_lb), upper=c(tau.1_ub, tau.2_ub), method="L-BFGS-B", hessian = FALSE)
  }
  tau_map   <- tau_optim$par
  tau.1_map <- tau_map[1]
  tau.2_map <- tau_map[2]
  logPO <- -tau_optim$value
  logPR <- dgamma(tau.1_map, shape = tau.1_pr[[1]], scale = tau.1_pr[[2]], log = TRUE) + dgamma(tau.2_map, shape = tau.2_pr[[1]], scale = tau.2_pr[[2]], log = TRUE)

  tau.2D_grid <- data.frame("tau.1" = tau.1_map, "tau.2" = tau.2_map, "logML" = logPO-logPR, "logPR" = logPR, "logPO" = logPO, "weight" = 1)

  # Grid-specific posterior store objects for j-th combination
  po_grid.list <- vector("list", M); names(po_grid.list) <- var.slct
  logML.mm <- logDet_V.pr.mm <- c()
  for (mm in 1:M){
    # Set-up for m-th equation
    own.slct.mm       <- own.slct[,mm]
    other.slct.mm     <- other.slct[,mm]
    K.mm              <- K + (mm - 1) # Get dimension of m-th equation 
    y.mm <- Y[,mm]
    
    # Prior variance-covariance  elements of m-th equation that change with respect to tau
    psi_B.mm <- lam_B.mm <- lam_B[,mm] 
    psi_B.mm[own.slct.mm]   <- tau.1_map*lam_B.mm[own.slct.mm]
    psi_B.mm[other.slct.mm] <- tau.2_map*lam_B.mm[other.slct.mm] 
    
    if(mm == 1){ 
      x.mm <- X    # Design matrices for the 1-st equation
      bl_pr.mm <- c(B_pr[,mm]) # Prior mean for 1-st equation
      psi_L.mm   <- NULL      # For 1-st equation we do not have contemp. relations 
    }else{
      x.mm <- cbind(X, -Y[,1:(mm-1)])  # Design matrices for the m-th equation
      bl_pr.mm <- c(B_pr[,mm], L_pr[mm,1:(mm-1)]) # Prior mean for m-th equation
      psi_L.mm   <- psi_L[mm,1:(mm-1)] # Prior variances for (m-1)-th contemp. relations 
    }
    
    v_pr.mm <- c(psi_B.mm, psi_L.mm) # Prior variances for m-th equation
    v_pr.inv.mm  <- 1/v_pr.mm
    # Prior DoF and scaling for m-th equation 
    s_pr.mm <-  s_pr + mm - M
    S_pr.mm  <- S_pr[mm,mm]
    
    pr_mm.list <- list(
      "bl_pr.mm"     = bl_pr.mm, 
      "v_pr.mm"     = v_pr.mm, 
      "v_pr.inv.mm" = v_pr.inv.mm,
      "s_pr.mm"     = s_pr.mm, 
      "S_pr.mm"     = S_pr.mm)
    
    po_mm.list         <- get.post.mom.mm(y.mm = y.mm, x.mm = x.mm, N = N, pr_mm.list = pr_mm.list)
    po_grid.list[[mm]] <- po_mm.list
    
    logML.mm[mm] <- get.logML.mm(po_mm.list = po_mm.list, pr_mm.list = pr_mm.list)
    logDet_V.pr.mm[mm] <- sum(log(v_pr.mm))
    
  }
  po_grid.list <- list(po_grid.list)
  
  tau.2D_grid$logDet_V.pr <- mean(logDet_V.pr.mm) 
  tau.2D_grid$logML <- sum(logML.mm)
  tau.2D_grid$logPO <- tau.2D_grid$logPR + tau.2D_grid$logML 
}
  
###--------------------------------------------------------------------------###
###---------------------------- Store objects -------------------------------###
###--------------------------------------------------------------------------###
B_draw  <- array(0, c(K, M))
L_draw  <- diag(M)
H_draw  <- c(1:M)

A_store               <- array(NA, c(nsave,K,M))
B_store               <- array(NA, c(nsave,K,M))
dimnames(A_store)     <- dimnames(B_store) <- list(1:nsave,Xlbl, Ylbl)
L_store               <- array(NA, c(nsave,M,M))
L.inv_store           <- array(NA, c(nsave,M,M))
L_store               <- array(NA, c(nsave,M,M))
dimnames(L_store)     <- list(1:nsave,Ylbl, Ylbl)
L.inv_store           <- array(NA, c(nsave,M,M))
dimnames(L.inv_store) <- list(1:nsave,Ylbl, Ylbl)
H_store               <- array(NA, c(nsave,M))
dimnames(H_store)     <- list(1:nsave,Ylbl)
SIG_store             <- array(NA, c(nsave,M,M))
dimnames(SIG_store)   <- list(1:nsave,Ylbl, Ylbl)
tau_pr_store           <- array(NA, c(nsave,2))
dimnames(tau_pr_store) <- list(1:nsave,c("tau.1", "tau.2"))
logDet_pr_store           <- array(NA, c(nsave, 1))
dimnames(logDet_pr_store) <- list(1:nsave, "logDet")

pb <- txtProgressBar(min = 0, max = ntot, style = 3)
irep <- 1
start <- Sys.time()

for(irep in 1:ntot){
###------------------------------ Step 1: -----------------------------------###
###---- Select global shrinkage parameters according to posterior weights ---### 
###--------------------------------------------------------------------------###
if(tau_gs == 0) tau_ind  <- 1 else tau_ind <- sample(1:tau_gs, size=1, prob=tau.2D_grid[,"weight"])
tau_pr_draw      <- tau.2D_grid[tau_ind,c("tau.1","tau.2")]  
logDet_V.pr_draw <- tau.2D_grid[tau_ind,"logDet_V.pr"]  
logDet_pr_draw   <- 0

po_draw     <- po_grid.list[[tau_ind]]
###----------------- Draw structural VAR parameters equation by equation -------------------### 
for (mm in 1:M){
  K.mm <- K + (mm - 1)
  po.mm_draw <- po_draw[[mm]]
    
###------------------------------ Step 2: -----------------------------------###
###------------ Draw sigma2_m marginally from a inverse Gamma --------------### 
###--------------------------------------------------------------------------###
  s_po.mm     <- po.mm_draw[["s_po.mm"]]
  S_po.mm     <- po.mm_draw[["S_po.mm"]]
  h.mm <- 1/rgamma(1, s_po.mm,S_po.mm)
  H_draw[mm] <- h.mm
  h.mm <- (diag(SIG.ols)^2)[mm]
###------------------------------ Step 3: -----------------------------------###
###---- Draw equation-specific coefficients and covariances ----###
###-------  from the conditional posterior (a_m,l_m)|sigma2_m,DATA ----------###
###--------------------------------------------------------------------------###
  Vchol_po.mm <- po.mm_draw[["Vchol_po.mm"]]
  bl_po.mm     <- po.mm_draw[["bl_po.mm"]]
  
  bl.mm <- bl_po.mm + sqrt(h.mm)*Vchol_po.mm%*%rnorm(K.mm)
  
  if(mm == 1){
      B_draw[,mm] <- bl.mm
    }else{
      B_draw[,mm] <- bl.mm[1:K]
      L_draw[mm,1:(mm-1)] <- bl.mm[(K+1):K.mm]
    }
  
  logDet_pr_draw <- logDet_pr_draw + (logDet_V.pr_draw + K.mm*log(h.mm))
}
 
###----------------- Obtain reduced-form VAR coefficients -------------------### 
L.inv_draw <- solve(L_draw)
A_draw     <- B_draw%*%t(L.inv_draw)
SIG_draw   <- L.inv_draw%*%diag(H_draw)%*%t(L.inv_draw)

###------------------------------ Step 4: -----------------------------------###
###----------------------------- Storage  -----------------------------------###
###--------------------------------------------------------------------------###
if(irep %in% save.set)
{
  save.ind <- save.ind + 1
  A_store[save.ind,,]     <- A_draw
  B_store[save.ind,,]     <- B_draw
  H_store[save.ind,]      <- H_draw
  L_store[save.ind,,]     <- L_draw
  L.inv_store[save.ind,,] <- L.inv_draw
  SIG_store[save.ind,,]   <- SIG_draw
  
  tau_pr_store[save.ind,]    <- as.numeric(tau_pr_draw)
  logDet_pr_store[save.ind,] <- logDet_pr_draw
} 
setTxtProgressBar(pb, irep)
}
end <- Sys.time()
# Time in minutes
est.time        <- c(as.numeric((ts(end)-ts(start))/60), ntot)
names(est.time) <- c("Minutes", "Draws")

est.mod <- list("A"        = A_store,
                "B"        = B_store,
                "H"        = H_store,
                "L"        = L_store, 
                "L.inv"    = L.inv_store, 
                "SIG"      = SIG_store, 
                
                "tau_pr"     = tau_pr_store,
                "logDet_pr"  = logDet_pr_store, 
                
                "OLS.est"  = OLS.est, 
                "est.time" = est.time, 
                
                "Ylbl"     = Ylbl,
                "Xlbl"     = Xlbl, 
                "X"        = X, 
                "Y"        = Y,
                "N"        = N,
                "Mp"       = M*p,
                "K"        = K)
}
