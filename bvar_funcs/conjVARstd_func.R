###############################################################################
###-------------------------- Auxiliary functions --------------------------### 
###--------------------------- for conjugate VAR ---------------------------### 
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

###----------- Get standard posterior moments of a conjugate VAR ------------###
get.post.mom <- function(X = X, Y = Y, N = N, pr_list = pr_list){
  A_pr      <- pr_list[["A_pr"]] 
  v_pr      <- pr_list[["v_pr"]] 
  v_pr.inv  <- pr_list[["v_pr.inv"]] 
  s_pr      <- pr_list[["s_pr"]] 
  S_pr      <- pr_list[["S_pr"]] 
  
  # Posterior moments
  V_po.inv     <- (crossprod(X) + diag(v_pr.inv)) # Posterior precision
  Vchol_po.inv <- t(chol(V_po.inv))               # Cholesky of posterior precision
  V_po         <- solve(V_po.inv)                 # Posterior variance-covariance
  V_po.chol    <- t(chol(V_po))                   # Cholesky of posterior variance-covariance
  A_po         <- V_po%*%(crossprod(X, Y) + diag(v_pr.inv)%*%A_pr) # Posterior mean
  a_po.vec     <- as.vector(A_po)                 # Vectorized coefficients
  
  S_po         <- S_pr + crossprod(Y) + t(A_pr)%*%diag(v_pr.inv)%*%A_pr -  t(A_po)%*%V_po.inv%*%A_po # Posterior scaling x DoF
  S_po.inv     <- solve(S_po)                     # Inverse of posterior scaling x DoF 
  s_po         <- N + s_pr                        # Posterior DoF
  
  return(list("A_po"         = A_po, 
              "a_po.vec"     = a_po.vec, 
              "V_po"         = V_po, 
              "V_po.chol"    = V_po.chol,
              "V_po.inv"     = V_po.inv, 
              "Vchol_po.inv" = Vchol_po.inv,
              "s_po"         = s_po, 
              "S_po"         = S_po,
              "S_po.inv"     = S_po.inv))
}    

###------- Marginal likelihood (simplyfied: without constant terms) ---------###
get.logML <- function(po_list = po_list, pr_list = pr_list, M = M){
  V_pr     <- diag(pr_list[["v_pr"]])
  
  Vchol_po.inv <- po_list[["Vchol_po.inv"]]
  S_po     <- po_list[["S_po"]]
  s_po     <- po_list[["s_po"]]
  
  # Ingredients for marginal likelihood
  S_po.det <- 2*sum(log(diag(chol(S_po))))
  V_po.det <- 2*sum(log(diag(Vchol_po.inv)))
  V_pr.det <- 2*sum(log(diag(chol(V_pr))))
  
  # Evaluate marginal likelihood
  logML <- -M/2*(V_pr.det + V_po.det) + (-s_po/2)*S_po.det
  
  if (is.infinite(logML)) logML <- -10^10
  return(logML)
}  
  
###---------------- Objective function for ML maximization ------------------###
obj_func <- function(Y = Y, X = X, own.slct = own.slct, other.slct = other.slct, N = N, M = M, Mp = Mp,
                     lam_A = lam_A, A_pr = A_pr, s_pr = s_pr, S_pr = S_pr, 
                     tau = tau, tau_pr = tau_pr){
  v_pr <- lam_A
  v_pr[1:Mp] <- tau*v_pr[1:Mp] 
  v_pr.inv  <- 1/v_pr
  
  pr_list <- list(
    "A_pr"      = A_pr, 
    "v_pr"      = v_pr, 
    "v_pr.inv"  = v_pr.inv,
    "s_pr"      = s_pr, 
    "S_pr"      = S_pr)
  
  # Get posterior moments
  po_list <- get.post.mom(X = X, Y = Y, N = N, pr_list = pr_list)
  # Compute marginal likelihood
  logML <- get.logML(po_list = po_list, pr_list = pr_list, M = M)
  logPR <- dgamma(tau, shape = tau_pr[[1]], scale = tau_pr[[2]], log = TRUE)
  
  obj_value <- -(logML + logPR)
  return(obj_value)
}

###--------------------------------------------------------------------------###
###----------------- Standard conjugate BVAR function -----------------------###
###--------------------------------------------------------------------------###
conjVARstd.func <- function(data.setup = data.setup, 
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

M  <- ncol(Y)   # No. of endogenous variables
N  <- nrow(Y)   # No. of observations
K  <- ncol(X)   # No. of coefficients per equation
Mp <- M*p       # No. of coefficients related to lags
v  <- (M+1)*M/2 # No. of free elements in the variance-covariance matrix
k  <- M*K       # k = (M*K) = (M*(M*p+1)) # No. of all coefficients of the system 
  
###--------------------------------------------------------------------------###
###----------------- Estimate a VAR(p) model with OLS -----------------------###
###--------------------------------------------------------------------------###
if(K > N){
  tXX <- crossprod(X) + 1e-3*diag(1,K)
  df <- N
}else{
  tXX <- crossprod(X) + 1e-6*diag(1,K)
  df <- (N-K+M+2)
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
###------------- consistent with symmetric Minnesota prior ------------------###
###--------------------------------------------------------------------------###

###------------- Prior mean of structural VAR coefficients ------------------###
prmean              <- rep(pr.mean,M) 
A_pr                <- matrix(0,K,M) 
diag(A_pr[1:M,1:M]) <- prmean

###------------- Local scalings according to Minnesota moments --------------### 
###------------- Elements of prior variance-covariance matrix ---------------### 
lam_A          <- rep(1, K)
lam_A[1:(M*p)] <- 1/rep((1:p)^2, each = M)
if(cons) lam_A <- lam_A/c(rep(s2.ARp, p),0.01^2) else lam_A <- lam_A/rep(s2.ARp, p)

###----------- Prior moments of structural error variances ------------------###
s_pr <-  M + 2
S_pr <- diag(as.numeric(s2.ARp))

###------- Global scalings (differentiate between own and other lags) -------###
tau_pr <- gammacoef(mode = 1, sd = 0.4) 
tau_lb <- 1/k 
tau_ub <- 10/M

###--------------------------------------------------------------------------###
###----------------------- Hyperparameter set-up ----------------------------###
###------------------------- Posterior set-up -------------------------------###
###------ for changing tau's update the dummies based on a Griddy-Gibbs -----###
###--------------------------------------------------------------------------###

if(tau_gs > 1){
  tau_grid <- data.frame("tau" = seq(tau_lb, tau_ub, length.out = tau_gs))
  tau_grid$logPR <- tau_grid$logML <- tau_grid$logPO <- tau_grid$logDet_pr <- NA
  # Posterior set-up
  po_grid.list <- vector("list", tau_gs)
  
  for(jj in 1:tau_gs){
    # Define prior moments
    tau.jj  <- tau_grid[jj,"tau"]
    v_pr.jj <- lam_A
    v_pr.jj[1:Mp] <- tau.jj*v_pr.jj[1:Mp] 
    v_pr.inv.jj  <- 1/v_pr.jj
    
    pr_jj.list <- list(
      "A_pr"      = A_pr, 
      "v_pr"      = v_pr.jj, 
      "v_pr.inv"  = v_pr.inv.jj,
      "s_pr"      = s_pr, 
      "S_pr"      = S_pr)
    
    # Get posterior moments
    po_jj.list <- get.post.mom(X = X, Y = Y, N = N, pr_list = pr_jj.list)
    po_grid.list[[jj]] <- po_jj.list
    # Compute marginal likelihood
    logML.jj <- get.logML(po_list = po_jj.list, pr_list = pr_jj.list, M = M)
    logPR.jj <- dgamma(tau.jj, shape = tau_pr[[1]], scale = tau_pr[[2]], log = TRUE)
    tau_grid[jj,"logML"] <- logML.jj 
    tau_grid[jj,"logPR"] <- logPR.jj 
    tau_grid[jj,"logPO"] <- logML.jj + logPR.jj
    tau_grid[jj,"logDet_pr"] <- sum(log(v_pr.jj))
    
  }
  
  tau_grid[,"logPOresc"] <- tau_grid[,"logPO"] - max(tau_grid[,"logPO"])
  tau_grid[,"weight"] <- exp(tau_grid[,"logPOresc"])/sum(exp(tau_grid[,"logPOresc"])) 
  tau_map <- tau_grid[which(tau_grid[,"weight"] == max(tau_grid[,"weight"])),"tau"]

}else{
  if(tau_gs == 0){ # FLAT prior
    tau_optim <- list("par" = 100, "value" = NA)
  }else{           # Optimize ML 
    tau_optim <- optim(par=0.1, fn = obj_func, Y = Y, X = X, own.slct = own.slct, other.slct = other.slct, N = N, M = M, Mp = Mp,
                          lam_A = lam_A, A_pr = A_pr, s_pr = s_pr, S_pr = S_pr, 
                          tau_pr = tau_pr, lower=tau_lb, upper=tau_ub, method="L-BFGS-B", hessian = TRUE)
  }
  tau_map  <- tau_optim$par
  tau_grid <- data.frame("tau" = tau_map)
  tau_grid$logML <- -tau_optim$value
  tau_grid$weight <- 1
  # Define prior moments
  v_pr <- lam_A
  v_pr[1:(M*p)] <- tau_map*v_pr[1:(M*p)] 
  v_pr.inv  <- 1/v_pr
  tau_grid$logDet_pr <- sum(log(v_pr))
  
  pr_list <- list(
    "A_pr"      = A_pr, 
    "v_pr"      = v_pr, 
    "v_pr.inv"  = v_pr.inv,
    "s_pr"      = s_pr, 
    "S_pr"      = S_pr)
  
  # Get posterior moments
  po_list <- get.post.mom(X = X, Y = Y, N = N, pr_list = pr_list)
  po_grid.list <- list(po_list)
}

###--------------------------------------------------------------------------###
###---------------------------- Store objects -------------------------------###
###--------------------------------------------------------------------------###
A_store                <- array(NA, c(nsave,K,M))
dimnames(A_store)      <- list(1:nsave,Xlbl, Ylbl)
L_store                <- array(NA, c(nsave,M,M))
dimnames(L_store)      <- list(1:nsave,Ylbl, Ylbl)
L.inv_store            <- array(NA, c(nsave,M,M))
dimnames(L.inv_store)  <- list(1:nsave,Ylbl, Ylbl)
SIG_store              <- array(NA, c(nsave,M,M))
dimnames(SIG_store)    <- list(1:nsave,Ylbl, Ylbl)
tau_pr_store              <- array(NA, c(nsave,1))
dimnames(tau_pr_store)    <- list(1:nsave,"tau")
logDet_pr_store           <- array(NA, c(nsave, 1))
dimnames(logDet_pr_store) <- list(1:nsave, "logDet")

pb <- txtProgressBar(min = 0, max = ntot, style = 3)
irep <- 1
start <- Sys.time()

for(irep in 1:ntot){
###------------------------------ Step 1: -----------------------------------###
###---- Select global shrinkage parameters according to posterior weights ---### 
###--------------------------------------------------------------------------###
if(tau_gs == 0) tau_ind  <- 1 else tau_ind <- sample(1:tau_gs, size=1, prob=tau_grid[,"weight"])

tau_pr_draw     <- tau_grid[tau_ind,"tau"]  
logDet_pr_draw  <- tau_grid[tau_ind,"logDet_pr"]  

po_draw     <- po_grid.list[[tau_ind]]

###------------------------------ Step 2: -----------------------------------###
###---------------- Draw SIG marginally from a inverse Wishart --------------### 
###--------------------------------------------------------------------------###
s_po     <- po_draw[["s_po"]]
S_po.inv <- po_draw[["S_po.inv"]]

SIG.inv_draw <- matrix(rWishart(1,s_po,S_po.inv),M,M)
SIG_draw     <- solve(SIG.inv_draw)
L.inv_draw   <- t(chol(SIG_draw))
L_draw       <- solve(L.inv_draw)

###------------------------------ Step 3: -----------------------------------###
###-- Draw coefficients vec(A) from the conditional posterior vec(A)|SIG,Y --###
###--------------------------------------------------------------------------###
V_po.chol <- po_draw[["V_po.chol"]]
a_po.vec  <- po_draw[["a_po.vec"]]

a_shks <- matrix(rnorm(k),K,M) # Define a K x M matrix of standard normal errors
#a_shks <- as.vector(V1.chol%*%a.shks%*%chol(SIG))
a_shks <- as.vector(V_po.chol%*%a_shks%*%t(L.inv_draw))
a.vec_draw <- a_po.vec + a_shks
A_draw     <- matrix(a.vec_draw,K,M) # Reshape in matrix form

###------------------------------ Step 4: -----------------------------------###
###----------------------------- Storage  -----------------------------------###
###--------------------------------------------------------------------------###
if(irep %in% save.set){
    save.ind <- save.ind + 1
    A_store[save.ind,,]        <- A_draw
    SIG_store[save.ind,,]      <- SIG_draw
    L.inv_store[save.ind,,]    <- L.inv_draw
    L_store[save.ind,,]        <- L_draw
    tau_pr_store[save.ind,]    <- tau_pr_draw 
    logDet_pr_store[save.ind,] <- M*logDet_pr_draw + K*(2*sum(log(diag(L.inv_draw))))
}
setTxtProgressBar(pb, irep)
}
end <- Sys.time()
# Time in minutes
est.time        <- c(as.numeric((ts(end)-ts(start))/60), ntot)
names(est.time) <- c("Minutes", "Draws")

est.mod <- list("A"     = A_store,
                "L"     = L_store, 
                "L.inv" = L.inv_store, 
                "SIG"   = SIG_store, 
                
                "tau_pr"    = tau_pr_store,
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
