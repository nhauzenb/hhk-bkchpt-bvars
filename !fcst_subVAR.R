###--------------------------------------------------------------------------###
###------------------ Forecast code for subspace VARs -----------------------###
###--------------------------------------------------------------------------###
rm(list = ls())

w.dir <- "HHK_bookchpt/"; .libPaths("RLIB_loc")
w.dir <- ""
error.dir <-  paste0(w.dir, "errors/"); dir.create(error.dir, showWarnings =FALSE)
est.dir   <-  paste0(w.dir, "est/") ; dir.create(est.dir, showWarnings = FALSE)
fcst.dir  <-  paste0(w.dir,"fcsts/"); dir.create(fcst.dir, showWarnings = FALSE)

require(scoringRules)
require(MASS)
require(stochvol)
require(RcppArmadillo)
require(Rcpp)
require(dbarts)
require(zoo)
require(shrinkTVP)
require(dplyr)
require(mvtnorm)

func.dir  <-  paste0(w.dir, "bvar_funcs/")
source(paste0(func.dir, "conjVARsub_func.R"))

###--------------------------------------------------------------------------###
###----------------------------- Preliminaries ------------------------------###
###---------------------- same for all specifications -----------------------###
###--------------------------------------------------------------------------###
trans         <- "I0" # Data transformation 
p             <- 4 # Number of lags
str.smp       <- 1970 - p/4 # Sample start
cons          <- TRUE # Intercept
fcst.hor      <- 8 # Forecast horizon
# Target variables for forecast exercise
fcst.vars     <- c("GDPC1", "UNRATE", "CPIAUCSL", "FEDFUNDS", "GS10") 
no.fcstvars   <- length(fcst.vars)
nsave         <- 4000 # Number of stored draws
fcst          <- TRUE
estim         <- !fcst
###--------------------------------------------------------------------------###
###----------------------------- Grid set-up --------------------------------###
###--------------------------------------------------------------------------###
str.qrt     <- 1998
end.qrt     <- 2022+3/4
end.in      <- seq(str.qrt, end.qrt, 1/4)   # end of in-sample
grid.sizes  <- c("S", "SM", "M") #Model sizes

###----------------------- Conjugate VARs -----------------------------------###
grid.conj   <- "conjVAR-SUBMIN"

grid.full <- expand.grid(
  "endsmp"  = end.in,
  "size"    = grid.sizes,
  "model"   = grid.conj,
  stringsAsFactors = FALSE)

run <- as.integer(Sys.getenv("SGE_TASK_ID")) # Assign run to cluster
run <- 10

###--------------------------------------------------------------------------###
###-------------------------- Extract specification -------------------------###
###--------------------------------------------------------------------------###
run.slct    <- grid.full[run,]          # Selected run
run.ismp    <- grid.full[run, "endsmp"] # Final insample period
insmp.dates <- seq(str.smp + p/4, run.ismp, 1/4)
run.endsmp  <-gsub(as.character(zoo::yearqtr(run.ismp)), pattern = " ", replacement = "")
run.h1      <- run.ismp + 1/4          # 1-step ahead
run.hfh     <- run.ismp + fcst.hor/4   # h-step ahead
run.fcstsmp <- gsub(as.character(zoo::yearqtr(seq(run.h1, run.hfh, 1/4))), pattern = " ", replacement = "")

run.size    <- grid.full[run,"size"]  # Information set 
run.mdl     <- grid.full[run,"model"] # Model/prior for conditional mean
run.cm      <- unlist(strsplit(run.mdl, "-"))
run.pr      <- run.cm[[2]] # Prior type
run.cm      <- run.cm[[1]] # Model type 

# Set grid according to specification 
if(run.pr == "FLAT"){
  tau_gs <- 0
}else if(run.pr %in% c("MINo", "ASYMo")){
  tau_gs <- 1
}else{
  tau_gs <- 900 
}

###--------------------------------------------------------------------------###
###----------------------- Model and MCMC set-up ----------------------------###
###--------------------------------------------------------------------------###
mcmc.setup <- list( 
  "nsave" = nsave,                
  "nburn" = 2000,     # burn-in
  "nthin" = 2        # store every nthin-th draw
) 

mdl.setup <- list(
  "prior"         = run.pr,       
  "pr.mean"       = 0,       # Prior mean 
  "tau_gs"        = tau_gs   # Grid size of Griddy Gibbs, which is used to select shrinkage parameter based on optimizing ML
)  

###--------------------------------------------------------------------------###
###----------------------------- Data set-up --------------------------------###
###--------------------------------------------------------------------------###
# Create Y and Yho
load(paste0(w.dir, "fred_data/fred_QD.rda"))
Yho  <- window(Xraw.stat, start = run.h1 ,  end = run.hfh, frequency = 4, extend = TRUE)
Yho  <- Yho[,fcst.vars]

var.slct <- info.sets[[run.size]]
YY       <- window(Xraw.stat, start = str.smp, end = run.ismp, frequency = 4)
YY       <- YY[,var.slct]
# Normalize data
Ymn <- apply(YY, 2, mean)
Ysd <- apply(YY, 2, sd)
YY <- apply(YY, 2, function(x) (x-mean(x))/sd(x))
M <- ncol(YY)
NN <- nrow(YY)
  
data.setup <- list(
  "Y"           = YY,
  "Yho"         = Yho,
  "Ymn"         = Ymn,
  "Ysd"         = Ysd,
  "M"           = M,            # No. of endogenous variables
  "NN"          = NN,           # No. of total observations
  "p"           = p,            # No. of lags
  "cons"        = cons,         # constant
  "fcst.hor"    = fcst.hor,     # Forecast horizon
  "insmp.dates" = insmp.dates
) 

###--------------------------------------------------------------------------###
###----------------------------- Estimation ---------------------------------###
###--------------------------------------------------------------------------###
# Directories 
if(fcst) file.dir <- paste0(fcst.dir, run.size, "/") else file.dir <- paste0(est.dir, run.size, "/")
dir.create(file.dir, showWarnings = FALSE)
folder.name <- paste0(file.dir, run.mdl,"-homo/")
dir.create(folder.name, showWarnings = FALSE)
file.name <- paste0(folder.name, "endsmp-", run.endsmp)
file.check <- paste0(file.name, ".rda")
#file.check <- paste0("endsmp-", run.endsmp, ".rda")
folder.files <- list.files(folder.name)

if (!(file.check %in% folder.files)){
  outputFile <- file(paste0(error.dir,  paste0(run.slct, collapse = "_"), ".txt"))
  tryCatch({
    if(run.mdl %in% c("conjVAR-FLAT", "conjVAR-MINg", "conjVAR-MINo")){
      est.mod <- conjVARstd.func(data.setup, mdl.setup, mcmc.setup)
    }else if(run.mdl %in% c("conjVAR-ASYMg", "conjVAR-ASYMo")){
      est.mod <- conjVARasym.func(data.setup, mdl.setup, mcmc.setup)
    }else if(run.mdl == "conjVAR-SUBMIN"){
      est.mod <- conjVARsub.func(data.setup, mdl.setup, mcmc.setup)
    } 
  }, error = function(e) {
    writeLines(as.character(e), outputFile)
  })
}else{
  stop("\n\n File already exists!")
}

list2env(est.mod, .GlobalEnv)

###--------------------------------------------------------------------------###
###---------------- Ex-post prediction step ---------------------------------###
###--------------------------------------------------------------------------###
# Reject draws with eigen values > 1.05
eig_crit   <- 1.05

# Initialization 
J_th                <- matrix(0,Mp,M)
diag(J_th[1:M,1:M]) <- 1

# Store objects
fcst_store <- fcst.mu_store <- array(NA, c(nsave,M,fcst.hor))
fcst.var_store              <- array(NA, c(nsave,M,M,fcst.hor))
dimnames(fcst_store)        <- dimnames(fcst.mu_store) <- list(1:nsave, Ylbl, run.fcstsmp)
dimnames(fcst.var_store)    <- list(1:nsave, Ylbl, Ylbl, run.fcstsmp)
stab_store <- eig_store     <-  c()

start <- Sys.time()
pb    <- txtProgressBar(min = 0, max = nsave, style = 3) #start progress bar
jrep  <- 1

for(jrep in 1:nsave){
  X_th <- matrix(c(Y[N,],X[N,1:(M*(p-1))]),Mp,1) 
  
  # Specify conditional mean and create companion
  cmu_th <- matrix(0,Mp,1)
  F_th <- matrix(0,Mp,Mp)
  diag(F_th[(M+1):(Mp),1:(Mp-M)]) <- 1
  
  A_th       <- est.mod$A[jrep,,]
  F_th[1:M,] <- t(A_th[1:Mp,])
  if(cons) cmu_th[1:M,1] <- A_th[K,] else  cmu_th[1:M,1] <- 0
  eig_F <- max(Re(eigen(F_th)$values))
  
  # Specify conditional variances 
  O_th       <- matrix(0,Mp,Mp)
  SIG_th     <- est.mod$SIG[jrep,,]

  # Initialization
  fcst_var  <- array( NA,c(M,M,fcst.hor))
  fcst_mu   <- matrix(NA,    M,fcst.hor)
  fcst_draw <- matrix(NA,    M,fcst.hor)
  rownames(fcst_draw) <- rownames(fcst_mu) <- Ylbl
  dimnames(fcst_var)  <- list(Ylbl,Ylbl,run.fcstsmp)
  
  for (hh in 1:fcst.hor){
    X_th     <- F_th%*%X_th + cmu_th
    O_th     <- F_th%*%O_th%*%t(F_th) + J_th%*%SIG_th%*%t(J_th)
    cholO_th <- try(t(chol(O_th[1:M,1:M])),silent=TRUE)
    
    if (is(cholO_th,"try-error")){
      yfc_th <- as.numeric(rmvnorm(1,mean= X_th[1:M],sigma= O_th[1:M,1:M]))
    }else{
      yfc_th <- as.numeric(X_th[1:M]+cholO_th%*%rnorm(M,0,1))
    }  
    
    fcst_draw[,hh] <- (yfc_th)*Ysd + Ymn
    fcst_var[,,hh] <- diag(Ysd)%*%O_th[1:M,1:M]%*%diag(Ysd)
    fcst_mu[,hh]   <- X_th[1:M] + Ymn
  }
  
  eig_store[jrep]         <- eig_F
  stab_store[jrep]        <- (eig_F < eig_crit)
  fcst_store[jrep,,]      <- fcst_draw
  fcst.mu_store[jrep,,]   <- fcst_mu
  fcst.var_store[jrep,,,] <- fcst_var
  
  setTxtProgressBar(pb, jrep)
}

###--------------------------------------------------------------------------###
###----------------------------- Evaluation ---------------------------------###
###--------------------------------------------------------------------------###
fcst_store         <- fcst_store[stab_store,fcst.vars,]
fcst.cov           <- array(0, c(no.fcstvars, no.fcstvars,fcst.hor))
dimnames(fcst.cov) <- list(fcst.vars, fcst.vars, run.fcstsmp)
for(pp in 1:fcst.hor) fcst.cov[,,pp] <- cov(fcst_store[,,pp])

fcst.smry <- apply(fcst_store, c(2,3), quantile, c(0.05, 0.5, 0.95), na.rm = T)
fcst.med  <- fcst.smry["50%",,]
fcst.sd   <- apply(fcst_store, c(2,3), sd, na.rm = T)

Yho           <- t(Yho)
colnames(Yho) <- run.fcstsmp

afeM <- abs(Yho - fcst.med)
sfeM <- (Yho - fcst.med)^2
# Density est. evaluation 
lplT  <- crpsT <- rep(NA, fcst.hor)
lplM  <- crpsM <- matrix(NA, no.fcstvars, fcst.hor)
colnames(crpsM) <- names(crpsT) <- colnames(lplM) <- names(lplT) <- run.fcstsmp
rownames(crpsM) <- rownames(lplM) <- fcst.vars

for(hh in 1:fcst.hor){
  lplT[[hh]]  <- mvtnorm::dmvnorm(Yho[,hh], fcst.med[,hh], fcst.cov[,,hh], log = TRUE)
  crpsT[[hh]] <- es_sample(Yho[,hh], t(fcst_store[,,hh]))
  lplM[,hh]   <- dnorm(Yho[,hh], fcst.med[,hh], fcst.sd[,hh], log = TRUE)
  if(!all(is.na(Yho[,hh]))) for(jj in 1:no.fcstvars) crpsM[jj,hh] <- crps_sample(y = Yho[jj,hh], dat = fcst_store[,jj,hh])
}

###--------------------------------------------------------------------------###
###-------------------------- Posterior means -------------------------------###
###--------------------------------------------------------------------------###
post.smry <- list(
  "A"     = apply(est.mod$A   , c(2,3), mean), 
  "SIG"   = apply(est.mod$SIG, c(2,3) , mean), 
  
  "logDet_pr" = est.mod$logDet_pr,
  "tau_pr"    = apply(est.mod$tau_pr, 2, mean),
  "q_pr"      = table(est.mod$q_pr))

OLS.est   <- est.mod$OLS.est
  
###--------------------------------------------------------------------------###
###----------------------------- Storage ------------------------------------###
###--------------------------------------------------------------------------###
ret.obj = list(
  "fcst_smry" = fcst.smry,
  "fcst_cov"  = fcst.cov, 
  "fcst_sd"   = fcst.sd, 
  "Yho"       = Yho, 
  
  "afeM"      = afeM, 
  "sfeM"      = sfeM, 
  "crpsM"     = crpsM, 
  "lplM"      = lplM, 
  "crpsT"     = crpsT, 
  "lplT"      = lplT, 
  
  "insmp_prd" = run.endsmp, 
  "ho_prds"   = run.fcstsmp,
  "est_time"  = est.time, 
  
  "post_smry" = post.smry,
  "OLS_est"   = OLS.est
)

save(file = paste0(file.name, ".rda"), list = c("ret.obj"))

