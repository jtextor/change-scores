###############################
## Peter's Simulation R-code ##
###############################

# Running this code generates CSV files in the current directory.
# These CSV files contain the results of all simulations reported
# in the paper. The file whose name starts with "final-" contains
# the numbers reported in Table 2 in the manuscript. The numbers 
# in the manuscript were generated using 100000 replicates. Change
# the line below to increase the replicates accordingly.

# setup simulation parameters
N       <- 1000
Nreps   <- 10

if( packageVersion('dagitty')<"0.2.3" ){
	warning("Please install at least version 0.2.3 of the dagitty package!")
	stop("Use this command: devtools::install_github('jtextor/dagitty/r')")
}

require(dagitty)
require(MASS)
require(rpsychi)

# function to execute multiple simulations and summarise findings
#Mu=Nmu; Sigma=Sig; Nobs=N; Nsims=Nreps; Seed=13; itn=1
runSims <- function(Mu,Sigma,Nobs,Nsims,Seed) {
  start <- Sys.time()
  Sum   <- NULL
  for (itn in 1:Nsims) {
    seed     <- Seed*Nobs*itn
    dat      <- data.frame(mvrnorm(Nobs,Mu,Sigma,empirical=FALSE)); names(dat) <-c("Y0","Y1","U","X")
    dat$DY   <- dat$Y1 - dat$Y0
    Means    <- apply(dat,2,mean)
    SDs      <- apply(dat,2,sd)
    CorMat   <- cor(dat)
    Corrs    <- CorMat[lower.tri(CorMat)]
    names(Corrs) <- c("Y0_Y1","Y0_U","Y0_X","Y0_DY","Y1_U","Y1_X","Y1_DY","U_X","U_DY","X_DY")
    mod1     <- lm(DY~X,data=dat)
    mod2     <- lm(Y1~X,data=dat)
    mod3     <- lm(Y1~X+Y0,data=dat)
    Coeffs   <- c(mod1$coefficients[2],mod2$coefficients[2],mod3$coefficients[2])
    names(Coeffs) <- c("Beta1","Beta2","Beta3")
    DataSum  <- list(Means,SDs,Corrs,Coeffs)
    Sum      <- rbind(Sum,unlist(DataSum)) 
    }
  end  <- Sys.time(); print(end-start)
  Sims <- apply(Sum,2,function(x) {quantile(x,c(0.025,0.5,0.975))})
  return(Sims)
}

DAG_base_confounder <- function(pU_WC0,pU_IC0,pU_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1) {
  dag <- dagitty(paste0("dag{  U->WC0 [beta=",pU_WC0,"]     U->IC0 [beta=",pU_IC0,"]     U->IC1 [beta=",pU_IC1,"] 
                             IC0->WC0 [beta=",pWC0_IC0,"] WC0->IC1 [beta=",pWC0_IC1,"] IC0->IC1 [beta=",pIC0_IC1,"]}"))
  return(dag) }

DAG_base_mediator <- function(pU_WC0,pU_IC0,pU_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1) {
  dag <- dagitty(paste0("dag{  U->WC0 [beta=",pU_WC0,"]     U->IC0 [beta=",pU_IC0,"]     U->IC1 [beta=",pU_IC1,"] 
                             IC0<-WC0 [beta=",pWC0_IC0,"] WC0->IC1 [beta=",pWC0_IC1,"] IC0->IC1 [beta=",pIC0_IC1,"]}"))
  return(dag) }

DAG_base_mediator2 <- function(pU1_WC0,pU1_IC0,pU1_IC1,pU2_IC0,pU2_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1) {
  dag <- dagitty(paste0("dag{ U1->WC0 [beta=",pU1_WC0,"]   U1->IC0 [beta=",pU1_IC0,"]   U1->IC1 [beta=",pU1_IC1,"] 
                             IC0<-WC0 [beta=",pWC0_IC0,"]  U2->IC0 [beta=",pU2_IC0,"]   U2->IC1 [beta=",pU2_IC1,"] 
                                                          WC0->IC1 [beta=",pWC0_IC1,"] IC0->IC1 [beta=",pIC0_IC1,"]}"))
  return(dag) }

##########################################
## Setup initial values for simulations ##
##########################################

# set means and variances
WCmu    <- 9.5
WCvar   <- 1.6^2
IC0mu   <- 4.0
IC0var  <- 0.74^2
IC1mu   <- 4.2
IC1var  <- 0.74^2

# set vectors for simulations
Nmu     <- c(IC0mu, IC1mu, 0,WCmu) 
Nvar    <- c(IC0var,IC1var,1,WCvar)

# set consistent path coefficient
pIC0_IC1 <- 0.65
pWC0_IC1 <- EffSize <- 0.433

# final output file identifier
Name    <- paste0("final-",(Nreps/1000),"k","-",substr(Sys.time(),1,10),".csv")

################
# Scenario 1A ##
################

pU_IC0   <- 0.0; pU_WC0 <- 0.0; pU_IC1 <- 0.0; pWC0_IC0 <- 0.0; 
dag      <- DAG_base_confounder(pU_WC0,pU_IC0,pU_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1)
Cor      <- impliedCovarianceMatrix(dag)
Sig      <- r2cov(sqrt(Nvar),Cor)
Sim1a    <- data.frame(runSims(Nmu,Sig,N,Nreps,13))
Filename <- paste0("1A-full-",Name); write.csv(Sim1a,file=Filename,row.names=FALSE)

################
# Scenario 1B ##
################

pU_IC0   <- 0.2; pU_IC1 <- 0.05; pU_WC0 <- 0.5; pWC0_IC0 <- 0.0
dag      <- DAG_base_confounder(pU_WC0,pU_IC0,pU_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1)
Cor      <- impliedCovarianceMatrix(dag)
Sig      <- r2cov(sqrt(Nvar),Cor)
Sim1b    <- data.frame(runSims(Nmu,Sig,N,Nreps,17))
Filename <- paste0("1B-full-",Name); write.csv(Sim1b,file=Filename,row.names=FALSE)

# set consistent path coefficient
pWC0_IC0 <- 0.5

#################
## Scenario 2A ##
#################

pU_IC0   <- 0.0; pU_IC1 <- 0.0; pU_WC0 <- 0.0
dag      <- DAG_base_confounder(pU_WC0,pU_IC0,pU_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1)
Cor      <- impliedCovarianceMatrix(dag)
Sig      <- r2cov(sqrt(Nvar),Cor)
Sim2a    <- data.frame(runSims(Nmu,Sig,N,Nreps,19))
Filename <- paste0("2A-full-",Name); write.csv(Sim2a,file=Filename,row.names=FALSE)

################
# Scenario 2B ##
################

pU_IC0   <- 0.2; pU_IC1 <- 0.05; pU_WC0 <- 0.5
dag      <- DAG_base_confounder(pU_WC0,pU_IC0,pU_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1)
Cor      <- impliedCovarianceMatrix(dag)
Sig      <- r2cov(sqrt(Nvar),Cor)
Sim2b    <- data.frame(runSims(Nmu,Sig,N,Nreps,23))
Filename <- paste0("2B-full-",Name); write.csv(Sim2b,file=Filename,row.names=FALSE)

#################
## Scenario 3A ##
#################

pU_IC0   <- 0.0; pU_IC1 <- 0.0; pU_WC0 <- 0.0
pWC0_IC1 <- EffSize-(pIC0_IC1*pWC0_IC0)
dag      <- DAG_base_mediator(pU_WC0,pU_IC0,pU_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1)
Cor      <- impliedCovarianceMatrix(dag)
Sig      <- r2cov(sqrt(Nvar),Cor)
Sim3a    <- data.frame(runSims(Nmu,Sig,N,Nreps,29))
Filename <- paste0("3A-full-",Name); write.csv(Sim3a,file=Filename,row.names=FALSE)

################
# Scenario 3B ##
################

pU_IC0   <- 0.2; pU_IC1 <- 0.05; pU_WC0 <- 0.5
pWC0_IC1 <- EffSize-(pIC0_IC1*pWC0_IC0)
dag      <- DAG_base_mediator(pU_WC0,pU_IC0,pU_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1)
Cor      <- impliedCovarianceMatrix(dag)
Sig      <- r2cov(sqrt(Nvar),Cor)
Sim3b    <- data.frame(runSims(Nmu,Sig,N,Nreps,31))
Filename <- paste0("3B-full-",Name); write.csv(Sim3b,file=Filename,row.names=FALSE)

# reset vetors for simulations
Nmu     <- c(IC0mu, IC1mu, 0,0,WCmu) 
Nvar    <- c(IC0var,IC1var,1,1,WCvar)

runSims2 <- function(Mu,Sigma,Nobs,Nsims,Seed) {
  start  <- Sys.time()
  Sum    <- NULL
  for (itn in 1:Nsims) {
    seed     <- Seed*Nobs*itn
    dat      <- data.frame(mvrnorm(Nobs,Mu,Sigma,empirical=FALSE)); names(dat) <-c("Y0","Y1","U1","U2","X")
    dat$DY   <- dat$Y1 - dat$Y0
    Means    <- apply(dat,2,mean)
    SDs      <- apply(dat,2,sd)
    CorMat   <- cor(dat)
    Corrs    <- CorMat[lower.tri(CorMat)]
    names(Corrs) <- c("Y0_Y1","Y0_U1","Y0_U2","Y0_X","Y0_DY","Y1_U1","Y1_U2","Y1_X","Y1_DY","U1_U2","U1_X","U1_DY","U2_X","U2_DY","X_DY")
    mod1     <- lm(DY~X,data=dat)
    mod2     <- lm(Y1~X,data=dat)
    mod3     <- lm(Y1~X+Y0,data=dat)
    Coeffs   <- c(mod1$coefficients[2],mod2$coefficients[2],mod3$coefficients[2])
    names(Coeffs) <- c("Beta1","Beta2","Beta3")
    DataSum  <- list(Means,SDs,Corrs,Coeffs)
    Sum      <- rbind(Sum,unlist(DataSum)) 
  }
  end  <- Sys.time(); print(end-start)
  Sims <- apply(Sum,2,function(x) {quantile(x,c(0.025,0.5,0.975))})
  return(Sims)
}

##################
## Scenario 3Ap ##
##################

pU1_IC0  <- 0.0; pU1_IC1 <- 0.0; pU1_WC0  <- 0.0
pU2_IC0  <- 0.5; pU2_IC1 <- 0.2
pWC0_IC1 <- EffSize-(pIC0_IC1*pWC0_IC0)
dag      <- DAG_base_mediator2(pU1_WC0,pU1_IC0,pU1_IC1,pU2_IC0,pU2_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1)
Cor      <- impliedCovarianceMatrix(dag)
Sig      <- r2cov(sqrt(Nvar),Cor)
Sim4a    <- data.frame(runSims2(Nmu,Sig,N,Nreps,37))
Filename <- paste0("3Ap-full-",Name); write.csv(Sim4a,file=Filename,row.names=FALSE)

#################
# Scenario 3Bp ##
#################

pU1_IC0  <- 0.2; pU1_IC1 <- 0.05; pU1_WC0  <- 0.5
pU2_IC0  <- 0.5; pU2_IC1 <- 0.2
pWC0_IC1 <- EffSize-(pIC0_IC1*pWC0_IC0)
dag      <- DAG_base_mediator2(pU1_WC0,pU1_IC0,pU1_IC1,pU2_IC0,pU2_IC1,pWC0_IC0,pWC0_IC1pIC0_IC1)
Cor      <- impliedCovarianceMatrix(dag)
Sig      <- r2cov(sqrt(Nvar),Cor)
Sim4b    <- data.frame(runSims2(Nmu,Sig,N,Nreps,41))
Filename <- paste0("3Bp-full-",Name); write.csv(Sim4b,file=Filename,row.names=FALSE)

#########################
## Process all results ##
#########################

Summ     <- NULL
for (itn in c("1A","1B","2A","2B","3A","3B","3Ap","3Bp")) {
  Label    <- paste0(itn,"-Full")
  Filename <- paste0(Label,"-",Name)
  File     <- read.csv(Filename)[,c("Beta1","Beta2","Beta3")]
  rownames(File) <- paste(Label,c("2.5%","50%","97.5%"))
  Summ     <- rbind(Summ,File) }

write.table(signif(Summ), Name, sep=",",row.names=TRUE,col.names=NA)


