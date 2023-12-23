####################################################################################
####################### Missing CpG Imputation ####################################
####################################################################################
## Golden standard for missing CpG for GrimAge2 paper
CpGmedian = read.csv("datMiniAnnotation3_Gold.csv")

## Pheno data
datPheno = read.csv("Master.csv")
order = datPheno$sampleid

## methylation CpGS
betavalue = readRDS("beta_filt.RDS")

## Rows are samples, columns are CpGs; Order matches with pheno data.
beta_value = as.data.frame(t(betavalue))
beta_value = beta_value[order,] 

## Start Imputation with Golden standard median value 
missing = CpGmedian$CpG [which(!CpGmedian$CpG %in% colnames(dat.meth))]
dat.meth[,missing] <- 0

for(i in 1:length(missing)){
  dat.meth[,missing[i]] <- CpGmedian$gold[which(CpGmedian$CpG==missing[i])]
}

saveRDS(dat.meth, file = "Imputed_Betadf.RDS")

####################################################################################
####################### Epigenetic Clock Derivation ################################
####################################################################################
library(methylCIPHER)
library(dnaMethyAge)

datPheno = read.csv("Master.csv")

## Load DNAm data
dat.meth = readRDS("Imputed_Betadf.RDS")
dat.meth_T = t(dat.meth) 

## Hannum, Horvath1, Horvath2, PEDBE, PhenoAge
datPheno$Hannum = calcHannum(dat.meth) 
datPheno$Horvath1 = calcHorvath1(dat.meth)      
datPheno$Horvath2 = calcHorvath2(dat.meth)
datPheno$PedBE = calcPEDBE(dat.meth)  
datPheno$PhenoAge = calcPhenoAge(dat.meth)


## Wu's clock: coefficient downloaded from Wu's paper.
load("coefWu.rda")
Wu.CpG = dat.meth[,coefWu$CpGmarker[-1]]
myWu = function(vec){
  x = c(1,vec)
  wu = as.numeric(x) %*% as.numeric(coefWu$CoefficientTraining)
  a = ifelse (wu < 0, 3 * exp(wu) - 1, 3 * wu + 2)
  return(a)
}

datPheno$Wu = apply(Wu.CpG,1,myWu)

## DNAmTL
datPheno$DNAmTL = methyAge(dat.meth_T,"LuA2019")[,"mAge"]

## DunedinPoAm38, DunedinPACE
datPheno$DunedinPoAm38 = calcDunedinPoAm38(dat.meth, imputation = F)
datPheno$DunedinPACE = methyAge(dat.meth_T,clock = "DunedinPACE")[,"mAge"] 

write.csv(datPheno,file = "ClocksMaster.csv",row.names = FALSE)

####################################################################################
####################### Age Accleration  Derivation ################################
####################################################################################
ClockMaster = read.csv("ClocksMaster.csv")

fitPedBE = lm(ClockMaster$PedBE~ClockMaster$ageinyear)
ClockMaster$AccAgePedBE = fitPedBE$residuals

fitWu = lm(ClockMaster$Wu~ClockMaster$ageinyear)
ClockMaster$AccAgeWu = fitWu$residuals

fitHorvath = lm(ClockMaster$Horvath1~ClockMaster$ageinyear)
ClockMaster$AccAgeHorvath = fitHorvath$residuals

fitskinHorvath= lm(ClockMaster$Horvath2~ClockMaster$ageinyear)
ClockMaster$AccAgeskinHorvath = fitskinHorvath$residuals

fitHannum = lm(ClockMaster$Hannum~ClockMaster$ageinyear)
ClockMaster$AccAgeHannum = fitHannum$residuals

fitPhenoAge = lm(ClockMaster$PhenoAge~ClockMaster$ageinyear)
ClockMaster$AccAgePhenoAge = fitPhenoAge$residuals

fitGrimAge2 = lm(ClockMaster$GrimAge2~ClockMaster$ageinyear)
ClockMaster$AccAgeGrimAge2 = fitGrimAge2$residuals

fitDNAmTL = lm(ClockMaster$DNAmTL~ClockMaster$ageinyear)
ClockMaster$AccAgeDNAmTL = fitDNAmTL$residuals

write.csv(ClockMaster,"AccAgeMaster.csv")

######################################################################################
###################### Cell Type Proportion Estimation ###############################
######################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EpiDISH")

library(EpiDISH)

beta_filt=readRDS("beta_filt.RDS")

## Convert to data frame
beta_filt_dataframe = as.data.frame(beta_filt)

## Obtain Epi, Fib, IC
out = epidish(beta.m = beta_filt_dataframe, ref.m = centEpiFibIC.m, method = "RPC")
outEst = out$estF
outEstF = as.data.frame(outEst)
outEstF$sampleid = rownames(outEstF)

## Merge cell type with Master
library(dplyr)
read.csv("AccAgeMaster.csv")
AccAgeMaster = left_join(AccAgeMaster,outEstF, by = "sampleid" )

write.csv(AccAgeMaster,file = "AccAgeMaster.csv", row.names = FALSE)

####################################################################################
####################### BayesianEAA--Multivariate Analysis  ########################
####################################################################################
library(MCMCglmm)
set.seed(2023)

ClockMaster = read.csv("AccAgeMaster.csv")

## Set prior for G matrix and R matrix
prior = list(
  G = list(G1 = list(V = diag(rep(1,10)), nu = 10)),
  R = list(R1 = list(V = diag(rep(1,10)), nu = 10))
)


## Multivariate Linear Mixed Model
fit = MCMCglmm(cbind(AccAgePedBE,AccAgeWu,AccAgeHorvath, AccAgeskinHorvath,
                     AccAgeHannum,AccAgePhenoAge,AccAgeGrimAge2,AccAgeDNAmTL,
                     DunedinPACE,DunedinPoAm38) ~ 
                 trait + trait : (drought + Female + Epi + Fib + gravida+ Birth_Season1) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",10),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 20000,burnin = 5000,thin = 1)

summary(fit)

## Posterior diagnostic plots
plot.MCMCglmm(fitdrought1)


