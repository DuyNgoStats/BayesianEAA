####################################################################################
####################### Epigenetic Clock Derivation ################################
####################################################################################
library(methylclock)

setwd("/Users/WMU-Guest/Desktop/Xi Qiao project")

################################### M value Clock ##################################
M = readRDS("dataset/M_df_213.RDS")
CpGMissing = checkClocks(M)
## By default the package computes the different clocks 
## when there are more than 80% of the required CpGs of each method. 
## All clock meet the requirement

Clock = DNAmAge(M,toBetas = TRUE)
write.csv(Clock,file = "EpigeneticClock/Clock_by_M.csv",row.names = FALSE)
## Bayesian Neural Network (BNN) Bayesian method cannot be estimated

################################## Beta value Clock ################################
CpGMissing = checkClocks(beta)
## By default the package computes the different clocks 
## when there are more than 80% of the required CpGs of each method. 
## All clock meet the requirement

Clock = DNAmAge(beta)
clock = Clock[order(Clock$id),]
write.csv(clock,file = "EpigeneticClock/Clock_by_Beta.csv",row.names = FALSE)
## Bayesian Neural Network (BNN) Bayesian method cannot be estimated




####################################################################################
############################ Accelerated Age Estimation  ###########################
####################################################################################
library(MCMCglmm)
set.seed(2022)

ClockMaster = read.csv("ClocksMaster.csv") ## 213  50

#################### Variables re-level and standardize ##########################
ClockMaster$mother = factor(ClockMaster$mother)
ClockMaster$drought = relevel(factor(ClockMaster$drought), ref = "0")
ClockMaster$region = relevel(factor(ClockMaster$region), ref = "1")

ClockMaster$sex = factor(ClockMaster$sex)

ClockMaster$Cows = scale(ClockMaster$Cows,center = T, scale = T)
ClockMaster$Shoats = scale(ClockMaster$Shoats,center = T, scale = T)
ClockMaster$Camels = scale(ClockMaster$Camels,center = T, scale = T)

ClockMaster$WaterInsecurity = factor(ClockMaster$WaterInsecurity)
ClockMaster$Voluntarynotforcedhazardousherdingortreeclimbingduringpregnancy =
  factor(ClockMaster$Voluntarynotforcedhazardousherdingortreeclimbingduringpregnancy)

ClockMaster$Birth_Season1 = ifelse(ClockMaster$Birth_Season=="Dry","Dry","MediumOrWet" )

############## AccAge estimated from regression residuals ########################
fitPedBE = lm(ClockMaster$PedBE~ClockMaster$ageinyear)
ClockMaster$AccAgePedBE = fitPedBE$residuals

fitHorvath = lm(ClockMaster$Horvath~ClockMaster$ageinyear)
ClockMaster$AccAgeHorvath = fitHorvath$residuals

fitHannum = lm(ClockMaster$Hannum~ClockMaster$ageinyear)
ClockMaster$AccAgeHannum = fitHannum$residuals

fitskinHorvath= lm(ClockMaster$skinHorvath~ClockMaster$ageinyear)
ClockMaster$AccAgeskinHorvath = fitskinHorvath$residuals

write.csv(ClockMaster,"AccAgeMaster.csv")

################# Clocks and AccAge Correlation matrix #############################
library(ggcorrplot)
cc = ClockMaster[,c("ageinyear","PedBE","Horvath","Hannum","skinHorvath")] 
ggcorrplot(cor(cc),lab = TRUE,title = "Correlation Matrix of Epigenetic clocks") + 
  scale_fill_gradient2(breaks = c(0,1),limit = c(0,1),low = "white", high = "red")

aa = ClockMaster[,c("ageinyear","AccAgePedBE","AccAgeHorvath","AccAgeHannum","AccAgeskinHorvath")] 
ggcorrplot(cor(aa),lab = TRUE,title = "Correlation Matrix of Accelerated Age") + 
  scale_fill_gradient2(breaks = c(0,1),limit = c(0,1),low = "white", high = "red")





####################################################################################
############################## Multivariate Analysis  ##############################
####################################################################################

ClockMaster = read.csv("AccAgeMaster.csv")

## Set prior for G matrix and R matrix
prior = list(
  G = list(G1 = list(V = diag(rep(1,4)), nu = 4)),
  R = list(R1 = list(V = diag(rep(1,4)), nu = 4))
)

## Region reference level: highland
ClockMaster$region = relevel(factor(ClockMaster$region),ref = "1")

#################### Multivariate Model for Association  Analysis ################
############################## Exposure: drought #################################
fit = MCMCglmm(cbind(AccAgePedBE, AccAgeHorvath, AccAgeHannum, AccAgeskinHorvath) ~ 
                 trait + trait : ( drought +sex + Epi + Fib + gravida) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",4),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 10000,burnin = 2000,thin = 1)
summary(fit)

############################## Exposure: region #################################
fit = MCMCglmm(cbind(AccAgePedBE, AccAgeHorvath, AccAgeHannum, AccAgeskinHorvath) ~ 
                 trait + trait : ( region +sex + Epi + Fib + gravida + WaterInsecurity +  
                                     Voluntarynotforcedhazardousherdingortreeclimbingduringpregnancy+ 
                                     Long.term_Prengancy_Stress_Standardized_Instrument) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",4),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 10000,burnin = 2000,thin = 1)
summary(fit)

############################## Exposure: FirstTrimesterMeanTemp #################################
fit = MCMCglmm(cbind(AccAgePedBE, AccAgeHorvath, AccAgeHannum, AccAgeskinHorvath) ~ 
                 trait + trait : ( FirstTrimesterMeanTemp +sex + Epi + Fib + gravida + WaterInsecurity +  
                                     PregestationRainfallZ+ 
                                     Long.term_Prengancy_Stress_Standardized_Instrument) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",4),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 10000,burnin = 2000,thin = 1)
summary(fit)

############################## Exposure: SecondTrimesterMeanTemp #################################
fit = MCMCglmm(cbind(AccAgePedBE, AccAgeHorvath, AccAgeHannum, AccAgeskinHorvath) ~ 
                 trait + trait : ( SecondTrimesterMeanTemp +sex + Epi + Fib + gravida + WaterInsecurity +  
                                     Voluntarynotforcedhazardousherdingortreeclimbingduringpregnancy+ 
                                     Long.term_Prengancy_Stress_Standardized_Instrument) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",4),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 10000,burnin = 2000,thin = 1)
summary(fit)


############################## Exposure: ThirdTrimesterMeanTemp #################################
fit = MCMCglmm(cbind(AccAgePedBE, AccAgeHorvath, AccAgeHannum, AccAgeskinHorvath) ~ 
                 trait + trait : ( ThirdTrimesterMeanTemp +sex + Epi + Fib + gravida + WaterInsecurity +  
                                     Voluntarynotforcedhazardousherdingortreeclimbingduringpregnancy+ 
                                     Long.term_Prengancy_Stress_Standardized_Instrument) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",4),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 10000,burnin = 2000,thin = 1)
summary(fit)


############################## Exposure: TwelveMonthsMeanTemp  #################################
fit = MCMCglmm(cbind(AccAgePedBE, AccAgeHorvath, AccAgeHannum, AccAgeskinHorvath) ~ 
                 trait + trait : ( TwelveMonthsMeanTemp +sex + Epi + Fib + gravida + WaterInsecurity +  
                                     Voluntarynotforcedhazardousherdingortreeclimbingduringpregnancy+ 
                                     Long.term_Prengancy_Stress_Standardized_Instrument) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",4),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 10000,burnin = 2000,thin = 1)
summary(fit)

################################### Sensitivity Analysis ########################################
fit = MCMCglmm(cbind(AccAgePedBE, AccAgeHorvath, AccAgeHannum, AccAgeskinHorvath) ~ 
                 trait + trait : ( FirstTrimesterMeanTemp + drought +sex + Epi + Fib + gravida + WaterInsecurity +  
                                     PregestationRainfallZ+ 
                                     Long.term_Prengancy_Stress_Standardized_Instrument) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",4),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 10000,burnin = 2000,thin = 1)

summary(fit)

fit = MCMCglmm(cbind(AccAgePedBE, AccAgeHorvath, AccAgeHannum, AccAgeskinHorvath) ~ 
                 trait + trait : ( SecondTrimesterMeanTemp + drought +sex + Epi + Fib + gravida + WaterInsecurity +  
                                     Voluntarynotforcedhazardousherdingortreeclimbingduringpregnancy+ 
                                     Long.term_Prengancy_Stress_Standardized_Instrument) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",4),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 10000,burnin = 2000,thin = 1)

summary(fit)

fit = MCMCglmm(cbind(AccAgePedBE, AccAgeHorvath, AccAgeHannum, AccAgeskinHorvath) ~ 
                 trait + trait : ( ThirdTrimesterMeanTemp + drought +sex + Epi + Fib + gravida + WaterInsecurity +  
                                     Voluntarynotforcedhazardousherdingortreeclimbingduringpregnancy+ 
                                     Long.term_Prengancy_Stress_Standardized_Instrument) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",4),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 10000,burnin = 2000,thin = 1)

summary(fit)

fit = MCMCglmm(cbind(AccAgePedBE, AccAgeHorvath, AccAgeHannum, AccAgeskinHorvath) ~ 
                 trait + trait : ( TwelveMonthsMeanTemp + drought +sex + Epi + Fib + gravida + WaterInsecurity +  
                                     Voluntarynotforcedhazardousherdingortreeclimbingduringpregnancy+ 
                                     Long.term_Prengancy_Stress_Standardized_Instrument) - 1, 
               random = ~us(trait):mother, rcov = ~ us(trait):units,
               family = rep("gaussian",4),
               data = ClockMaster,prior = prior,
               verbose = FALSE,nitt = 10000,burnin = 2000,thin = 1)

summary(fit)

######################### Variable correlation checking ###############################

tmp = ClockMaster[,c("drought","region", "ageinyear","ThirdTrimesterMeanTemp","sex" , "Epi", "Fib", "gravida", 
                     "Long.term_Prengancy_Stress_Standardized_Instrument", "PregestationRainfallZ",
                     "WaterInsecurity", "Voluntarynotforcedhazardousherdingortreeclimbingduringpregnancy")] 
cor(tmp)














