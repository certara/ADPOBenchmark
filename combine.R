rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggplot2)

home_dir <- getwd()   
# nlme has one extra THETA, for residual error
NLME_Std_all <- read.csv(file.path("NLme","NLMEStdResults.csv")) 
NLME_Dual_all <- read.csv(file.path("NLme_dual","NLMEDualResults.csv")) 
NONMEM_all <- read.csv(file.path("NONMEM","NONMEMResults.csv"))
# sometime, if .rdata file save an certain time, model gets repeated??, so, only use first
NONMEM_all <- NONMEM_all %>% 
  distinct(Model_num,.keep_all = TRUE)
NLME_Std_times <- NLME_Std_all %>% 
  select(Model_num,Est_time,Cov_time) %>% 
  filter(Est_time > 0) %>% 
  filter(Cov_time > 0)

NLME_Dual_times <- NLME_Dual_all %>%  
  select(Model_num,Est_time,Cov_time) %>% 
  filter(Est_time > 0) %>% 
  filter(Cov_time > 0)
NONMEM_times <- NONMEM_all %>%  
select(Model_num,Est_time,Cov_time) %>% 
  filter(Est_time > 0) %>% 
  filter(Cov_time > 0)
summary(NONMEM_times)
geoMeanNONMENEst  <- round(exp(mean(log(NONMEM_times$Est_time))),1)
geoMeanNONMENCov  <- round(exp(mean(log(NONMEM_times$Cov_time))),1)
geoMeanNLMEStdEst  <- round(exp(mean(log(NLME_Std_times$Est_time))),1)
geoMeanNLMEStdCov  <- round(exp(mean(log(NLME_Std_times$Cov_time))),1)
geoMeanNLMEDualEst  <-round(exp(mean(log(NLME_Dual_times$Est_time))),1)
geoMeanNLMEDualCov  <- round(exp(mean(log(NLME_Dual_times$Cov_time))),1)
summary(NLME_Std_times)
summary(NLME_Dual_times)

# set with both NLME and NONMEM
times <- inner_join(x=NLME_Dual_times,y=inner_join(x=NONMEM_times, y=NLME_Std_times,by="Model_num") )
colnames(times) <- c("Model_Num","NLME_Dual_Est","NLME_Dual_Cov","NONMEM_Est","NONMEM_Cov","NLME_Std_Est","NLME_Std_Cov")
Sum_times <- c("GMean_NONMEM_EST" = round(exp(mean(log(times$NONMEM_Est))),0),
               "GMean_NLME_Std_EST" = round(exp(mean(log(times$NLME_Std_Est))),0),
               "GMean_NLME_Dual_EST" = round(exp(mean(log(times$NLME_Dual_Est))),0),
               "GMean_NONMEM_Cov" = round(exp(mean(log(times$NONMEM_Cov))),1),
               "GMean_NLME_Std_Cov" = round(exp(mean(log(times$NLME_Std_Cov))),1),
               "GMean_NLME_Dual_Cov" = round(exp(mean(log(times$NLME_Dual_Cov))),1)
               )
 
Sum_times 
Sum_success <- c("FSuccess_NONMEM" = mean(NONMEM_all$Success),
                 "FSuccess_NLME_std" = mean(NLME_Std_all$Success),
                 "FSuccess_NLME_Dual" = mean(NLME_Dual_all$Success),
                 "FCov_NONMEM" = mean(NONMEM_all$Covar),
                 "FCov_NLME_std" = mean(NLME_Std_all$Covar),
                 "FCov_NLME_Dual" = mean(NLME_Dual_all$Covar))
Sum_success
NONMEM_all <- NONMEM_all %>% 
  mutate(absNPDECmaxMean = if_else(NPDECmaxMean > -990, abs(NPDECmaxMean), NA))%>% 
  mutate(absNPDECminMean = if_else(NPDECminMean > -990, abs(NPDECminMean), NA)) 
NLME_Std_all <- NLME_Std_all %>% 
  mutate(absNPDECmaxMean = if_else(NPDECmaxMean > -990, abs(NPDECmaxMean), NA))%>% 
  mutate(absNPDECminMean = if_else(NPDECminMean > -990, abs(NPDECminMean), NA)) 
NLME_Dual_all <- NLME_Dual_all %>% 
  mutate(absNPDECmaxMean = if_else(NPDECmaxMean > -990, abs(NPDECmaxMean), NA))%>% 
  mutate(absNPDECminMean = if_else(NPDECminMean > -990, abs(NPDECminMean), NA))
summary(NLME_Dual_all)
summary(NONMEM_all)
summary(NLME_Std_all) 
SumNPDE <- c("Cmax_NONMEM" = mean(NONMEM_all$absNPDECmaxMean,na.rm = TRUE),
             "CmaxNLME_std" = mean(NLME_Std_all$absNPDECmaxMean,na.rm = TRUE),
             "CmaxNLME_Dual" = mean(NLME_Dual_all$absNPDECmaxMean,na.rm = TRUE),
             "Cmin_NONMEM" = mean(NONMEM_all$absNPDECminMean,na.rm = TRUE),
             "Cmin_NLME_std" = mean(NLME_Std_all$absNPDECminMean,na.rm = TRUE),
             "Cmin_NLME_Dual" = mean(NLME_Dual_all$absNPDECminMean,na.rm = TRUE))

SumNPDE
# GmeansEst <- c(NLME_std = exp(mean(log(NLME_Std_all %>% filter(Est_time > 0)))),
#                NLME_dual = exp(mean(log(NLME_Dual_all %>% filter(Est_time > 0)))),
#             NONMEM = exp(mean(log(NONMEM_all$Est_time%>% filter(Est_time > 0)))))  
# # GmeansCov <- c(NLME = exp(mean(log(NLME_Standard$Cov_time))),
#                NLME = exp(mean(log(NONMEM$Cov_time)))) 
AllResults <- rbind(NLME_Dual_all,NLME_Std_all,NONMEM_all)
ggplot(AllResults,
       aes(x=n_theta,y=Est_time,color = Algorithm)) +
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_log10() +
  xlab("Number of estimated THETA parameters")+
  ylab("Estimation time (seconds)") 
ggsave("Estimation Run time vs NTheta by algorithm.jpeg",device="jpeg",height = 7, width=10)
ggplot(AllResults,
       aes(x=n_omega,y=Est_time,color = Algorithm))+
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_log10()  +
  xlab("Number of estimated OMEGA parameters")+
  ylab("Estimation time (seconds)") 
ggsave("Estimation Run time vs NOmegas by algorithm.jpeg",device="jpeg",height = 7, width=10)
ggplot(AllResults ,
       aes(x=nparms,y=Est_time,color = Algorithm)) +
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_log10() +
  xlab("Number of estimated parameters (NTHETA+NOMEGA)")+
  ylab("Estimation time (seconds)")  
ggsave("Estimation Run time vs NParms by algorithm.jpeg",device="jpeg",height = 7, width=10)
ggplot(AllResults,
       aes(x=n_theta,y=NPDECmaxMean,color = Algorithm)) +
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_log10() +
  xlab("Number of estimated parameters (NTHETA+NOMEGA)")+
  ylab("Cmax NPDE Mean") 
ggsave("Cmax NPDE Mean vs NParms by algorithm.jpeg",device="jpeg",height = 7, width=10)

ggplot(AllResults,
       aes(x=n_theta,y=NPDECminMean,color = Algorithm)) +
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_log10() +
  xlab("Number of estimated parameters (NTHETA+NOMEGA)")+
  ylab("Cmin NPDE Mean") 
ggsave("Cmin NPDE Mean vs NParms by algorithm.jpeg",device="jpeg",height = 7, width=10)


ggplot(AllResults,
       aes(x=n_theta,y=NPDECmaxSD,color = Algorithm)) +
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_log10() +
  xlab("Number of estimated parameters (NTHETA+NOMEGA)")+
  ylab("Cmax NPDE SD") 
ggsave("Cmax NPDE SD vs NParms by algorithm.jpeg",device="jpeg",height = 7, width=10)

ggplot(AllResults,
       aes(x=n_theta,y=NPDECminSD,color = Algorithm)) +
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_log10() +
  xlab("Number of estimated parameters (NTHETA+NOMEGA)")+
  ylab("Cmin NPDE SD") 
ggsave("Cmin NPDE SD vs NParms by algorithm.jpeg",device="jpeg",height = 7, width=10)



ggplot(AllResults,
       aes(x=n_theta,y=Iterations,color = Algorithm)) +
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_log10() +
  xlab("Number of estimated parameters (NTHETA+NOMEGA)")+
  ylab("Number of iterations") 
ggsave("Number of iterations vs NParms by algorithm.jpeg",device="jpeg",height = 7, width=10)
# 
cnames <- colnames(NLME_Dual_all)  
colnames(NLME_Dual_all) <- paste0(cnames,"_Dual")
colnames(NLME_Std_all) <- paste0(cnames,"_Std")

colnames(NONMEM_all) <- paste0(cnames,"_NONMEM")
wide_data <- cbind(NONMEM_all, NLME_Std_all ,  NLME_Dual_all) %>%
  mutate(Check_ModelNum = (Model_num_Dual== Model_num_Std)) %>%
  mutate(Check_nparms = ((n_theta_Dual+n_omega_Dual)== (n_theta_Std+n_omega_Std))) %>%
  mutate(Ratio_Dual_Std = Est_time_Dual/Est_time_Std)  %>%  
  mutate(NParms=n_theta_Dual+n_omega_Dual)
wide_data %>% filter(Ratio_Dual_Std < 1)
nmtest <- AllResults %>% filter(Algorithm !="NLME_dual")
dualtest <- AllResults %>% filter(Algorithm !="NONMEM")
#stats
t.test(formula = Est_time ~ Algorithm,  # Formula
       data = nmtest) # Dataframe containing the variables
t.test(formula = Cov_time ~ Algorithm,  # Formula
       data = nmtest) # Dataframe containing the variables


t.test(formula = Est_time ~ Algorithm,  # Formula
       data = dualtest) # Dataframe containing the variables
t.test(formula = Cov_time ~ Algorithm,  # Formula
       data = dualtest) # Dataframe containing the variables

table_NONMEM <- table(nmtest$Success,
                      nmtest$Algorithm)
chisq.test(table_NONMEM)
table_NONMEM <- table(nmtest$Covar,
                      nmtest$Algorithm)
table_NONMEM
chisq.test(table_NONMEM)

table_DUAL <- table(dualtest$Success,
                    dualtest$Algorithm)
table_DUAL
chisq.test(table_DUAL)
table_DUAL <- table(dualtest$Covar,
                    dualtest$Algorithm)
table_DUAL
chisq.test(table_DUAL)
any(wide_data$Check_nparms==FALSE)
any(wide_data$Check_ModelNum==FALSE)
ggplot(wide_data,aes(x=factor(NParms),y=Ratio_Dual_Std)) +
  geom_boxplot()+
  ylab("Ratio of Estimation time, Dual Numbers/Standard")+
  xlab("Number of Estimated Parameters (THETA+OMEGA)") +
  scale_y_log10()
ggsave("Est ratio vs NParms.jpeg",device="jpeg",width=9,height = 6)
ggplot(wide_data,aes(x=factor(n_theta_Dual),y=Ratio_Dual_Std)) +
  geom_boxplot()+
  ylab("Ratio of Estimation time, Dual Numbers/Standard")+
  xlab("Number of Estimated THETA") +
  scale_y_log10()
ggsave("Est ratio vs NTheta.jpeg",device="jpeg",width=9,height = 6)

 
 