library(dplyr)
library(ggplot2)
NONMEM_data <- read.csv("NONMEMresults.csv")
# calculate status = 7 if crash/didn't complete, 1 is success, 4 if errors
NONMEM_data <- NONMEM_data %>%
  mutate(Status = if_else(messages == 37,1,4))  %>%
  mutate(Status = if_else(messages == -99999,7,Status))  %>%
  relocate(Status, .after = Success)

NONMEM_data_filtered <- NONMEM_data %>%
  filter(Status < 7)
NLME_standard_data <- read.csv("NLMEResults_Standard.csv")
# update algorithm name
NLME_standard_data$Algorithm <- "NLME Standard"
NLME_standard_data_filtered <- NLME_standard_data

NMIDs <- NONMEM_data_filtered %>%
  select(Model_num)
NLMEIDs <- NLME_standard_data_filtered %>%
  select(Model_num)
all_ids <- inner_join(NMIDs, NLMEIDs)

plot_data <- rbind(NONMEM_data_filtered,NLME_standard_data_filtered)

plot_data <- plot_data %>% select(Model_num,n_theta,n_omega,Algorithm,
                                 Est_time,Cov_time)
plot_data <- plot_data %>%
  filter(Model_num %in% all_ids$Model_num)

# one additional THETA for NONMEM, for residual error
plot_data <- plot_data %>%
  mutate(n_theta = if_else(Algorithm == "NONMEM",n_theta + 1, n_theta))

colnames(plot_data) <- c("Model_num","Number of Thetas","Number of  Omegas","Algorithm","Estimation","Covariance" )

library(reshape)
# need to collapse in to just Time, nparms for Est_time/Cov_time and n_theta, neta
new_a <- melt(plot_data, variable_name = c("Parameter"),id = c("Model_num","Algorithm"),
          measure.vars = c("Estimation","Covariance"))
new_b <- melt(plot_data, variable_name = c("n_parms"),id = c("Model_num","Algorithm"),
              measure.vars = c("Number of Thetas","Number of  Omegas"))
new_data <- inner_join(new_a,new_b,by=c("Model_num","Algorithm"),relationship = "many-to-many")
colnames(new_data) <- c("Model_num" ,"Algorithm" ,"Parameter", "Time","which_parm","n_parms")
# convert from seconds to  minutes
new_data <- new_data %>%
  mutate(Time = Time/60)
# remove crashes
new_data <-  new_data %>%
  filter(Time > 0)
TimeNParm <- ggplot(new_data) +
  geom_point(aes(x=n_parms, y = Time,colour = Algorithm, fill=Algorithm, shape=Algorithm), size = 1.2) +
  geom_smooth(aes(n_parms,Time,col=Algorithm),method="lm",se=FALSE, linewidth=0.8)+
  facet_grid(Parameter~which_parm, scales = "free") +
  ylab("Time (minutes)") +
  theme(strip.text = element_text(size = 10),
        axis.title.x= element_blank(),
        legend.position="bottom")+
  scale_y_log10(breaks = c(0.01, 0.1, 1,10,100),labels = c(0.01, 0.1, 1,10,100))

TimeNParm
ggsave("TimeNParm.png",
       plot= TimeNParm,
       device="png",
       width=6,
       height=4)

######################## fraction converged
all_data <- rbind(NONMEM_data,
                  NLME_standard_data)

table_data <- all_data %>%
  dplyr::select(n_theta, n_omega, Success,Covar,Algorithm) %>%
  group_by(Algorithm) %>%
  summarise(Converge = 100*mean(Success),Covar= 100*mean(Covar))
print(table_data)
write.csv(table_data,"Table.csv")
