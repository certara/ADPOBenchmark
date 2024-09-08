# places where path is hardcoded:
# pydarwin_makedata

rm(list=ls())
home_dir <- getwd()
library(Metrics)
source(file.path(home_dir,"make_data.R"))
source(file.path(home_dir,"check_data.R"))
source(file.path(home_dir,"ReadNM_xml.R"))
source(file.path(home_dir,"CompileResultsNONMEM.R"))
source(file.path(home_dir,"CleanUp.R"))
source(file.path(home_dir,"CompileResultsNLME.r"))
source(file.path(home_dir,"GetTrueParms.r"))
source(file.path(home_dir,"DATA","make_data.r"))
source(file.path(home_dir,"CopyNONMEMControls.R"))
source(file.path(home_dir,"CopyNLMEControls.R"))
# no need to recreate models
# this freezes when run from here
# need to step through r file, run pydarwin from command line
#make_data(home_dir)
# make NONMEM control files
# note that the run directory in options.json must be an absolute path
# this must be set manually, current value is c:\git\adpoBenchmark\nonmem\run
#setwd(file.path(home_dir,"pydarwin_nonmem"))
#shell("python -m darwin.run_search_in_folder .")
# uncomment $EST and $COV copy to home_dir/NONMEM/run?.mod
CopyNONMEMControls(home_dir)
CopyNLMEControls(home_dir)
for(i in 1:72){
  check_data(home_dir,i)
  }
run_nlme(home_dir,"standard")
run_nlme(home_dir,"ADPO")
run_NONMEM(home_dir)
NONMEM_data <- read.csv(file.path(home_dir,"NONMEMResults.csv"))
# Append NM to column names
colnames(NONMEM_data) <- paste0("NM",colnames(NONMEM_data))
NLME_standard_data <- read.csv(file.path(home_dir,"NLMEResults_standard.csv"))

colnames(NLME_standard_data) <- paste0("NLME_ST",colnames(NLME_standard_data))
NLME_standard_data <- NLME_standard_data %>%
  distinct(NLME_STModel_num, .keep_all = TRUE) %>%
  filter(NLME_STModel_num<=8)
# append NLME_ST to column names
NLME_ADPO_data <- read.csv(file.path(home_dir,"NLMEResults_ADPO.csv")) %>%
  dplyr::distinct(Model_num,.keep_all=TRUE)
# append NLME_ADPO to column names
# data set for run time ratio/RMSE ratio/MAE ratio
all_data <- cbind(NONMEM_data,
                  NLME_standard_data)
# check that genome is the same
check <- all_data %>%
  mutate(This_compOK = (NMthis_comp == NLME_STthis_comp),
         This_etaOK = (NMthis_eta == NLME_STthis_eta),
         This_thisvwtOK = (NMthis_vwt == NLME_STthis_vwt),
         This_gammaOK = (NMthis_gamma == NLME_STthis_gamma)) %>%
  select(This_compOK,This_etaOK,This_thisvwtOK,This_gammaOK)
print(check)
ratios <-  all_data %>%
  mutate(NM_NLMEStEst = NMEst_time/NLME_STEst_time,
         NM_NLMEStCOV = NMCov_time/NLME_STCov_time,
         NM_NLMERMSE = NMRMSE/NLME_STRMSE,
         NM_NLMEMAE = NMMAE/NLME_STMAE) %>%
  select(NM_NLMEStEst,NM_NLMEStCOV,NM_NLMERMSE,NM_NLMEMAE)
library(ggplot2)
NM_NLMEST_ESTTime <- ggplot(all_data) +
  geom_histogram(aes(x=NMEst_time,colour = NMAlgorithm,fill=NMAlgorithm),alpha = 0.4) +
  geom_histogram(aes(x=NLME_STEst_time,colour = NLME_STAlgorithm,fill=NLME_STAlgorithm),alpha = 0.4) +
  xlab("Estimation time (seconds)")
ggsave("NMvsNLMESTEstTime.jpeg",
       plot= NM_NLMEST_ESTTime,
       device="jpeg",
       width=8,
       height=5)
print(NM_NLMEST_ESTTime)
NM_NLMEST_COVTime <- ggplot(all_data) +
  geom_histogram(aes(x=NMCov_time,colour = NMAlgorithm,fill=NMAlgorithm),alpha = 0.4) +
  geom_histogram(aes(x=NLME_STCov_time,colour = NLME_STAlgorithm,fill=NLME_STAlgorithm),alpha = 0.4) +
  xlab("Estimation time (seconds)")
ggsave("NMvsNLMESTCovTime.jpeg",
       plot= NM_NLMEST_ESTTime,
       device="jpeg",
       width=8,
       height=5)
print(NM_NLMEST_COVTime)
# ratio of each "true" parameters

TrueParms <- GetTrueParms(home_dir) %>%
  dplyr::summarise(VmaxGeoMean = exp(mean(log(VMAX))),
                   KMGeoMean = exp(mean(log(KM))),
                   VGeoMean = exp(mean(log(V2))),
                   KAGeoMean = exp(mean(log(KA))),
                   VmaxCV = 100*sd(log(VMAX)),
                   KMCV = 100*sd(log(KM)),
                   VCV = 100*sd(log(V2)),
                   KACV = 100*sd(log(KA)))
fraction <- all_data %>%
  mutate(NM_NLMEStEst = NMEst_time/NLME_STEst_time,
TrueParms
