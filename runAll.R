rm(list=ls())
home_dir <- getwd()

# used by NLME
if (Sys.info()['sysname'] == "Linux") {
  pyDarwinInterpreter <- "/home/user/venv/bin/python3"
  nmfe_path <- "/opt/nm751/util/nmfe75"

  gcc_dir <- dirname(system("which gcc", intern = TRUE))
  nlme_dir <- "/home/user/InstallDirNLME/"
  if (grepl("Ubuntu", Sys.info()["version"])) {
    Sys.setenv("PML_BIN_DIR" = "UBUNTU2204")
  }

  INSTALLDIRStandard <- nlme_dir
  INSTALLDIRADPO <- "home/user/NLME_Engine_ADPO"
  INSTALLDIRnoHessian <- "home/user/NLME_Engine_noHessian"

} else {
  pyDarwinInterpreter <- "C:/python/venv/bin/python"
  nmfe_path <- "C:/nm751/util/nmfe75"

  # for Windows NLME engine
  gcc_dir <- "C:\\Program Files\\Certara\\mingw64"
  nlme_dir <- "C:\\Program Files\\Certara\\NLME_Engine"

  INSTALLDIRStandard <- nlme_dir
  INSTALLDIRADPO <- "D:/NLME_Engine_ADPO"
  INSTALLDIRnoHessian <- "D:/NLME_Engine_noHessian"
}

nlme_dirs <- c(standard = INSTALLDIRStandard,
               ADPO = INSTALLDIRADPO,
               NoHessian = INSTALLDIRnoHessian)

library(Metrics)
library(ggplot2)
library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(lubridate)
library(readtext)
library(Certara.RsNLME)
library(Certara.RDarwin)
library(xml2)

source(file.path(home_dir,"make_data.R"))
source(file.path(home_dir,"check_data.R"))
source(file.path(home_dir,"ReadNM_xml.R"))
source(file.path(home_dir,"CompileResultsNONMEM.R"))
source(file.path(home_dir,"CleanUp.R"))
source(file.path(home_dir,"CompileResultsNLME.R"))
source(file.path(home_dir,"GetTrueParms.R"))

make_data(home_dir, pyDarwinInterpreter, nmfe_path, nlme_dir, gcc_dir)
# make NONMEM control files
for(i in 1:72){
  check_data(home_dir,i)
}

run_NONMEM(home_dir, nmfe_path)
NONMEM_data <- read.csv(file.path(home_dir,"NONMEMResults.csv"))
# Append NM to column names
colnames(NONMEM_data) <- paste0("NM",colnames(NONMEM_data))
# filter not finished
NONMEM_data <- filter(NONMEM_data, NONMEM_data$NMEst_time > 0)

run_nlme(home_dir,"standard", nlme_dirs)
NLME_standard_data <- read.csv(file.path(home_dir,"NLMEResults_standard.csv"))
colnames(NLME_standard_data) <- paste0("NLME_ST",colnames(NLME_standard_data))
NLME_standard_data <- NLME_standard_data %>%
  distinct(NLME_STModel_num, .keep_all = TRUE)

run_nlme(home_dir,"NoHessian", nlme_dirs)
run_nlme(home_dir,"ADPO", nlme_dirs)

# append NLME_ST to column names
NLME_ADPO_data <- read.csv(file.path(home_dir,"NLMEResults_ADPO.csv")) %>%
  dplyr::distinct(Model_num,.keep_all=TRUE)
# append NLME_ADPO to column names
# data set for run time ratio/RMSE ratio/MAE ratio
NLME_NOHessian_data <- NLME_standard_data %>%
  distinct(NLME_STModel_num, .keep_all = TRUE) %>%
  filter(NLME_STModel_num<=8)
# append NLME_ST to column names

all_data <- inner_join(NONMEM_data,
                       NLME_standard_data,
                       by = c("NMModel_num" = "NLME_STModel_num"))

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
