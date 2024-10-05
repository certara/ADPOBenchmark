NONMEM_data <- read.csv(file.path(home_dir,"NONMEMResults.csv"))
# Append NM to column names
colnames(NONMEM_data) <- paste0("NM",colnames(NONMEM_data))
NONMEM_data <- filter(NONMEM_data, NONMEM_data$NMEst_time > 0)
NLME_standard_data <- read.csv(file.path(home_dir,"NLMEResults_standard.csv"))

colnames(NLME_standard_data) <- paste0("NLME_ST",colnames(NLME_standard_data))
NLME_standard_data <- NLME_standard_data %>%
  distinct(NLME_STModel_num, .keep_all = TRUE)

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
