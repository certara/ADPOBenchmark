GetNMParms <- function(xml_file){
  data <- read_xml(xml_file, encoding = "ASCII")
  control_node <- data %>%
    xml_find_all("//nm:control_stream") %>%
    xml_text()  %>%
    str_split("\n")
  data_file <- control_node[[1]][grep("^\\$DATA",control_node[[1]])] %>%
    str_replace("\\$DATA ","") %>%
    str_replace("IGNORE\\=@","") %>%
    str_replace("REWIND\\=@","") %>%
    str_replace("ACCEPT\\=@","") %>%
    str_trim()

  iterations_node <- xml_find_all(data, "//nm:monitor")
  if(length(iterations_node) == 0){
    # something failed
    iterations <- -999
  }else{
   iterations_children <- xml_children(iterations_node)[[2]]
   iterations <- as.integer(xml_attr(iterations_children,"iteration"))
  }
messages_node <- xml_find_all(data, "//nm:termination_txtmsg")
messages_children <- xml_children(messages_node)
message_contents <- xml_contents(messages_node)
messages <- as.numeric(xml_text(message_contents))
if(length(xml_contents(messages_node)) == 0){
  finished <- FALSE
}else{
  finished <- TRUE
}
if(37 %in% messages){
  success <- TRUE
}else{
  success <- FALSE
}
message_char <- paste(messages, collapse = "|")
status <- data %>%
  xml_find_all("//nm:termination_status") %>%
  xml_text() %>%
  as.numeric()
Est_time <- data %>%
  xml_find_all("//nm:estimation_elapsed_time") %>%
  xml_text() %>%
  as.numeric()

Cov_time <- data %>%
  xml_find_all("//nm:covariance_elapsed_time") %>%
  xml_text() %>%
  as.numeric()

covariance_node <- xml_find_all(data, "//nm:covariance")
if (length(covariance_node) == 0) {
  Covar <- FALSE
}else{
  Covar <- TRUE
}

theta <- xml_find_all(data, "//nm:theta") %>%
  xml_children() %>%
  xml_text()

omega <- xml_find_all(data, "//nm:omega") %>%
  xml_children()

omega_dim <- length(omega)
omega_vec <-   xml_children(omega) %>%
  xml_text()
# omega filled down rows first
omega <- matrix(0,nrow=omega_dim,ncol=omega_dim)
cur_omega <- 0
for(thiscol in 1:omega_dim){
  for(thisrow in 1:thiscol){
    cur_omega <- cur_omega + 1
    omega[thiscol,thisrow] <- as.numeric(omega_vec[cur_omega])
  }
}

sigma <- xml_find_all(data, "//nm:sigma") %>%
  xml_children()

sigma_dim <- length(sigma)
sigma_vec <-   xml_children(sigma) %>%
  xml_text()
# sigma filled down rows first
sigma <- matrix(0,nrow=sigma_dim,ncol=sigma_dim)
cur_sigma <- 0
for(thiscol in 1:sigma_dim){
  for(thisrow in 1:thiscol){
    cur_sigma <- cur_sigma + 1
    sigma[thiscol,thisrow] <- as.numeric(sigma_vec[cur_sigma])
  }
}

if(Covar){
  cov_node <- xml2::xml_find_all(data, "//nm:covariance")
  children <- xml_children(cov_node)
  dim <- xml_length(cov_node, only_elements = TRUE)
  xml_text(children[[2]])
  cov <- matrix(-999, nrow=dim,ncol=dim)
  cur_row <- 0
  for(this_row in children){
    row_children <- xml_children(this_row)
    cur_row <-  cur_row + 1
    cur_col <- 0
    for(this_col in row_children){
      cur_col <- cur_col + 1
      cov[cur_row,cur_col] <- cov[cur_col,cur_row] <- as.numeric(xml_text(this_col))
    }
  }
  # remove 0 cols and rows
  zero_rows = apply(cov, 1, function(row) any(row != 0 ))

  CovarMat <-  cov[zero_rows,zero_rows]

  theta_se <- xml_find_all(data, "//nm:thetase") %>%
    xml_children() %>%
    xml_text()

  omegase <- xml_find_all(data, "//nm:omegase") %>%
    xml_children()

  omegase_vec <-   xml_children(omegase) %>%
    xml_text()
  # omega filled down rows first
  omega_se <- matrix(0,nrow=omega_dim,ncol=omega_dim)
  cur_omega <- 0
  for(thiscol in 1:omega_dim){
    for(thisrow in 1:thiscol){
      cur_omega <- cur_omega + 1
      omega_se[thiscol,thisrow] <- as.numeric(omegase_vec[cur_omega])
    }
  }

  sigma_se <- xml_find_all(data, "//nm:sigmase") %>%
    xml_children()

  sigmase_vec <-   xml_children(sigma_se) %>%
    xml_text()
  # sigma filled down rows first
  sigma_se <- matrix(0,nrow=sigma_dim,ncol=sigma_dim)
  cur_sigma <- 0
  for(thiscol in 1:sigma_dim){
    for(thisrow in 1:thiscol){
      cur_sigma <- cur_sigma + 1
      sigma_se[thiscol,thisrow] <- as.numeric(sigmase_vec[cur_sigma])
    }
  }
}else{
  theta_se <- NULL
  omega_se <- NULL
  sigma_se <- NULL
  CovarMat <-  -999
}
OFV <-   xml_find_all(data, "//nm:final_objective_function")  %>%
  xml_text() %>%
  as.numeric()
Max_Eigen <- xml_find_all(data, "//nm:covariance_status")
return(list(
  OFV = OFV,
  theta = theta,
  omega = omega,
  sigma = sigma,
  theta_se = theta_se,
  omega_se = omega_se,
  sigma_se = sigma_se,
  dataFile = data_file,
  iterations = iterations,
  messages = messages,
  success = success,
  covar= Covar,
  CovStatus = Max_Eigen,
  CovarMat = CovarMat,
  EstTime = Est_time,
  CovTime = Cov_time
))
}
