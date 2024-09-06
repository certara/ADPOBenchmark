CleanUp <- function(dir){
  if(file.exists(file.path(dir,"temp_dir"))) unlink(file.path(dir,"temp_dir"),recursive = TRUE, force = TRUE)
  files <- list.files(dir,
                      full.names = TRUE,
                      recursive = TRUE)

  Savedfiles <- files[grep("Run[0-9].*",files)]
  files <-  files[! files %in% Savedfiles]
  for(file in files){
    if(file.exists(file)) file.remove(file)
  }

}
