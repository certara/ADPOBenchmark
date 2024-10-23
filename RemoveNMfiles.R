RemoveNMfiles <- function(which_dir){
  # only keep .mod and .xml
  if(file.exists(file.path(which_dir/"temp_dir"))){
    unlink(file.path(which_dir/"temp_dir"), recursive = TRUE)
  }
   FileList <- list.files(path = mydir, pattern = pattern, full.names = TRUE, all.files = TRUE )

}

RemoveNLMEfiles <- function(which_dir){
  # only keep .mod and .xml
  if(file.exists(file.path(which_dir/"temp_dir"))){
    unlink(file.path(which_dir/"temp_dir"), recursive = TRUE)
  }
  FileList <- list.files(path = mydir, pattern = pattern, full.names = TRUE, all.files = TRUE )
  
}