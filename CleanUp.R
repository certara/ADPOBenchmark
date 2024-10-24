CleanUpNM <- function(wd) {
  if (file.exists(file.path(wd, "temp_dir")))
    unlink(file.path(wd, "temp_dir"),
           recursive = TRUE,
           force = TRUE)
  files <- grep(
    list.files(
      path = wd,
      full.names = TRUE,
      recursive = TRUE
    ),
    pattern = "Run[0-9].*",
    invert = TRUE,
    value = TRUE
  )
  cfiles <-
    list.files(
      path = wd,
      pattern = "co.$",
      full.names = TRUE,
      recursive = TRUE  )
  sfiles <-
    list.files(
      path = wd,
      pattern = "sh.$",
      full.names = TRUE,
      recursive = TRUE  )
  file.remove(c(files,cfiles,sfiles))
}
CleanUpNLME <- function(wd) {
  # just the big stuff
  # keep err? and .log file

  files <- grep(
    list.files(
      path = wd,
      full.names = TRUE,
      recursive = TRUE
    ),
    pattern = "csv",
    value = TRUE
  )

  file.remove(files)
  if(file.exists(file.path(wd,"dmp.txt"))) file.remove(file.path(wd,"dmp.txt"))
  if(file.exists(file.path(wd,"out.txt"))) file.remove(file.path(wd,"out.txt"))
}
