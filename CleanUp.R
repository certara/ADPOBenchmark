CleanUp <- function(wd) {
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

  file.remove(files)
}
