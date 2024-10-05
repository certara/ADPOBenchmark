 jobType <-"GENERIC"
 num_samples <- 1
 exeFileExtension <- ".ps1"
 parallelMethod <- "none"
 num_processes <- 1
 localWorkingDir <- "D:/git/ADPOBenchmark/NLME/1/Run1"
 SharedWorkingDir <- "D:/git/ADPOBenchmark/GEN241005094744.81/Shared"
 files_to_copy <- "jobArgsCombined.txt data1.txt cols1.txt test.mdl"
 model_file <- "test.mdl"
 extra_args_file <- "jobArgsCombined.txt"
 MpiExecutable <- "NLME7.exe"
 MpiArgument <- "MPINO"
 MpiLocal <- "NO"
 MpiNumCores <- 1
 assign("MpiExecutable", MpiExecutable, envir = nlmeEnv)
 assign("MpiArgument", MpiArgument, envir = nlmeEnv)
 assign("MpiLocal", MpiLocal, envir = nlmeEnv)
 assign("MpiNumCores", MpiNumCores, envir = nlmeEnv)

 control_lines <- c('WorkFlow,jobArgsCombined.txt:1,*,out.txt,err1.txt err2.txt ps1err.txt fort.27 integration_errors.txt nlme7engine.log out.txt dmp.txt,  ')
runNLMESample <- function(indx, eArgsFile, ofn, extraArgs="", seed=-1, max_tries=1, exePostfix = "", SharedWorkingDir = SharedWorkingDir)
{
    if (SharedWorkingDir == "") {
        SharedWorkingDir <- get("SharedWorkingDir", envir = nlmeEnv)
    }
    source(file.path(SharedWorkingDir, "myglobaldefs.r"), local = TRUE)
    if (missing(eArgsFile)) {
        eArgsFile <- extra_args_file
    }
    if (dirname(eArgsFile) == ".") {
        eArgsFile <- file.path(SharedWorkingDir, eArgsFile)
    }
    if (jobType == "BOOTSTRAP" && seed != -1) {
        if (.Platform$OS.type == "unix") {
            extraArgs <- sprintf("/boot ${bootSeed} /bootsamp %d -boottry ${numTries}", 
                indx)
        }
        else {
            extraArgs <- sprintf("/boot $bootSeed /bootsamp %d -boottry $numTries", 
                indx)
        }
    }
    baseDirectory <- SharedWorkingDir
    if (tolower(parallelMethod) %in% c("multicore", "none", "local_mpi", 
        "mpi")) {
        baseJobDirectory <- baseDirectory
    }
    else {
        baseJobDirectory <- file.path(baseDirectory, "registry")
    }
    indx02d <- sprintf("%02d", indx%%100)
    indx03d <- sprintf("%03d", indx)
    if (indx == 0) {
        workingDir <- baseJobDirectory
        newFilePath <- file.path(baseJobDirectory, paste0("exNLME", 
            exeFileExtension))
    }
    else {
        workingDir <- file.path(baseJobDirectory, "jobs", indx02d, 
            indx)
        dir.create(workingDir, recursive = TRUE, showWarnings = FALSE)
        newFilePath <- file.path(workingDir, paste0("exNLME", 
            exeFileExtension))
        files_to_copy <- unlist(strsplit(files_to_copy, split = " "))
        copy_filesWarnLong(file.path(SharedWorkingDir, files_to_copy), 
            workingDir, overwrite = TRUE)
        copy_filesWarnLong(eArgsFile, workingDir, overwrite = TRUE)
    }
    statusFile <- file.path(baseDirectory, paste0("S_", indx03d, 
        ".status"))
    statusBackupFile <- file.path(baseDirectory, paste0("S_", 
        indx03d, ".status.bak"))
    logFile <- file.path(workingDir, ".status")
    outFile <- file.path(workingDir, ofn)
    if (nchar(eArgsFile) > 0) {
        extraArgsFile <- paste0("@", basename(eArgsFile))
    }
    MpiExecutable <- get("MpiExecutable", envir = nlmeEnv)
    MpiArgument <- get("MpiArgument", envir = nlmeEnv)
    MpiNumCores <- get("MpiNumCores", envir = nlmeEnv)
    MpiLocal <- get("MpiLocal", envir = nlmeEnv)
    MpiExecutable <- gsub(".exe", paste0(exePostfix, ".exe"), 
        MpiExecutable, fixed = TRUE)
    if (.Platform$OS.type == "unix") {
        argsToPaste <- shQuote(paste(extraArgs, extraArgsFile), 
            type = "cmd")
        commandString <- paste(shQuote(file.path(baseDirectory, 
            paste0("execNLMECmd", exeFileExtension))), "RUN", 
            shQuote(file.path(baseDirectory, model_file), type = "cmd"), 
            shQuote(workingDir, type = "cmd"), MpiArgument, MpiLocal, 
            MpiNumCores, "\"\"", "NLME_DIR", "\"\"", argsToPaste, 
            "\"\"", exePostfix)
        NLMEexeLinking <- ifelse(baseDirectory != workingDir, 
            paste0("\nln -s ", shQuote(file.path(baseDirectory, 
                MpiExecutable), type = "cmd"), " ", shQuote(workingDir, 
                type = "cmd")), "")
        cat(paste0("#!/bin/bash", "\nset -x", "\ndeclare -i numTries", 
            "\ndeclare -i bootSeed", "\ncd ", shQuote(workingDir), 
            "\nrm -f ", shQuote(outFile), NLMEexeLinking, "\nchmod 777 ", 
            MpiExecutable, "\necho 'RUNNING' >", shQuote(logFile, 
                type = "cmd"), "\necho 'RUNNING' >", shQuote(statusFile, 
                type = "cmd"), "\necho '", format(Sys.time(), 
                "%Y-%m-%d %H:%M:%S %Z"), "' >>", shQuote(statusFile, 
                type = "cmd"), "\nmaxTries=", max_tries, "\nnumTries=1", 
            "\nbootSeed=", seed, "\nstartTime=`date +\"%x %T\"`", 
            "\nwhile [ ${numTries} -le ${maxTries} ]; do", "\n  ", 
            commandString, "\n  stopTime=`date +\"%x %T\"`", 
            "\n  ((numTries++))", "\n  bootSeed=${bootSeed}+1", 
            "\n  sleep 1", "\n  if [ -s ", shQuote(outFile, type = "cmd"), 
            " ]", "\n  then", "\n    echo 'SUCCESS' >", shQuote(logFile, 
                type = "cmd"), "\n    echo 'SUCCESS' >", shQuote(statusFile, 
                type = "cmd"), "\n    echo ${startTime} >>", 
            shQuote(statusFile, type = "cmd"), "\n    echo ${stopTime} >>", 
            shQuote(statusFile, type = "cmd"), "\n    exit", 
            "\n  fi", "\ndone", "\necho 'FAILED' >", shQuote(logFile, 
                type = "cmd"), "\necho 'FAILED' >", shQuote(statusFile, 
                type = "cmd")), file = newFilePath, sep = "\n", 
            append = FALSE)
        Sys.chmod(newFilePath, mode = "0777", use_umask = TRUE)
        Sys.chmod(list.files(path = Sys.getenv("INSTALLDIR"), 
            pattern = "*.sh|TDL5$", full.names = TRUE, recursive = TRUE), 
            mode = "777")
        if (indx != 0 && (tolower(parallelMethod) %in% c("multicore", 
            "none", "local_mpi", "mpi"))) {
            system(paste("nohup", shQuote(newFilePath, type = "cmd"), 
                " > ", shQuote(file.path(workingDir, "log"), 
                  type = "cmd"), " 2>&1  &"))
        }
        else {
            system(paste0(shQuote(newFilePath, type = "cmd"), 
                " > ", shQuote(file.path(workingDir, "log"), 
                  type = "cmd"), " 2>&1 ", indx, " "))
        }
    }
    else {
        argsToPaste <- shQuote(paste(extraArgs, extraArgsFile))
        argsToPaste <- paste0("\"\"\"", argsToPaste, "\"\"\"")
        commandString <- paste("RUN", shQuote(file.path(baseDirectory, 
            model_file)), shQuote(workingDir), MpiArgument, MpiLocal, 
            MpiNumCores, "\"\"\"\"\"\"", "NLME_DIR", "\"\"\"\"\"\"", 
            argsToPaste, "\"\"\"\"\"\"", exePostfix)
        copy_filesWarnLong(file.path(Sys.getenv("INSTALLDIR"), 
            "NLMESample.ps1"), newFilePath, overwrite = TRUE)
        if (indx != 0) {
            copy_filesWarnLong(file.path(baseDirectory, MpiExecutable), 
                workingDir, overwrite = TRUE)
        }
        ps_args <- c("-noninteractive", "-executionpolicy", "remotesigned", 
            "-File", shQuote(newFilePath), shQuote(gsub("/", 
                "\\", logFile, fixed = TRUE)), shQuote(gsub("/", 
                "\\", statusFile, fixed = TRUE)), shQuote(gsub("/", 
                "\\", statusBackupFile, fixed = TRUE)), shQuote(gsub("/", 
                "\\", outFile, fixed = TRUE)), paste("-bootseed ", 
                seed), max_tries, shQuote(baseDirectory), shQuote(commandString))
        ps_wait <- (indx == 0) || !(tolower(parallelMethod) %in% 
            c("multicore", "none", "local_mpi", "mpi"))
        cat("powershell ", paste0(ps_args, " ", collapse = " "), 
            file = file.path(workingDir, "mpilog.txt"))
        system2("powershell", args = ps_args, stdout = file.path(workingDir, 
            "ps1log.txt"), stderr = file.path(workingDir, "ps1err.txt"), 
            wait = ps_wait)
    }
}
getExtraArgumentFilename<- function(line)
{
    fileRec <- unlist(strsplit(line, split = ","))[2]
    tokens <- unlist(strsplit(fileRec, split = ":"))
    if (length(tokens) > 2) {
        file <- paste0(tokens[1], ":", tokens[2])
    }
    else {
        file <- tokens[1]
    }
    file
}
getRunSuccessFilename<- function(control_lines_given)
{
    files <- sapply(strsplit(c(control_lines_given), split = ","), 
        function(x) x[4])
    files
}
getExtraArgumentFilenameIndex<- function(line)
{
    fileRec <- unlist(strsplit(line, split = ","))[2]
    tokens <- unlist(strsplit(fileRec, split = ":"))
    if (length(tokens) > 2) {
        fileIndex <- tokens[3]
    }
    else {
        fileIndex <- tokens[2]
    }
    fileIndex
}
getExePostfix<- function(line)
{
    exePostfix <- gsub("^\\s+|\\s+$", "", (unlist(strsplit(line, 
        split = ","))[7]))
    if (is.na(exePostfix)) {
        exePostfix <- ""
    }
    exePostfix
}
