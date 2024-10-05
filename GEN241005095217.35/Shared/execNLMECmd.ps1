#
# execNLMECmd.ps1 : Runs NLME on windows from command line.
#
# Args :
#        $RUN_MODE   = COMPILE_AND_RUN COMPILE RUN
#        $MODELFILE  = PML file to use for the model
#        $LCWD       = Full path Working directory to start on local host
#        $MPIFLAG    = MPIYES | MPINO
#        $LOCAL_HOST =  YES | NO
#        $NUM_NODES  = Number of mpi nodes
#        $SHARED_DRIVE = Location of shared drive for remote MPI
#        $RCWD         = Name of the working directory on remote/shared drive
#        
#
#                          This directory will be created on either :
#                          --  Shared Drive OR
#                          --  %USERPROFILE%
#                              on remote nodes and will be cleaned after the
#                              run is finished.
# 
#        $FILES        = List of files to copy to remote node or shared directory
#        $NLME_ARGS    = Arguments passed on to xxx.exe
# 

param(
    [string] $RUN_MODE,
    [string] $MODELFILE,
    [string] $WORKING_DIR,
    [string] $MPIFLAG,
    [string] $LOCAL_HOST,
    [Int32]  $NUM_NODES,
    [string] $SHARED_DRIVE = "",
    [string] $RCWD,
    [string] $FILES,
    [string] $NLME_ARGS,
    [string] $CMD_HASHCODE = "",
    [string] $NLME_EXE_POSTFIX = "",
    [string] $MPI_ARGS
)

# $PSBoundParameters

$ScriptDir = Split-Path -parent $MyInvocation.MyCommand.Path
. "$ScriptDir\common.ps1"


$RCWD += Get-Random
$FILES = $FILES.Trim()

"WORKING_DIR=$WORKING_DIR, MPIFLAG=$MPIFLAG, LOCAL_HOST=$LOCAL_HOST, NUM_NODES=$NUM_NODES, SHARED_DRIVE=$SHARED_DRIVE"

[string] $HASH_TO_USE = if ($CMD_HASHCODE -ine "") {$CMD_HASHCODE} else {$env:NLME_HASH}
$HASH_TO_USE = if ($HASH_TO_USE) {"/hash $HASH_TO_USE"}

[string] $INSTALLDIR = if ($env:INSTALLDIR) {$env:INSTALLDIR -replace '\\$', ''} else {"$env:PhoenixDir\application\lib\NLME\Executables"}

[string] $NLMEGCCDir64 = if ($env:NLMEGCCDir64) {$env:NLMEGCCDir64 -replace '\\$', ''} else {"C:\PHSTMinGW64"}
[string] $PhoenixMPIDir = if ($env:PhoenixMSMPIDir) {$env:PhoenixMSMPIDir -replace '\\$', ''} else {$env:MSMPI_BIN}

[string] $DEBUG_FLAG = if ($env:NLME_BUILD_DEBUG) {"-O2"} else {"-g"}
[string] $BUILD_FLAGS="-std=c++11 -m64 -I""$INSTALLDIR"" -malign-double -fverbose-asm -c"

[string] $CCOMPILER="$NLMEGCCDir64\bin\gcc.exe"
[string] $FCOMPILER="$NLMEGCCDir64\bin\gfortran.exe"

[string] $TDL5 = "$INSTALLDIR\TDL5.exe"

# Define any undefined args - note if OUTERRESULTFILE is not defined (normally the case in the examples)
# it is hard coded to out.txt below.

[string] $OUTERRESULTFILE = if ($env:OUTERRESULTFILE) {$env:OUTERRESULTFILE} else {"out.txt"}
[Int32] $SUBMODELS = if ($env:SUBMODELS) {$env:SUBMODELS} else {0}

[string] $gccpath = "$NLMEGCCDir64\bin;$NLMEGCCDir64\libexec\gcc\x86_64-w64-mingw32\8.4.0"

$env:path="$gccpath;$PhoenixMPIDir;$env:path"

Set-Location (Get-Item -LiteralPath $WORKING_DIR).FullName

if ($RUN_MODE -ine "RUN")
{
    "model=$MODELFILE, nlmeDir=$INSTALLDIR"

    # ------------------------------------------------------
    if (Test-Path .\Work -PathType container)
    {
        "Deleting files..."

        run_cmd "rmdir /S /Q .\Work"
    }

    run_cmd "md .\Work"

    die "Unable to create Work directory"

    # ------------------------------------------------------
    "Translating..."

    run_cmd """$TDL5"" $HASH_TO_USE /L .\$MODELFILE .\Work 1> log.txt 2>&1"

    die "ERROR in model translation"

    if (-not (Test-Path .\Work\Model.cpp -PathType leaf))
    {
        "ERROR in generating Model.cpp"

        exit 22
    }

    # ------------------------------------------------------
    "Compiling..."

    run_app $CCOMPILER "$BUILD_FLAGS $DEBUG_FLAG .\Work\Model.cpp"

    die "ERROR compiling Model.cpp"

    run_cmd "move /y Model.o .\Work"

    # ------------------------------------------------------
    "Linking..."

    if ($MPIFLAG -ieq "MPINO")
    {
        $mpi_lib = """$INSTALLDIR\libMPI_STUB.a"""
    }
    else
    {
        $mpi_lib = """$INSTALLDIR\libmsmpi.a"""
        $exe_prefix = "mpi"
    }

    run_app $FCOMPILER "$DEBUG_FLAG --enable-stdcall-fixup -static .\Work\Model.o ""$INSTALLDIR\libNLME7.a"" ""$INSTALLDIR\libcrlibm.a"" ""$INSTALLDIR\libNLME7_FORT.a"" ""$INSTALLDIR\libLAPACK.a"" ""$INSTALLDIR\libBLAS.a"" $mpi_lib -lstdc++ -o $($exe_prefix)NLME7$NLME_EXE_POSTFIX.exe"

    die "ERROR linking"

    if ($RUN_MODE -ieq "COMPILE")
    {
        exit 0
    }
}

if ($MPIFLAG -ieq "MPINO")
{
    $NLME_ARGS = $NLME_ARGS -replace "EXECUTION_DIR\\", ""
    $NLME_ARGS = $NLME_ARGS -replace '"', ''

    run_cmd ".\NLME7$NLME_EXE_POSTFIX.exe $NLME_ARGS 1> err1.txt 2> err2.txt"

    bye
}

if ($LOCAL_HOST -ieq "YES")
{
    $EXECUTION_DIR = "."
}
else
{
    $FILES = $FILES -replace '"', ''
    $FILES_TO_COPY = "mpiNLME7$NLME_EXE_POSTFIX.exe $FILES"

    [string[]] $filesToCopy = $FILES_TO_COPY -split " "

    $SHARED_DRIVE = $SHARED_DRIVE -replace '"', ''

    $EXECUTION_DIR = "$SHARED_DRIVE\$RCWD"
    $WD = "-wdir $EXECUTION_DIR"

    CopyFiles $filesToCopy $EXECUTION_DIR
}

$EXECUTABLE = "$EXECUTION_DIR\mpiNLME7$NLME_EXE_POSTFIX.exe"

$NLME_ARGS = $NLME_ARGS -replace "EXECUTION_DIR\\", $(if ($EXECUTION_DIR -eq ".") {""} else {"$EXECUTION_DIR\"})
$NLME_ARGS = $NLME_ARGS -replace '"', ''

# 
# Run the executable
# 

if (-not $MPI_ARGS)
{
    $MPI_ARGS = "-n $NUM_NODES $WD"
}

run_cmd """$PhoenixMPIDir\mpiexec"" $MPI_ARGS $EXECUTABLE $NLME_ARGS 1> err1.txt 2> err2.txt"

if ($EXECUTION_DIR -eq ".")
{
    bye
}

$res = $LASTEXITCODE

#
# If we created a directory to work in, then clean it up.
# First copy all files from temporary directory to the run directoy
#
"$EXECUTION_DIR"
copy $EXECUTION_DIR\* .

$cleanupCommand = "rmdir /S /Q $EXECUTION_DIR"

run_cmd $cleanupCommand

if ($LASTEXITCODE)
{
    "WARNING: unable to cleanup local host"
}

bye $res

# SIG # Begin signature block
# MIIn7wYJKoZIhvcNAQcCoIIn4DCCJ9wCAQExDTALBglghkgBZQMEAgEweQYKKwYB
# BAGCNwIBBKBrMGkwNAYKKwYBBAGCNwIBHjAmAgMBAAAEEB/MO2BZSwhOtyTSxil+
# 81ECAQACAQACAQACAQACAQAwMTANBglghkgBZQMEAgEFAAQgF7fWENzhHhWMwb4q
# sTGa2BRMkuK7XsxlyWyWtiZu7O2ggiDuMIIF3zCCBMegAwIBAgIQTkDkN1Tt5owA
# AAAAUdOUfzANBgkqhkiG9w0BAQsFADCBvjELMAkGA1UEBhMCVVMxFjAUBgNVBAoT
# DUVudHJ1c3QsIEluYy4xKDAmBgNVBAsTH1NlZSB3d3cuZW50cnVzdC5uZXQvbGVn
# YWwtdGVybXMxOTA3BgNVBAsTMChjKSAyMDA5IEVudHJ1c3QsIEluYy4gLSBmb3Ig
# YXV0aG9yaXplZCB1c2Ugb25seTEyMDAGA1UEAxMpRW50cnVzdCBSb290IENlcnRp
# ZmljYXRpb24gQXV0aG9yaXR5IC0gRzIwHhcNMjEwNTA3MTU0MzQ1WhcNMzAxMTA3
# MTYxMzQ1WjBpMQswCQYDVQQGEwJVUzEWMBQGA1UECgwNRW50cnVzdCwgSW5jLjFC
# MEAGA1UEAww5RW50cnVzdCBDb2RlIFNpZ25pbmcgUm9vdCBDZXJ0aWZpY2F0aW9u
# IEF1dGhvcml0eSAtIENTQlIxMIICIjANBgkqhkiG9w0BAQEFAAOCAg8AMIICCgKC
# AgEAp4GP9xRFtmJD8tiu0yVeSE9Rv8V9n1AcNdHWfmEqlBltJ0akphpd91RRaoAi
# xqhmdU1Ug8leaBur9ltksK2tIL1U70ZrbQLnRa519o6KuTIui7h3HFJNeYhfpToY
# yVAslyctv9oAfWN/7zLsRodj25qfw1ohNnv5m9XKoG5yLPzh8Z5wTQhWFW+Qq/tI
# urnXwYJ4hWUuf7XJwOIUtzcRQQbiiuCo9uV+tngFAcNg7U8HQS4KE0njkJt/3b36
# rL9kUdFcm7T1XOdc/zubpaAa130JssK3/24cvMh95ukr/HKzFOlKVRKEnEQldR32
# KvBPpSA9aCXrYZd8D+W2PfOuw8ERvBuOzOBHMF5CAIZx41isBsplH3uUpktXZwx+
# Xq14Z1tV417rx9jsTG6Gy/Pc+J+HqnJYEg99pvj4Qjk7PCzkMk1JjODhAMI4oJz6
# hD5B3G5WrsYaW/RnaAUBzRu/roe8nVP2Lui2a+SZ3sVPh1io0mUeyB/Vcm7uWRxX
# OwlyndfKt5DGzXtFkpFCA0x9P8ryqrjCDobzEJ9GLqRmhmhaaBhwKTgRgGBrikOj
# c2zjs2s3/+adZwGSht8vSNH7UGDVXP4h0wFCY/7vcLQXwI+o7tPBS18S6v39Lg6H
# RGDjqfTCGKPj/c4MhCIN86d42pPz2zjPuS8zxv8HPF6+RdMCAwEAAaOCASswggEn
# MA4GA1UdDwEB/wQEAwIBhjASBgNVHRMBAf8ECDAGAQH/AgEBMB0GA1UdJQQWMBQG
# CCsGAQUFBwMDBggrBgEFBQcDCDA7BgNVHSAENDAyMDAGBFUdIAAwKDAmBggrBgEF
# BQcCARYaaHR0cDovL3d3dy5lbnRydXN0Lm5ldC9ycGEwMwYIKwYBBQUHAQEEJzAl
# MCMGCCsGAQUFBzABhhdodHRwOi8vb2NzcC5lbnRydXN0Lm5ldDAwBgNVHR8EKTAn
# MCWgI6Ahhh9odHRwOi8vY3JsLmVudHJ1c3QubmV0L2cyY2EuY3JsMB0GA1UdDgQW
# BBSCutY9l86fz3Hokjev/bO1aTVXzzAfBgNVHSMEGDAWgBRqciZ60B7vfec7aVHU
# bI2fkBJmqzANBgkqhkiG9w0BAQsFAAOCAQEAH15BBLaDcCRTLFVzHWU6wOy0ewSY
# Xlk4EwmkWZRCXlC/T2xuJSCQk1hADfUZtGLuJF7CAVgVAh0QCW+o1PuSfjc4Pi8U
# fY8dQzZks2YTXxTMpXH3WyFLxpe+3JX8cH0RHNMh3dAkOSnF/goapc97ee46b97c
# v+kR3RaDCNMsjX9NqBR5LwVhUjjrYPMUaH3LsoqtwJRc5CYOLIrdRsPO5FZRxVbj
# hbhNm0VyiwfxivtJuF/R8paBXWlSJPEII9LWIw/ri9d+i8GTa/rxYntY6VCbl24X
# iA3hxkOY14FhtoWdR+yxnq4/IDtDndiiHODUfAjCr3YG+GJmerb3+sivNTCCBnAw
# ggRYoAMCAQICEHHvVXSvNVTDWixp9m9La80wDQYJKoZIhvcNAQENBQAwaTELMAkG
# A1UEBhMCVVMxFjAUBgNVBAoMDUVudHJ1c3QsIEluYy4xQjBABgNVBAMMOUVudHJ1
# c3QgQ29kZSBTaWduaW5nIFJvb3QgQ2VydGlmaWNhdGlvbiBBdXRob3JpdHkgLSBD
# U0JSMTAeFw0yMTA1MDcxOTIwNDVaFw00MDEyMjkyMzU5MDBaME8xCzAJBgNVBAYT
# AlVTMRYwFAYDVQQKEw1FbnRydXN0LCBJbmMuMSgwJgYDVQQDEx9FbnRydXN0IENv
# ZGUgU2lnbmluZyBDQSAtIE9WQ1MyMIICIjANBgkqhkiG9w0BAQEFAAOCAg8AMIIC
# CgKCAgEAnpl2FxqeVhcIYyyTTYNhhDM0ArbZYqDewg65IEzIV50P3VRbDQzWAd0v
# SOGRCeHbyBUEgrZ78NjWMDsZcXD7qKaX9ildpAyp9FM+V9sMTm78dttfJOmqX0Pj
# k+cOz8olvMRMMAtaD+YG9OVuDJlmWE+DYcJzfFwibwFFxQ/3QE9kS9AXCqkOHgIv
# oY9M8mdQ2z7kn8JPP3TrMaTQlNCZvDCSCWrLJM2i2HZS0E51mE9kWtJeg/RYwF1q
# dcTYP2Q6ixQN2Hbh6rlr5xFwSRE4YxNu8cb6vRBFNQfmdhXQdRaqwkNX/qv+Y3NG
# IqC48+THcEYJ+ak3QZqzS2wfcHKjB/Y1knQRZG75AtXAkpXxl1l+De6iJfJxVbib
# jb/N7q7d+wznrjJOUI2h39Fzv8HOf3Xaq7/QrYI4xeeI7aJtOoYRt9ew4aiLOwxB
# F5pf5FuYyJ0An/dz0sPpnwWHeSGD1gvt0cwIn+DxxclYulNf1Iexi1mo0l7NadA+
# +sQ5Ca+0te3nPPoih9Zz+ReVasMc9VV4X9T6C8BbP4x4FQ5aTDpu5SaY0CfMIN/A
# hjt6jWVGftlhXqn0rj7U/K9FxzqzhQRKi8gJXbN7AihZ44Z9gKJYQGZi4DhVg6uf
# KUEmurvp2GT4trsoc80VSteec+NmTLFRnYEji8iGd7K2LDcicCECAwEAAaOCASww
# ggEoMBIGA1UdEwEB/wQIMAYBAf8CAQAwHQYDVR0OBBYEFO+funmwc/IlHnicA1Kc
# G1OE3o3tMB8GA1UdIwQYMBaAFIK61j2Xzp/PceiSN6/9s7VpNVfPMDMGCCsGAQUF
# BwEBBCcwJTAjBggrBgEFBQcwAYYXaHR0cDovL29jc3AuZW50cnVzdC5uZXQwMQYD
# VR0fBCowKDAmoCSgIoYgaHR0cDovL2NybC5lbnRydXN0Lm5ldC9jc2JyMS5jcmww
# DgYDVR0PAQH/BAQDAgGGMBMGA1UdJQQMMAoGCCsGAQUFBwMDMEUGA1UdIAQ+MDww
# MAYEVR0gADAoMCYGCCsGAQUFBwIBFhpodHRwOi8vd3d3LmVudHJ1c3QubmV0L3Jw
# YTAIBgZngQwBBAEwDQYJKoZIhvcNAQENBQADggIBAF7zhpk1wXnu8BCv15JD0oXQ
# W+CYoOBxUckUy1CaYA6wBCZJsc/pupsIodpDXdRvI6K5+EDHR/5NAQtIkD/I3Gq0
# PlM1KL5ASkeFah53QMyAf2G0PE95qOajpqB+RIZxvxIblYFy9w2l0N5nn8aiuPFq
# +fz+dGbGZOZ5PWoDYU5LH8wgYssCGOxj7X5xP5a6C15oImfsH8DSBRZmsbKk6vzF
# laONEqX1je8bIM2Z9+cy81lxH92U5nnlUiMQVir8WTi/v3klkmrH/atnd3GxBH01
# rRTBPqj8IxdWCBh813oia5FqzDVFbU87nUOdBbid8/w0IVwEGDJXODTByoMjRqaI
# IyHGfhSAq7HvuwusCT/uU5Exs+JURKq1fTA8LCOc6D+jWOpACBejIF96wAzbqv8D
# FgMNdGQimpReMDV2E/XT4ePgB8rZ6kWIRpxU1RDi8zIJQLbnXBcy/syv623PYDx1
# 8+5cYEBVG7VZr3IjaE2cdAQMEMmvUFunDWYPluWaleAgohrQsO44SZ4qZ56RlmyY
# 28QQbWB8Hm5I57Z+rzMHEnHvvZU7vqmD1EJ9t6c011+GkbWvVljaVX0Xvdu8zWRB
# FY0xUQZPtC6yiz2c803jWANUzKyI+FI8TktGCSUZ/xXnp5hGLn266uPjfP/5uRmV
# vna5DXmyAlEaSsifiMJDMIIGqjCCBJKgAwIBAgIQSmrxbxietarWh8M9boaJ7zAN
# BgkqhkiG9w0BAQ0FADBPMQswCQYDVQQGEwJVUzEWMBQGA1UEChMNRW50cnVzdCwg
# SW5jLjEoMCYGA1UEAxMfRW50cnVzdCBDb2RlIFNpZ25pbmcgQ0EgLSBPVkNTMjAe
# Fw0yMjEwMjAxNjIxMDVaFw0yMzEwMjAxNjIxMDVaMG4xCzAJBgNVBAYTAlVTMREw
# DwYDVQQIEwhNaXNzb3VyaTEUMBIGA1UEBxMLU2FpbnQgTG91aXMxGjAYBgNVBAoT
# EUNlcnRhcmEgVVNBLCBJbmMuMRowGAYDVQQDExFDZXJ0YXJhIFVTQSwgSW5jLjCC
# AiIwDQYJKoZIhvcNAQEBBQADggIPADCCAgoCggIBANvZ4fp6onI1BfP5ptCU4O+P
# aEoGd3wBj7saF6AUYDzKZmnfAMNJBCnKIanmgLlPHXRSUD82RsJqcyze7+KA2a6k
# ZlALgO0zyKSP2F08fIMC5KhOpL58cdFh4xPKeCSVP2Ht9eolWeMh7fGmw6B0+SPz
# Gbg9cQBhwMx5yClz1zk/lDFWHlOyec+OOkImKWJQ+Sjj39eaV8nuSDB0owLOknd0
# Q25luUojUOE4PR4SHSR+DbF7qVRtrvqTHeICwMzC5/jjc94rTrXokff7+CoMNtZQ
# q+gWlm+ihhXG3ydQaiKHk42Hd2+tyWOFp5wywNzOZOOVH6RZPtNjdb2Skiwq4p13
# 9mRPHiupxDbisiRAsekqI9IFVd9K/Sac5wQA57JtyLvxw1BHRjzrqEmqFpubhPW/
# I9LooEdHzL1VXZQgkTr0f5M7Os8Bhi8uTLZbOssfvnL57IZkrsLCnrdkDGMsKOn5
# 94phD0Hz0FlA5IhzjORSqY+5niMa6VhjX/100RhjNyYxA4BV8S8bw5Ux5bi8ejzS
# 5FQfkdi+0/R99YtOPVIibiGrdQmpQxv8vtImYcfr3a5mDKFmS8/9yOj8XnGLP+OV
# /qNEDdBEDggyHyTqxpy+DjJVFr6sxIH7gVJWbztHU1XvgeDgVz6rX3ysAcyrc1Tg
# iBK3aNiu0fVbSpEUyD9jAgMBAAGjggFhMIIBXTAMBgNVHRMBAf8EAjAAMB0GA1Ud
# DgQWBBS/fct/2kM1i0BZlal6HORrrvmAzjAfBgNVHSMEGDAWgBTvn7p5sHPyJR54
# nANSnBtThN6N7TBnBggrBgEFBQcBAQRbMFkwIwYIKwYBBQUHMAGGF2h0dHA6Ly9v
# Y3NwLmVudHJ1c3QubmV0MDIGCCsGAQUFBzAChiZodHRwOi8vYWlhLmVudHJ1c3Qu
# bmV0L292Y3MyLWNoYWluLnA3YzAxBgNVHR8EKjAoMCagJKAihiBodHRwOi8vY3Js
# LmVudHJ1c3QubmV0L292Y3MyLmNybDAOBgNVHQ8BAf8EBAMCB4AwEwYDVR0lBAww
# CgYIKwYBBQUHAwMwTAYDVR0gBEUwQzA3BgpghkgBhvpsCgEDMCkwJwYIKwYBBQUH
# AgEWG2h0dHBzOi8vd3d3LmVudHJ1c3QubmV0L3JwYTAIBgZngQwBBAEwDQYJKoZI
# hvcNAQENBQADggIBADdZUY2ZsUkPQjSV8B+psXvxa36HpU8SplM4zglL1eOBhaiV
# mS81K0GAEd4jcwWwCh4gU2UXRuARvE5BQDl0IC7q4TzDBjbWMPREy/XPOL0zs3sq
# edfjSW6QyDyvqxvNjqlfqWMt33ao3EvvqjVrjtzG4MVnsPM3LszINaMEz0TAm/Ie
# 9SU3pb88Y1Xa3fS/gPHjrSqy2okWjOf+xTIK+GmRgKdBHdTcA3WL+akstSoPQ6xw
# EXcj7js4chyHRvYlLM4eKt9ypX1JyMfZ2SSxSASyeS+l+08sFY5dfE6xSPFmb/oD
# FKHMjjOzNjMiOTqfdCbAWVJhciA8BddK0Ipz8CmMbiO+OuyZf6g1nGEQXOYWIPBT
# aBixxBnVS+tn7XaZaohf/AN/7Atfv0txjQ+jKq6IhWjre5G+QgIkUli+Ew/DamDn
# 3g8PZV9DaIvoJtivapJgbgUAUD/LcQinM+pRVltB9NTZwumapv1TaOSYJllRmt3s
# OM07JLyLBJxeDftKdZYi4Q+GPkUBEBtyz8P+z4qGqc1wW2+uOP70cQncat8CsZJT
# Kb7oVf3MnyBma80ZCge9UatLf6U4AH1mUNPmKm586D8JBZLwOe3cpFbDLyTlQJsq
# Rvx9m0mBscs9iTo/tFl4KUFtkmzXXXo6ChFMFLsZ7lU8zxS60YQ7woBZfl7aMIIG
# 7DCCBNSgAwIBAgIQMA9vrN1mmHR8qUY2p3gtuTANBgkqhkiG9w0BAQwFADCBiDEL
# MAkGA1UEBhMCVVMxEzARBgNVBAgTCk5ldyBKZXJzZXkxFDASBgNVBAcTC0plcnNl
# eSBDaXR5MR4wHAYDVQQKExVUaGUgVVNFUlRSVVNUIE5ldHdvcmsxLjAsBgNVBAMT
# JVVTRVJUcnVzdCBSU0EgQ2VydGlmaWNhdGlvbiBBdXRob3JpdHkwHhcNMTkwNTAy
# MDAwMDAwWhcNMzgwMTE4MjM1OTU5WjB9MQswCQYDVQQGEwJHQjEbMBkGA1UECBMS
# R3JlYXRlciBNYW5jaGVzdGVyMRAwDgYDVQQHEwdTYWxmb3JkMRgwFgYDVQQKEw9T
# ZWN0aWdvIExpbWl0ZWQxJTAjBgNVBAMTHFNlY3RpZ28gUlNBIFRpbWUgU3RhbXBp
# bmcgQ0EwggIiMA0GCSqGSIb3DQEBAQUAA4ICDwAwggIKAoICAQDIGwGv2Sx+iJl9
# AZg/IJC9nIAhVJO5z6A+U++zWsB21hoEpc5Hg7XrxMxJNMvzRWW5+adkFiYJ+9Uy
# UnkuyWPCE5u2hj8BBZJmbyGr1XEQeYf0RirNxFrJ29ddSU1yVg/cyeNTmDoqHvzO
# WEnTv/M5u7mkI0Ks0BXDf56iXNc48RaycNOjxN+zxXKsLgp3/A2UUrf8H5VzJD0B
# KLwPDU+zkQGObp0ndVXRFzs0IXuXAZSvf4DP0REKV4TJf1bgvUacgr6Unb+0ILBg
# frhN9Q0/29DqhYyKVnHRLZRMyIw80xSinL0m/9NTIMdgaZtYClT0Bef9Maz5yIUX
# x7gpGaQpL0bj3duRX58/Nj4OMGcrRrc1r5a+2kxgzKi7nw0U1BjEMJh0giHPYla1
# IXMSHv2qyghYh3ekFesZVf/QOVQtJu5FGjpvzdeE8NfwKMVPZIMC1Pvi3vG8Aij0
# bdonigbSlofe6GsO8Ft96XZpkyAcSpcsdxkrk5WYnJee647BeFbGRCXfBhKaBi2f
# A179g6JTZ8qx+o2hZMmIklnLqEbAyfKm/31X2xJ2+opBJNQb/HKlFKLUrUMcpEmL
# QTkUAx4p+hulIq6lw02C0I3aa7fb9xhAV3PwcaP7Sn1FNsH3jYL6uckNU4B9+rY5
# WDLvbxhQiddPnTO9GrWdod6VQXqngwIDAQABo4IBWjCCAVYwHwYDVR0jBBgwFoAU
# U3m/WqorSs9UgOHYm8Cd8rIDZsswHQYDVR0OBBYEFBqh+GEZIA/DQXdFKI7RNV8G
# EgRVMA4GA1UdDwEB/wQEAwIBhjASBgNVHRMBAf8ECDAGAQH/AgEAMBMGA1UdJQQM
# MAoGCCsGAQUFBwMIMBEGA1UdIAQKMAgwBgYEVR0gADBQBgNVHR8ESTBHMEWgQ6BB
# hj9odHRwOi8vY3JsLnVzZXJ0cnVzdC5jb20vVVNFUlRydXN0UlNBQ2VydGlmaWNh
# dGlvbkF1dGhvcml0eS5jcmwwdgYIKwYBBQUHAQEEajBoMD8GCCsGAQUFBzAChjNo
# dHRwOi8vY3J0LnVzZXJ0cnVzdC5jb20vVVNFUlRydXN0UlNBQWRkVHJ1c3RDQS5j
# cnQwJQYIKwYBBQUHMAGGGWh0dHA6Ly9vY3NwLnVzZXJ0cnVzdC5jb20wDQYJKoZI
# hvcNAQEMBQADggIBAG1UgaUzXRbhtVOBkXXfA3oyCy0lhBGysNsqfSoF9bw7J/Ra
# oLlJWZApbGHLtVDb4n35nwDvQMOt0+LkVvlYQc/xQuUQff+wdB+PxlwJ+TNe6qAc
# Jlhc87QRD9XVw+K81Vh4v0h24URnbY+wQxAPjeT5OGK/EwHFhaNMxcyyUzCVpNb0
# llYIuM1cfwGWvnJSajtCN3wWeDmTk5SbsdyybUFtZ83Jb5A9f0VywRsj1sJVhGbk
# s8VmBvbz1kteraMrQoohkv6ob1olcGKBc2NeoLvY3NdK0z2vgwY4Eh0khy3k/ALW
# PncEvAQ2ted3y5wujSMYuaPCRx3wXdahc1cFaJqnyTdlHb7qvNhCg0MFpYumCf/R
# oZSmTqo9CfUFbLfSZFrYKiLCS53xOV5M3kg9mzSWmglfjv33sVKRzj+J9hyhtal1
# H3G/W0NdZT1QgW6r8NDT/LKzH7aZlib0PHmLXGTMze4nmuWgwAxyh8FuTVrTHurw
# ROYybxzrF06Uw3hlIDsPQaof6aFBnf6xuKBlKjTg3qj5PObBMLvAoGMs/FwWAKjQ
# xH/qEZ0eBsambTJdtDgJK0kHqv3sMNrxpy/Pt/360KOE2See+wFmd7lWEOEgbsau
# sfm2usg1XTN2jvF8IAwqd661ogKGuinutFoAsYyr4/kKyVRd1LlqdJ69SK6YMIIG
# 9TCCBN2gAwIBAgIQOUwl4XygbSeoZeI72R0i1DANBgkqhkiG9w0BAQwFADB9MQsw
# CQYDVQQGEwJHQjEbMBkGA1UECBMSR3JlYXRlciBNYW5jaGVzdGVyMRAwDgYDVQQH
# EwdTYWxmb3JkMRgwFgYDVQQKEw9TZWN0aWdvIExpbWl0ZWQxJTAjBgNVBAMTHFNl
# Y3RpZ28gUlNBIFRpbWUgU3RhbXBpbmcgQ0EwHhcNMjMwNTAzMDAwMDAwWhcNMzQw
# ODAyMjM1OTU5WjBqMQswCQYDVQQGEwJHQjETMBEGA1UECBMKTWFuY2hlc3RlcjEY
# MBYGA1UEChMPU2VjdGlnbyBMaW1pdGVkMSwwKgYDVQQDDCNTZWN0aWdvIFJTQSBU
# aW1lIFN0YW1waW5nIFNpZ25lciAjNDCCAiIwDQYJKoZIhvcNAQEBBQADggIPADCC
# AgoCggIBAKSTKFJLzyeHdqQpHJk4wOcO1NEc7GjLAWTkis13sHFlgryf/Iu7u5WY
# +yURjlqICWYRFFiyuiJb5vYy8V0twHqiDuDgVmTtoeWBIHIgZEFsx8MI+vN9Xe8h
# msJ+1yzDuhGYHvzTIAhCs1+/f4hYMqsws9iMepZKGRNcrPznq+kcFi6wsDiVSs+F
# UKtnAyWhuzjpD2+pWpqRKBM1uR/zPeEkyGuxmegN77tN5T2MVAOR0Pwtz1UzOHoJ
# HAfRIuBjhqe+/dKDcxIUm5pMCUa9NLzhS1B7cuBb/Rm7HzxqGXtuuy1EKr48TMys
# igSTxleGoHM2K4GX+hubfoiH2FJ5if5udzfXu1Cf+hglTxPyXnypsSBaKaujQod3
# 4PRMAkjdWKVTpqOg7RmWZRUpxe0zMCXmloOBmvZgZpBYB4DNQnWs+7SR0MXdAUBq
# tqgQ7vaNereeda/TpUsYoQyfV7BeJUeRdM11EtGcb+ReDZvsdSbu/tP1ki9Sheja
# RFEqoswAyodmQ6MbAO+itZadYq0nC/IbSsnDlEI3iCCEqIeuw7ojcnv4VO/4ayew
# hfWnQ4XYKzl021p3AtGk+vXNnD3MH65R0Hts2B0tEUJTcXTC5TWqLVIS2SXP8NPQ
# kUMS1zJ9mGzjd0HI/x8kVO9urcY+VXvxXIc6ZPFgSwVP77kv7AkTAgMBAAGjggGC
# MIIBfjAfBgNVHSMEGDAWgBQaofhhGSAPw0F3RSiO0TVfBhIEVTAdBgNVHQ4EFgQU
# Aw8xyJEqk71j89FdTaQ0D9KVARgwDgYDVR0PAQH/BAQDAgbAMAwGA1UdEwEB/wQC
# MAAwFgYDVR0lAQH/BAwwCgYIKwYBBQUHAwgwSgYDVR0gBEMwQTA1BgwrBgEEAbIx
# AQIBAwgwJTAjBggrBgEFBQcCARYXaHR0cHM6Ly9zZWN0aWdvLmNvbS9DUFMwCAYG
# Z4EMAQQCMEQGA1UdHwQ9MDswOaA3oDWGM2h0dHA6Ly9jcmwuc2VjdGlnby5jb20v
# U2VjdGlnb1JTQVRpbWVTdGFtcGluZ0NBLmNybDB0BggrBgEFBQcBAQRoMGYwPwYI
# KwYBBQUHMAKGM2h0dHA6Ly9jcnQuc2VjdGlnby5jb20vU2VjdGlnb1JTQVRpbWVT
# dGFtcGluZ0NBLmNydDAjBggrBgEFBQcwAYYXaHR0cDovL29jc3Auc2VjdGlnby5j
# b20wDQYJKoZIhvcNAQEMBQADggIBAEybZVj64HnP7xXDMm3eM5Hrd1ji673LSjx1
# 3n6UbcMixwSV32VpYRMM9gye9YkgXsGHxwMkysel8Cbf+PgxZQ3g621RV6aMhFII
# Rhwqwt7y2opF87739i7Efu347Wi/elZI6WHlmjl3vL66kWSIdf9dhRY0J9Ipy//t
# Ldr/vpMM7G2iDczD8W69IZEaIwBSrZfUYngqhHmo1z2sIY9wwyR5OpfxDaOjW1PY
# qwC6WPs1gE9fKHFsGV7Cg3KQruDG2PKZ++q0kmV8B3w1RB2tWBhrYvvebMQKqWzT
# IUZw3C+NdUwjwkHQepY7w0vdzZImdHZcN6CaJJ5OX07Tjw/lE09ZRGVLQ2TPSPhn
# Z7lNv8wNsTow0KE9SK16ZeTs3+AB8LMqSjmswaT5qX010DJAoLEZKhghssh9BXEa
# Syc2quCYHIN158d+S4RDzUP7kJd2KhKsQMFwW5kKQPqAbZRhe8huuchnZyRcUI0B
# IN4H9wHU+C4RzZ2D5fjKJRxEPSflsIZHKgsbhHZ9e2hPjbf3E7TtoC3ucw/ZELqd
# mSx813UfjxDElOZ+JOWVSoiMJ9aFZh35rmR2kehI/shVCu0pwx/eOKbAFPsyPfip
# g2I2yMO+AIccq/pKQhyJA9z1XHxw2V14Tu6fXiDmCWp8KwijSPUV/ARP380hHHrl
# 9Y4a1LlAMYIGWTCCBlUCAQEwYzBPMQswCQYDVQQGEwJVUzEWMBQGA1UEChMNRW50
# cnVzdCwgSW5jLjEoMCYGA1UEAxMfRW50cnVzdCBDb2RlIFNpZ25pbmcgQ0EgLSBP
# VkNTMgIQSmrxbxietarWh8M9boaJ7zALBglghkgBZQMEAgGgfDAQBgorBgEEAYI3
# AgEMMQIwADAZBgkqhkiG9w0BCQMxDAYKKwYBBAGCNwIBBDAcBgorBgEEAYI3AgEL
# MQ4wDAYKKwYBBAGCNwIBFTAvBgkqhkiG9w0BCQQxIgQg/NzIBEgvqXSTI8mqOiW8
# qAWLtUFX5xo2VrtvK9Rf7bYwDQYJKoZIhvcNAQEBBQAEggIAo6dN1a2qmPB6MHk1
# N/Oxf5eh+0D73gXfQjU8Bjh6pyNB7m1BTA9Mt/LEBc00CT4MXRFDc/+uZ63Uidf/
# S/AHR220FvUwuR9vIV4wGh7SWKW2qZc3C4JDVAReaOvhmb51flNQbOdEIbLQPpiH
# xV8tSPBV3KrJELI3WclCC24SDdQymymuC1Wgp9NHfUMGQZ0+pN6/bAhUa7XCXu1b
# FdWdHqCeio+7uANqGVlb7fBljOmWTDT/nag4M6CP40tulaMMiO2dlgpuhy37TG1A
# Vc7PI4fplzWoCaLBh3W95Ks0DJ8iElbBW5UbHdfyPb2X6EtplXciGgc5a0nPwoL6
# sTZNUv7aE2uoocjimNmh8VeXNKQiVRYW+/bBBRjMKP4hpHzaYz2Kv0llA6gDqKqk
# SsI//pan2/yEd8L4fep1XQxsxsNjfgRxwlcmMUCggujKe7JmYGLJn1L5Tx6NIj9T
# hqb19eC7Ey1plTYOk0vPzqIZPBPrVSx7dNyXw2biJcIIM+Mr7b0wFkMmn5T6aSui
# RQIFE+SuUtpJ2Q6CIF+2l3lKAb+492ALjRlvSpf0pE2TmSz8R38EeTyVjD/J+yn1
# NDUmPFXVzlF+3KB3HIRVm4h45PteapQ/l5nwquw6EufIcy0Xa8Ga81tRTbx2Er6l
# BfOXGzoQ+/lYGcMW86SpE5koDFChggNLMIIDRwYJKoZIhvcNAQkGMYIDODCCAzQC
# AQEwgZEwfTELMAkGA1UEBhMCR0IxGzAZBgNVBAgTEkdyZWF0ZXIgTWFuY2hlc3Rl
# cjEQMA4GA1UEBxMHU2FsZm9yZDEYMBYGA1UEChMPU2VjdGlnbyBMaW1pdGVkMSUw
# IwYDVQQDExxTZWN0aWdvIFJTQSBUaW1lIFN0YW1waW5nIENBAhA5TCXhfKBtJ6hl
# 4jvZHSLUMA0GCWCGSAFlAwQCAgUAoHkwGAYJKoZIhvcNAQkDMQsGCSqGSIb3DQEH
# ATAcBgkqhkiG9w0BCQUxDxcNMjMxMDE5MjI0ODU5WjA/BgkqhkiG9w0BCQQxMgQw
# DQxkrgsKZ9KP4PTIh6TlAn4bJJmdKWRrTWbIzKlRb9nJl4C7Q+EWeOhj4lbVrLMq
# MA0GCSqGSIb3DQEBAQUABIICAJJevRg7VIjUvLPWMF1WxkmFrZhhBr0TuYJGghSO
# kXzdcFOB0rUFby9+XylT4QPL1qsUdFYYlU/1lqcb+mnK4t//cKqX8Qt+jdC6tf4x
# OeZQBaD4n9fZhOeyH2fX6auTXvYDuE5Pfz23LBSYpFQo12c/54afeffwqQBOHWCq
# nS+hS2Z2J4iNo8q8FmtAtSfPepXE4sQIbxtig73ni3nGpMqvs+Ra1MX+ZGHK4Tdi
# gT70zeuNhEpkcda6d5EQ6/bfTOYjFAX1GC56POtt25QukMgOQ5MRmC335aT2eVi5
# MB6tvQko9uHNbSNYsfvd50AlSWZwOxrlCYaZ6vpOap/5B/ylVg6aXQu14uTQypi4
# Lin23EnjZ8Jm5ghg0rs0xTkAMKXpTkeRkO3w5ijKVTUyXMzvwAb/uPZRi0tDmK4n
# nEBKau4x8oq/uRSOau9UtB0S+9vAV3vkPR7Cg+pIrx83qwP2aTd1OZ2CIdOpNebC
# p1w1P6Ux/pMzBZa7KBHeV+Wn6aA1ARhPB8MjIM5xPItpNVlNXDz8/FpphHqw5aLJ
# z95i88M+gMPOjrvkABpQ8N+dxr3jxCLtvO8Afpd9UlnzSjdsvqvmhZ1Le/1rRKSN
# Z+WmaFOesCMs7msL3qXhYDJYZzBaQBjYVk4YiGOers8evt96ko+1Fejz8EXH6axI
# OQYr
# SIG # End signature block
