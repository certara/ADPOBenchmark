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

[string] $DEBUG_FLAG = if ($env:NLME_BUILD_DEBUG) {"-g"} else {"-O3"}
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
# MIIu5QYJKoZIhvcNAQcCoIIu1jCCLtICAQExDzANBglghkgBZQMEAgEFADB5Bgor
# BgEEAYI3AgEEoGswaTA0BgorBgEEAYI3AgEeMCYCAwEAAAQQH8w7YFlLCE63JNLG
# KX7zUQIBAAIBAAIBAAIBAAIBADAxMA0GCWCGSAFlAwQCAQUABCBgpWWbDLTAowru
# 2lwLbWImjQPbHBWrGxS3wKObE7Tg9KCCEswwggXfMIIEx6ADAgECAhBOQOQ3VO3m
# jAAAAABR05R/MA0GCSqGSIb3DQEBCwUAMIG+MQswCQYDVQQGEwJVUzEWMBQGA1UE
# ChMNRW50cnVzdCwgSW5jLjEoMCYGA1UECxMfU2VlIHd3dy5lbnRydXN0Lm5ldC9s
# ZWdhbC10ZXJtczE5MDcGA1UECxMwKGMpIDIwMDkgRW50cnVzdCwgSW5jLiAtIGZv
# ciBhdXRob3JpemVkIHVzZSBvbmx5MTIwMAYDVQQDEylFbnRydXN0IFJvb3QgQ2Vy
# dGlmaWNhdGlvbiBBdXRob3JpdHkgLSBHMjAeFw0yMTA1MDcxNTQzNDVaFw0zMDEx
# MDcxNjEzNDVaMGkxCzAJBgNVBAYTAlVTMRYwFAYDVQQKDA1FbnRydXN0LCBJbmMu
# MUIwQAYDVQQDDDlFbnRydXN0IENvZGUgU2lnbmluZyBSb290IENlcnRpZmljYXRp
# b24gQXV0aG9yaXR5IC0gQ1NCUjEwggIiMA0GCSqGSIb3DQEBAQUAA4ICDwAwggIK
# AoICAQCngY/3FEW2YkPy2K7TJV5IT1G/xX2fUBw10dZ+YSqUGW0nRqSmGl33VFFq
# gCLGqGZ1TVSDyV5oG6v2W2Swra0gvVTvRmttAudFrnX2joq5Mi6LuHccUk15iF+l
# OhjJUCyXJy2/2gB9Y3/vMuxGh2Pbmp/DWiE2e/mb1cqgbnIs/OHxnnBNCFYVb5Cr
# +0i6udfBgniFZS5/tcnA4hS3NxFBBuKK4Kj25X62eAUBw2DtTwdBLgoTSeOQm3/d
# vfqsv2RR0VybtPVc51z/O5uloBrXfQmywrf/bhy8yH3m6Sv8crMU6UpVEoScRCV1
# HfYq8E+lID1oJethl3wP5bY9867DwRG8G47M4EcwXkIAhnHjWKwGymUfe5SmS1dn
# DH5erXhnW1XjXuvH2OxMbobL89z4n4eqclgSD32m+PhCOTs8LOQyTUmM4OEAwjig
# nPqEPkHcblauxhpb9GdoBQHNG7+uh7ydU/Yu6LZr5JnexU+HWKjSZR7IH9Vybu5Z
# HFc7CXKd18q3kMbNe0WSkUIDTH0/yvKquMIOhvMQn0YupGaGaFpoGHApOBGAYGuK
# Q6NzbOOzazf/5p1nAZKG3y9I0ftQYNVc/iHTAUJj/u9wtBfAj6ju08FLXxLq/f0u
# DodEYOOp9MIYo+P9zgyEIg3zp3jak/PbOM+5LzPG/wc8Xr5F0wIDAQABo4IBKzCC
# AScwDgYDVR0PAQH/BAQDAgGGMBIGA1UdEwEB/wQIMAYBAf8CAQEwHQYDVR0lBBYw
# FAYIKwYBBQUHAwMGCCsGAQUFBwMIMDsGA1UdIAQ0MDIwMAYEVR0gADAoMCYGCCsG
# AQUFBwIBFhpodHRwOi8vd3d3LmVudHJ1c3QubmV0L3JwYTAzBggrBgEFBQcBAQQn
# MCUwIwYIKwYBBQUHMAGGF2h0dHA6Ly9vY3NwLmVudHJ1c3QubmV0MDAGA1UdHwQp
# MCcwJaAjoCGGH2h0dHA6Ly9jcmwuZW50cnVzdC5uZXQvZzJjYS5jcmwwHQYDVR0O
# BBYEFIK61j2Xzp/PceiSN6/9s7VpNVfPMB8GA1UdIwQYMBaAFGpyJnrQHu995ztp
# UdRsjZ+QEmarMA0GCSqGSIb3DQEBCwUAA4IBAQAfXkEEtoNwJFMsVXMdZTrA7LR7
# BJheWTgTCaRZlEJeUL9PbG4lIJCTWEAN9Rm0Yu4kXsIBWBUCHRAJb6jU+5J+Nzg+
# LxR9jx1DNmSzZhNfFMylcfdbIUvGl77clfxwfREc0yHd0CQ5KcX+Chqlz3t57jpv
# 3ty/6RHdFoMI0yyNf02oFHkvBWFSOOtg8xRofcuyiq3AlFzkJg4sit1Gw87kVlHF
# VuOFuE2bRXKLB/GK+0m4X9HyloFdaVIk8Qgj0tYjD+uL136LwZNr+vFie1jpUJuX
# bheIDeHGQ5jXgWG2hZ1H7LGerj8gO0Od2KIc4NR8CMKvdgb4YmZ6tvf6yK81MIIG
# cDCCBFigAwIBAgIQce9VdK81VMNaLGn2b0trzTANBgkqhkiG9w0BAQ0FADBpMQsw
# CQYDVQQGEwJVUzEWMBQGA1UECgwNRW50cnVzdCwgSW5jLjFCMEAGA1UEAww5RW50
# cnVzdCBDb2RlIFNpZ25pbmcgUm9vdCBDZXJ0aWZpY2F0aW9uIEF1dGhvcml0eSAt
# IENTQlIxMB4XDTIxMDUwNzE5MjA0NVoXDTQwMTIyOTIzNTkwMFowTzELMAkGA1UE
# BhMCVVMxFjAUBgNVBAoTDUVudHJ1c3QsIEluYy4xKDAmBgNVBAMTH0VudHJ1c3Qg
# Q29kZSBTaWduaW5nIENBIC0gT1ZDUzIwggIiMA0GCSqGSIb3DQEBAQUAA4ICDwAw
# ggIKAoICAQCemXYXGp5WFwhjLJNNg2GEMzQCttlioN7CDrkgTMhXnQ/dVFsNDNYB
# 3S9I4ZEJ4dvIFQSCtnvw2NYwOxlxcPuoppf2KV2kDKn0Uz5X2wxObvx2218k6apf
# Q+OT5w7PyiW8xEwwC1oP5gb05W4MmWZYT4NhwnN8XCJvAUXFD/dAT2RL0BcKqQ4e
# Ai+hj0zyZ1DbPuSfwk8/dOsxpNCU0Jm8MJIJasskzaLYdlLQTnWYT2Ra0l6D9FjA
# XWp1xNg/ZDqLFA3YduHquWvnEXBJEThjE27xxvq9EEU1B+Z2FdB1FqrCQ1f+q/5j
# c0YioLjz5MdwRgn5qTdBmrNLbB9wcqMH9jWSdBFkbvkC1cCSlfGXWX4N7qIl8nFV
# uJuNv83urt37DOeuMk5QjaHf0XO/wc5/ddqrv9CtgjjF54jtom06hhG317DhqIs7
# DEEXml/kW5jInQCf93PSw+mfBYd5IYPWC+3RzAif4PHFyVi6U1/Uh7GLWajSXs1p
# 0D76xDkJr7S17ec8+iKH1nP5F5Vqwxz1VXhf1PoLwFs/jHgVDlpMOm7lJpjQJ8wg
# 38CGO3qNZUZ+2WFeqfSuPtT8r0XHOrOFBEqLyAlds3sCKFnjhn2AolhAZmLgOFWD
# q58pQSa6u+nYZPi2uyhzzRVK155z42ZMsVGdgSOLyIZ3srYsNyJwIQIDAQABo4IB
# LDCCASgwEgYDVR0TAQH/BAgwBgEB/wIBADAdBgNVHQ4EFgQU75+6ebBz8iUeeJwD
# UpwbU4Teje0wHwYDVR0jBBgwFoAUgrrWPZfOn89x6JI3r/2ztWk1V88wMwYIKwYB
# BQUHAQEEJzAlMCMGCCsGAQUFBzABhhdodHRwOi8vb2NzcC5lbnRydXN0Lm5ldDAx
# BgNVHR8EKjAoMCagJKAihiBodHRwOi8vY3JsLmVudHJ1c3QubmV0L2NzYnIxLmNy
# bDAOBgNVHQ8BAf8EBAMCAYYwEwYDVR0lBAwwCgYIKwYBBQUHAwMwRQYDVR0gBD4w
# PDAwBgRVHSAAMCgwJgYIKwYBBQUHAgEWGmh0dHA6Ly93d3cuZW50cnVzdC5uZXQv
# cnBhMAgGBmeBDAEEATANBgkqhkiG9w0BAQ0FAAOCAgEAXvOGmTXBee7wEK/XkkPS
# hdBb4Jig4HFRyRTLUJpgDrAEJkmxz+m6mwih2kNd1G8jorn4QMdH/k0BC0iQP8jc
# arQ+UzUovkBKR4VqHndAzIB/YbQ8T3mo5qOmoH5EhnG/EhuVgXL3DaXQ3mefxqK4
# 8Wr5/P50ZsZk5nk9agNhTksfzCBiywIY7GPtfnE/lroLXmgiZ+wfwNIFFmaxsqTq
# /MWVo40SpfWN7xsgzZn35zLzWXEf3ZTmeeVSIxBWKvxZOL+/eSWSasf9q2d3cbEE
# fTWtFME+qPwjF1YIGHzXeiJrkWrMNUVtTzudQ50FuJ3z/DQhXAQYMlc4NMHKgyNG
# pogjIcZ+FICrse+7C6wJP+5TkTGz4lREqrV9MDwsI5zoP6NY6kAIF6MgX3rADNuq
# /wMWAw10ZCKalF4wNXYT9dPh4+AHytnqRYhGnFTVEOLzMglAtudcFzL+zK/rbc9g
# PHXz7lxgQFUbtVmvciNoTZx0BAwQya9QW6cNZg+W5ZqV4CCiGtCw7jhJnipnnpGW
# bJjbxBBtYHwebkjntn6vMwcSce+9lTu+qYPUQn23pzTXX4aRta9WWNpVfRe927zN
# ZEEVjTFRBk+0LrKLPZzzTeNYA1TMrIj4UjxOS0YJJRn/FeenmEYufbrq4+N8//m5
# GZW+drkNebICURpKyJ+IwkMwggZxMIIEWaADAgECAhALBv9gFrEjvi/8h6aykuQa
# MA0GCSqGSIb3DQEBDQUAME8xCzAJBgNVBAYTAlVTMRYwFAYDVQQKEw1FbnRydXN0
# LCBJbmMuMSgwJgYDVQQDEx9FbnRydXN0IENvZGUgU2lnbmluZyBDQSAtIE9WQ1My
# MB4XDTI0MDExMjE3MDIxM1oXDTI2MDExMjE3MDIxMlowbjELMAkGA1UEBhMCVVMx
# ETAPBgNVBAgTCE1pc3NvdXJpMRQwEgYDVQQHEwtTYWludCBMb3VpczEaMBgGA1UE
# ChMRQ2VydGFyYSBVU0EsIEluYy4xGjAYBgNVBAMTEUNlcnRhcmEgVVNBLCBJbmMu
# MIICIjANBgkqhkiG9w0BAQEFAAOCAg8AMIICCgKCAgEA4XrJJ+d7N1KI8BD4Abzp
# 8Q7XVvPMnd3BsdsYO4zu5FSqICUyMI+2VNpH65QTSapJ53/pbCDY3kAe9t4xrey7
# hmrubtmq0sbb14BBajHqyZYd95Vw+q3ZetofE5kjwGBQ8NbN63fefiFXgQLovFzI
# G7WXWROKZxg612oQ07LTxnlmH6ihgLRD+4fa+Szuea+lB0tUmOySgX6JnLfhxddI
# 6RQtEL5341bsi3ZdlGmGtYRuJ//hO4rl/2NWk3Q1dNIThlNWL7OmoLK6pdzL6K5C
# fBevozUStKWlJNYDODfYpVHcB1Ur0PyhE+xATP7hvEaYlCwP5Lx78dYkq7t8caMh
# ZTcETjIxhlBNVEnK/WKB6fm2YuDDIr7KxGsxJS238gFZFmTgsNWgIolEJbtp6Eij
# iqg+sMdH6bHlRs1T3KKPufvGK4d4IOZs5gRZ8V2C3HNz4MqReVWJBqV8fwDgMaGW
# d7Co39Xge/BpsohxNISvdKbtXSsk4rl9iYnobcJkQoElD+Vg8EE4J3/YADdhzdRU
# D/lrt4J8hkYIvxAoKxfrDVRbg3OVbfaQgVFgSXL6Bes9iHOz3BxF6vpl8FCFmewy
# sir2mbIjonGGFNLMIwD7ifgkzqb9fRVyWrOrzQLZMHURCuQNaBNBWMtvIRjP5If4
# tvCe0UeRrt1JW3in8y41sAsCAwEAAaOCASgwggEkMAwGA1UdEwEB/wQCMAAwHQYD
# VR0OBBYEFBACsSii96bGwQe7pimxIEXYV+KpMB8GA1UdIwQYMBaAFO+funmwc/Il
# HnicA1KcG1OE3o3tMGcGCCsGAQUFBwEBBFswWTAjBggrBgEFBQcwAYYXaHR0cDov
# L29jc3AuZW50cnVzdC5uZXQwMgYIKwYBBQUHMAKGJmh0dHA6Ly9haWEuZW50cnVz
# dC5uZXQvb3ZjczItY2hhaW4ucDdjMDEGA1UdHwQqMCgwJqAkoCKGIGh0dHA6Ly9j
# cmwuZW50cnVzdC5uZXQvb3ZjczIuY3JsMA4GA1UdDwEB/wQEAwIHgDATBgNVHSUE
# DDAKBggrBgEFBQcDAzATBgNVHSAEDDAKMAgGBmeBDAEEATANBgkqhkiG9w0BAQ0F
# AAOCAgEAZveTUaMJxht6HqF7iR/Oi7ZG4x3HX4ieN+vgp9TQmbQzAf+Sygub6P3L
# F7xNijdEZCTdaU2TSDT1GaLPxV7iJSGO9yNZn/FW3caxDnIWGs/dUjUhNOdO0SZB
# pjaKx2CEesRzkDNBxidGlpuxEPipLbgpKH9yMvg9Ww2bPVosm+Q2frj5t/Hxa4AM
# HRXCLdHsLJ/O6vBw0Gm0gPh/0xty4Y8dWN0zL9du5yHxjv8tZM5DIacX5rGiRBnP
# pcZ4JCCTUheYB38lT3A4R68mDDAZCLQTLU1ICUSZiik/eVjLvKIzG5GBmAzInYMx
# PvhrJ/6AIWJsejSV0KjZqRny8ugY2EXqFXvDVi/+duXpmBiGsgkyjumqX8sl2Cqh
# WYuBSg+1Fvf6vfFBZFckHCx8iLsjS0ftQDlXUU2MVgbMniu8pxnD5bbeQaOZc5Oi
# XxgfHtQe4KlVmn3DM0+DYpAaFwMeIlDp7fSmTkjdNC5UpsP1pEUqesD/gc/Ub78A
# NuXjnyvZsvQAs7IpvLU2gd/hVOJmw/lcBD29mMBgfDm/hOz4eQJEzdgvJu1kRpEe
# hALhn8fdSxf1nVaLWpf7zx7pXLdCC8BSdZw0Jf+jXYJ0/U4kr3i/FEONEbxLd9Cf
# /zvyZplF/N7wkGtVWACObiQmgiPr6uM7yoK7n3JZ8+HIdv3FAPgxghtvMIIbawIB
# ATBjME8xCzAJBgNVBAYTAlVTMRYwFAYDVQQKEw1FbnRydXN0LCBJbmMuMSgwJgYD
# VQQDEx9FbnRydXN0IENvZGUgU2lnbmluZyBDQSAtIE9WQ1MyAhALBv9gFrEjvi/8
# h6aykuQaMA0GCWCGSAFlAwQCAQUAoHwwEAYKKwYBBAGCNwIBDDECMAAwGQYJKoZI
# hvcNAQkDMQwGCisGAQQBgjcCAQQwHAYKKwYBBAGCNwIBCzEOMAwGCisGAQQBgjcC
# ARUwLwYJKoZIhvcNAQkEMSIEIOXGTdxi+yu8Ga9vX5pv8Y9K4xPFttp5TJSmov4I
# 1l/mMA0GCSqGSIb3DQEBAQUABIICAM15jD0tyISJOE1wvMNF3CemtLxOzSzs7vOe
# zOCIHn0kCadZKpggYjYIEk3WrMfy8mBqMLQmtVM/rHBqpFKGcnqkEWfMlAlkX8r1
# 2r/kqQzvXyeXg3RxXfE8IouBaJBTg1JOpwSA6vViBQCWc/xcb1Wkta/JR0+n7Ub2
# NYqRJugSwJQEVw6AdqNuErJItvfIy7eKs7A9gYxrECSKL+FXWhpAakkEkRkKBsj2
# b9CFH5D9wOIP2FQ5obwzUWloaAM+RTgAlPOYRj3XxPTcA2mYGQ6slUpEs41a+UGZ
# vAbambCXrNfDOJFKmUHDuv+CsIbQrF6ZrVO4i5JJxOh6bMweEC2bO7EWoHvsJywO
# C0KRlV1aPcdhsv3gFqnutp6a5TMOULEWTDjxoyyQBtl+Z3Il4jcat65rMoXUewXo
# mSuSZ5SSnyopYPYO6iOpmb4yfUezIBmbNWAp3aU8vk9rQkNh312u/A51gvxhTu3D
# x/ArH8tj41eVnHtftOxmDQref/fFf158HIuBJ8GaZLwFYZ3u1fDbCf3NJw80KaKM
# W60R0RM+/3RIMoJXbJ/dXs4tBGRkBuNVV0fVeNT5Z9uaaSDUJUkY4DBDVaQuH47o
# ERbw/baPSu4NJOr58+xTeJVsmm4HLng5FSOja6OnAmdK61VXZrKfidlJiJDc0HdC
# 4WKGhnFWoYIYXzCCGFsGCisGAQQBgjcDAwExghhLMIIYRwYJKoZIhvcNAQcCoIIY
# ODCCGDQCAQMxDTALBglghkgBZQMEAgMwgfMGCyqGSIb3DQEJEAEEoIHjBIHgMIHd
# AgEBBgpghkgBhvpsCgMFMDEwDQYJYIZIAWUDBAIBBQAEIJ9iueABJCuHVD7Pk6We
# TE4uQSC+iG7nSgJsikmU764oAgh60OB5OfGyvhgPMjAyNDA0MTcyMTE1MTFaMAMC
# AQGgeaR3MHUxCzAJBgNVBAYTAkNBMRAwDgYDVQQIEwdPbnRhcmlvMQ8wDQYDVQQH
# EwZPdHRhd2ExFjAUBgNVBAoTDUVudHJ1c3QsIEluYy4xKzApBgNVBAMTIkVudHJ1
# c3QgVGltZXN0YW1wIEF1dGhvcml0eSAtIFRTQTKgghMOMIIF3zCCBMegAwIBAgIQ
# TkDkN1Tt5owAAAAAUdOUfzANBgkqhkiG9w0BAQsFADCBvjELMAkGA1UEBhMCVVMx
# FjAUBgNVBAoTDUVudHJ1c3QsIEluYy4xKDAmBgNVBAsTH1NlZSB3d3cuZW50cnVz
# dC5uZXQvbGVnYWwtdGVybXMxOTA3BgNVBAsTMChjKSAyMDA5IEVudHJ1c3QsIElu
# Yy4gLSBmb3IgYXV0aG9yaXplZCB1c2Ugb25seTEyMDAGA1UEAxMpRW50cnVzdCBS
# b290IENlcnRpZmljYXRpb24gQXV0aG9yaXR5IC0gRzIwHhcNMjEwNTA3MTU0MzQ1
# WhcNMzAxMTA3MTYxMzQ1WjBpMQswCQYDVQQGEwJVUzEWMBQGA1UECgwNRW50cnVz
# dCwgSW5jLjFCMEAGA1UEAww5RW50cnVzdCBDb2RlIFNpZ25pbmcgUm9vdCBDZXJ0
# aWZpY2F0aW9uIEF1dGhvcml0eSAtIENTQlIxMIICIjANBgkqhkiG9w0BAQEFAAOC
# Ag8AMIICCgKCAgEAp4GP9xRFtmJD8tiu0yVeSE9Rv8V9n1AcNdHWfmEqlBltJ0ak
# phpd91RRaoAixqhmdU1Ug8leaBur9ltksK2tIL1U70ZrbQLnRa519o6KuTIui7h3
# HFJNeYhfpToYyVAslyctv9oAfWN/7zLsRodj25qfw1ohNnv5m9XKoG5yLPzh8Z5w
# TQhWFW+Qq/tIurnXwYJ4hWUuf7XJwOIUtzcRQQbiiuCo9uV+tngFAcNg7U8HQS4K
# E0njkJt/3b36rL9kUdFcm7T1XOdc/zubpaAa130JssK3/24cvMh95ukr/HKzFOlK
# VRKEnEQldR32KvBPpSA9aCXrYZd8D+W2PfOuw8ERvBuOzOBHMF5CAIZx41isBspl
# H3uUpktXZwx+Xq14Z1tV417rx9jsTG6Gy/Pc+J+HqnJYEg99pvj4Qjk7PCzkMk1J
# jODhAMI4oJz6hD5B3G5WrsYaW/RnaAUBzRu/roe8nVP2Lui2a+SZ3sVPh1io0mUe
# yB/Vcm7uWRxXOwlyndfKt5DGzXtFkpFCA0x9P8ryqrjCDobzEJ9GLqRmhmhaaBhw
# KTgRgGBrikOjc2zjs2s3/+adZwGSht8vSNH7UGDVXP4h0wFCY/7vcLQXwI+o7tPB
# S18S6v39Lg6HRGDjqfTCGKPj/c4MhCIN86d42pPz2zjPuS8zxv8HPF6+RdMCAwEA
# AaOCASswggEnMA4GA1UdDwEB/wQEAwIBhjASBgNVHRMBAf8ECDAGAQH/AgEBMB0G
# A1UdJQQWMBQGCCsGAQUFBwMDBggrBgEFBQcDCDA7BgNVHSAENDAyMDAGBFUdIAAw
# KDAmBggrBgEFBQcCARYaaHR0cDovL3d3dy5lbnRydXN0Lm5ldC9ycGEwMwYIKwYB
# BQUHAQEEJzAlMCMGCCsGAQUFBzABhhdodHRwOi8vb2NzcC5lbnRydXN0Lm5ldDAw
# BgNVHR8EKTAnMCWgI6Ahhh9odHRwOi8vY3JsLmVudHJ1c3QubmV0L2cyY2EuY3Js
# MB0GA1UdDgQWBBSCutY9l86fz3Hokjev/bO1aTVXzzAfBgNVHSMEGDAWgBRqciZ6
# 0B7vfec7aVHUbI2fkBJmqzANBgkqhkiG9w0BAQsFAAOCAQEAH15BBLaDcCRTLFVz
# HWU6wOy0ewSYXlk4EwmkWZRCXlC/T2xuJSCQk1hADfUZtGLuJF7CAVgVAh0QCW+o
# 1PuSfjc4Pi8UfY8dQzZks2YTXxTMpXH3WyFLxpe+3JX8cH0RHNMh3dAkOSnF/goa
# pc97ee46b97cv+kR3RaDCNMsjX9NqBR5LwVhUjjrYPMUaH3LsoqtwJRc5CYOLIrd
# RsPO5FZRxVbjhbhNm0VyiwfxivtJuF/R8paBXWlSJPEII9LWIw/ri9d+i8GTa/rx
# YntY6VCbl24XiA3hxkOY14FhtoWdR+yxnq4/IDtDndiiHODUfAjCr3YG+GJmerb3
# +sivNTCCBm8wggRXoAMCAQICECW8K/MpyhB/Hqm6iIXUnTswDQYJKoZIhvcNAQEN
# BQAwaTELMAkGA1UEBhMCVVMxFjAUBgNVBAoMDUVudHJ1c3QsIEluYy4xQjBABgNV
# BAMMOUVudHJ1c3QgQ29kZSBTaWduaW5nIFJvb3QgQ2VydGlmaWNhdGlvbiBBdXRo
# b3JpdHkgLSBDU0JSMTAeFw0yMTA1MDcxOTIyMTRaFw00MDEyMjkyMzU5MDBaME4x
# CzAJBgNVBAYTAlVTMRYwFAYDVQQKEw1FbnRydXN0LCBJbmMuMScwJQYDVQQDEx5F
# bnRydXN0IFRpbWUgU3RhbXBpbmcgQ0EgLSBUUzIwggIiMA0GCSqGSIb3DQEBAQUA
# A4ICDwAwggIKAoICAQC1AyoGtoRPNMyeMb7qjsZ7biAkDwPXvYE2M+Zv0j67xJ6q
# oMxmXUJgNFHiLWGDujyeaLhLw2aOpd4rupstQaXe0MtXBS2I2cBGiG08NQ0ZkKy4
# DBnwTMXbRVvcO8K8jUQA4Dj//13IzwiaPdSy63uVw8SlAOBiAWRZX4zje4up+UW3
# xrCiCjdDuEaBq4Z+fy/e8F/rzSDMpS0x46gumZvgeN30212CY30wOYh+JAbmfGCE
# eMhcKeWVy/V7T89Y3JDPp6J7FFTE4DeYMMGbtq6cKfZrJUPnEmo+GYu+wOeB10ow
# CH58jd8880iTId6Bg2qdAD7XYLrRs2IIlum2SQA49Fx2Ddp3aj2gld4eocxZel6f
# z+l2XUDytRW1YGgs81rJI4PY9RpraSikttSuYgbeJkW93ulWd6rcZLBBzcwT8V1x
# dLKUCEtPMm5+cLh36dUyN8J63kIS6HEc4thiv6prQYYGW+ZpviYJ9JfC/kz0gHKE
# btvexQepjhWibeEb4AkP9aAHoLvEd3MJPAeTjQG1EmctTRm1uMXJEKtwz0L/pScd
# 1hLW5BhEYPs5XYS7ZrVTEp0DFIJlKbTsSXL9s0PlwwIpJLof+Li+XaO3Lqn8z2LZ
# +pfEE3jjVblaeoTr/7vPaYjAtvmLYIVBEFDHBRDSXnadPjXs9k+K+RJ7P68LNwID
# AQABo4IBLDCCASgwEgYDVR0TAQH/BAgwBgEB/wIBADAdBgNVHQ4EFgQUJg/wxEgI
# G83dkfVUVLazs/yZ8QgwHwYDVR0jBBgwFoAUgrrWPZfOn89x6JI3r/2ztWk1V88w
# MwYIKwYBBQUHAQEEJzAlMCMGCCsGAQUFBzABhhdodHRwOi8vb2NzcC5lbnRydXN0
# Lm5ldDAxBgNVHR8EKjAoMCagJKAihiBodHRwOi8vY3JsLmVudHJ1c3QubmV0L2Nz
# YnIxLmNybDAOBgNVHQ8BAf8EBAMCAYYwEwYDVR0lBAwwCgYIKwYBBQUHAwgwRQYD
# VR0gBD4wPDAwBgRVHSAAMCgwJgYIKwYBBQUHAgEWGmh0dHA6Ly93d3cuZW50cnVz
# dC5uZXQvcnBhMAgGBmeBDAEEAjANBgkqhkiG9w0BAQ0FAAOCAgEAdj1GaIVfCcDO
# yfjHuNd+p1w7C0ZzziJTizj2Ebp3xMKHIY8n2QyV6+hL5VzXkBVvqCosimrgIhE0
# efq9lnnIdhbNsUTqcVEPm1XJGHzVgnmc86a3k6kFOHICBpehqLJ5fl4I4m5seZqo
# h5TOf49VNkAPnz9R1Wa+e6uG5m6Huk5jXbHYjh/LZ8MNcNp665OyFITSPn2TPxYM
# NqBceQCfC27lhCrYiMFtBLc385KacOA7A/3NuyeCzi/8jeSyyr74JYXG7XTIPTVf
# OAk9eU/rG+BBXqV0gT9RFcD4SYiPursF1K1FgjN5wSWNX1Q9keS4nxeYAF2tKOVP
# Xxv7+FS1pcQk/PB2O/gNXsxHsMqqu25R31O1SRrxYIe3+f1pBnVfc9YRkPKAWI7l
# ww8DmIwEU7Mph98/97DpTFeBJER5aP4bNgfWZT3sb9bCtaphfGYG7NLlaYD4cZIu
# XOIRRhhFS9b6BWTvu94GykMlvd+NyQF0YYjb8MemPeMMcbx/S+fI4G7g2oD5AJ7A
# ayXVo7pcK/7EYCAUSgcjMeUay5FEspp7Q/FbmLUhS7gxOyJU7nlh95qUG2YnKsbf
# 4WVd73E55lAl/Yc0ua5dfCc752WT+CiEsW+GkyyTk7Zwr6HuyKRhqYQ7+wq3+Lht
# Ju5HTvVeBfqcDxF918uRrkMg9xVZY7wwgga0MIIEnKADAgECAhBbcCbMlvZ4GruF
# 9hH1bbtuMA0GCSqGSIb3DQEBDQUAME4xCzAJBgNVBAYTAlVTMRYwFAYDVQQKEw1F
# bnRydXN0LCBJbmMuMScwJQYDVQQDEx5FbnRydXN0IFRpbWUgU3RhbXBpbmcgQ0Eg
# LSBUUzIwHhcNMjQwMTE5MTY0NzQ3WhcNMzUwNDE4MDAwMDAwWjB1MQswCQYDVQQG
# EwJDQTEQMA4GA1UECBMHT250YXJpbzEPMA0GA1UEBxMGT3R0YXdhMRYwFAYDVQQK
# Ew1FbnRydXN0LCBJbmMuMSswKQYDVQQDEyJFbnRydXN0IFRpbWVzdGFtcCBBdXRo
# b3JpdHkgLSBUU0EyMIICIjANBgkqhkiG9w0BAQEFAAOCAg8AMIICCgKCAgEAqoYE
# OF6PaL+D9Vr9VJvFfTp1ncSnLU9t6dAFH1HjM7svXzqxllSK6Qh8NK2Jg1WknwLM
# IwvYG3pApMyfQuoTf3y44LdKAgXig0kEbwaGzXNBqYPUmGf69FIZeuNKWSiHVhdd
# SPGGkQu4ImTbQfldVLU1pG443AgNGlYYQMN+mDxCM4QNxaVhUc4gbU8Ay0LwqHUb
# 20b+Kdwbntf4GAVRdjCbdL2VHxlTZRVHLFZja+m6SKwKOLbBcv0gCqN0GmsHf9Hd
# rBfOtRzHeokM7G0cMI0F8K89l8w1tLUFA2a6nnb8OdrImtYSEuRBwoUiQPDLuojp
# 0ofCq8Y0O+WrDQAGDga1i3vRCyLaPKjJVnvwNQSW6llGjI/UoLWpg7DOhPtLROVB
# qBbzr9rRoCdw3wfvN/Oukc7UIX+GmNxe7o/A2kfbacoQuZGVgBVj8SsawpahH8L3
# PNT2fSQHJahUlG8KVdvbJENuLjuie0m7tdYYj9kEs77qx7VkmkvOUmEeKwUeYzdG
# nbHJ1V6HpOrXNLIhQhe4Oig6XqXtPv03F39jIPJ71l/K8xQ/4c7/ineUZm2JweDs
# fwRwOGQn9acXfU3KDIEbxeXxNsV6rn0ppEc1OPoN9FMDKQX8b6GLyc3xuIhA09Lb
# niUxrdfmWtgEtIS7BEZhZv9dMt780z58Thjvft8CAwEAAaOCAWUwggFhMAwGA1Ud
# EwEB/wQCMAAwHQYDVR0OBBYEFPV2GvgQmJKhG3epACzxlWICC3knMB8GA1UdIwQY
# MBaAFCYP8MRICBvN3ZH1VFS2s7P8mfEIMGgGCCsGAQUFBwEBBFwwWjAjBggrBgEF
# BQcwAYYXaHR0cDovL29jc3AuZW50cnVzdC5uZXQwMwYIKwYBBQUHMAKGJ2h0dHA6
# Ly9haWEuZW50cnVzdC5uZXQvdHMyLWNoYWluMjU2LnA3YzAxBgNVHR8EKjAoMCag
# JKAihiBodHRwOi8vY3JsLmVudHJ1c3QubmV0L3RzMmNhLmNybDAOBgNVHQ8BAf8E
# BAMCB4AwFgYDVR0lAQH/BAwwCgYIKwYBBQUHAwgwTAYDVR0gBEUwQzA3BgpghkgB
# hvpsCgEHMCkwJwYIKwYBBQUHAgEWG2h0dHBzOi8vd3d3LmVudHJ1c3QubmV0L3Jw
# YTAIBgZngQwBBAIwDQYJKoZIhvcNAQENBQADggIBAKmrfb8aAIVb3O1xJl6Ugq9c
# gkv6HDnFU7XDBt0DYH75YZpBIMKuQRcupUIIkQlelzCYgUXWsrWEPYvphwfaAT/g
# CFhnESCUHsAWjmN3vZtsBY09tcuaMalKXGIyohPOkJwNx5BPZ8PgqH+HhEvX8sEh
# DxDnF7n7vQnMvoqmAf5Ngk9pIJp1a+QN91AmU/wz9/4brqdqwKjrHq8i0z1gFZ+6
# 5NUppLVXn7Fl9rFMYdXSyNq3rKoYHyAYiqb49Qf5civ2Y9glnBb++5TfhnSiILTy
# CN8W7zmAdjqSsdCWg2rafFOJWRsNXPG7KfIhT2EsJIn4dgl/2WiQjlcMZNV2AHFZ
# 89SEyDyhiH+ob/O9bn+wqI7mk2zpFMV1HAwrzvIH+7Wu1EExv8HMaZgYrlsIj6tc
# ZLmEar1cOKHfT0K3S1tS0973O8ufb8JZQiJOCxi3Isgv/GoJhe1QKVF6xJRLtnFl
# ikqGmkt4S1aKod4vi5NbMsyhue+ptgzYBgsXML8Nb4+TrMsR9fHHAJ7QGdecX45U
# fGupQztj3MFEq72MOkPwcj8klc2EkV0hAA14aw1cIySfTK80yxRa3rHkRVD9r2+n
# BYKnc8/P6ZLqcyqx4d2iA+YgvB1nGlbCLvasX8pOgbDmWh1zz9IU81B4KAVOFW6F
# JPgzqIivdG30Us6MqISeMYIEFjCCBBICAQEwYjBOMQswCQYDVQQGEwJVUzEWMBQG
# A1UEChMNRW50cnVzdCwgSW5jLjEnMCUGA1UEAxMeRW50cnVzdCBUaW1lIFN0YW1w
# aW5nIENBIC0gVFMyAhBbcCbMlvZ4GruF9hH1bbtuMAsGCWCGSAFlAwQCA6CCAYkw
# GgYJKoZIhvcNAQkDMQ0GCyqGSIb3DQEJEAEEMBwGCSqGSIb3DQEJBTEPFw0yNDA0
# MTcyMTE1MTFaMCkGCSqGSIb3DQEJNDEcMBowCwYJYIZIAWUDBAIDoQsGCSqGSIb3
# DQEBDTBPBgkqhkiG9w0BCQQxQgRAQdPIY7ZVlMXiVb6v6MXRq1YEEA39l756Ga8t
# uaLojwWZmfiVC37Mgp21uqQ89VIzKS24XtE89RC25SwoIBWhTzCB0AYLKoZIhvcN
# AQkQAi8xgcAwgb0wgbowgbcwCwYJYIZIAWUDBAIDBEA5EUIuFwI+qpkkmXQODsjo
# 0nLTVfxc9mz5EVavl1U05ICv07x8TFtX79H/vNt1FGXg1AVahU6bETnZ9+xV1f4k
# MGYwUqRQME4xCzAJBgNVBAYTAlVTMRYwFAYDVQQKEw1FbnRydXN0LCBJbmMuMScw
# JQYDVQQDEx5FbnRydXN0IFRpbWUgU3RhbXBpbmcgQ0EgLSBUUzICEFtwJsyW9nga
# u4X2EfVtu24wCwYJKoZIhvcNAQENBIICAElef3xqX5k9GwzSXty7pLl5lBIkhZFh
# 6U+GHa+hKY4HA0EOmwacQESuJmxbm4vPq5/5rWWXNt94aYz1cYxe5Azlr3zz6bfA
# MTNxZQ2MOM4bCX50G2Eyd4LNZ+FEDUN3xHJHo22T5mhIDmflnVN4Lg4sJ/Yi2pDq
# u2aNwV9sEkCqJmFZxrVaQDsKEvYgDTOfmNWw+IaLT16dHyJ4Sc9rbAi5t2GWssBS
# JHmfMl72zu+u53zDXEItjVP3pdJJ+E4Bo8jKsAUBoJHgqH/30k+0AvXBUdEgfVCB
# ViISyfM3Gq3ryd0EU3M9WdCQJ/jVG95nPoy37N+AB5/Ijx+0NCPGnPHebos2RAyb
# Sl0EW6Al2HcLq305NB37EYfnw15WzcQmtyUFbD8xn4He6LE8a+Hq1yHON9c0SLSZ
# AGjdeGHsaS4SKO5E7reOS5ICuc02MzgnTdphnn71IorqlziXpsHuUiphmmR/M3zY
# 98F4ZPAz+NCkdYbCXyMwW/ob+h5YEbPVrKUp+3unOANqp1Et2EOmwqZNg3xv5D5F
# hxMUe93CE/4co3/bwWIZXzPnU8YYAtIv+EeRbcBoSlxv67VBCSxjXjbsswMBKSk1
# ANqUiBnUfb7QmtFSlcxseqCdO+SqdddAEoxxNyYsfZLXiHMgX1RpWzH+1PXEQfNB
# kCuU5HeoBKUj
# SIG # End signature block
