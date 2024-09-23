#

param(
    [parameter(mandatory=$true)]
    [string] 
    $logFile,

    [parameter(mandatory=$true)]
    [string] 
    $statusFile,

    [parameter(mandatory=$true)]
    [string] 
    $statusBackupFile,

    [parameter(mandatory=$true)]
    [string] 
    $outFile,

    [parameter(mandatory=$true)]
    [Int32]
    $bootSeed,

    [parameter(mandatory=$true)]
    [Int32]
    $maxTries,

    [parameter(mandatory=$true)]
    [string] 
    $baseDir,

    [parameter(mandatory=$true)]
    [string] 
    $commandString
)

. "$baseDir\common.ps1"

$sysTime = Get-Date -Format 'yyyy/MM/dd HH:mm:ss K'

"Sample command line: $commandString"

$statusRetries = 1
$numTries = 1

"RUNNING" | Out-File -Encoding ascii -FilePath $logFile
"RUNNING" | Out-File -Encoding ascii -FilePath $statusFile

"' $sysTime '" | Out-File -Encoding ascii -FilePath $statusFile -Append

$startTime = Get-Date -Format 'yyyy/MM/dd HH:mm:ss'

function report_status([string] $status)
{
    $moved = $false

    while ($statusRetries -le 20)
    {
        $statusRetries += 1

        $status | Out-File -Encoding ascii -FilePath $logFile

        if (-not $moved)
        {
            $statusFile = (Get-Item -LiteralPath $statusFile).FullName

            move -force $statusFile $statusBackupFile

            if (-not $?)
            {
                "File is locked on move"
                timeout /t 1
                continue
            }

            $moved = $true
        }

        $status | Out-File -Encoding ascii -FilePath $statusFile

        if ($?)
        {
            return
        }

        "File is locked on write"
        timeout /t 1
    }
}

function report_time([string] $time_str)
{
    while ($statusRetries -le 20)
    {
        $statusRetries += 1

        $time_str | Out-File -Encoding ascii -FilePath $statusFile -Append

        if ($?)
        {
            return
        }

        "File is locked on write"
        timeout /t 1
    }
}

while ($numTries -le $maxTries)
{
    $cs = $commandString -replace '\$bootSeed', $bootSeed -replace '\$numTries', $numTries

    run_app "powershell.exe" "-noninteractive -executionpolicy remotesigned -file ""$baseDir\execNLMECmd.ps1"" $cs" -silent $true

    $stopTime = Get-Date -Format 'yyyy/MM/dd HH:mm:ss'

    if (Test-Path $outFile -PathType leaf)
    {
        report_status "SUCCESS"

        report_time $startTime
        report_time $stopTime

        exit 0
    }

    $numTries += 1
    $bootSeed += 1
}

report_status "FAILED"

exit 1

# SIG # Begin signature block
# MIIu5gYJKoZIhvcNAQcCoIIu1zCCLtMCAQExDzANBglghkgBZQMEAgEFADB5Bgor
# BgEEAYI3AgEEoGswaTA0BgorBgEEAYI3AgEeMCYCAwEAAAQQH8w7YFlLCE63JNLG
# KX7zUQIBAAIBAAIBAAIBAAIBADAxMA0GCWCGSAFlAwQCAQUABCAmwwqwyTsP5BFp
# RyOMPwejaa5hqrjewrpXNCEPMJObeaCCEswwggXfMIIEx6ADAgECAhBOQOQ3VO3m
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
# /zvyZplF/N7wkGtVWACObiQmgiPr6uM7yoK7n3JZ8+HIdv3FAPgxghtwMIIbbAIB
# ATBjME8xCzAJBgNVBAYTAlVTMRYwFAYDVQQKEw1FbnRydXN0LCBJbmMuMSgwJgYD
# VQQDEx9FbnRydXN0IENvZGUgU2lnbmluZyBDQSAtIE9WQ1MyAhALBv9gFrEjvi/8
# h6aykuQaMA0GCWCGSAFlAwQCAQUAoHwwEAYKKwYBBAGCNwIBDDECMAAwGQYJKoZI
# hvcNAQkDMQwGCisGAQQBgjcCAQQwHAYKKwYBBAGCNwIBCzEOMAwGCisGAQQBgjcC
# ARUwLwYJKoZIhvcNAQkEMSIEIIZIN9bcSKgj2uZ3HlQyuEuAcow/u3D+MVWSFdfK
# 07/+MA0GCSqGSIb3DQEBAQUABIICAGzH1OjyAgh2aGsXfexXck1XZbaCHgtnNPdI
# ucTM57Vwx062p8dbm6aJ5q+OdOUiic47SnNHTWF3QpGaAGu0yBIrzyItFKHlemTq
# Kkc5iiWIlrdtrwPMGDtJo/ju48No/z88OXQMRWqlIDKauytXvh3pIxPMnw2mPIB8
# rLW9J10ocBkJXINgEd+1/66lBtG8HJO6BCNrsB8nTahd1BheE2ZnL4kShDejMYOq
# p/kjyZvvkWBwwZf5eolwQw6cbPcvaRBURmS5fnW9ZGUiUYWj85/A9y/1pyMJqd5M
# aGX6H2VUEMuejZGiOFvOVPxm2kr3S/QtOErH0UkLWZ6qJnxmX9eOa5gTQpG2o+Gn
# ixQBASLHIfWE1PVrQZE5Luc91ZNg4HvnTMJWL4Kcci4CXUgx2nNrfms6b9EV8UOT
# A689DZIm6q49xWA79XxTA+b6lx8lHonGqlN6EmZ7NJtBxGYNMZmr8iIdB6nhB84Y
# +XKq1IbLYYPUX1qkEHT+rYv+3hBuUuO4DEVmbccJzxAz8iMfbrqgHjBP7Rs/K9Nx
# QTz6B76/j5eWfIWwkYFaIXwE8OhKsmQCFMbWEVYYXxmV6HZvh8GqR0dy9zb3I4yv
# BqO/FQF5x2mZM9JwihjfuyWxQ9IbriXPx1MSjUFsVXXz2oskrMce3ORmwJkUvhBF
# wVqSh95goYIYYDCCGFwGCisGAQQBgjcDAwExghhMMIIYSAYJKoZIhvcNAQcCoIIY
# OTCCGDUCAQMxDTALBglghkgBZQMEAgMwgfQGCyqGSIb3DQEJEAEEoIHkBIHhMIHe
# AgEBBgpghkgBhvpsCgMFMDEwDQYJYIZIAWUDBAIBBQAEIEN7Iv47IRZfqk7faoEb
# 1Fd++0zx5SRgN5cT+d7UCtfQAgkA64HfG/zzS+sYDzIwMjQwNDE3MjExNTEyWjAD
# AgEBoHmkdzB1MQswCQYDVQQGEwJDQTEQMA4GA1UECBMHT250YXJpbzEPMA0GA1UE
# BxMGT3R0YXdhMRYwFAYDVQQKEw1FbnRydXN0LCBJbmMuMSswKQYDVQQDEyJFbnRy
# dXN0IFRpbWVzdGFtcCBBdXRob3JpdHkgLSBUU0EyoIITDjCCBd8wggTHoAMCAQIC
# EE5A5DdU7eaMAAAAAFHTlH8wDQYJKoZIhvcNAQELBQAwgb4xCzAJBgNVBAYTAlVT
# MRYwFAYDVQQKEw1FbnRydXN0LCBJbmMuMSgwJgYDVQQLEx9TZWUgd3d3LmVudHJ1
# c3QubmV0L2xlZ2FsLXRlcm1zMTkwNwYDVQQLEzAoYykgMjAwOSBFbnRydXN0LCBJ
# bmMuIC0gZm9yIGF1dGhvcml6ZWQgdXNlIG9ubHkxMjAwBgNVBAMTKUVudHJ1c3Qg
# Um9vdCBDZXJ0aWZpY2F0aW9uIEF1dGhvcml0eSAtIEcyMB4XDTIxMDUwNzE1NDM0
# NVoXDTMwMTEwNzE2MTM0NVowaTELMAkGA1UEBhMCVVMxFjAUBgNVBAoMDUVudHJ1
# c3QsIEluYy4xQjBABgNVBAMMOUVudHJ1c3QgQ29kZSBTaWduaW5nIFJvb3QgQ2Vy
# dGlmaWNhdGlvbiBBdXRob3JpdHkgLSBDU0JSMTCCAiIwDQYJKoZIhvcNAQEBBQAD
# ggIPADCCAgoCggIBAKeBj/cURbZiQ/LYrtMlXkhPUb/FfZ9QHDXR1n5hKpQZbSdG
# pKYaXfdUUWqAIsaoZnVNVIPJXmgbq/ZbZLCtrSC9VO9Ga20C50WudfaOirkyLou4
# dxxSTXmIX6U6GMlQLJcnLb/aAH1jf+8y7EaHY9uan8NaITZ7+ZvVyqBuciz84fGe
# cE0IVhVvkKv7SLq518GCeIVlLn+1ycDiFLc3EUEG4orgqPblfrZ4BQHDYO1PB0Eu
# ChNJ45Cbf929+qy/ZFHRXJu09VznXP87m6WgGtd9CbLCt/9uHLzIfebpK/xysxTp
# SlUShJxEJXUd9irwT6UgPWgl62GXfA/ltj3zrsPBEbwbjszgRzBeQgCGceNYrAbK
# ZR97lKZLV2cMfl6teGdbVeNe68fY7Exuhsvz3Pifh6pyWBIPfab4+EI5Ozws5DJN
# SYzg4QDCOKCc+oQ+QdxuVq7GGlv0Z2gFAc0bv66HvJ1T9i7otmvkmd7FT4dYqNJl
# Hsgf1XJu7lkcVzsJcp3XyreQxs17RZKRQgNMfT/K8qq4wg6G8xCfRi6kZoZoWmgY
# cCk4EYBga4pDo3Ns47NrN//mnWcBkobfL0jR+1Bg1Vz+IdMBQmP+73C0F8CPqO7T
# wUtfEur9/S4Oh0Rg46n0whij4/3ODIQiDfOneNqT89s4z7kvM8b/BzxevkXTAgMB
# AAGjggErMIIBJzAOBgNVHQ8BAf8EBAMCAYYwEgYDVR0TAQH/BAgwBgEB/wIBATAd
# BgNVHSUEFjAUBggrBgEFBQcDAwYIKwYBBQUHAwgwOwYDVR0gBDQwMjAwBgRVHSAA
# MCgwJgYIKwYBBQUHAgEWGmh0dHA6Ly93d3cuZW50cnVzdC5uZXQvcnBhMDMGCCsG
# AQUFBwEBBCcwJTAjBggrBgEFBQcwAYYXaHR0cDovL29jc3AuZW50cnVzdC5uZXQw
# MAYDVR0fBCkwJzAloCOgIYYfaHR0cDovL2NybC5lbnRydXN0Lm5ldC9nMmNhLmNy
# bDAdBgNVHQ4EFgQUgrrWPZfOn89x6JI3r/2ztWk1V88wHwYDVR0jBBgwFoAUanIm
# etAe733nO2lR1GyNn5ASZqswDQYJKoZIhvcNAQELBQADggEBAB9eQQS2g3AkUyxV
# cx1lOsDstHsEmF5ZOBMJpFmUQl5Qv09sbiUgkJNYQA31GbRi7iRewgFYFQIdEAlv
# qNT7kn43OD4vFH2PHUM2ZLNmE18UzKVx91shS8aXvtyV/HB9ERzTId3QJDkpxf4K
# GqXPe3nuOm/e3L/pEd0WgwjTLI1/TagUeS8FYVI462DzFGh9y7KKrcCUXOQmDiyK
# 3UbDzuRWUcVW44W4TZtFcosH8Yr7Sbhf0fKWgV1pUiTxCCPS1iMP64vXfovBk2v6
# 8WJ7WOlQm5duF4gN4cZDmNeBYbaFnUfssZ6uPyA7Q53Yohzg1HwIwq92BvhiZnq2
# 9/rIrzUwggZvMIIEV6ADAgECAhAlvCvzKcoQfx6puoiF1J07MA0GCSqGSIb3DQEB
# DQUAMGkxCzAJBgNVBAYTAlVTMRYwFAYDVQQKDA1FbnRydXN0LCBJbmMuMUIwQAYD
# VQQDDDlFbnRydXN0IENvZGUgU2lnbmluZyBSb290IENlcnRpZmljYXRpb24gQXV0
# aG9yaXR5IC0gQ1NCUjEwHhcNMjEwNTA3MTkyMjE0WhcNNDAxMjI5MjM1OTAwWjBO
# MQswCQYDVQQGEwJVUzEWMBQGA1UEChMNRW50cnVzdCwgSW5jLjEnMCUGA1UEAxMe
# RW50cnVzdCBUaW1lIFN0YW1waW5nIENBIC0gVFMyMIICIjANBgkqhkiG9w0BAQEF
# AAOCAg8AMIICCgKCAgEAtQMqBraETzTMnjG+6o7Ge24gJA8D172BNjPmb9I+u8Se
# qqDMZl1CYDRR4i1hg7o8nmi4S8NmjqXeK7qbLUGl3tDLVwUtiNnARohtPDUNGZCs
# uAwZ8EzF20Vb3DvCvI1EAOA4//9dyM8Imj3Usut7lcPEpQDgYgFkWV+M43uLqflF
# t8awogo3Q7hGgauGfn8v3vBf680gzKUtMeOoLpmb4Hjd9NtdgmN9MDmIfiQG5nxg
# hHjIXCnllcv1e0/PWNyQz6eiexRUxOA3mDDBm7aunCn2ayVD5xJqPhmLvsDngddK
# MAh+fI3fPPNIkyHegYNqnQA+12C60bNiCJbptkkAOPRcdg3ad2o9oJXeHqHMWXpe
# n8/pdl1A8rUVtWBoLPNaySOD2PUaa2kopLbUrmIG3iZFvd7pVneq3GSwQc3ME/Fd
# cXSylAhLTzJufnC4d+nVMjfCet5CEuhxHOLYYr+qa0GGBlvmab4mCfSXwv5M9IBy
# hG7b3sUHqY4Vom3hG+AJD/WgB6C7xHdzCTwHk40BtRJnLU0ZtbjFyRCrcM9C/6Un
# HdYS1uQYRGD7OV2Eu2a1UxKdAxSCZSm07Ely/bND5cMCKSS6H/i4vl2jty6p/M9i
# 2fqXxBN441W5WnqE6/+7z2mIwLb5i2CFQRBQxwUQ0l52nT417PZPivkSez+vCzcC
# AwEAAaOCASwwggEoMBIGA1UdEwEB/wQIMAYBAf8CAQAwHQYDVR0OBBYEFCYP8MRI
# CBvN3ZH1VFS2s7P8mfEIMB8GA1UdIwQYMBaAFIK61j2Xzp/PceiSN6/9s7VpNVfP
# MDMGCCsGAQUFBwEBBCcwJTAjBggrBgEFBQcwAYYXaHR0cDovL29jc3AuZW50cnVz
# dC5uZXQwMQYDVR0fBCowKDAmoCSgIoYgaHR0cDovL2NybC5lbnRydXN0Lm5ldC9j
# c2JyMS5jcmwwDgYDVR0PAQH/BAQDAgGGMBMGA1UdJQQMMAoGCCsGAQUFBwMIMEUG
# A1UdIAQ+MDwwMAYEVR0gADAoMCYGCCsGAQUFBwIBFhpodHRwOi8vd3d3LmVudHJ1
# c3QubmV0L3JwYTAIBgZngQwBBAIwDQYJKoZIhvcNAQENBQADggIBAHY9RmiFXwnA
# zsn4x7jXfqdcOwtGc84iU4s49hG6d8TChyGPJ9kMlevoS+Vc15AVb6gqLIpq4CIR
# NHn6vZZ5yHYWzbFE6nFRD5tVyRh81YJ5nPOmt5OpBThyAgaXoaiyeX5eCOJubHma
# qIeUzn+PVTZAD58/UdVmvnurhuZuh7pOY12x2I4fy2fDDXDaeuuTshSE0j59kz8W
# DDagXHkAnwtu5YQq2IjBbQS3N/OSmnDgOwP9zbsngs4v/I3kssq++CWFxu10yD01
# XzgJPXlP6xvgQV6ldIE/URXA+EmIj7q7BdStRYIzecEljV9UPZHkuJ8XmABdrSjl
# T18b+/hUtaXEJPzwdjv4DV7MR7DKqrtuUd9TtUka8WCHt/n9aQZ1X3PWEZDygFiO
# 5cMPA5iMBFOzKYffP/ew6UxXgSREeWj+GzYH1mU97G/WwrWqYXxmBuzS5WmA+HGS
# LlziEUYYRUvW+gVk77veBspDJb3fjckBdGGI2/DHpj3jDHG8f0vnyOBu4NqA+QCe
# wGsl1aO6XCv+xGAgFEoHIzHlGsuRRLKae0PxW5i1IUu4MTsiVO55YfealBtmJyrG
# 3+FlXe9xOeZQJf2HNLmuXXwnO+dlk/gohLFvhpMsk5O2cK+h7sikYamEO/sKt/i4
# bSbuR071XgX6nA8RfdfLka5DIPcVWWO8MIIGtDCCBJygAwIBAgIQW3AmzJb2eBq7
# hfYR9W27bjANBgkqhkiG9w0BAQ0FADBOMQswCQYDVQQGEwJVUzEWMBQGA1UEChMN
# RW50cnVzdCwgSW5jLjEnMCUGA1UEAxMeRW50cnVzdCBUaW1lIFN0YW1waW5nIENB
# IC0gVFMyMB4XDTI0MDExOTE2NDc0N1oXDTM1MDQxODAwMDAwMFowdTELMAkGA1UE
# BhMCQ0ExEDAOBgNVBAgTB09udGFyaW8xDzANBgNVBAcTBk90dGF3YTEWMBQGA1UE
# ChMNRW50cnVzdCwgSW5jLjErMCkGA1UEAxMiRW50cnVzdCBUaW1lc3RhbXAgQXV0
# aG9yaXR5IC0gVFNBMjCCAiIwDQYJKoZIhvcNAQEBBQADggIPADCCAgoCggIBAKqG
# BDhej2i/g/Va/VSbxX06dZ3Epy1PbenQBR9R4zO7L186sZZUiukIfDStiYNVpJ8C
# zCML2Bt6QKTMn0LqE398uOC3SgIF4oNJBG8Ghs1zQamD1Jhn+vRSGXrjSlkoh1YX
# XUjxhpELuCJk20H5XVS1NaRuONwIDRpWGEDDfpg8QjOEDcWlYVHOIG1PAMtC8Kh1
# G9tG/incG57X+BgFUXYwm3S9lR8ZU2UVRyxWY2vpukisCji2wXL9IAqjdBprB3/R
# 3awXzrUcx3qJDOxtHDCNBfCvPZfMNbS1BQNmup52/DnayJrWEhLkQcKFIkDwy7qI
# 6dKHwqvGNDvlqw0ABg4GtYt70Qsi2jyoyVZ78DUElupZRoyP1KC1qYOwzoT7S0Tl
# QagW86/a0aAncN8H7zfzrpHO1CF/hpjcXu6PwNpH22nKELmRlYAVY/ErGsKWoR/C
# 9zzU9n0kByWoVJRvClXb2yRDbi47ontJu7XWGI/ZBLO+6se1ZJpLzlJhHisFHmM3
# Rp2xydVeh6Tq1zSyIUIXuDooOl6l7T79Nxd/YyDye9ZfyvMUP+HO/4p3lGZticHg
# 7H8EcDhkJ/WnF31NygyBG8Xl8TbFeq59KaRHNTj6DfRTAykF/G+hi8nN8biIQNPS
# 254lMa3X5lrYBLSEuwRGYWb/XTLe/NM+fE4Y737fAgMBAAGjggFlMIIBYTAMBgNV
# HRMBAf8EAjAAMB0GA1UdDgQWBBT1dhr4EJiSoRt3qQAs8ZViAgt5JzAfBgNVHSME
# GDAWgBQmD/DESAgbzd2R9VRUtrOz/JnxCDBoBggrBgEFBQcBAQRcMFowIwYIKwYB
# BQUHMAGGF2h0dHA6Ly9vY3NwLmVudHJ1c3QubmV0MDMGCCsGAQUFBzAChidodHRw
# Oi8vYWlhLmVudHJ1c3QubmV0L3RzMi1jaGFpbjI1Ni5wN2MwMQYDVR0fBCowKDAm
# oCSgIoYgaHR0cDovL2NybC5lbnRydXN0Lm5ldC90czJjYS5jcmwwDgYDVR0PAQH/
# BAQDAgeAMBYGA1UdJQEB/wQMMAoGCCsGAQUFBwMIMEwGA1UdIARFMEMwNwYKYIZI
# AYb6bAoBBzApMCcGCCsGAQUFBwIBFhtodHRwczovL3d3dy5lbnRydXN0Lm5ldC9y
# cGEwCAYGZ4EMAQQCMA0GCSqGSIb3DQEBDQUAA4ICAQCpq32/GgCFW9ztcSZelIKv
# XIJL+hw5xVO1wwbdA2B++WGaQSDCrkEXLqVCCJEJXpcwmIFF1rK1hD2L6YcH2gE/
# 4AhYZxEglB7AFo5jd72bbAWNPbXLmjGpSlxiMqITzpCcDceQT2fD4Kh/h4RL1/LB
# IQ8Q5xe5+70JzL6KpgH+TYJPaSCadWvkDfdQJlP8M/f+G66nasCo6x6vItM9YBWf
# uuTVKaS1V5+xZfaxTGHV0sjat6yqGB8gGIqm+PUH+XIr9mPYJZwW/vuU34Z0oiC0
# 8gjfFu85gHY6krHQloNq2nxTiVkbDVzxuynyIU9hLCSJ+HYJf9lokI5XDGTVdgBx
# WfPUhMg8oYh/qG/zvW5/sKiO5pNs6RTFdRwMK87yB/u1rtRBMb/BzGmYGK5bCI+r
# XGS5hGq9XDih309Ct0tbUtPe9zvLn2/CWUIiTgsYtyLIL/xqCYXtUClResSUS7Zx
# ZYpKhppLeEtWiqHeL4uTWzLMobnvqbYM2AYLFzC/DW+Pk6zLEfXxxwCe0BnXnF+O
# VHxrqUM7Y9zBRKu9jDpD8HI/JJXNhJFdIQANeGsNXCMkn0yvNMsUWt6x5EVQ/a9v
# pwWCp3PPz+mS6nMqseHdogPmILwdZxpWwi72rF/KToGw5lodc8/SFPNQeCgFThVu
# hST4M6iIr3Rt9FLOjKiEnjGCBBYwggQSAgEBMGIwTjELMAkGA1UEBhMCVVMxFjAU
# BgNVBAoTDUVudHJ1c3QsIEluYy4xJzAlBgNVBAMTHkVudHJ1c3QgVGltZSBTdGFt
# cGluZyBDQSAtIFRTMgIQW3AmzJb2eBq7hfYR9W27bjALBglghkgBZQMEAgOgggGJ
# MBoGCSqGSIb3DQEJAzENBgsqhkiG9w0BCRABBDAcBgkqhkiG9w0BCQUxDxcNMjQw
# NDE3MjExNTEyWjApBgkqhkiG9w0BCTQxHDAaMAsGCWCGSAFlAwQCA6ELBgkqhkiG
# 9w0BAQ0wTwYJKoZIhvcNAQkEMUIEQKaGjygpvqCXC6hDVGKg5jfmD7JCiqmdfBFF
# 5tXlHujy9vbUSPAL6cknnPIorbuVudaOlHkixmXm2wWPSTRGFNEwgdAGCyqGSIb3
# DQEJEAIvMYHAMIG9MIG6MIG3MAsGCWCGSAFlAwQCAwRAORFCLhcCPqqZJJl0Dg7I
# 6NJy01X8XPZs+RFWr5dVNOSAr9O8fExbV+/R/7zbdRRl4NQFWoVOmxE52ffsVdX+
# JDBmMFKkUDBOMQswCQYDVQQGEwJVUzEWMBQGA1UEChMNRW50cnVzdCwgSW5jLjEn
# MCUGA1UEAxMeRW50cnVzdCBUaW1lIFN0YW1waW5nIENBIC0gVFMyAhBbcCbMlvZ4
# GruF9hH1bbtuMAsGCSqGSIb3DQEBDQSCAgA/izNIPyritZZB3LmQBx09aE8HOz35
# TXMEbS/u51guIAO+Z3Vao+3C//SkTsCUdPywAglJcc/2bskH/ndciaPrky/HnZgC
# +6yWVZC37j9RdThJLA7UPYun6ulMJY39gSzDM4SlGQI+5LFHLBNR43clgLfhnWuy
# cEdL7VHs2MMke4QiCEjvXnsb5Tz/stK8k6WG9vIZMCNJ31soo0zS0iIfNgXHe3c/
# XVSMbYx2nbxMjRcBfvDsl/OqYZDMljy90rQr73iCK2ZhVDrfvNbFhEgZ6qg6bPCY
# DqXlFiomh6sSd+hrXjv9G5PSZ5td6M81euRzjSqzicck3MK/TpLXSkHxGlQfrPTS
# EwUV9tczYEUmIF0ZLdcJNrlfIpfjQJUt5v6+G7KPpAmq/42jRFXXP9wbeT22q8Ci
# 1OvuTlgbR7vQqYstCa5l0fnHApY4z83P5Ro27pgLFhUywiq84BHNg7c+gwj1uhJ+
# rKgn+ht23H71/QJC/Q2CgAWb9U8w3LI9tw3D7ZL6BSmyB7gJTNFJl6YZhN2VKy/L
# uBJ0+P2D2Q9SvyZzZexs4m6bXH5+His8lOBDEEMoJNW9OO6JbbfxZZmuy4agsou3
# h6i1eom4NBLB266X2uyLI3BY8J7yqVMZ3hJYW8qvPfAp2zig1N8K6IKtvRgTntcb
# 1GrFk036ZDE0Pw==
# SIG # End signature block
