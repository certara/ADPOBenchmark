\n$env:NLME_HASH = 1770978959
$INSTALLDIR="D:\nlme_459"
if($INSTALLDIR -eq "" -or $INSTALLDIR -eq $null)
{
  "Installation directory is not specified"
  exit 1
}
powershell -noninteractive -executionpolicy remotesigned -File $INSTALLDIR\generic_run.ps1 none $INSTALLDIR $shared_directory D:\git\ADPOBenchmark\NLME\62\run62 D:\git\ADPOBenchmark\NLME\62\run62\jobControlFile.txt 1 WorkFlow