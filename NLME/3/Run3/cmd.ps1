\n$env:NLME_HASH = 1770978959
$INSTALLDIR="D:\NLME_Engine_ADPO"
if($INSTALLDIR -eq "" -or $INSTALLDIR -eq $null)
{
  "Installation directory is not specified"
  exit 1
}
powershell -noninteractive -executionpolicy remotesigned -File $INSTALLDIR\generic_run.ps1 none $INSTALLDIR $shared_directory D:\git\ADPOBenchmark\NLME\3\Run3 D:\git\ADPOBenchmark\NLME\3\Run3\jobControlFile.txt 1 WorkFlow
