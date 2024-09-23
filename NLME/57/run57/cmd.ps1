\n$env:NLME_HASH = 1770978959
$INSTALLDIR="D:\NLME_Engine_noHessian"
if($INSTALLDIR -eq "" -or $INSTALLDIR -eq $null)
{
  "Installation directory is not specified"
  exit 1
}
powershell -noninteractive -executionpolicy remotesigned -File $INSTALLDIR\generic_run.ps1 none $INSTALLDIR $shared_directory D:\git\ADPOBenchmark\NLME\57\run57 D:\git\ADPOBenchmark\NLME\57\run57\jobControlFile.txt 1 WorkFlow
