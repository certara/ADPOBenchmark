\n$env:NLME_HASH = 1770978959
$INSTALLDIR="c:\program files\Certara\nnlme_engine_dualmay21_old"
if($INSTALLDIR -eq "" -or $INSTALLDIR -eq $null)
{
  "Installation directory is not specified"
  exit 1
}
powershell -noninteractive -executionpolicy remotesigned -File $INSTALLDIR\generic_run.ps1 none $INSTALLDIR $shared_directory C:\git\adpoBenchmark\NLME\8\Run8 C:\git\adpoBenchmark\NLME\8\Run8\jobControlFile.txt 1 WorkFlow
