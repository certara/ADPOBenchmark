
$INSTALLDIR="c:\program files\Certara\nlme_engine"
if($INSTALLDIR -eq "" -or $INSTALLDIR -eq $null)
{
  "Installation directory is not specified"
  exit 1
}
powershell -noninteractive -executionpolicy remotesigned -File $INSTALLDIR\generic_run.ps1 none $INSTALLDIR $shared_directory C:\git\adpo_speed\CompileResults\NLME\19\Run19 C:\git\adpo_speed\CompileResults\NLME\19\Run19\jobControlFile.txt 1 WorkFlow
