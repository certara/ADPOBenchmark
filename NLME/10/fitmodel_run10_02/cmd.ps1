
$INSTALLDIR="C:\Program Files\Certara\NLME_Engine"
if($INSTALLDIR -eq "" -or $INSTALLDIR -eq $null)
{
  "Installation directory is not specified"
  exit 1
}
$shared_directory="C:\git\adpoBenchmark\NLME\10\fitmodel_run10_02"
powershell -noninteractive -executionpolicy remotesigned -File $INSTALLDIR\generic_run.ps1 local_mpi $INSTALLDIR $shared_directory C:\git\adpoBenchmark\NLME\10\fitmodel_run10_02 C:\git\adpoBenchmark\NLME\10\fitmodel_run10_02\jobControlFile.txt 4 WorkFlow
