$PROBLEM    Dual Numbers benchmark - NONMEM
$ABBR DERIV2=NO 
$INPUT       ID TIME AMT DV WT
$DATA      c:/git/adpoBenchmark/data/sim_0_2_1_0.csv IGNORE=@
  ;; data file name must match out file name in the template file for simulations in data folder
  ;;from data\template.txt
  ;; $TABLE ID TIME AMT IOBS EVID WT FILE=c:\git\dual_numbers\data\OUT_0_2_1_0.DAT NOPRINT NOHEADER NOAPPEND
  ;; and file read/write in make_data.r
  ;; file_name <- paste0(this_comp,"_",this_eta,"_",this_vwt,"_",this_gamma)
$SUBROUTINE ADVAN6 TOL=7
$MODEL
  COMP=(DEPOT,DEFDOSE)
  COMP=(CENTRAL,NODOSE,DEFOBS)

$PK      
  CWT = WT/70 

  TVVMAX= THETA(1)   
  VMAX=TVVMAX*EXP(ETA(1)) 
  TVKM = THETA(2)
  KM = TVKM *EXP(ETA(2))  
  TVV2=THETA(3) *CWT**THETA(5)
  V2=TVV2 *EXP(ETA(3))    
  TVKA=THETA(4) 
  KA=TVKA   

  SC = V2
$ERROR 
 IPRED = F   

  Y=F*EXP(EPS(2)) + EPS(1)
$DES
  CONC = A(2)/SC
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1)-VMAX*CONC/(KM+CONC) 


$THETA   
  (100,4000)	; THETA(1) VMAX UNITS =  mass/time
  (100,4000)	; THETA(2) KM UNITS = mass/volume
  (1,50) 		; THETA(3) V  UNITS = volume
  (0.01,1.2,3) 		; THETA(4) KA UNITS = 1/time    

  (0,1.1,10) 	;; THETA(5) V~WT

  ;; Start OMEGA
$OMEGA  
  0.1		; ETA(1) ETA ON VMAX
  0.1		; ETA(2) ETA ON KM
  0.1		; ETA(3) ETA ON V2
$SIGMA     
  (1) ;; EPS(1) ADDITIVE
  (0.3) ;; EPS(2) PROPORTIONAL
$EST METHOD=COND INTER MAXEVALS=9999 NOHABORT NOOMEGABOUNDTEST NOSIGMABOUNDTEST NOTHETABOUNDTEST
$COV UNCOND PRECOND=2
  ;;; Model Identifier =  0,2,1,0

;; Phenotype: ([('COMP', 0), ('ETAs', 2), ('V~WT', 1), ('GAMMA', 0)])
;; Genotype: [0, 2, 1, 0]
;; Num non-influential tokens: 0

$TABLE ID TIME IPRED DV NOAPPEND ONEHEADER NOPRINT FILE= Run1Preds.dat