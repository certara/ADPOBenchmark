$PROBLEM    ADPO benchmark - NONMEM
$ABBR DERIV2=NO 
$INPUT   ID TIME AMT DV EVID BQL WT
$DATA      c:/git/adpoBenchMark/data/pydarwin_makedata/data_sim.csv IGNORE=@
$SUBROUTINE ADVAN6 TOL=7
$MODEL
  COMP=(DEPOT,DEFDOSE)
  COMP=(CENTRAL,NODOSE,DEFOBS)
  COMP=(PERI,NODOSE)
$PK      
  CWT = WT/70 

  TVVMAX= THETA(1)   
  VMAX=TVVMAX*EXP(ETA(1)) 
  TVKM = THETA(2)
  KM = TVKM  
  TVV2=THETA(3) 
  V2=TVV2 *EXP(ETA(2))    
  TVKA=THETA(4) 
  KA=TVKA   
  K23=THETA(5)
  K32=THETA(6)
  SC = V2
$ERROR    
  IPRED = F
  IOBS = F*EXP(EPS(2)) + EPS(1)
  Y=F*EXP(EPS(2)) + EPS(1)
$DES
  CONC = A(2)/SC
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1)-VMAX*CONC/(KM+CONC) -K23*A(2) + K32*A(3)
  DADT(3) = K23*A(2)-K32*A(3)

$THETA   
  (100,4000,10000)	; THETA(1) VMAX UNITS =  mass/time
  (100,1000,3000)	; THETA(2) KM UNITS = mass/volume
  (1,50) 		; THETA(3) V  UNITS = volume
  (0.01,1.2,3) 		; THETA(4) KA UNITS = 1/time    


  (0.0001,2,10)	 ;; THETA(5) K23
  (0.0001,3,10)	 ;; THETA(6) K32 
  ;; Start OMEGA5
$OMEGA  
  0.1		; ETA(1) ETA ON VMAX
  0.1		; ETA(2) ETA ON V2

$SIGMA     
  (1) ;; EPS(1) ADDITIVE
  (0.3) ;; EPS(2) PROPORTIONAL
$SIM (2345) ONLYSIM 
$TABLE ID TIME AMT IOBS EVID WT FILE=c:\git\adpoBenchmark\data\OUT_1_1_0_0.DAT NOPRINT NOHEADER NOAPPEND
  ;;; Model Identifier =  1,1,0,0
$TABLE ID VMAX KM V2 KA FIRSTONLY NOAPPEND NOPRINT FILE= c:\git\adpoBenchmark\data\PARMS_1_1_0_0.DAT

;; Phenotype: ([('COMP', 1), ('ETAs', 1), ('V~WT', 0), ('GAMMA', 0)])
;; Genotype: [1, 1, 0, 0]
;; Num non-influential tokens: 0
