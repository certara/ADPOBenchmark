$PROBLEM    Dual Numbers benchmark
$INPUT       ID TIME AMT DV WT
  ;; $DATA      c:/git/adpoBenchmark/data/sim_0_3_0_1.csv IGNORE=@
$DATA      D:/git/adpoBenchmark/data/sim_0_3_0_1.csv IGNORE=@

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
  TVV2=THETA(3) 
  V2=TVV2 *EXP(ETA(3))    
  TVKA=THETA(4) 
  KA=TVKA *EXP(ETA(4))   

  SC = V2
$ERROR    
  IPRED = F
  Y=F*EXP(EPS(2)) + EPS(1)
$DES
  CONC = A(2)/SC
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1)-VMAX*CONC**THETA(5)/(KM**THETA(5)+CONC**THETA(5)) 


$THETA   
  (100,4000)	; THETA(1) VMAX UNITS =  MG/HR
  (100,4000)	; THETA(2) KM UNITS =  CONC
  (0,50) 	; THETA(3) V  UNITS = L
  (0,1.2) 	; THETA(4) KA UNITS = 1/HR    
  (0,1.7) 	;; THETA(5) GAMMA
; empty $OMEGA    
$OMEGA 
  0.3		; ETA(1) ETA ON VMAX
  0.3		; ETA(2) ETA ON KM
  0.3		; ETA(3) ETA ON V2
  0.3		; ETA(4) ETA ON KA
$SIGMA     
  (1) ;; EPS(1) ADDITIVE
  (0.3) ;; EPS(2) PROPORTIONAL
  ;;$EST METHOD=COND INTER MAXEVALS=9999 NOHABORT NOOMEGABOUNDTEST NOSIGMABOUNDTEST NOTHETABOUNDTEST
  ;;$COV UNCOND PRECOND=2
  ;;; Model Identifier =  0,3,0,1

$TABLE ID TIME IPRED DV NOAPPEND ONEHEADER NOPRINT FILE= Run1Preds.dat
;; Phenotype: ([('COMP', 0), ('ETAs', 3), ('V~WT', 0), ('GAMMA', 1)])
;; Genotype: [0, 3, 0, 1]
;; Num non-influential tokens: 0
