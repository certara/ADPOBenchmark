$SIZES MAXFCN = 5000000
$PROBLEM    Dual Numbers benchmark
$INPUT       ID TIME AMT DV WT
  ;; $DATA      c:/git/adpoBenchmark/data/sim_2_2_0_0.csv IGNORE=@
$DATA      ..\\..\\data\\sim_2_2_0_0.csv IGNORE=@

$SUBROUTINE ADVAN6 TOL=7
$MODEL
  COMP=(DEPOT,DEFDOSE)
  COMP=(CENTRAL,NODOSE,DEFOBS)
  COMP=(PERI1,NODOSE)
  COMP=(PERI2,NODOSE)
$PK
  CWT = WT/70

  TVVMAX= THETA(1)
  VMAX=TVVMAX*EXP(ETA(1)) 
  TVKM = THETA(2)
  KM = TVKM *EXP(ETA(2)) 
  TVV2=THETA(3) 
  V2=TVV2 *EXP(ETA(3)) 
  TVKA=THETA(4)
  KA=TVKA 
  K23=THETA(5)
  K32=THETA(6)
  K24=THETA(7)
  K42=THETA(8)
  SC = V2
$ERROR
  IPRED = F
  Y=F*EXP(EPS(2)) + EPS(1)
$DES
  CONC = A(2)/SC
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1)-VMAX*CONC/(KM+CONC) -K23*A(2)+K32*A(3) -K24*A(2)+K42*A(4)
  DADT(3) = K23*A(2)-K32*A(3)
  DADT(4) = K24*A(2)-K42*A(4)

$THETA
  (100,4000)	; THETA(1) VMAX UNITS =  MG/HR
  (100,4000)	; THETA(2) KM UNITS =  CONC
  (0,50) 	; THETA(3) V  UNITS = L
  (0,1.2) 	; THETA(4) KA UNITS = 1/HR


  (0,0.7)	 ;; THETA(5) K23
  (0,1)	 ;; THETA(6) K32
  (0.0001,0.1) 	 ;; THETA(7) K24
  (0.0001,0.05) 	 ;; THETA(8) K42
; empty $OMEGA
$OMEGA 
  0.1		; ETA(1) ETA ON VMAX
  0.1		; ETA(2) ETA ON KM
  0.1		; ETA(3) ETA ON V2
$SIGMA
  (1) ;; EPS(1) ADDITIVE
  (0.3) ;; EPS(2) PROPORTIONAL
  ;;$EST METHOD=COND INTER MAXEVALS=9999 NOHABORT NOOMEGABOUNDTEST NOSIGMABOUNDTEST NOTHETABOUNDTEST
  ;;$COV UNCOND PRECOND=2
  ;;; Model Identifier =  2,2,0,0

$TABLE ID TIME IPRED DV NOAPPEND ONEHEADER NOPRINT FILE= Run1Preds.dat

;; Phenotype: ([('COMP', 2), ('ETAs', 2), ('V~WT', 0), ('GAMMA', 0)])
;; Genotype: [2, 2, 0, 0]
;; Num non-influential tokens: 0
