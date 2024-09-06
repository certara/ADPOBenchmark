      MODULE NMPRD4P
      USE SIZES, ONLY: DPSIZE
      USE NMPRD4,ONLY: VRBL
      IMPLICIT NONE
      SAVE
      REAL(KIND=DPSIZE), DIMENSION (:),POINTER ::COM
      REAL(KIND=DPSIZE), POINTER ::CWT,TVVMAX,VMAX,TVKM,KM,TVV2,V2
      REAL(KIND=DPSIZE), POINTER ::TVKA,KA,SC,IPRED,Y,CONC,A000032
      REAL(KIND=DPSIZE), POINTER ::A000034,A000036,A000037,D000001
      REAL(KIND=DPSIZE), POINTER ::D000002,D000003,D000039,D000038
      REAL(KIND=DPSIZE), POINTER ::D000037,C000033,C000032,D000042
      REAL(KIND=DPSIZE), POINTER ::D000041,D000040,A000038,A000039
      REAL(KIND=DPSIZE), POINTER ::A000040,E000005,F000082,E000013
      REAL(KIND=DPSIZE), POINTER ::E000015,F000097,F000094,F000098
      REAL(KIND=DPSIZE), POINTER ::F000095
      CONTAINS
      SUBROUTINE ASSOCNMPRD4
      COM=>VRBL
      CWT=>COM(000001);TVVMAX=>COM(000002);VMAX=>COM(000003)
      TVKM=>COM(000004);KM=>COM(000005);TVV2=>COM(000006)
      V2=>COM(000007);TVKA=>COM(000008);KA=>COM(000009)
      SC=>COM(000010);IPRED=>COM(000011);Y=>COM(000012)
      CONC=>COM(000013);A000032=>COM(000014);A000034=>COM(000015)
      A000036=>COM(000016);A000037=>COM(000017);D000001=>COM(000018)
      D000002=>COM(000019);D000003=>COM(000020);D000039=>COM(000021)
      D000038=>COM(000022);D000037=>COM(000023);C000033=>COM(000024)
      C000032=>COM(000025);D000042=>COM(000026);D000041=>COM(000027)
      D000040=>COM(000028);A000038=>COM(000029);A000039=>COM(000030)
      A000040=>COM(000031);E000005=>COM(000032);F000082=>COM(000033)
      E000013=>COM(000034);E000015=>COM(000035);F000097=>COM(000036)
      F000094=>COM(000037);F000098=>COM(000038);F000095=>COM(000039)
      END SUBROUTINE ASSOCNMPRD4
      END MODULE NMPRD4P
      SUBROUTINE MODEL (IDNO,NCM,NPAR,IR,IATT,LINK)                           
      USE PRMOD_CHAR, ONLY: NAME                                              
      USE SIZES,     ONLY: DPSIZE,ISIZE,SD
      USE PRDIMS,    ONLY: GPRD,HPRD,GERD,HERD,GPKD
      INTEGER(KIND=ISIZE) :: IDNO,NCM,NPAR,IR,IATT,LINK,I,J                   
      DIMENSION :: IATT(IR,*),LINK(IR,*)                                      
      SAVE
      INTEGER(KIND=ISIZE), DIMENSION (2,7) :: MOD
      CHARACTER(LEN=SD), DIMENSION(2) :: CMOD
      DATA (MOD(I,  1),I=  1,  2)/&
      1,1 /
      DATA (MOD(I,  2),I=  1,  2)/&
      1,1 /
      DATA (MOD(I,  3),I=  1,  2)/&
      1,0 /
      DATA (MOD(I,  4),I=  1,  2)/&
      0,1 /
      DATA (MOD(I,  5),I=  1,  2)/&
      1,0 /
      DATA (MOD(I,  6),I=  1,  2)/&
      0,0 /
      DATA (MOD(I,  7),I=  1,  2)/&
      0,0 /
      DATA (CMOD(I),I=  1,  2) &
      /'DEPOT','CENTRAL'/
      FORALL (I=1:2) NAME(I)=CMOD(I)
      FORALL (I=1:2,J=1:7) IATT(I,J)=MOD(I,J)
      IDNO=9999                                                               
      NCM=  2
      NPAR=004
      RETURN
      END
      SUBROUTINE PK(ICALL,IDEF,THETA,IREV,EVTREC,NVNT,INDXS,IRGG,GG,NETAS)      
      USE NMPRD4P
      USE SIZES,     ONLY: DPSIZE,ISIZE
      USE PRDIMS,    ONLY: GPRD,HPRD,GERD,HERD,GPKD
      USE NMBAYES_REAL,    ONLY: PRIORINFO
      USE NMPRD_REAL,ONLY: ETA,EPS                                            
      USE NMPRD_INT, ONLY: MSEC=>ISECDER,MFIRST=>IFRSTDER,COMACT,COMSAV,IFIRSTEM
      USE NMPRD_INT, ONLY: MDVRES,ETASXI,NPDE_MODE,NOFIRSTDERCODE
      USE NMPRD_REAL, ONLY: DV_LOQ,CDF_L,DV_LAQ,CDF_LA
      USE NMPRD_INT, ONLY: IQUIT
      USE PROCM_INT, ONLY: NEWIND=>PNEWIF                                       
      USE NMBAYES_REAL, ONLY: LDF                                             
      IMPLICIT REAL(KIND=DPSIZE) (A-Z)                                          
      REAL(KIND=DPSIZE) :: EVTREC                                               
      SAVE
      INTEGER(KIND=ISIZE) :: FIRSTEM
      INTEGER(KIND=ISIZE) :: ICALL,IDEF,IREV,NVNT,INDXS,IRGG,NETAS              
      DIMENSION :: IDEF(7,*),THETA(*),EVTREC(IREV,*),INDXS(*),GG(IRGG,GPKD+1,*) 
      FIRSTEM=IFIRSTEM
      IF (ICALL <= 1) THEN                                                      
      CALL ASSOCNMPRD4
      IDEF(   1,0001)=  -9
      IDEF(   1,0002)=  -1
      IDEF(   1,0003)=   0
      IDEF(   1,0004)=   0
      IDEF(   2,0003)=   0
      IDEF(   2,0004)=   0
      IDEF(   3,0002)=   5
      CALL GETETA(ETA)                                                          
      IF (IQUIT == 1) RETURN                                                    
      RETURN                                                                    
      ENDIF                                                                     
      IF (NEWIND /= 2) THEN
      IF (ICALL == 4) THEN
      CALL SIMETA(ETA)
      ELSE
      CALL GETETA(ETA)
      ENDIF
      IF (IQUIT == 1) RETURN
      ENDIF
 !  level            0
      WT=EVTREC(NVNT,005)
      CWT=WT/70.D0 
      TVVMAX=THETA(001) 
      B000001=DEXP(ETA(001)) 
      VMAX=TVVMAX*B000001 
      TVKM=THETA(002) 
      B000003=DEXP(ETA(002)) 
      KM=TVKM*B000003 
      TVV2=THETA(003) 
      B000005=DEXP(ETA(003)) 
      V2=TVV2*B000005 
      TVKA=THETA(004) 
      KA=TVKA 
      SC=V2 
      P000001=SC 
      P000002=KA 
      P000003=KM 
      P000004=VMAX 
      IF (FIRSTEM == 1) THEN
!                      A000032 = DERIVATIVE OF VMAX W.R.T. ETA(001)
      A000032=TVVMAX*B000001 
!                      A000034 = DERIVATIVE OF KM W.R.T. ETA(002)
      A000034=TVKM*B000003 
!                      A000036 = DERIVATIVE OF V2 W.R.T. ETA(003)
      A000036=TVV2*B000005 
!                      A000037 = DERIVATIVE OF SC W.R.T. ETA(003)
      A000037=A000036 
!                      A000038 = DERIVATIVE OF P000001 W.R.T. ETA(003)
      A000038=A000037 
!                      A000039 = DERIVATIVE OF P000003 W.R.T. ETA(002)
      A000039=A000034 
!                      A000040 = DERIVATIVE OF P000004 W.R.T. ETA(001)
      A000040=A000032 
      GG(0001,1,1)=P000001
      GG(0001,0004,1)=A000038
      GG(0002,1,1)=P000002
      GG(0003,1,1)=P000003
      GG(0003,0003,1)=A000039
      GG(0004,1,1)=P000004
      GG(0004,0002,1)=A000040
      GG(0005,1,1)=SC
      GG(0005,0004,1)=A000037
      ELSE
      GG(0001,1,1)=P000001
      GG(0002,1,1)=P000002
      GG(0003,1,1)=P000003
      GG(0004,1,1)=P000004
      GG(0005,1,1)=SC
      ENDIF
      RETURN
      END
      SUBROUTINE ERROR (ICALL,IDEF,THETA,IREV,EVTREC,NVNT,INDXS,F,G,HH)       
      USE NMPRD4P
      USE SIZES,     ONLY: DPSIZE,ISIZE
      USE PRDIMS,    ONLY: GPRD,HPRD,GERD,HERD,GPKD
      USE NMPRD_REAL,ONLY: ETA,EPS                                            
      USE NMPRD_INT, ONLY: MSEC=>ISECDER,MFIRST=>IFRSTDER,IQUIT,IFIRSTEM
      USE NMPRD_INT, ONLY: MDVRES,ETASXI,NPDE_MODE,NOFIRSTDERCODE
      USE NMPRD_REAL, ONLY: DV_LOQ,CDF_L,DV_LAQ,CDF_LA
      USE NMPRD_INT, ONLY: NEWL2
      USE PROCM_INT, ONLY: NEWIND=>PNEWIF                                       
      IMPLICIT REAL(KIND=DPSIZE) (A-Z)                                        
      REAL(KIND=DPSIZE) :: EVTREC                                             
      SAVE
      INTEGER(KIND=ISIZE) :: ICALL,IDEF,IREV,NVNT,INDXS                       
      DIMENSION :: IDEF(*),THETA(*),EVTREC(IREV,*),INDXS(*)                   
      REAL(KIND=DPSIZE) :: G(GERD,*),HH(HERD,*)                               
      INTEGER(KIND=ISIZE) :: FIRSTEM
      FIRSTEM=IFIRSTEM
      IF (ICALL <= 1) THEN                                                    
      CALL ASSOCNMPRD4
      IDEF(2)=-1
      IDEF(3)=000
      RETURN
      ENDIF
      IF (ICALL == 4) THEN
      IF (NEWL2 == 1) THEN
      CALL SIMEPS(EPS)
      IF (IQUIT == 1) RETURN
      ENDIF
      ENDIF
      D000001=G(001,1)
      D000002=G(002,1)
      D000003=G(003,1)
 !  level            0
      IPRED=F 
      B000001=DEXP(EPS(002)) 
      Y=F*B000001+EPS(001) 
!                      C000032 = DERIVATIVE OF Y W.R.T. EPS(002)
      C000032=F*B000001 
!                      C000033 = DERIVATIVE OF Y W.R.T. EPS(001)
      C000033=1.D0 
      IF (FIRSTEM == 1) THEN !1
!                      D000037 = DERIVATIVE OF Y W.R.T. ETA(003)
      D000037=B000001*D000003 
!                      D000038 = DERIVATIVE OF Y W.R.T. ETA(002)
      D000038=B000001*D000002 
!                      D000039 = DERIVATIVE OF Y W.R.T. ETA(001)
      D000039=B000001*D000001 
!                      D000040 = DERIVATIVE OF C000032 W.R.T. ETA(003)
      D000040=B000001*D000003 
!                      D000041 = DERIVATIVE OF C000032 W.R.T. ETA(002)
      D000041=B000001*D000002 
!                      D000042 = DERIVATIVE OF C000032 W.R.T. ETA(001)
      D000042=B000001*D000001 
      G(001,1)=D000039
      G(002,1)=D000038
      G(003,1)=D000037
      ENDIF !1
      HH(001,1)=C000033
      HH(002,1)=C000032
      IF (FIRSTEM == 1) THEN !2
      HH(002,002)=D000042
      HH(002,003)=D000041
      HH(002,004)=D000040
      ENDIF !2
      F=Y
      RETURN
      END
      SUBROUTINE TOL(NRD,ANRD,NRDC,ANRDC)
      USE SIZES,     ONLY: ISIZE
      INTEGER(KIND=ISIZE) :: NRD(0:*), ANRD(0:*), NRDC(0:*), ANRDC(0:*)
      NRD(1)=7 
      RETURN
      END
      SUBROUTINE DES (A,P,T,DADT,IR,DA,DP,DT)                                 
      USE NMPRD4P
      USE SIZES,     ONLY: DPSIZE,ISIZE
      USE PRDIMS,    ONLY: GPRD,HPRD,GERD,HERD,GPKD
      USE NMPRD_INT, ONLY: IERPRD,IERPRDU,NETEXT,IQUIT                        
      USE NMPRD_CHAR,ONLY: ETEXT                                              
      USE NMPRD_INT, ONLY: MSEC=>ISECDER,MFIRST=>IFRSTDER,IFIRSTEM,IFIRSTEMJAC
      USE PRCOM_INT, ONLY: MITER
      USE NMPRD_INT, ONLY: MDVRES,ETASXI,NPDE_MODE,NOFIRSTDERCODE
      USE NMPRD_REAL, ONLY: DV_LOQ,CDF_L,DV_LAQ,CDF_LA
      USE PROCM_REAL,ONLY: THETA=>THETAS
      USE PRMOD_INT, ONLY: ICALL=>ICALLD,IDEFD,IDEFA
      IMPLICIT REAL(KIND=DPSIZE) (A-Z)                                        
      SAVE
      INTEGER(KIND=ISIZE) :: IR                                               
      DIMENSION :: A(*),P(*),DADT(*),DA(IR,*),DP(IR,*),DT(*)                  
      INTEGER(KIND=ISIZE) :: FIRSTEM,IFIRSTEMJACIN
      IF(MITER==1.OR.MITER==4) IFIRSTEM=1
      FIRSTEM=IFIRSTEM
      IFIRSTEMJACIN=IFIRSTEMJAC
      IF(NOFIRSTDERCODE/=1) THEN
      IFIRSTEMJAC=FIRSTEM
      ELSE
      IFIRSTEMJAC=0
      ENDIF
      IF(IFIRSTEMJACIN==-2) RETURN
      IF (ICALL == 1) THEN
      CALL ASSOCNMPRD4
      IDEFD(1)=  5
      IDEFD(2)=0
      DA(   1,1)=0000014280
      DA(   2,1)=0000028441
      DA(   3,1)=0000028560
      DA(   4,1)=0000000000
      DP(   1,1)=0000014399
      DP(   2,1)=0000028441
      DP(   3,1)=0000028560
      DP(   4,1)=0000028679
      DP(   5,1)=0000028798
      DP(   6,1)=0000000000
      DT(   1)=0000000000
      RETURN
      ENDIF
 !  level            0
 !  level            0
      CONC=A(2)/P(001) 
      DADT(1)=-P(002)*A(1) 
      IF(CONC == 0.D0.AND.THETA(005) <= 0.D0)THEN 
      IERPRD=1
      ETEXT(1)='DES SUBROUTINE: ERROR IN COMPUTATION'
      ETEXT(2)='ATTEMPT TO COMPUTE 0**POWER WITH POWER<=0.'
      RETURN
      ENDIF 
      IF(CONC <  0.D0)THEN 
      IERPRD=1
      ETEXT(1)='DES SUBROUTINE: ERROR IN COMPUTATION'
      ETEXT(2)='ATTEMPT TO COMPUTE BASE**POWER WITH BASE<0.'
      RETURN
      ENDIF 
      B000003=0.D0 
      IF(CONC == 0.D0)THEN 
      B000003=1.D0 
      ENDIF 
      B000004=1.D0-B000003 
      B000005=CONC+B000003 
      B000006=B000005**THETA(005) 
      B000007=B000004*B000006 
      IF(P(003) == 0.D0.AND.THETA(005) <= 0.D0)THEN 
      IERPRD=1
      ETEXT(1)='DES SUBROUTINE: ERROR IN COMPUTATION'
      ETEXT(2)='ATTEMPT TO COMPUTE 0**POWER WITH POWER<=0.'
      RETURN
      ENDIF 
      IF(P(003) <  0.D0)THEN 
      IERPRD=1
      ETEXT(1)='DES SUBROUTINE: ERROR IN COMPUTATION'
      ETEXT(2)='ATTEMPT TO COMPUTE BASE**POWER WITH BASE<0.'
      RETURN
      ENDIF 
      B000008=0.D0 
      IF(P(003) == 0.D0)THEN 
      B000008=1.D0 
      ENDIF 
      B000009=1.D0-B000008 
      B000010=P(003)+B000008 
      B000011=B000010**THETA(005) 
      B000012=B000009*B000011 
      IF(CONC == 0.D0.AND.THETA(005) <= 0.D0)THEN 
      IERPRD=1
      ETEXT(1)='DES SUBROUTINE: ERROR IN COMPUTATION'
      ETEXT(2)='ATTEMPT TO COMPUTE 0**POWER WITH POWER<=0.'
      RETURN
      ENDIF 
      IF(CONC <  0.D0)THEN 
      IERPRD=1
      ETEXT(1)='DES SUBROUTINE: ERROR IN COMPUTATION'
      ETEXT(2)='ATTEMPT TO COMPUTE BASE**POWER WITH BASE<0.'
      RETURN
      ENDIF 
      B000013=0.D0 
      IF(CONC == 0.D0)THEN 
      B000013=1.D0 
      ENDIF 
      B000014=1.D0-B000013 
      B000015=CONC+B000013 
      B000016=B000015**THETA(005) 
      B000017=B000014*B000016 
      B000018=B000012+B000017 
      DADT(2)=P(002)*A(1)-P(004)*B000007/B000018 
      IF (FIRSTEM == 1) THEN ! 1
      B000001=1.D0/P(001) 
      B000002=-A(2)/P(001)/P(001) 
!                      E000005 = DERIVATIVE OF DADT(1) W.R.T. A(001)
      E000005=-P(002) 
!                      F000082 = DERIVATIVE OF DADT(1) W.R.T. P(002)
      F000082=-A(1) 
      B000019=THETA(005)-1.D0 
      B000020=B000005**B000019 
      B000021=B000020*THETA(005) 
!                      E000007 = DERIVATIVE OF B000006 W.R.T. A(002)
      E000007=B000021*B000001 
!                      E000008 = DERIVATIVE OF B000007 W.R.T. A(002)
      E000008=B000004*E000007 
      B000022=THETA(005)-1.D0 
      B000023=B000015**B000022 
      B000024=B000023*THETA(005) 
!                      E000010 = DERIVATIVE OF B000016 W.R.T. A(002)
      E000010=B000024*B000001 
!                      E000011 = DERIVATIVE OF B000017 W.R.T. A(002)
      E000011=B000014*E000010 
!                      E000013 = DERIVATIVE OF DADT(2) W.R.T. A(001)
      E000013=P(002) 
      B000025=-P(004)/B000018 
!                      E000014 = DERIVATIVE OF DADT(2) W.R.T. A(002)
      E000014=B000025*E000008 
      B000026=P(004)*B000007/B000018/B000018 
!                      E000015 = DERIVATIVE OF DADT(2) W.R.T. A(002)
      E000015=B000026*E000011+E000014 
      B000027=THETA(005)-1.D0 
      B000028=B000005**B000027 
      B000029=B000028*THETA(005) 
!                      F000084 = DERIVATIVE OF B000006 W.R.T. P(001)
      F000084=B000029*B000002 
!                      F000085 = DERIVATIVE OF B000007 W.R.T. P(001)
      F000085=B000004*F000084 
      B000030=THETA(005)-1.D0 
      B000031=B000010**B000030 
      B000032=B000031*THETA(005) 
!                      F000088 = DERIVATIVE OF B000012 W.R.T. P(003)
      F000088=B000009*B000032 
      B000033=THETA(005)-1.D0 
      B000034=B000015**B000033 
      B000035=B000034*THETA(005) 
!                      F000090 = DERIVATIVE OF B000016 W.R.T. P(001)
      F000090=B000035*B000002 
!                      F000091 = DERIVATIVE OF B000017 W.R.T. P(001)
      F000091=B000014*F000090 
!                      F000094 = DERIVATIVE OF DADT(2) W.R.T. P(002)
      F000094=A(1) 
      B000036=-B000007/B000018 
!                      F000095 = DERIVATIVE OF DADT(2) W.R.T. P(004)
      F000095=B000036 
      B000037=-P(004)/B000018 
!                      F000096 = DERIVATIVE OF DADT(2) W.R.T. P(001)
      F000096=B000037*F000085 
      B000038=P(004)*B000007/B000018/B000018 
!                      F000097 = DERIVATIVE OF DADT(2) W.R.T. P(001)
      F000097=B000038*F000091+F000096 
!                      F000098 = DERIVATIVE OF DADT(2) W.R.T. P(003)
      F000098=B000038*F000088 
      ENDIF !1
      IF (FIRSTEM == 1) THEN !2
      DA(   1,1)=E000005
      DA(   2,1)=E000013
      DA(   3,1)=E000015
      DP(   1,1)=F000082
      DP(   2,1)=F000097
      DP(   3,1)=F000094
      DP(   4,1)=F000098
      DP(   5,1)=F000095
      ENDIF !2
      RETURN
      END
      SUBROUTINE FSIZESR(NAME_FSIZES,F_SIZES)
      USE SIZES, ONLY: ISIZE
      INTEGER(KIND=ISIZE), DIMENSION(*) :: F_SIZES
      CHARACTER(LEN=*),    DIMENSION(*) :: NAME_FSIZES
      NAME_FSIZES(01)='LTH'; F_SIZES(01)=5
      NAME_FSIZES(02)='LVR'; F_SIZES(02)=5
      NAME_FSIZES(03)='LVR2'; F_SIZES(03)=0
      NAME_FSIZES(04)='LPAR'; F_SIZES(04)=14
      NAME_FSIZES(05)='LPAR3'; F_SIZES(05)=0
      NAME_FSIZES(06)='NO'; F_SIZES(06)=0
      NAME_FSIZES(07)='MMX'; F_SIZES(07)=1
      NAME_FSIZES(08)='LNP4'; F_SIZES(08)=0
      NAME_FSIZES(09)='LSUPP'; F_SIZES(09)=1
      NAME_FSIZES(10)='LIM7'; F_SIZES(10)=0
      NAME_FSIZES(11)='LWS3'; F_SIZES(11)=0
      NAME_FSIZES(12)='MAXIDS'; F_SIZES(12)=60
      NAME_FSIZES(13)='LIM1'; F_SIZES(13)=0
      NAME_FSIZES(14)='LIM2'; F_SIZES(14)=0
      NAME_FSIZES(15)='LIM3'; F_SIZES(15)=0
      NAME_FSIZES(16)='LIM4'; F_SIZES(16)=0
      NAME_FSIZES(17)='LIM5'; F_SIZES(17)=0
      NAME_FSIZES(18)='LIM6'; F_SIZES(18)=0
      NAME_FSIZES(19)='LIM8'; F_SIZES(19)=0
      NAME_FSIZES(20)='LIM10'; F_SIZES(20)=0
      NAME_FSIZES(21)='LIM11'; F_SIZES(21)=0
      NAME_FSIZES(22)='LIM13'; F_SIZES(22)=0
      NAME_FSIZES(23)='LIM15'; F_SIZES(23)=0
      NAME_FSIZES(24)='LIM16'; F_SIZES(24)=0
      NAME_FSIZES(25)='MAXRECID'; F_SIZES(25)=0
      NAME_FSIZES(26)='PC'; F_SIZES(26)=0
      NAME_FSIZES(27)='PCT'; F_SIZES(27)=1
      NAME_FSIZES(28)='PIR'; F_SIZES(28)=6
      NAME_FSIZES(29)='PD'; F_SIZES(29)=7
      NAME_FSIZES(30)='PAL'; F_SIZES(30)=0
      NAME_FSIZES(31)='MAXFCN'; F_SIZES(31)=0
      NAME_FSIZES(32)='MAXIC'; F_SIZES(32)=0
      NAME_FSIZES(33)='PG'; F_SIZES(33)=0
      NAME_FSIZES(34)='NPOPMIXMAX'; F_SIZES(34)=0
      NAME_FSIZES(35)='MAXOMEG'; F_SIZES(35)=3
      NAME_FSIZES(36)='MAXPTHETA'; F_SIZES(36)=8
      NAME_FSIZES(37)='MAXITER'; F_SIZES(37)=20
      NAME_FSIZES(38)='ISAMPLEMAX'; F_SIZES(38)=0
      NAME_FSIZES(39)='DIMTMP'; F_SIZES(39)=0
      NAME_FSIZES(40)='DIMCNS'; F_SIZES(40)=0
      NAME_FSIZES(41)='DIMNEW'; F_SIZES(41)=0
      NAME_FSIZES(42)='PDT'; F_SIZES(42)=4
      NAME_FSIZES(43)='LADD_MAX'; F_SIZES(43)=0
      NAME_FSIZES(44)='MAXSIDL'; F_SIZES(44)=0
      NAME_FSIZES(45)='NTT'; F_SIZES(45)=5
      NAME_FSIZES(46)='NOMEG'; F_SIZES(46)=3
      NAME_FSIZES(47)='NSIGM'; F_SIZES(47)=2
      NAME_FSIZES(48)='PPDT'; F_SIZES(48)=1
      RETURN
      END SUBROUTINE FSIZESR
      SUBROUTINE MUMODEL2(THETA,MU_,ICALL,IDEF,NEWIND,&
      EVTREC,DATREC,IREV,NVNT,INDXS,F,G,H,IRGG,GG,NETAS)
      USE NMPRD4P
      USE SIZES,     ONLY: DPSIZE,ISIZE
      USE PRDIMS,    ONLY: GPRD,HPRD,GERD,HERD,GPKD
      USE NMBAYES_REAL,    ONLY: PRIORINFO
      USE NMPRD_REAL,ONLY: ETA,EPS
      USE NMPRD_INT, ONLY: MSEC=>ISECDER,MFIRST=>IFRSTDER,COMACT,COMSAV,IFIRSTEM
      USE NMPRD_INT, ONLY: MDVRES,ETASXI,NPDE_MODE,NOFIRSTDERCODE
      USE NMPRD_REAL, ONLY: DV_LOQ,CDF_L,DV_LAQ,CDF_LA
      USE NMPRD_INT, ONLY: IQUIT
      USE NMBAYES_REAL, ONLY: LDF
      IMPLICIT REAL(KIND=DPSIZE) (A-Z)
      REAL(KIND=DPSIZE)   :: MU_(*)
      INTEGER NEWIND
      REAL(KIND=DPSIZE) :: EVTREC
      SAVE
      INTEGER(KIND=ISIZE) :: FIRSTEM
      INTEGER(KIND=ISIZE) :: ICALL,IDEF,IREV,NVNT,INDXS,IRGG,NETAS
      DIMENSION :: IDEF(7,*),THETA(*),EVTREC(IREV,*),INDXS(*),GG(IRGG,GPKD+1,*)
      RETURN
      END
