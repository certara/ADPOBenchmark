      MODULE NMPRD4P
      USE SIZES, ONLY: DPSIZE
      USE NMPRD4,ONLY: VRBL
      IMPLICIT NONE
      SAVE
      REAL(KIND=DPSIZE), DIMENSION (:),POINTER ::COM
      REAL(KIND=DPSIZE), POINTER ::CWT,TVVMAX,VMAX,TVKM,KM,TVV2,V2
      REAL(KIND=DPSIZE), POINTER ::TVKA,KA,SC,IPRED,Y,CONC,A000032
      REAL(KIND=DPSIZE), POINTER ::A000034,A000036,A000038,A000039
      REAL(KIND=DPSIZE), POINTER ::A000040,D000001,D000002,D000003
      REAL(KIND=DPSIZE), POINTER ::D000004,D000005,D000042,D000045
      REAL(KIND=DPSIZE), POINTER ::D000044,D000041,D000043,C000033
      REAL(KIND=DPSIZE), POINTER ::C000032,D000047,D000046,A000041
      REAL(KIND=DPSIZE), POINTER ::A000042,A000043,A000044,E000007
      REAL(KIND=DPSIZE), POINTER ::F000085,E000010,E000012,F000091
      REAL(KIND=DPSIZE), POINTER ::F000088,F000092,F000089
      CONTAINS
      SUBROUTINE ASSOCNMPRD4
      COM=>VRBL
      CWT=>COM(000001);TVVMAX=>COM(000002);VMAX=>COM(000003)
      TVKM=>COM(000004);KM=>COM(000005);TVV2=>COM(000006)
      V2=>COM(000007);TVKA=>COM(000008);KA=>COM(000009)
      SC=>COM(000010);IPRED=>COM(000011);Y=>COM(000012)
      CONC=>COM(000013);A000032=>COM(000014);A000034=>COM(000015)
      A000036=>COM(000016);A000038=>COM(000017);A000039=>COM(000018)
      A000040=>COM(000019);D000001=>COM(000020);D000002=>COM(000021)
      D000003=>COM(000022);D000004=>COM(000023);D000005=>COM(000024)
      D000042=>COM(000025);D000045=>COM(000026);D000044=>COM(000027)
      D000041=>COM(000028);D000043=>COM(000029);C000033=>COM(000030)
      C000032=>COM(000031);D000047=>COM(000032);D000046=>COM(000033)
      A000041=>COM(000034);A000042=>COM(000035);A000043=>COM(000036)
      A000044=>COM(000037);E000007=>COM(000038);F000085=>COM(000039)
      E000010=>COM(000040);E000012=>COM(000041);F000091=>COM(000042)
      F000088=>COM(000043);F000092=>COM(000044);F000089=>COM(000045)
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
      KM=TVKM 
      TVV2=THETA(003) 
      B000003=DEXP(ETA(002)) 
      V2=TVV2*B000003 
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
!                      A000036 = DERIVATIVE OF V2 W.R.T. ETA(002)
      A000036=TVV2*B000003 
!                      A000039 = DERIVATIVE OF SC W.R.T. ETA(002)
      A000039=A000036 
!                      A000041 = DERIVATIVE OF P000001 W.R.T. ETA(002)
      A000041=A000039 
!                      A000043 = DERIVATIVE OF P000004 W.R.T. ETA(001)
      A000043=A000032 
      GG(0001,1,1)=P000001
      GG(0001,0003,1)=A000041
      GG(0002,1,1)=P000002
      GG(0003,1,1)=P000003
      GG(0004,1,1)=P000004
      GG(0004,0002,1)=A000043
      GG(0005,1,1)=SC
      GG(0005,0003,1)=A000039
      ELSE
      GG(0001,1,1)=P000001
      GG(0002,1,1)=P000002
      GG(0003,1,1)=P000003
      GG(0004,1,1)=P000004
      GG(0005,1,1)=SC
      ENDIF
      IF (MSEC == 1) THEN
!                      A000034 = DERIVATIVE OF A000032 W.R.T. ETA(001)
      A000034=TVVMAX*B000001 
!                      A000038 = DERIVATIVE OF A000036 W.R.T. ETA(002)
      A000038=TVV2*B000003 
!                      A000040 = DERIVATIVE OF A000039 W.R.T. ETA(002)
      A000040=A000038 
!                      A000042 = DERIVATIVE OF A000041 W.R.T. ETA(002)
      A000042=A000040 
!                      A000044 = DERIVATIVE OF A000043 W.R.T. ETA(001)
      A000044=A000034 
      GG(0001,0003,0003)=A000042
      GG(0004,0002,0002)=A000044
      GG(0005,0003,0003)=A000040
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
 !  level            0
      IPRED=F 
      B000001=DEXP(EPS(002)) 
      Y=F*B000001+EPS(001) 
!                      C000032 = DERIVATIVE OF Y W.R.T. EPS(002)
      C000032=F*B000001 
!                      C000033 = DERIVATIVE OF Y W.R.T. EPS(001)
      C000033=1.D0 
      IF (FIRSTEM == 1) THEN !1
!                      D000041 = DERIVATIVE OF Y W.R.T. ETA(002)
      D000041=B000001*D000002 
!                      D000042 = DERIVATIVE OF Y W.R.T. ETA(001)
      D000042=B000001*D000001 
!                      D000046 = DERIVATIVE OF C000032 W.R.T. ETA(002)
      D000046=B000001*D000002 
!                      D000047 = DERIVATIVE OF C000032 W.R.T. ETA(001)
      D000047=B000001*D000001 
      G(001,1)=D000042
      G(002,1)=D000041
      ENDIF !1
      HH(001,1)=C000033
      HH(002,1)=C000032
      IF (FIRSTEM == 1) THEN !2
      HH(002,002)=D000047
      HH(002,003)=D000046
      ENDIF !2
      IF (MSEC == 1) THEN
      D000003=G(001,002)
      D000004=G(002,002)
      D000005=G(002,003)
!                      D000043 = DERIVATIVE OF D000041 W.R.T. ETA(002)
      D000043=B000001*D000005 
!                      D000044 = DERIVATIVE OF D000042 W.R.T. ETA(002)
      D000044=B000001*D000004 
!                      D000045 = DERIVATIVE OF D000042 W.R.T. ETA(001)
      D000045=B000001*D000003 
      G(001,002)=D000045
      G(002,002)=D000044
      G(002,003)=D000043
      ENDIF
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
      IDEFD(1)=  0
      IDEFD(2)=0
      DA(   1,1)=0000014280
      DA(   2,1)=0000028441
      DA(   3,1)=0000028560
      DA(   4,1)=0000028562
      DA(   5,1)=0000014312
      DA(   6,1)=0000028473
      DA(   7,1)=0000028591
      DA(   8,1)=0000028593
      DA(   9,1)=0000028594
      DA(  10,1)=0000000000
      DP(   1,1)=0000014399
      DP(   2,1)=0000028441
      DP(   3,1)=0000028560
      DP(   4,1)=0000028679
      DP(   5,1)=0000028798
      DP(   6,1)=0000028442
      DP(   7,1)=0000028680
      DP(   8,1)=0000028799
      DP(   9,1)=0000028444
      DP(  10,1)=0000028682
      DP(  11,1)=0000028801
      DP(  12,1)=0000028445
      DP(  13,1)=0000028683
      DP(  14,1)=0000000000
      DT(   1)=0000000000
      RETURN
      ENDIF
 !  level            0
 !  level            0
      CONC=A(2)/P(001) 
      DADT(1)=-P(002)*A(1) 
      B000006=P(003)+CONC 
      DADT(2)=P(002)*A(1)-P(004)*CONC/B000006 
      IF (FIRSTEM == 1) THEN ! 1
      B000001=1.D0/P(001) 
      B000002=-A(2)/P(001)/P(001) 
!                      E000007 = DERIVATIVE OF DADT(1) W.R.T. A(001)
      E000007=-P(002) 
!                      F000085 = DERIVATIVE OF DADT(1) W.R.T. P(002)
      F000085=-A(1) 
!                      E000010 = DERIVATIVE OF DADT(2) W.R.T. A(001)
      E000010=P(002) 
      B000007=-P(004)/B000006 
!                      E000011 = DERIVATIVE OF DADT(2) W.R.T. A(002)
      E000011=B000007*B000001 
      B000008=P(004)*CONC/B000006/B000006 
!                      E000012 = DERIVATIVE OF DADT(2) W.R.T. A(002)
      E000012=B000008*B000001+E000011 
!                      F000088 = DERIVATIVE OF DADT(2) W.R.T. P(002)
      F000088=A(1) 
      B000014=-CONC/B000006 
!                      F000089 = DERIVATIVE OF DADT(2) W.R.T. P(004)
      F000089=B000014 
      B000015=-P(004)/B000006 
!                      F000090 = DERIVATIVE OF DADT(2) W.R.T. P(001)
      F000090=B000015*B000002 
      B000016=P(004)*CONC/B000006/B000006 
!                      F000091 = DERIVATIVE OF DADT(2) W.R.T. P(001)
      F000091=B000016*B000002+F000090 
!                      F000092 = DERIVATIVE OF DADT(2) W.R.T. P(003)
      F000092=B000016 
      ENDIF !1
      IF (MSEC == 1) THEN 
      B000003=A(2)/P(001)/P(001)/P(001) 
      B000004=A(2)/P(001)/P(001)/P(001) 
!                      F000083 = DERIVATIVE OF B000002 W.R.T. P(001)
      F000083=B000004+B000003 
!                      F000084 = DERIVATIVE OF F000081 W.R.T. P(001)
      F000084=F000083 
      B000005=-1.D0/P(001)/P(001) 
!                      E000006 = DERIVATIVE OF F000081 W.R.T. A(002)
      E000006=B000005 
!                      E000008 = DERIVATIVE OF F000085 W.R.T. A(001)
      E000008=-1.D0 
      B000010=P(004)/B000006/B000006 
!                      E000014 = DERIVATIVE OF B000007 W.R.T. A(002)
      E000014=B000010*B000001 
!                      E000015 = DERIVATIVE OF E000011 W.R.T. A(002)
      E000015=B000001*E000014 
      B000011=P(004)/B000006/B000006 
!                      E000016 = DERIVATIVE OF B000008 W.R.T. A(002)
      E000016=B000011*B000001 
      B000012=-P(004)*CONC/B000006/B000006/B000006 
!                      E000017 = DERIVATIVE OF B000008 W.R.T. A(002)
      E000017=B000012*B000001+E000016 
      B000013=-P(004)*CONC/B000006/B000006/B000006 
!                      E000018 = DERIVATIVE OF B000008 W.R.T. A(002)
      E000018=B000013*B000001+E000017 
!                      E000019 = DERIVATIVE OF E000012 W.R.T. A(002)
      E000019=B000001*E000018 
!                      E000020 = DERIVATIVE OF E000012 W.R.T. A(002)
      E000020=E000015+E000019 
!                      F000094 = DERIVATIVE OF F000087 W.R.T. P(001)
      F000094=F000084 
      B000018=-1.D0/B000006 
!                      F000095 = DERIVATIVE OF B000014 W.R.T. P(001)
      F000095=B000018*B000002 
      B000019=CONC/B000006/B000006 
!                      F000096 = DERIVATIVE OF B000014 W.R.T. P(001)
      F000096=B000019*B000002+F000095 
!                      F000098 = DERIVATIVE OF F000089 W.R.T. P(003)
      F000098=B000019 
!                      F000099 = DERIVATIVE OF F000089 W.R.T. P(001)
      F000099=F000096 
      B000020=-1.D0/B000006 
      B000021=P(004)/B000006/B000006 
!                      F000101 = DERIVATIVE OF B000015 W.R.T. P(001)
      F000101=B000021*B000002 
!                      F000103 = DERIVATIVE OF F000090 W.R.T. P(003)
      F000103=B000002*B000021 
!                      F000104 = DERIVATIVE OF F000090 W.R.T. P(001)
      F000104=B000002*F000101 
!                      F000105 = DERIVATIVE OF F000090 W.R.T. P(004)
      F000105=B000002*B000020 
!                      F000106 = DERIVATIVE OF F000090 W.R.T. P(001)
      F000106=B000015*F000084+F000104 
      B000022=CONC/B000006/B000006 
      B000023=P(004)/B000006/B000006 
!                      F000108 = DERIVATIVE OF B000016 W.R.T. P(001)
      F000108=B000023*B000002 
      B000024=-P(004)*CONC/B000006/B000006/B000006 
!                      F000109 = DERIVATIVE OF B000016 W.R.T. P(001)
      F000109=B000024*B000002+F000108 
      B000025=-P(004)*CONC/B000006/B000006/B000006 
!                      F000111 = DERIVATIVE OF B000016 W.R.T. P(001)
      F000111=B000025*B000002+F000109 
!                      F000112 = DERIVATIVE OF B000016 W.R.T. P(003)
      F000112=B000025+B000024 
!                      F000113 = DERIVATIVE OF F000091 W.R.T. P(003)
      F000113=B000002*F000112 
!                      F000114 = DERIVATIVE OF F000091 W.R.T. P(001)
      F000114=B000002*F000111 
!                      F000115 = DERIVATIVE OF F000091 W.R.T. P(004)
      F000115=B000002*B000022 
!                      F000116 = DERIVATIVE OF F000091 W.R.T. P(001)
      F000116=B000016*F000094+F000114 
!                      F000117 = DERIVATIVE OF F000091 W.R.T. P(001)
      F000117=F000106+F000116 
!                      F000118 = DERIVATIVE OF F000091 W.R.T. P(004)
      F000118=F000105+F000115 
!                      F000119 = DERIVATIVE OF F000091 W.R.T. P(003)
      F000119=F000103+F000113 
!                      F000120 = DERIVATIVE OF F000092 W.R.T. P(003)
      F000120=F000112 
!                      F000121 = DERIVATIVE OF F000092 W.R.T. P(001)
      F000121=F000111 
!                      F000122 = DERIVATIVE OF F000092 W.R.T. P(004)
      F000122=B000022 
!                      F000124 = DERIVATIVE OF F000093 W.R.T. P(001)
      F000124=F000084 
!                      E000022 = DERIVATIVE OF F000087 W.R.T. A(002)
      E000022=E000006 
!                      E000023 = DERIVATIVE OF F000088 W.R.T. A(001)
      E000023=1.D0 
      B000026=CONC/B000006/B000006 
!                      E000024 = DERIVATIVE OF B000014 W.R.T. A(002)
      E000024=B000026*B000001 
!                      E000025 = DERIVATIVE OF F000089 W.R.T. A(002)
      E000025=E000024 
      B000027=P(004)/B000006/B000006 
!                      E000026 = DERIVATIVE OF B000015 W.R.T. A(002)
      E000026=B000027*B000001 
!                      E000027 = DERIVATIVE OF F000090 W.R.T. A(002)
      E000027=B000002*E000026 
!                      E000028 = DERIVATIVE OF F000090 W.R.T. A(002)
      E000028=B000015*E000006+E000027 
      B000028=-P(004)*CONC/B000006/B000006/B000006 
!                      E000029 = DERIVATIVE OF B000016 W.R.T. A(002)
      E000029=B000028*B000001 
      B000029=-P(004)*CONC/B000006/B000006/B000006 
!                      E000030 = DERIVATIVE OF B000016 W.R.T. A(002)
      E000030=B000029*B000001+E000029 
!                      E000031 = DERIVATIVE OF F000091 W.R.T. A(002)
      E000031=B000002*E000030 
!                      E000032 = DERIVATIVE OF F000091 W.R.T. A(002)
      E000032=B000016*E000022+E000031 
!                      E000033 = DERIVATIVE OF F000091 W.R.T. A(002)
      E000033=E000028+E000032 
!                      E000034 = DERIVATIVE OF F000092 W.R.T. A(002)
      E000034=E000030 
!                      E000036 = DERIVATIVE OF F000093 W.R.T. A(002)
      E000036=E000006 
      ENDIF !msec
      IF (FIRSTEM == 1) THEN !2
      DA(   1,1)=E000007
      DA(   2,1)=E000010
      DA(   3,1)=E000012
      DP(   1,1)=F000085
      DP(   2,1)=F000091
      DP(   3,1)=F000088
      DP(   4,1)=F000092
      DP(   5,1)=F000089
      ENDIF !2
      IF (MSEC == 1) THEN
      DA(   4,1)=E000020
      DA(   5,1)=E000008
      DA(   6,1)=E000023
      DA(   7,1)=E000033
      DA(   8,1)=E000034
      DA(   9,1)=E000025
      DP(   6,1)=F000117
      DP(   7,1)=F000121
      DP(   8,1)=F000099
      DP(   9,1)=F000119
      DP(  10,1)=F000120
      DP(  11,1)=F000098
      DP(  12,1)=F000118
      DP(  13,1)=F000122
      ENDIF
      RETURN
      END
      SUBROUTINE FSIZESR(NAME_FSIZES,F_SIZES)
      USE SIZES, ONLY: ISIZE
      INTEGER(KIND=ISIZE), DIMENSION(*) :: F_SIZES
      CHARACTER(LEN=*),    DIMENSION(*) :: NAME_FSIZES
      NAME_FSIZES(01)='LTH'; F_SIZES(01)=4
      NAME_FSIZES(02)='LVR'; F_SIZES(02)=4
      NAME_FSIZES(03)='LVR2'; F_SIZES(03)=0
      NAME_FSIZES(04)='LPAR'; F_SIZES(04)=10
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
      NAME_FSIZES(28)='PIR'; F_SIZES(28)=14
      NAME_FSIZES(29)='PD'; F_SIZES(29)=7
      NAME_FSIZES(30)='PAL'; F_SIZES(30)=0
      NAME_FSIZES(31)='MAXFCN'; F_SIZES(31)=0
      NAME_FSIZES(32)='MAXIC'; F_SIZES(32)=0
      NAME_FSIZES(33)='PG'; F_SIZES(33)=0
      NAME_FSIZES(34)='NPOPMIXMAX'; F_SIZES(34)=0
      NAME_FSIZES(35)='MAXOMEG'; F_SIZES(35)=2
      NAME_FSIZES(36)='MAXPTHETA'; F_SIZES(36)=7
      NAME_FSIZES(37)='MAXITER'; F_SIZES(37)=20
      NAME_FSIZES(38)='ISAMPLEMAX'; F_SIZES(38)=0
      NAME_FSIZES(39)='DIMTMP'; F_SIZES(39)=0
      NAME_FSIZES(40)='DIMCNS'; F_SIZES(40)=0
      NAME_FSIZES(41)='DIMNEW'; F_SIZES(41)=0
      NAME_FSIZES(42)='PDT'; F_SIZES(42)=4
      NAME_FSIZES(43)='LADD_MAX'; F_SIZES(43)=0
      NAME_FSIZES(44)='MAXSIDL'; F_SIZES(44)=0
      NAME_FSIZES(45)='NTT'; F_SIZES(45)=4
      NAME_FSIZES(46)='NOMEG'; F_SIZES(46)=2
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