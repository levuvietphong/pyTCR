        PROGRAM WRTFAST
C
C  ***  This program selects random genesis locations and generates random  ***
C  ***    tracks from monthly 850 and 250 mb mean flow and variance         ***
C  ***       It then runs a toy model to produce an intensity estimate      ***
C  ***       It can run off of models or NCEP reanalysis data, and can      ***
C  ***       alternatively draw from best track genesis pdfs.               ***
C
C  ***               Created 12/2016  by K. Emanuel                         ***
C  ***                Copyright WindRiskTech, 2016                          ***
C  ***                     Revised April, 2022                              ***
C
C  ***   Modified May 2015 to ingest monthly mixed layer depths and gammas  ***
C  ***               on 181 x 360 grid, if they exist                       ***
C
        USE, intrinsic :: ieee_arithmetic
C
        PARAMETER(NTG=600)
C
c
C  ***      Flow-related files  ***
C
        REAL*8, ALLOCATABLE :: U850(:,:,:), U250(:,:,:),V850(:,:,:)
        REAL*8, ALLOCATABLE :: V250(:,:,:)
        REAL*8, ALLOCATABLE :: UVAR250(:,:,:), VVAR250(:,:,:)
        REAL*8, ALLOCATABLE :: UVAR850(:,:,:),VVAR850(:,:,:)
        REAL*8, ALLOCATABLE :: U250V250(:,:,:),U250U850(:,:,:)
        REAL*8, ALLOCATABLE :: V250V850(:,:,:), U850V850(:,:,:)
        REAL*8, ALLOCATABLE :: U250V850(:,:,:), V250U850(:,:,:)
        REAL*8, ALLOCATABLE :: VOR(:,:,:)
        REAL*8, ALLOCATABLE :: VTEM(:,:,:), RH(:,:,:), T600(:,:,:)
        REAL*8, ALLOCATABLE :: VPOTC(:,:,:)
        REAL, ALLOCATABLE :: X1VERT(:), Y1VERT(:)
        REAL, ALLOCATABLE :: X2VERT(:), Y2VERT(:)
c
        REAL AMP(15), PHASEU850(15), PHASEV850(15)
        REAL PHASEU250(15),PHASEV250(15)
        REAL*8 A(4,4), PCA(4), A1(4,4), PCA1(4), B1(4,4), PCB1(4)
        REAL*8 A2(4,4),B2(4,4),PCA2(4),PCB2(4)
        REAL*8 A3(4,4),B3(4,4),PCA3(4),PCB3(4)
        REAL*8 A4(4,4),B4(4,4),PCA4(4),PCB4(4)
        REAL LONGINC, LATINC, LONGSTART, LATSTART
        REAL BASINFAC(7)
C
C  ***   Print quantities   **
C
        INTEGER MONTH, DATE
C
C  ***       Other quantities    ***
C
        REAL*8 TT,DTT,TIME, ALAT, ALOG,PRINTIME,PIAMP
        REAL*8 ALATOLD
        REAL LATMAX, LOGMAX, LOGMIN, LATMING, LATMAXG
        REAL LATMAXT,LOGMAXT,LOGMINT
        INTEGER DAY(12)
        INTEGER ISEED, year, rerandom
        REAL, ALLOCATABLE :: router(:),vouter(:)
C
C  ***     Dimension ocean variables   ***
C
        REAL DISTANCE(600),UTIME(600)
        REAL HMIX0(600),BATHTRACK(600),VPTRACK(600)
        REAL LANDFRAC(600)
        REAL STRAT(600),RHTRACK(600),T600TRACK(600)
        REAL CDTRACK(600)
C
C  ***     Climatological ocean and potential intensity files  ***
C
        REAL BATH(1440,721), BATHM(1440,721)
        REAL*8 CDRAG(1440,721)
        REAL*8, ALLOCATABLE :: THERMAL(:,:,:),MIXDEPTH(:,:,:)
        REAL, ALLOCATABLE :: THERMALd(:,:,:),MIXDEPTHd(:,:,:)
        REAL LATTE(600),LOGGTE(600)
        INTEGER MONTE(600),DATETE(600),HOURTE(600),WRITETIME(600)
C
C  ***      Shear-related files  ***
C
        REAL SHEARTRACK(600), U850TR(600), V850TR(600)
C
C  ***      Dimension graphics arrays       ***
C
        REAL RMG(NTG), RTMG(NTG), LATTIME(NTG), LOGTIME(NTG)
        REAL HMG(NTG), UTRANS(NTG), VTRANS(NTG)
        REAL SSTT(NTG),  HMG0(NTG)
        REAL VMG(NTG), T600G(NTG),RHG(NTG)
        REAL VMGSE(NTG), RMGSE(NTG)
        REAL VSHG(NTG), VPOTG(NTG), U850G(NTG), V850G(NTG)
        REAL MUG(NTG), MDG(NTG), PCG(NTG,2)
        REAL*8 TIMEG(NTG)
        REAL HATM, HATMI
C
C  ***      Character arrays    ***
C
        CHARACTER*2 AYEAR,AMON,ADAY,AHOUR,BASIN,ISTNUM,DATEC,MOND
        CHARACTER*2 ABAS,ABASU,ASNUM,ABAS2,SHOUR
        CHARACTER*3 MONC,DATETEMP
        CHARACTER*1 LATN,LOGGB,ILATI,ILONGI,AK
        CHARACTER*80 JFILE,FNAMEA,FNAMEB,PATHBUF
        CHARACTER*95 NAM
        CHARACTER*30 NAM2
        CHARACTER*20 ANAM,MODEL,RUN,FLAV
        CHARACTER*11 SNAME
        CHARACTER*15 SDATE
        CHARACTER*80 DIREC3
        CHARACTER*120 DIRECTORY
        CHARACTER*4 IDENT, FCST, FN2, FCSTO,AYEAR4
        CHARACTER*4 FCST2, gmeth, shape
        CHARACTER*3 lev
        CHARACTER*50 SFILENAME
        CHARACTER*70 polyfile
        CHARACTER GENPOINTS*8
        CHARACTER AF1*5, CITYNAME*90, AMON1*1, BAS*2
        CHARACTER ch1*1,ch2*2,ch3*3,ch4*4,ch5*5
C
	character(len=:), allocatable :: filenum
	character*4 ng
C
C  ***   Models parameter arrays  ***
C
        Real latincr(200),longincr(200),latstartr(200),piampr(200)
        Real shearfacr(200),basinfacr(200,7),longstartr(200)
        Integer nlor(200),nlar(200)
        Character*8 mname(200), mdirec(200), MDIR
C
C  ***       Other quantities    ***
C
        REAL H_A, MEQ, MF, MT, MUMAX, MDMIN
        REAL*8 DT,DT1,DTG,TT1,DT0
        integer*8 adate(3)
        integer restart
        logical ex
        real xtemp(2), ytemp(2)
c
        read(*,*)idum   ! File number, read in from shell script
        if(idum.lt.10)then
         allocate(character(1) :: filenum)
         ng="(I1)"
        else
         allocate(character(2) :: filenum)
         ng="(I2)"
        end if
        write(filenum,ng)idum
c
        data day/31,28,31,30,31,30,31,31,30,31,30,31/
C
C  ***   Read in city parameters  ***
C
        OPEN(UNIT=11,FILE='wrtfast'//filenum//'.in',STATUS='OLD')
        READ(11,*)BAS,MODEL,FLAV,RUN,gmeth,shape,CITYNAME,
     1   CLAT,CLOG,CRAD,NUMTRACKS,VMCRIT,restart,rerandom,adate(1),
     2   polyfile
        CLOSE(11)
        ADATE(2)=1
        ADATE(3)=1                
C
C  ***  Value of ocean data latitude points
C
        NLAOC=181
        IF(MODEL(1:2).EQ.'dk')NLAOC=180
c
C  VMCRIT is the critical maximum wind (kts) within city radius below which track is ignored
C
        IF(CLOG.LT.0.0)CLOG=CLOG+360.0
c
        clog1=clog
        clat1=clat
        jint=0
c
        IBAS=6
        IF(BAS.EQ.'AL')IBAS=1
        IF(BAS.EQ.'EP')IBAS=2
        IF(BAS.EQ.'WP')IBAS=3
        IF(BAS.EQ.'IO')IBAS=4
        IF(BAS.EQ.'SH')IBAS=5
        IF(BAS.EQ.'MT')IBAS=7
c
c
c  ***    The follwing are attributes of particular models....also see PARAMETER statment   ***
c
        DIREC3='data/'//trim(model)//'/'//trim(flav)//'/'//trim(run)
C
c       Grid and other specifications
c
c  ***  The parameter PIAMP is adjusted to bring the 20th century potential intensity more in line with obs  ***
c  ***      Default values  ***
c
        lev='250'
        shearcoef=1.0
        piamp=1.0
c
        DO I=1,7
          BASINFAC(I)=1.0
        END DO
c
c  ***   Read in Model Parameters   ***
c
        open(unit=11,file='Models.txt',status='old')
        read(11,*)
        read(11,*)
        nmods=0
        do i=1,100
          read(11,*,end=9)mname(i),mdirec(i),nlor(i),nlar(i),
     1     longincr(i),latincr(i),longstartr(i),latstartr(i),
     2     piampr(i),shearfacr(i),(basinfacr(i,j),j=1,7)
           nmods=nmods+1
        end do
    9   continue
        close(11)
c
        nstate=0
        do i=1,nmods
          if(trim(model).eq.trim(mname(i)))then
              nstate=1
              imod=i
          endif
        end do
        if(nstate.eq.0)then
           print*, 'Model not found'
           stop
        endif
c
        MDIR=mdirec(imod)
        NLO=nlor(imod)
        NLA=nlar(imod)
        LONGINC=longincr(imod)
        LATINC=latincr(imod)
        LONGSTART=longstartr(imod)
        LATSTART=latstartr(imod)
        PIAMP=piampr(imod)
        shearcoef=shearfacr(imod)
        basinfac(1:7)=basinfacr(imod,1:7)
c
        IF(TRIM(MODEL).EQ.'ccsm3')THEN
           lev='200'
        ENDIF
C
        IF(TRIM(MDIR).EQ.'archive')THEN
          DIREC3='/archive1/emanuel/riskproject/'//
     1     TRIM(DIREC3)
        ELSEIF(TRIM(MDIR).EQ.'wrt2')THEN
          DIREC3='/opt/wrt2/'//TRIM(DIREC3)
        ENDIF
c
C
C  ***  Value of ocean data dimensions
C
        NLAOC=181
        NLOOC=360
        IF(MODEL(1:2).EQ.'dk')THEN
           NLAOC=NLA
           NLOOC=NLO
        END IF
C
        ALLOCATE(U850(NLO,NLA,12),U250(NLO,NLA,12),V850(NLO,NLA,12))
        ALLOCATE(V250(NLO,NLA,12),UVAR250(NLO,NLA,12))
        ALLOCATE(UVAR850(NLO,NLA,12),VVAR250(NLO,NLA,12))
        ALLOCATE(VVAR850(NLO,NLA,12),U250V250(NLO,NLA,12))
        ALLOCATE(U250U850(NLO,NLA,12),V250V850(NLO,NLA,12))
        ALLOCATE(U850V850(NLO,NLA,12),VOR(NLO,NLA,12))
        ALLOCATE(U250V850(NLO,NLA,12),V250U850(NLO,NLA,12))
        ALLOCATE(VTEM(NLO,NLA,12),RH(NLO,NLA,12),VPOTC(NLO,NLA,12))
        ALLOCATE(T600(NLO,NLA,12))
        ALLOCATE(THERMAL(12,NLAOC,NLOOC),MIXDEPTH(12,NLAOC,NLOOC))
        ALLOCATE(THERMALd(12,NLAOC,NLOOC),MIXDEPTHd(12,NLAOC,NLOOC))
c
        CITYNAME=TRIM(MODEL)//'/'//TRIM(FLAV)//'/'//TRIM(RUN)//'/'//
     1   TRIM(CITYNAME)//'/'//TRIM(BAS)
c
c  ***   Open output file  ***
c
         OPEN(UNIT=15,FILE=TRIM(CITYNAME)//'/latlongs.out',
     1    status='unknown')
c
        mpoly=0
c
        if(shape.eq.'poly')then
c
         open(unit=11,file='polyfiles/'//trim(polyfile),status='old')
         nvert=0
 208     continue
         read(11,*,end=209)xdum,ydum,x2dum,y2dum
         nvert=nvert+1
         goto 208
 209     continue
         allocate(x1vert(nvert),y1vert(nvert),x2vert(nvert),
     1     y2vert(nvert))
         rewind(11)
         do i=1,nvert
          read(11,*)x1vert(i),y1vert(i),x2vert(i),y2vert(i)
          if(x1vert(i).lt.0.0)x1vert(i)=x1vert(i)+360.0
          if(x2vert(i).lt.0.0)x2vert(i)=x2vert(i)+360.0
         end do
         close(11)
c
        if(abs(x2vert(nvert)-x1vert(1)).lt.0.01.and.
     1   abs(y2vert(nvert)-y1vert(1)).lt.0.01)then
          mpoly=1
        end if
        end if
c
c  ***  Some constants ***
c
        RM=70.0  ! Initial radius of maximum winds (km)
        DTG=2.0D0 ! Graphics output interval (hours)
        DTG=DTG*3600.
C
C  ***   Beta drift   ***
C
        UDRIFT_nocov=-0.4
        UDRIFT_cov=-0.9
        VDRIFT_nocov=2.0
        VDRIFT_cov=1.4  !  Changed from 1.0 to 1.4 8/29/20
C
C  ***    Time scale in Fourier series wind fields, and weight given             ***
C  ***      to 250 and 850 hPa flow in setting storm motion                      ***
C 
        TAUMAX=15.0
        WEIGHT=0.8
c
c  ***    Global frequency factor (empirical, matched to present climate)       ***
c  ***             Used only for random seeding method                          ***
c
        globalfac=38.0   !   Changed from 55.0 to 38.0 on 12/14/2016
C
C  ***  Latitude and longitude boundaries on tracks (only ABS(LATMAX)) used      ***
C
c
c   The following bounds for genesis locations are used in all basins except global and Mediterranean
c
        LATMING=3.0
        LATMAXG=45.0
c
c  ***      Definitions of basin bounds and frequency factors for climo. These have been             ***
c  ***  calibrated to give the observed frequencies averaged between 1970 and 2005 in all basins     ***
c  ***          for all storms whose maximum winds exceed 40 knots. Driven by NCEP reanalysis.       ***
c
        IF(BAS.EQ.'AL')THEN
         LATMAX=40.0
         LOGMAX=355.0
         LOGMIN=260.0
         LATMAXT=65.0
         LOGMAXT=380.0
         LOGMINT=240.0
         climfac=46.23 
        ELSE IF(BAS.EQ.'EP')THEN
         LATMAX=35.0
         LOGMAX=285.0
         LOGMIN=180.0
         LATMAXT=60.0
         LOGMAXT=300.0
         LOGMINT=160.0
         climfac=65.93 
        ELSE IF(BAS.EQ.'AP')THEN
         LATMAX=40.0
         LOGMAX=355.0
         LOGMIN=180.0
         LATMAXT=65.0
         LOGMAXT=380.0
         LOGMINT=160.0
         climfac=57.00
        ELSE IF(BAS.EQ.'WP')THEN
         LATMAX=55.0
         LOGMAX=220.0
         LOGMIN=100.0 
         LATMAXT=65.0
         LOGMAXT=260.0
         LOGMINT=50.0
         climfac=53.5 
        ELSE IF(BAS.EQ.'IO')THEN
         LATMAX=32.0
         LOGMAX=100.0
         LOGMIN=40.0 
         LATMAXT=40.0
         LOGMAXT=150.0
         LOGMINT=20.0
         climfac=15.08
        ELSE IF(BAS.EQ.'SH'.OR.BAS.EQ.'BM')THEN
         LATMAX=50.0
         LOGMAX=290.0
         LOGMIN=15.0 
         LATMAXT=60.0
         LOGMAXT=300.0
         LOGMINT=5.0
         climfac=91.8
         if(BAS.EQ.'BM')climfac=39.93
        ELSE IF(BAS.EQ.'SA')THEN
         LATMAX=50.0
         LOGMAX=380.0
         LOGMIN=300.0
         LATMAXT=60.0
         LOGMAXT=380.0
         LOGMINT=290.0
        ELSE IF(BAS.EQ.'GB')THEN
         LATMAX=80.0
         LATMAXG=80.0
         LATMING=-65.0
         LOGMAX=359.0
         LOGMIN=0.0
         LATMAXT=87.0 
         LOGMAXT=380.0
         LOGMINT=0.0
        ELSE IF(BAS.EQ.'MT')THEN
         LATMAX=45.0
         LATMAXG=45.0
         LATMING=30.0
         LOGMAX=40.0 
         LOGMIN=0.0
         LATMAXT=80.0
         LOGMAXT=60.0
         LOGMINT=0.0
        END IF
c
c  ***   Revise latitude bound     ***
c
        latmaxg=min(latmaxg,latmax)
c
c  ***   Calculate area of basin    ***
c
        PI=ACOS(-1.0)
        ADEG=PI/180.
        areafac=(sin(adeg*latmaxg)-sin(adeg*latming))*(logmax-logmin)
        IF(BAS.EQ.'AL')then
         areafac=areafac-(sin(adeg*15.0)-sin(adeg*latming))*
     1     (275.0-logmin)
        END IF
        IF(BAS.EQ.'EP')THEN
         areafac=areafac-(sin(adeg*latmaxg)-sin(adeg*16.0))*
     1    (logmax-262.0)
        END IF
        IF(BAS.EQ.'AP')THEN
         areafac=areafac-(sin(adeg*latmaxg)-sin(adeg*16.0))*
     1    (285.0-262.0)-(sin(adeg*15.0)-sin(adeg*latming))*
     2    (275.0-260.0)
        END IF
C
        ENDTIME=30.*24.*3600.    
        DELPRINT=2.*3600.
c
        vsfac=1.2   ! Standard value 0.7
c
        flfac=PI/(6.0*3600.0)
        AFAC2=ADEG*6.38E6
        DISFACTOR=1./(60.*1852.)
        COEFS=2.*PI/(TAUMAX*24.*3600.)
        DTT=30.0D0*60.0D0
C
        IFLAG=0
c       ICOUNT=0
C
C   ***  Set default values of certain parameters
C
        TS=27.0
        TO=-70.0
        H_A=0.8
        PA=1012.0
        FC=5.0
        CECD=1.0
        CD=0.5
        CD1=3.0
        TRANSFAC=0.7
        CDM=3.0
        CDMAX=1.8
        CKMAX=2.0
        VCD=30.0
        IVMAXMIN=15
        CKEXP=2.0
        DT=240.0
        HATM=2500.0
        HATMI=1./HATM
C
        ES=6.112*EXP(17.67*TS/(243.5+TS))
        EA=H_A*ES
        QS=0.622*ES/PA
        TSA=TS+273.15
        TOA=TO+273.15
        EF=(TS-TO)/TSA
        CHI=2.5E6*EF*QS*(1.-H_A)
        SCHI=SQRT(CHI)
        FC=FC*1.0E-5
        CD=CD*0.001
        CD1=CD1*1.0E-5
        ATHA=1000.*LOG(TSA)+2.5E6*QS*H_A/TSA-QS*H_A*
     1    461.*LOG(H_A)
        THETA_E=EXP((ATHA-287.*LOG(0.001*PA))*0.001)
        PT=1000.0*(TOA/THETA_E)**3.484
        DELP=0.5*(PA-PT)
        TM=0.85*TS+0.15*TO
        ESM=6.112*EXP(17.67*TM/(243.5+TM))
        QSM=0.622*ESM/(PA-0.5*DELP)
        TM1=TM+3.0
        ESM1=6.112*EXP(17.67*TM1/(243.5+TM1))
        QSM1=0.622*ESM1/(PA-0.5*DELP)
        TMA=TM+273.15
        GRATB=(1.+2.5E6*QSM/(287.*TMA))/
     1   (1.+2.5E6*2.5E6*QSM/(461.*1005.*TMA*TMA))
        BETA=CHI/(287.*TSA)
        Q=0.5*H_A/(GRATB*(1.-H_A))
        RM=RM*1000.0
        RMINIT=RM
        RMBE=RM
        DT0=DT
        AKNOTFAC=1852.0/3600.0
        IF(gmeth.eq.'clim')THEN
         VM=15.*1852./3600
        ELSE
         VM=12.*1852./3600
        END IF
        VMINIT0=VM
	XM0=0.3
        CDV=CD1/CD
        RHAA=H_A
C
C   ***   Find Theoretical Maximum Wind and Minimum Pressure
C
        CALL THEORY(CECD,RHAA,TS,TO,PA,VMTH,PCTH,IFAIL)
C
C   ***   SET CERTAIN CONSTANTS    ***
C        
        TIME=0.0
        NT=TIME/DT
        DAMP=0.1
C
C   *** Specify initial intensification rate  ***
C
        IF(gmeth.eq.'clim')THEN
c         VINTEN=1.5/(6.*3600.)
         VINTEN=100.0
        ELSE
         VINTEN=100.0
        END IF
        VINTEN0=VINTEN
C
C   *** Read in monthly mean wind data  ***
C
        OPEN(11,file=TRIM(DIREC3)//'/u850',
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)U850
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/u'//lev,
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)U250
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/v850',
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)V850
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/v'//lev,
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)V250
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/u850u850',
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)UVAR850
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/u'//lev//'u'//lev,
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)UVAR250
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/v850v850',
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)VVAR850
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/v'//lev//'v'//lev,
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)VVAR250
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/u'//lev//'v'//lev,
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)U250V250
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/u'//lev//'u850',
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)U250U850
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/u850v850',
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)U850V850
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/v'//lev//'v850',
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)V250V850
        CLOSE(11)
C
        U250V850=0.0
        V250U850=0.0
        VDRIFT=VDRIFT_nocov
        UDRIFT=UDRIFT_nocov
C
        OPEN(11,file=TRIM(DIREC3)//'/u'//lev//'v850',
     1         status='old', access='direct',recl=8*NLO*NLA*12,
     2   err=31)
        READ(11,REC=1)U250V850
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/v'//lev//'u850',
     1         status='old', access='direct',recl=8*NLO*NLA*12,
     2   err=31)
        READ(11,REC=1)V250U850
        CLOSE(11)
C
        VDRIFT=VDRIFT_cov
        UDRIFT=UDRIFT_cov
C
   31   CONTINUE
C
C   ***   Calculate vorticity   ***
C
        do k=1,12
         do j=1,nla
          jp=j+1
          jp=min(jp,nla)
          jm=j-1
          jm=max(jm,1)
          coslat=cos((-90.+latinc*float(j-1))*adeg)
          coslatp=cos((-90.+latinc*float(jp-1))*adeg)
          coslatm=cos((-90.+latinc*float(jm-1))*adeg)
          do i=1,nlo
           ip=i+1
           if(i.eq.nlo)ip=1
           im=i-1
           if(i.eq.1)im=nlo
           vor(i,j,k)=((u850(i,jm,k)*coslatm-u850(i,jp,k)*coslatp)/
     1       (coslat*latinc)
     2      +(v850(ip,j,k)-v850(im,j,k))/longinc)/(2.0*111120.) 
          end do
         end do
        end do
C
C   ***  Read in potential intensity, humidity, and mixed layer files   ***
C
        OPEN(11,file=TRIM(DIREC3)//'/rh',
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)RH
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/T600',
     1   status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)T600
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/pi',
     1         status='old', access='direct',recl=8*NLO*NLA*12)
        READ(11,REC=1)VTEM
        CLOSE(11)
c
c  ***   Added factor of 1.1 on 2/20/2007 to bring histograms more in line with obs ***
c  ***          If statement added to deal with NaNs in VPOTC file                  ***
c
c       year=1000*(ichar(run(1:1))-48)+100*(ichar(run(2:2))-48)+
c     1   10*(ichar(run(3:3))-48)+ichar(run(4:4))-48
c
      DO 1380 K=1,12
        DO 1350 M=1,NLO
         DO 1350 N1=1,NLA
            IF(VTEM(M,N1,K).GT.0.0.AND.VTEM(M,N1,K).LE.200.0)THEN
c             VPOTC(M,N1,K)=PIAMP*0.94*MAX(VTEM(M,N1,K),0.0D0)
             VPOTC(M,N1,K)=PIAMP*0.9*MAX(VTEM(M,N1,K),0.0D0)
            ELSE
             VPOTC(M,N1,K)=0.0
            END IF
 1350   CONTINUE
C
 1380        CONTINUE
C
C ***  Find out whether time-evolving upper ocean data exist   ***
C
       inquire(file=TRIM(DIREC3)//'/hmix', exist=ex)
C
C ***  If time-evolving upper ocean data exist, use it, otherwise  ***
C ***                 use Levitus climatology                      ***
C
       IF(ex)then
C
        OPEN(11,file=TRIM(DIREC3)//'/hmix',
     1         status='old', access='direct',recl=8*NLAOC*NLOOC*12)
        READ(11,REC=1)MIXDEPTH
        CLOSE(11)
C
        OPEN(11,file=TRIM(DIREC3)//'/gamma',
     1         status='old', access='direct',recl=8*NLAOC*NLOOC*12)
        READ(11,REC=1)THERMAL
        CLOSE(11)
C
       ELSE
C
       DO K=1,12
         WRITE(AMON,1310)K
         OPEN(11,FILE='climdata/mix'//AMON, STATUS='OLD')
         OPEN(12,FILE='climdata/stratnorm'//AMON, STATUS='OLD')
         DO I=1,180
          READ(11,*)(MIXDEPTHd(K,I,J),J=1,360)
          READ(12,*)(THERMALd(K,I,J),J=1,360)
         END DO
         close(11)
         close(12)
         DO J=1,360
          DO I=1,180
            MIXDEPTH(K,I,J)=MIXDEPTHd(K,I,J)
            THERMAL(K,I,J)=THERMALd(K,I,J)
          END DO
          IF(MODEL(1:2).NE.'dk')THEN
            MIXDEPTH(K,181,J)=MIXDEPTH(K,180,J)
            THERMAL(K,181,J)=THERMAL(K,180,J) 
          END IF
c
c     Variants on mixed layer depth
c
         IF(FLAV(1:2).EQ.'eo'.OR.RUN(1:2).EQ.'pl'.OR.
     1      FLAV(1:2).EQ.'pe')THEN
          DO I=1,NLAOC
           MIXDEPTH(K,I,J)=100.0
          END DO
         END IF
C
        END DO
C
       END DO
C
       END IF
c
 1310 FORMAT(I2.2)
c
C   ***  Read in bathymetric and drag files   ***
C
        NAM='climdata/bathsmooth2.dat'
        OPEN(11,FILE=NAM,status='old',
     1    access='direct',recl=4*1440*721)
        READ(11,REC=1)BATH
        close(11)
c
        bathm=max(bath,-40.0)
C
        NAM='climdata/C_Drag500.dat'
        OPEN(11,FILE=NAM,status='old',
     1    access='direct',recl=8*1440*721)
        READ(11,REC=1)CDRAG
        close(11)
        CDRAG=CDRAG/(1.0+50.0*CDRAG) ! Accounts for application to gradient wind 12/2018
c        CDRAG=CDRAG/MINVAL(CDRAG) ! Normalizes to unity over water
C
C           ***  Re-initialize random numbers  ***
C
        if(rerandom.eq.1.or.gmeth.eq.'clim')then
          call random_seed ( PUT = adate(1:3))
        end if
c
        M=0
        NTOTAL=0
        NUMSTORMS=0
        IF(restart.eq.1)THEN
         OPEN(UNIT=12,FILE=TRIM(CITYNAME)//'/stats'//filenum//'.out',
     1    status='unknown')
         read(12,1223)NUMSTORMS,NTOTAL
         CLOSE(12)
        END IF
 1223   FORMAT(85X,I5,1X,I9)
c
c        IF(gmeth.eq.'clim')THEN
c          call random_seed ( PUT = adate(1:3))
c        END IF
C
        GENPOINTS='gen'//BAS//'.in'
        IF(BAS.NE.'SA'.AND.BAS.NE.'MT'.AND.BAS.NE.'AP')THEN
         OPEN(UNIT=18, FILE=GENPOINTS, status='old')
        END IF
C
    5   CONTINUE
C
    7   continue
C
        NTOTAL=NTOTAL+1
        IFLAG=0
C       IFLAGB=0
c
        IF(gmeth.ne.'clim')THEN
c
c       IMON=1+int(12.0*RAND(0))
        call random_number(arand)
        IMON=1+int(12.0*arand)
        IMON=MIN(IMON,12)
c        IDATE=1+INT(FLOAT(DAY(IMON))*RAND(0))
        call random_number(arand)
        IDATE=1+INT(FLOAT(DAY(IMON))*arand)
        IDATE=MIN(IDATE,DAY(IMON))
        call random_number(arand)
        IHOUR=6*int(4.0*arand)
c        IHOUR=6*int(4.0*rand(0))
        IF(IHOUR.EQ.24)IHOUR=0
        ymin=sin(adeg*latming)
        ymax=sin(adeg*latmaxg)
c        y=ymin+(ymax-ymin)*RAND(0)
        call random_number(arand)
        y=ymin+(ymax-ymin)*arand
        AILAT=ASIN(y)/adeg
        IF(BAS.EQ.'SH'.OR.BAS.EQ.'SA'.OR.BAS.EQ.'BM')AILAT=-AILAT
c        AILOG=LOGMIN+(LOGMAX-LOGMIN)*RAND(0)
        call random_number(arand)
        AILOG=LOGMIN+(LOGMAX-LOGMIN)*arand
c
        if(bas.eq.'AL'.and.((ailat.lt.15.8.and.ailog.lt.276.5).or. 
     1   (ailat.lt.9.0.and.ailog.lt.288.0)))then
         NTOTAL=NTOTAL-1
         call random_number(gdum)
         goto 7
        end if
c
        if(bas.eq.'EP'.and.ailat.gt.16.0.and.ailog.gt.262.0)then
          NTOTAL=NTOTAL-1
          call random_number(gdum)
         goto 7
        end if
c
        if(bas.eq.'AP'.and.ailat.gt.16.0.and.ailog.gt.262.0.and.
     1   ((ailat.lt.15.8.and.ailog.lt.276.5).or.
     2   (ailat.lt.9.0.and.ailog.lt.288.0)))then
          NTOTAL=NTOTAL-1
          call random_number(gdum)
         goto 7
        end if
c
        if(abs(ailat).lt.3.0)then
          NTOTAL=NTOTAL-1
          call random_number(gdum)
         goto 7
        end if
C
        ELSE
C
  121   CONTINUE
C
        READ(18,*,END=201)IMON, IDATE, IHOUR,AILAT,AILOG
        IF(AILOG.LT.0.0)AILOG=AILOG+360.0
         call random_number(gdum)
        GOTO 401
  201   CONTINUE
        REWIND(18)
        GOTO 121
  401   CONTINUE
C
        END IF
C
         JLAT=1+(AILAT-LATSTART)/LATINC
         IF(AILOG.GT.0.0)THEN
          JLOGG=1+(AILOG-LONGSTART)/LONGINC
         ELSE
          JLOGG=1+(360.+AILOG-LONGSTART)/LONGINC
         END IF
         JLATP=JLAT+1
         JLOGGP=JLOGG+1
         IF(JLOGGP.GT.NLO)JLOGGP=1
C
         IF(IDATE.GE.15)THEN
          IMONK=IMON
          IMONKP=IMON+1
          IF(IMONKP.EQ.13)IMONKP=1
          FAC=FLOAT(IDATE-15)/30.
         ELSE
          IMONK=IMON-1
          IF(IMONK.EQ.0)IMONK=12
          IMONKP=IMON
          FAC=FLOAT(15+IDATE)/30.
         END IF
         VOR1=(1.-FAC)*VOR(JLOGG,JLAT,IMONK)+
     1      FAC*VOR(JLOGG,JLAT,IMONKP)
         VOR2=(1.-FAC)*VOR(JLOGG,JLATP,IMONK)+FAC*
     1      VOR(JLOGG,JLATP,IMONKP)
         VOR3=(1.-FAC)*VOR(JLOGGP,JLAT,IMONK)+FAC*
     1      VOR(JLOGGP,JLAT,IMONKP)
         VOR4=(1.-FAC)*VOR(JLOGGP,JLATP,IMONK)+FAC*
     1      VOR(JLOGGP,JLATP,IMONKP)
         DELY=1.+(AILAT-LATSTART)/LATINC-FLOAT(JLAT)
         DELX=1.+(360.+AILOG-LONGSTART)/LONGINC-
     1      INT(1.+(360.+AILOG-LONGSTART)/LONGINC)
         D1=(1.-DELX)*(1.-DELY)
         D2=DELY*(1.-DELX)
         D3=DELX*(1.-DELY)
         D4=DELX*DELY
         VORTRACK=D1*VOR1+D2*VOR2+D3*VOR3+D4*VOR4
c
        IF(gmeth.ne.'clim')THEN
c
c   ***  Estimate probability of genesis  ***
c
        flocal=flfac*sin(adeg*ailat)
        wht=1.0e4*flocal+6.0e4*vortrack
        wht=wht*2.0e4*flocal
        wht=max(wht,0.03)
c
        prob=1.0-exp(-wht)
c
c   ***  Lines added 8/26/2020 as experiment to create better cut off at low latitudes
c
         dum=(abs(ailat)-1.0)/12.5
         dum=max(dum,0.0)
         dum=min(dum,1.0)
         plowlat=dum**3.5
         prob=prob*plowlat
c
c   ***  End of addition
c
c        prob=prob**3
c
        call random_number(arand)
        if(arand.ge.prob)goto 7
c        if(rand(0).ge.prob)goto 7
c
        END IF
c
        SIGN=1.0
        IF(AILAT.LT.0.0)SIGN=-1.0
        VDRIFTNEW=VDRIFT*SIGN
        NSTEPS=0
        IMONI=IMON
        IDATEI=IDATE
        IHOURI=IHOUR
C
        ALAT=AILAT
        ALOG=AILOG
        IMON=IMONI
        IDATE=IDATEI
        IHOUR=IHOURI
C
        TT=0.
        RMIN=1.0E14
        JTIME=0
        PRINTIME=-DTT
        TIME=3600.*FLOAT(IHOUR)
C
        SUMM=0.0
c
        DO I=1,15
         BI=FLOAT(I)
         SUMM=SUMM+1./(BI**3)
         AMP(I)=1./(BI**1.5)
         call random_number(PHPS)
         PHASEU850(I)=2.*PI*(PHPS-0.5)
         call random_number(PHPS)
         PHASEV850(I)=2.*PI*(PHPS-0.5)
         call random_number(PHPS)
         PHASEU250(I)=2.*PI*(PHPS-0.5)
         call random_number(PHPS)
         PHASEV250(I)=2.*PI*(PHPS-0.5)
        END DO
        SUMM=SQRT(2./SUMM)
C
C      ***  Time march storm position   ***
C
  101 CONTINUE
C
        TT=TT+DTT
        PRINTIME=PRINTIME+DTT
C
        TIME=TIME+DTT
        IF(TIME.GT.24.*3600.)THEN
         TIME=TIME-24.*3600.
         IDATE=IDATE+1
            IF((IMON.EQ.1.OR.IMON.EQ.3.OR.IMON.EQ.5.OR.IMON.EQ.7.OR.
     1     IMON.EQ.8.OR.IMON.EQ.10).AND.IDATE.GT.31)THEN
                  IDATE=1
              IMON=IMON+1
            END IF
            IF((IMON.EQ.4.OR.IMON.EQ.6.OR.IMON.EQ.9.OR.IMON.EQ.11)
     1                .AND.IDATE.GT.30)THEN
                  IDATE=1
              IMON=IMON+1
            END IF
            IF(IMON.EQ.2.AND.IDATE.GT.28)THEN
              IDATE=1
              IMON=IMON+1
            END IF
            IF(IMON.GT.12)IMON=1
        END IF
        IHOUR=TIME/3600.
        IMIN=TIME/60.
C
         JLAT=1+(ALAT-LATSTART)/LATINC
         IF(ALOG.GT.0.0)THEN
          JLOGG=1+(ALOG-LONGSTART)/LONGINC
         ELSE
          JLOGG=1+(360.+ALOG-LONGSTART)/LONGINC
         END IF
         JLATP=JLAT+1
         JLOGGP=JLOGG+1
         IF(JLOGG.GT.NLO)JLOGG=JLOGG-NLO
         IF(JLOGGP.GT.NLO)JLOGGP=JLOGGP-NLO
C
          IF(IDATE.GE.15)THEN
           IMONK=IMON
           IMONKP=IMON+1
           IF(IMONKP.EQ.13)IMONKP=1
           FAC=FLOAT(IDATE-15)/30.
          ELSE
           IMONK=IMON-1
           IF(IMONK.EQ.0)IMONK=12
           IMONKP=IMON
           FAC=FLOAT(15+IDATE)/30.
          END IF
          FAC=MIN(FAC,1.0)
          FAC=MAX(FAC,0.0)
C
           DELY=1.+(ALAT-LATSTART)/LATINC-FLOAT(JLAT)
           DELX=1.+(360.+ALOG-LONGSTART)/LONGINC-
     1      INT(1.+(360.+ALOG-LONGSTART)/LONGINC)
           D1=(1.-DELX)*(1.-DELY)
           D2=DELY*(1.-DELX)
           D3=DELX*(1.-DELY)
           D4=DELX*DELY
C
        do ii=1,4
          do jj=1,4
            A1(ii,jj)=0.0
            B1(ii,jj)=0.0
            A2(ii,jj)=0.0
            B2(ii,jj)=0.0
            A3(ii,jj)=0.0
            B3(ii,jj)=0.0
            A4(ii,jj)=0.0
            B4(ii,jj)=0.0
          end do
        end do
c
           SU1=(1.-FAC)*U850(JLOGG,JLAT,IMONK)+
     1            FAC*U850(JLOGG,JLAT,IMONKP)
           SU2=(1.-FAC)*U850(JLOGG,JLATP,IMONK)+FAC*
     1      U850(JLOGG,JLATP,IMONKP)
           SU3=(1.-FAC)*U850(JLOGGP,JLAT,IMONK)+FAC*
     1      U850(JLOGGP,JLAT,IMONKP)
           SU4=(1.-FAC)*U850(JLOGGP,JLATP,IMONK)+FAC*
     1      U850(JLOGGP,JLATP,IMONKP)
           SV1=(1.-FAC)*V850(JLOGG,JLAT,IMONK)+
     1            FAC*V850(JLOGG,JLAT,IMONKP)
           SV2=(1.-FAC)*V850(JLOGG,JLATP,IMONK)+FAC*
     1      V850(JLOGG,JLATP,IMONKP)
           SV3=(1.-FAC)*V850(JLOGGP,JLAT,IMONK)+FAC*
     1      V850(JLOGGP,JLAT,IMONKP)
           SV4=(1.-FAC)*V850(JLOGGP,JLATP,IMONK)+FAC*
     1      V850(JLOGGP,JLATP,IMONKP)
C
          U850TRACK=D1*SU1+D2*SU2+D3*SU3+D4*SU4
          V850TRACK=D1*SV1+D2*SV2+D3*SV3+D4*SV4
c
c           SU1=(1.-FAC)*UVAR850(JLOGG,JLAT,IMONK)+
c     1            FAC*UVAR850(JLOGG,JLAT,IMONKP)
c           SU2=(1.-FAC)*UVAR850(JLOGG,JLATP,IMONK)+FAC*
c     1      UVAR850(JLOGG,JLATP,IMONKP)
c           SU3=(1.-FAC)*UVAR850(JLOGGP,JLAT,IMONK)+FAC*
c     1      UVAR850(JLOGGP,JLAT,IMONKP)
c           SU4=(1.-FAC)*UVAR850(JLOGGP,JLATP,IMONK)+FAC*
c     1      UVAR850(JLOGGP,JLATP,IMONKP)
c           SV1=(1.-FAC)*VVAR850(JLOGG,JLAT,IMONK)+
c     1            FAC*VVAR850(JLOGG,JLAT,IMONKP)
c           SV2=(1.-FAC)*VVAR850(JLOGG,JLATP,IMONK)+FAC*
c     1      VVAR850(JLOGG,JLATP,IMONKP)
c           SV3=(1.-FAC)*VVAR850(JLOGGP,JLAT,IMONK)+FAC*
c     1      VVAR850(JLOGGP,JLAT,IMONKP)
c           SV4=(1.-FAC)*VVAR850(JLOGGP,JLATP,IMONK)+FAC*
c     1      VVAR850(JLOGGP,JLATP,IMONKP)
C
c          UVAR850TRACK=D1*SU1+D2*SU2+D3*SU3+D4*SU4
c          VVAR850TRACK=D1*SV1+D2*SV2+D3*SV3+D4*SV4
C
           SU1=(1.-FAC)*U250(JLOGG,JLAT,IMONK)+
     1            FAC*U250(JLOGG,JLAT,IMONKP)
           SU2=(1.-FAC)*U250(JLOGG,JLATP,IMONK)+FAC*
     1      U250(JLOGG,JLATP,IMONKP)
           SU3=(1.-FAC)*U250(JLOGGP,JLAT,IMONK)+FAC*
     1      U250(JLOGGP,JLAT,IMONKP)
           SU4=(1.-FAC)*U250(JLOGGP,JLATP,IMONK)+FAC*
     1      U250(JLOGGP,JLATP,IMONKP)
           SV1=(1.-FAC)*V250(JLOGG,JLAT,IMONK)+
     1            FAC*V250(JLOGG,JLAT,IMONKP)
           SV2=(1.-FAC)*V250(JLOGG,JLATP,IMONK)+FAC*
     1      V250(JLOGG,JLATP,IMONKP)
           SV3=(1.-FAC)*V250(JLOGGP,JLAT,IMONK)+FAC*
     1      V250(JLOGGP,JLAT,IMONKP)
           SV4=(1.-FAC)*V250(JLOGGP,JLATP,IMONK)+FAC*
     1      V250(JLOGGP,JLATP,IMONKP)
C
          U250TRACK=D1*SU1+D2*SU2+D3*SU3+D4*SU4
          V250TRACK=D1*SV1+D2*SV2+D3*SV3+D4*SV4
C
           SU1=(1.-FAC)*T600(JLOGG,JLAT,IMONK)+
     1            FAC*T600(JLOGG,JLAT,IMONKP)
           SU2=(1.-FAC)*T600(JLOGG,JLATP,IMONK)+FAC*
     1      T600(JLOGG,JLATP,IMONKP)
           SU3=(1.-FAC)*T600(JLOGGP,JLAT,IMONK)+FAC*
     1      T600(JLOGGP,JLAT,IMONKP)
           SU4=(1.-FAC)*T600(JLOGGP,JLATP,IMONK)+FAC*
     1      T600(JLOGGP,JLATP,IMONKP)
C
           T600TEMP=D1*SU1+D2*SU2+D3*SU3+D4*SU4
C
           SU1=(1.-FAC)*RH(JLOGG,JLAT,IMONK)+
     1            FAC*RH(JLOGG,JLAT,IMONKP)
           SU2=(1.-FAC)*RH(JLOGG,JLATP,IMONK)+FAC*
     1      RH(JLOGG,JLATP,IMONKP)
           SU3=(1.-FAC)*RH(JLOGGP,JLAT,IMONK)+FAC*
     1      RH(JLOGGP,JLAT,IMONKP)
           SU4=(1.-FAC)*RH(JLOGGP,JLATP,IMONK)+FAC*
     1      RH(JLOGGP,JLATP,IMONKP)
C
           RHTEMP=D1*SU1+D2*SU2+D3*SU3+D4*SU4
C
           A1(1,2)=U250V250(JLOGG,JLAT,IMONK)
           B1(1,2)=U250V250(JLOGG,JLAT,IMONKP)
           A2(1,2)=U250V250(JLOGG,JLATP,IMONK)
           B2(1,2)=U250V250(JLOGG,JLATP,IMONKP)
           A3(1,2)=U250V250(JLOGGP,JLAT,IMONK)
           B3(1,2)=U250V250(JLOGGP,JLAT,IMONKP)
           A4(1,2)=U250V250(JLOGGP,JLATP,IMONK)
           B4(1,2)=U250V250(JLOGGP,JLATP,IMONKP)
c
           A1(1,3)=U250U850(JLOGG,JLAT,IMONK)
           B1(1,3)=U250U850(JLOGG,JLAT,IMONKP)
           A2(1,3)=U250U850(JLOGG,JLATP,IMONK)
           B2(1,3)=U250U850(JLOGG,JLATP,IMONKP)
           A3(1,3)=U250U850(JLOGGP,JLAT,IMONK)
           B3(1,3)=U250U850(JLOGGP,JLAT,IMONKP)
           A4(1,3)=U250U850(JLOGGP,JLATP,IMONK)
           B4(1,3)=U250U850(JLOGGP,JLATP,IMONKP)
c
           A1(1,4)=U250V850(JLOGG,JLAT,IMONK)
           B1(1,4)=U250V850(JLOGG,JLAT,IMONKP)
           A2(1,4)=U250V850(JLOGG,JLATP,IMONK)
           B2(1,4)=U250V850(JLOGG,JLATP,IMONKP)
           A3(1,4)=U250V850(JLOGGP,JLAT,IMONK)
           B3(1,4)=U250V850(JLOGGP,JLAT,IMONKP)
           A4(1,4)=U250V850(JLOGGP,JLATP,IMONK)
           B4(1,4)=U250V850(JLOGGP,JLATP,IMONKP)
c
           A1(2,3)=V250U850(JLOGG,JLAT,IMONK)
           B1(2,3)=V250U850(JLOGG,JLAT,IMONKP)
           A2(2,3)=V250U850(JLOGG,JLATP,IMONK)
           B2(2,3)=V250U850(JLOGG,JLATP,IMONKP)
           A3(2,3)=V250U850(JLOGGP,JLAT,IMONK)
           B3(2,3)=V250U850(JLOGGP,JLAT,IMONKP)
           A4(2,3)=V250U850(JLOGGP,JLATP,IMONK)
           B4(2,3)=V250U850(JLOGGP,JLATP,IMONKP)
C    
           A1(2,4)=V250V850(JLOGG,JLAT,IMONK)
           B1(2,4)=V250V850(JLOGG,JLAT,IMONKP)
           A2(2,4)=V250V850(JLOGG,JLATP,IMONK)
           B2(2,4)=V250V850(JLOGG,JLATP,IMONKP)
           A3(2,4)=V250V850(JLOGGP,JLAT,IMONK)
           B3(2,4)=V250V850(JLOGGP,JLAT,IMONKP)
           A4(2,4)=V250V850(JLOGGP,JLATP,IMONK)
           B4(2,4)=V250V850(JLOGGP,JLATP,IMONKP)
C
           A1(3,4)=U850V850(JLOGG,JLAT,IMONK)
           B1(3,4)=U850V850(JLOGG,JLAT,IMONKP)
           A2(3,4)=U850V850(JLOGG,JLATP,IMONK)
           B2(3,4)=U850V850(JLOGG,JLATP,IMONKP)
           A3(3,4)=U850V850(JLOGGP,JLAT,IMONK)
           B3(3,4)=U850V850(JLOGGP,JLAT,IMONKP)
           A4(3,4)=U850V850(JLOGGP,JLATP,IMONK)
           B4(3,4)=U850V850(JLOGGP,JLATP,IMONKP)
C
           A1(3,3)=UVAR850(JLOGG,JLAT,IMONK)
                B1(3,3)=UVAR850(JLOGG,JLAT,IMONKP)
           A2(3,3)=UVAR850(JLOGG,JLATP,IMONK)
           B2(3,3)=UVAR850(JLOGG,JLATP,IMONKP)
           A3(3,3)=UVAR850(JLOGGP,JLAT,IMONK)
           B3(3,3)=UVAR850(JLOGGP,JLAT,IMONKP)
           A4(3,3)=UVAR850(JLOGGP,JLATP,IMONK)
           B4(3,3)=UVAR850(JLOGGP,JLATP,IMONKP)
c
           A1(4,4)=VVAR850(JLOGG,JLAT,IMONK)
                B1(4,4)=VVAR850(JLOGG,JLAT,IMONKP)
           A2(4,4)=VVAR850(JLOGG,JLATP,IMONK)
           B2(4,4)=VVAR850(JLOGG,JLATP,IMONKP)
           A3(4,4)=VVAR850(JLOGGP,JLAT,IMONK)
           B3(4,4)=VVAR850(JLOGGP,JLAT,IMONKP)
           A4(4,4)=VVAR850(JLOGGP,JLATP,IMONK)
           B4(4,4)=VVAR850(JLOGGP,JLATP,IMONKP)
C
           A1(1,1)=UVAR250(JLOGG,JLAT,IMONK)
                B1(1,1)=UVAR250(JLOGG,JLAT,IMONKP)
           A2(1,1)=UVAR250(JLOGG,JLATP,IMONK)
           B2(1,1)=UVAR250(JLOGG,JLATP,IMONKP)
           A3(1,1)=UVAR250(JLOGGP,JLAT,IMONK)
           B3(1,1)=UVAR250(JLOGGP,JLAT,IMONKP)
           A4(1,1)=UVAR250(JLOGGP,JLATP,IMONK)
           B4(1,1)=UVAR250(JLOGGP,JLATP,IMONKP)
C
           A1(2,2)=VVAR250(JLOGG,JLAT,IMONK)
                B1(2,2)=VVAR250(JLOGG,JLAT,IMONKP)
           A2(2,2)=VVAR250(JLOGG,JLATP,IMONK)
           B2(2,2)=VVAR250(JLOGG,JLATP,IMONKP)
           A3(2,2)=VVAR250(JLOGGP,JLAT,IMONK)
           B3(2,2)=VVAR250(JLOGGP,JLAT,IMONKP)
           A4(2,2)=VVAR250(JLOGGP,JLATP,IMONK)
           B4(2,2)=VVAR250(JLOGGP,JLATP,IMONKP)
C
           PIT1=(1.-FAC)*VPOTC(JLOGG,JLAT,IMONK)+
     1      FAC*VPOTC(JLOGG,JLAT,IMONKP)
           PIT2=(1.-FAC)*VPOTC(JLOGG,JLATP,IMONK)+FAC*
     1      VPOTC(JLOGG,JLATP,IMONKP)
           PIT3=(1.-FAC)*VPOTC(JLOGGP,JLAT,IMONK)+FAC*
     1      VPOTC(JLOGGP,JLAT,IMONKP)
           PIT4=(1.-FAC)*VPOTC(JLOGGP,JLATP,IMONK)+FAC*
     1      VPOTC(JLOGGP,JLATP,IMONKP)
C
C  ***   Added 4/25/2015 to set potential intensities to be zero over land ***
C
        flat=1+4*(latstart+90+latinc*(jlat-1))
        klat=flat
        flogg=1+4*(longstart+longinc*(jlogg-1))
        klogg=flogg
c
        if(klogg.gt.1440)klogg=klogg-1440
c
        klatp=klat+1
        kloggp=klogg+1
        if(kloggp.gt.1440)kloggp=kloggp-1440
        flat2=klat+4*latinc
        klat2=flat2
        flogg2=klogg+4*longinc
        klogg2=flogg2
        if(klogg2.gt.1440)klogg2=klogg2-1440
        klat2p=klat2+1
        klogg2p=klogg2+1
        if(klogg2p.gt.1440)klogg2p=klogg2p-1440
c
c  Reduce aliasing from very deep water
c
        batha=bathm(klogg,klat)
        bathb=bathm(klogg,klatp)
        bathc=bathm(kloggp,klat)
        bathd=bathm(kloggp,klatp)
        delya=flat-float(klat)
        delxa=flogg-float(int(flogg))
        D1a=(1.-DELXa)*(1.-DELYa)
        D2a=DELYa*(1.-DELXa)
        D3a=DELXa*(1.-DELYa)
        D4a=DELXa*DELYa
        bath1=D1a*batha+D2a*bathb+D3a*bathc+D4a*bathd
        bath1=min(bath1,1.0)
        bath1=max(bath1,0.0)
        bath1=1.0-bath1
c
        batha=bathm(klogg,klat2)
        bathb=bathm(klogg,klat2p)
        bathc=bathm(kloggp,klat2)
        bathd=bathm(kloggp,klat2p)
        delya=flat2-float(klat2)
        delxa=flogg-float(int(flogg))
        D1a=(1.-DELXa)*(1.-DELYa)
        D2a=DELYa*(1.-DELXa)
        D3a=DELXa*(1.-DELYa)
        D4a=DELXa*DELYa
        bath2=D1a*batha+D2a*bathb+D3a*bathc+D4a*bathd
        bath2=min(bath2,1.0)
        bath2=max(bath2,0.0)
        bath2=1.0-bath2
c
        batha=bathm(klogg2,klat)
        bathb=bathm(klogg2,klatp)
        bathc=bathm(klogg2p,klat)
        bathd=bathm(klogg2p,klatp)
        delya=flat-float(klat)
        delxa=flogg2-float(int(flogg2))
        D1a=(1.-DELXa)*(1.-DELYa)
        D2a=DELYa*(1.-DELXa)
        D3a=DELXa*(1.-DELYa)
        D4a=DELXa*DELYa
        bath3=D1a*batha+D2a*bathb+D3a*bathc+D4a*bathd
        bath3=min(bath3,1.0)
        bath3=max(bath3,0.0)
        bath3=1.0-bath3
c
        batha=bathm(klogg2,klat2)
        bathb=bathm(klogg2,klat2p)
        bathc=bathm(klogg2p,klat2)
        bathd=bathm(klogg2p,klat2p)
        delya=flat2-float(klat2)
        delxa=flogg2-float(int(flogg2))
        D1a=(1.-DELXa)*(1.-DELYa)
        D2a=DELYa*(1.-DELXa)
        D3a=DELXa*(1.-DELYa)
        D4a=DELXa*DELYa
        bath4=D1a*batha+D2a*bathb+D3a*bathc+D4a*bathd
        bath4=min(bath4,1.0)
        bath4=max(bath4,0.0)
        bath4=1.0-bath4
c
        pit1=pit1*bath1
        pit2=pit2*bath2
        pit3=pit3*bath3
        pit4=pit4*bath4
C   
        U850SUM=0.0
        V850SUM=0.0
        U250SUM=0.0
        V250SUM=0.0
        DO I=1,15
         U850SUM=U850SUM+AMP(I)*SIN(COEFS*FLOAT(I)*TT+PHASEU850(I))
         V850SUM=V850SUM+AMP(I)*SIN(COEFS*FLOAT(I)*TT+PHASEV850(I))
         U250SUM=U250SUM+AMP(I)*SIN(COEFS*FLOAT(I)*TT+PHASEU250(I))
         V250SUM=V250SUM+AMP(I)*SIN(COEFS*FLOAT(I)*TT+PHASEV250(I))
        END DO
C
C      Uses cholesky decomposition from numerical recipes
C
        CALL choldc(A1,4,4,PCA1)
        CALL choldc(A2,4,4,PCA2)
        CALL choldc(A3,4,4,PCA3)
        CALL choldc(A4,4,4,PCA4)
        CALL choldc(B1,4,4,PCB1)
        CALL choldc(B2,4,4,PCB2)
        CALL choldc(B3,4,4,PCB3)
        CALL choldc(B4,4,4,PCB4)
c
        do ii=1,4
           PCA(ii)=(1.-FAC)*(d1*PCA1(ii)+d2*PCA2(ii)+d3*PCA3(ii)+
     1      d4*PCA4(ii))+FAC*(d1*PCB1(ii)+d2*PCB2(ii)+d3*PCB3(ii)+
     2      d4*PCB4(ii))
           do jj=1,4
              A(ii,jj)=(1.-FAC)*(d1*A1(ii,jj)+d2*A2(ii,jj)+d3*
     1         A3(ii,jj)+d4*A4(ii,jj))+FAC*(d1*B1(ii,jj)+d2*
     2         B2(ii,jj)+d3*B3(ii,jj)+d4*B4(ii,jj))
           end do
        end do
C
        U250M=U250TRACK+SUMM*PCA(1)*U250SUM
        V250M=V250TRACK+SUMM*(A(2,1)*U250SUM+PCA(2)*V250SUM)
        U850M=U850TRACK+SUMM*(A(3,1)*U250SUM+A(3,2)*V250SUM+
     1    PCA(3)*U850SUM)
        V850M=V850TRACK+SUMM*(A(4,2)*V250SUM+A(4,3)*U850SUM+
     1    A(4,1)*U250SUM+PCA(4)*V850SUM)
C
c        PITMEAN=0.25*(PIT1+PIT2+PIT3+PIT4)
c        IF(PITMEAN.LE.20.0)THEN
c          WEIGHT=0.95
c        ELSE
c          WEIGHT=0.8
c        END IF
C
        UMEAN=WEIGHT*U850M+(1.-WEIGHT)*U250M+UDRIFT*COS(ALAT*ADEG)
        VMEAN=WEIGHT*V850M+(1.-WEIGHT)*V250M+VDRIFTNEW*COS(ALAT*ADEG)
        UNET=U250M-U850M
        VNET=V250M-V850M
c
c   These lines added 7/20/2011 to account for shear effects on storm motion
c
c        usdrift=-1.7*vnet/(3.0+abs(vnet))
c        vsdrift=1.7*unet/(3.0+abs(unet))
        usdrift=-0.07*vnet
        vsdrift=0.07*unet
        umean=umean+usdrift*sign
        vmean=vmean+vsdrift*sign
c
        VSHEAR=SQRT(UNET**2+VNET**2)
C
        ALAT=ALAT+DTT*VMEAN*DISFACTOR
        COSFAC=COS(ALAT*ADEG)
        ALOG=ALOG+DTT*UMEAN*DISFACTOR/COSFAC
        dlogg=alog
C
c   Added 7/9/2015: If track wanders equatorward of 2 degrees throw it away and start over
c
        IF(ABS(ALAT).LT.2.0)THEN
          NTOTAL=NTOTAL-1
          call random_number(gdum)
          GOTO 7
        END IF
C
        RADIUS=((ALAT-CLAT1)**2+(COSFAC*(dlogg-CLOG1))**2)
        RMIN=MIN(RADIUS,RMIN)
C
        IF(TT.LT.(1.5*DTT).OR.PRINTIME.GE.DELPRINT)THEN
         PRINTIME=0.0
         JTIME=JTIME+1
         MONTE(JTIME)=IMON
         DATETE(JTIME)=IDATE
         HOURTE(JTIME)=IHOUR
         LATTE(JTIME)=ALAT
         LOGGTE(JTIME)=ALOG
         SHEARTRACK(JTIME)=VSHEAR
         U850TR(JTIME)=U850M
         V850TR(JTIME)=V850M
c
c  For potential intensity, use maximum value at nearby grid points to avoid influence of land
c     zeros on storms near land (but still at sea)
c
         PIM=0.0
         IF(MAX(PIT1,PIT2,PIT2,PIT4).GT.0.0)THEN
          AD1=D1*MAX(MIN(PIT1,1.0),0.0)
          AD2=D2*MAX(MIN(PIT2,1.0),0.0)
          AD3=D3*MAX(MIN(PIT3,1.0),0.0)
          AD4=D4*MAX(MIN(PIT4,1.0),0.0)
          PIM=(AD1*PIT1+AD2*PIT2+AD3*PIT3+AD4*PIT4)/
     1      (AD1+AD2+AD3+AD4+1.0e-6)
         END IF
         PIM=MAX(PIM,0.0)
c
         RHM=RHTEMP
         RHM=MAX(RHM,0.0)
         RHM=MIN(RHM,100.0)
         T600M=T600TEMP
C
        VPTRACK(JTIME)=PIM
        RHTRACK(JTIME)=RHM
        T600TRACK(JTIME)=T600M
C
        END IF
C
        IF(TT.GT.ENDTIME.OR.ABS(ALAT).GT.LATMAXT.OR.ALOG.GT.LOGMAXT
     1    .OR.ALOG.LT.LOGMINT.OR.ABS(ALAT).LT.2.0)THEN
           GOTO 900
        ELSE
           GOTO 101
        END IF
c
  900 CONTINUE
C
        IF(shape .eq.'circ')then
c
        jint=1
        RMIN=0.001*SQRT(RMIN)/DISFACTOR
        IF(RMIN.LE.CRAD)THEN
         M=M+1
        ELSE
         GOTO 5
        END IF
c
        elseif(shape.eq.'poly')then
c
        CALL box(loggte,latte,jtime,x1vert,y1vert,x2vert,y2vert,
     1    nvert,nmap,xint,yint,jint)
c
        IF(nmap.eq.1)THEN
         M=M+1
         clog1=xint
         clat1=yint
        ELSE
         GOTO 5
        END IF
c
        end if
C
        CLOSE(11)
C
        DISTANCE(1)=0.0
        NST=0
        IMONI=MONTE(1)
        IDATEI=DATETE(1)
        ITIMEI=HOURTE(1)
c
        DO 28 I=1,JTIME
         IMON=MONTE(I)
         IDATE=DATETE(I)
         ITIME=HOURTE(I)
         ALAT=LATTE(I)
         ALOGG=LOGGTE(I)
         IF(I.NE.1.AND.ITIME.EQ.ITIMEOLD.AND.IDATE.EQ.IDATEOLD)GOTO 28
         NST=NST+1 
         IF(NST.GT.1)THEN
          DELX=AFAC2*COS(ALAT*ADEG)*(ALOGG-ALOGGOLD)
          DELY=AFAC2*(ALAT-ALATOLD)
          DELTIME=3600.*FLOAT(ITIME-ITIMEOLD)
          IF(IDATE.NE.IDATEOLD)DELTIME=DELTIME+3600*24
          DELDIS=SQRT(DELX*DELX+DELY*DELY)
          DISTANCE(NST)=DISTANCE(NST-1)+DELDIS
         END IF
         WRITETIME(NST)=ITIME
c
         ALATOLD=ALAT
         ALOGGOLD=ALOGG
         ITIMEOLD=ITIME
         ITIMEB=ITIME
         IMONOLD=IMON
         IDATEOLD=IDATE
   28   CONTINUE
        NST0=NST
C
         WRITE(MONC,1801)MONTE(NST0)
         WRITE(MOND,1801)MONTE(NST0)
         WRITE(DATEC,1801)DATETE(NST0)
         WRITE(DATETEMP,1805)WRITETIME(NST0)
 1801   FORMAT(I2.2)
 1805        FORMAT(I2.2,'Z')
C
        NSTORM=JTIME
        UTIME(1)=0.0
c
        DO J=1,NSTORM
         UTIME(J)=2.*3600.*FLOAT(J-1)
         ALAT=LATTE(J)
         ALOGG=LOGGTE(J)
         IMON=MONTE(J)
         IDATE=DATETE(J)
C
         IF(ex)THEN
          ILAT=ALAT+91.001
          IF(MODEL(1:2).EQ.'dk')THEN
            ILAT=(ALAT+90.0)*float(NLAOC)/181.0+0.501
          END IF
         ELSE
          ILAT=ALAT+90.501
         END IF
C
         ILOGG=ALOGG*float(NLOOC)/361.0+0.5001
         IF(ILOGG.LT.1)ILOGG=ILOGG+NLOOC
         IF(ILOGG.GT.NLOOC)ILOGG=ILOGG-NLOOC
         ILATP=ILAT+1
         ILATP=MIN(ILATP,NLAOC)
         ILOGGP=ILOGG+1
         IF(ILOGGP.GT.NLOOC)ILOGGP=ILOGGP-NLOOC
C
         IVLAT=1+(90.-ALAT)
         IVLATP=IVLAT+1
         IVLATP=MIN(IVLATP,181)
         IVLOGG=ALOGG+1.
         IF(IVLOGG.LT.1)IVLOGG=IVLOGG+360
         IF(IVLOGG.GT.360)IVLOGG=IVLOGG-360
         IVLOGGP=IVLOGG+1
         IF(IVLOGGP.GT.360)IVLOGGP=IVLOGGP-360
C
       KLAT=1+4.*(ALAT+90.)
       IF(ALOGG.GT.0.0)THEN
        KLOGG=1+4.*(ALOGG+0.01)
       ELSE
        KLOGG=1+4.*(360.01+ALOGG)
       END IF
       KLATP=KLAT+1
       KLOGGP=KLOGG+1
       IF(KLOGG.GT.1440)KLOGG=KLOGG-1440
       IF(KLOGGP.GT.1440)KLOGGP=KLOGGP-1440
C
          IF(IDATE.GE.15)THEN
           IMONK=IMON
           IMONKP=IMON+1
           IF(IMONKP.EQ.13)IMONKP=1
           FAC=FLOAT(IDATE-15)/30.
          ELSE
           IMONK=IMON-1
           IF(IMONK.EQ.0)IMONK=12
           IMONKP=IMON
           FAC=FLOAT(15+IDATE)/30.
          END IF
          DE1=(1.-FAC)*MIXDEPTH(IMONK,ILAT,ILOGG)+FAC*MIXDEPTH(IMONKP,
     1    ILAT,ILOGG)
          DE2=(1.-FAC)*MIXDEPTH(IMONK,ILATP,ILOGG)+FAC*MIXDEPTH(IMONKP,
     1    ILATP,ILOGG)
          DE3=(1.-FAC)*MIXDEPTH(IMONK,ILAT,ILOGGP)+FAC*MIXDEPTH(IMONKP,
     1    ILAT,ILOGGP)
          DE4=(1.-FAC)*MIXDEPTH(IMONK,ILATP,ILOGGP)+FAC*MIXDEPTH(
     1    IMONKP,ILATP,ILOGGP)
          DTH1=(1.-FAC)*THERMAL(IMONK,ILAT,ILOGG)+FAC*THERMAL(IMONKP,
     1    ILAT,ILOGG)
          DTH2=(1.-FAC)*THERMAL(IMONK,ILATP,ILOGG)+FAC*THERMAL(IMONKP,
     1    ILATP,ILOGG)
          DTH3=(1.-FAC)*THERMAL(IMONK,ILAT,ILOGGP)+FAC*THERMAL(IMONKP,
     1    ILAT,ILOGGP)
          DTH4=(1.-FAC)*THERMAL(IMONK,ILATP,ILOGGP)+FAC*THERMAL(
     1    IMONKP,ILATP,ILOGGP)
         BA1=BATH(KLOGG,KLAT)
         BA2=BATH(KLOGG,KLATP)
         BA3=BATH(KLOGGP,KLAT)
         BA4=BATH(KLOGGP,KLATP)
         CA1=CDRAG(KLOGG,KLAT)
         CA2=CDRAG(KLOGG,KLATP)
         CA3=CDRAG(KLOGGP,KLAT)
         CA4=CDRAG(KLOGGP,KLATP)
c
c   Calculate landfraction, linearly grading to about 3 meters elevation
c
         LANDFRAC(J)=0.25*(MIN(MAX(0.33*BA1,0.0),1.0)+
     1    MIN(MAX(0.33*BA2,0.0),1.0)+MIN(MAX(0.33*BA3,0.0),1.0)+
     2    MIN(MAX(0.33*BA4,0.0),1.0))
c
c  ***   These 8 lines added 2/20/07 after land mask removed from calculating VPOT  ***
c  ***      Now VPOT uses skin temperature over land. These lines zero out VPOT     ***
c  ***                                over land                                     ***
c
c  ***   Forgot that bathymetry points here need to coincide with potential intensity points.  ***
c  ***                    More lines added 12/11/2008 to accomplish this                       ***
c
          ALATI=ALAT+90.
          DELY=ALATI-INT(ALATI)
          COSFAC=COS(ALAT*ADEG)
          ALOGGI=ALOGG+180.
          delx=ALOGGI-int(ALOGGI)
          D1=(1.-DELX)*(1.-DELY)
          D2=DELY*(1.-DELX)
          D3=DELX*(1.-DELY)
          D4=DELX*DELY
          HMIX0(J)=D1*DE1+D2*DE2+D3*DE3+D4*DE4
          STRAT(J)=D1*DTH1+D2*DTH2+D3*DTH3+D4*DTH4
          STRAT(J)=MAX(STRAT(J),0.1)
C
          DELY=1.+4.*(ALAT+90.)-FLOAT(KLAT)
          DELX=1.+4.*(360.+ALOGG)-INT(1.+4.*(360.+ALOGG))
          D1=(1.-DELX)*(1.-DELY)
          D2=DELY*(1.-DELX)
          D3=DELX*(1.-DELY)
          D4=DELX*DELY
          BATHTRACK(J)=D1*BA1+D2*BA2+D3*BA3+D4*BA4
          CDTRACK(J)=D1*CA1+D2*CA2+D3*CA3+D4*CA4
c
          if(BATHTRACK(J).LT.0.0)CDTRACK(J)=8.0e-4 ! Added 1/19 to account for smoothing of ECMWF CD file
c
          IF(FLAV(1:2).EQ.'eo'.OR.FLAV(1:2).EQ.'pe')THEN
           BATHTRACK(J)=-1000.0
          END IF
c
          HMIX0(J)=MIN(HMIX0(J),-BATHTRACK(J))
          HMIX0(J)=MAX(HMIX0(J),20.0)
          UTEMP=UTIME(J)/(24.*3600.)
c
        END DO
C
        IF(gmeth.ne.'clim')THEN
C
C  Added these lines on August 2 2014 to do a partial calculation of GPI to use as a filter
C
        Tc=T600TRACK(1)-273.15
        est=6.112*exp(17.67*Tc/(243.5+Tc))
        qst=0.622*est/(600.0-est)
        pin=VPTRACK(1)
        chit=2.5e6*qst*(1.0-0.01*RHTRACK(1))/MAX(pin*pin,1.0)
        chit=MAX(chit,0.01)
        gpi=(chit**-1.333)*(max((pin-35.0),0.0))**2/
     1    ((1.0+0.04*SHEARTRACK(1))**4)
c
        IF(BATHTRACK(1).GE.0.0.OR.gpi.LT.100.0)THEN
c
c  End of 8/2/2014 additions/changes
c
         GOTO 5
c
        END IF
C
        END IF
C
C      *** Non-dimensionalize Parameters and Initial Conditions  ***
c
        TIME=UTIME(NSTORM)
C
        UT=(DISTANCE(2)-DISTANCE(1))/(UTIME(2)-UTIME(1))
c
        DO 33 I=1,NSTORM
         DISTANCE(I)=DISTANCE(I)
         UTIME(I)=UTIME(I)
   33   CONTINUE
C
C            ****   INITIALIZE FIELDS   ***
c
c  ***  Make random draw from a normal distribution
c
 9056       call random_number(arand)
            av1=2.*arand-1.
            call random_number(arand)
            av2=2.*arand-1.
            avr=av1**2+av2**2
            if(avr.ge.1.0)goto 9056
            avr=max(avr,1.0e-8)
            avfac=sqrt(-2.*log(avr)/avr)
            Gasdev=av2*avfac
            if(Gasdev.lt.0.0)goto 9056
c
        VMINIT=VMINIT0+5.0*AKNOTFAC*Gasdev
        VMINIT=0.5*VMINIT  !  Added May, 2014 to compensate for other changes made below. Changed to 0.5 10/2014
        VMINIT=MAX(VMINIT, 9.0*AKNOTFAC)
        VMINIT=MIN(VMINIT, 18.0*AKNOTFAC)
c
        VMINIT=VMINIT*1.2
c
c   Added 4/9/2007 to weight genesis points by variance of 850 hPa wind
c
c        VMINIT=VMINIT+1.0*SQRT(UVAR850TRACK+VVAR850TRACK)
c
c   Added May 2014 to account for changing deformation radii with climate
c
        TM=T600TRACK(1)-273.15
        ESM=6.112*EXP(17.67*TM/(243.5+TM))
        QSM=0.622*ESM/600.0
        reducfac=1.25*sqrt(QSM/0.01)   !  Deformation radius factor for initial radius of maximum winds
        VM=VMINIT   !   Added May, 2014. Omitting this earlier was an error. 
c
        TT=0.0D0
        TT1=0.0D0
        NTT=0
C
C    *** Initialize toy model   ***
C
	    V1=VM
	    V2=VM
            XM1=0.2
	    XM2=XM1
C
C             ***  SET TIME LOOPING PARAMETERS AND BEGIN TIME LOOP  ***
C
        ADIST=0.0
        DT=DT0
        DT1=DT
        TIMEPRINT=12.0
        NEXT=0
        NTCOUNT=0
        alatsold=latte(1)
        alogsold=loggte(1)
c
c  Switched to 8.0 on 5/31/05
c
        VTHRESH=10.0
        VTHRESH2=0.5
        TTHRESH=2.*24.*3600.
        VMAX3=2.*VTHRESH
        VMAX2=VMAX3
C
C           ***  PROGRAM RETURNS TO 77 AFTER EACH TIME STEP  ***
C
   77 CONTINUE
C
        TT=TT+DT
c
c   Switched from VMAX3 to VMAX2 in next statement, 5/20/05
c
        IF(TT.GT.TIME.OR.(V2.LT.VTHRESH.AND.TT.GT.TTHRESH))GOTO 705
C
C
C  **  Calculate certain along-track quantities
C
        KU=NSTORM
        DO 80 J=2,NSTORM
         IF(UTIME(J).GT.TT)KU=MIN(KU,J) 
   80   CONTINUE
        ADIST=DISTANCE(NSTORM)
        VPOTM=VPTRACK(NSTORM)
        RHM=RHTRACK(NSTORM)
        T600M=T600TRACK(NSTORM)
        VSHEAR=SHEARTRACK(NSTORM)
        U850T=U850TR(NSTORM)
        V850T=V850TR(NSTORM)
        CDY=CDTRACK(NSTORM)
        BATHY=BATHTRACK(NSTORM)
        BATHFRAC=LANDFRAC(NSTORM)
        ACCEL=0.0
        DUTINV=1./(UTIME(KU)-UTIME(KU-1))
        IF(KU.LE.NSTORM)THEN
         ADIST=(DISTANCE(KU)*(TT-UTIME(KU-1))+DISTANCE(KU-1)*(UTIME(KU)-
     1    TT))*DUTINV
         VPOTM=(VPTRACK(KU)*(TT-UTIME(KU-1))+VPTRACK(KU-1)*(UTIME(KU)-
     1    TT))*DUTINV
         RHM=(RHTRACK(KU)*(TT-UTIME(KU-1))+RHTRACK(KU-1)*(UTIME(KU)-
     1    TT))*DUTINV
         T600M=(T600TRACK(KU)*(TT-UTIME(KU-1))+T600TRACK(KU-1)*
     1    (UTIME(KU)-TT))*DUTINV
         VSHEAR=(SHEARTRACK(KU)*(TT-UTIME(KU-1))+SHEARTRACK(KU-1)*
     1    (UTIME(KU)-TT))*DUTINV
         U850T=(U850TR(KU)*(TT-UTIME(KU-1))+U850TR(KU-1)*
     1    (UTIME(KU)-TT))*DUTINV
         V850T=(V850TR(KU)*(TT-UTIME(KU-1))+V850TR(KU-1)*
     1    (UTIME(KU)-TT))*DUTINV
c
         vstar=vsfac*sqrt(U850T**2+V850T**2)
c
         BATHY=(BATHTRACK(KU)*(TT-UTIME(KU-1))+BATHTRACK(KU-1)*
     1    (UTIME(KU)-TT))/(UTIME(KU)-UTIME(KU-1))
         BATHFRAC=(LANDFRAC(KU)*(TT-UTIME(KU-1))+LANDFRAC(KU-1)*
     1    (UTIME(KU)-TT))/(UTIME(KU)-UTIME(KU-1))
         CDY=(CDTRACK(KU)*(TT-UTIME(KU-1))+CDTRACK(KU-1)*
     1    (UTIME(KU)-TT))/(UTIME(KU)-UTIME(KU-1))
         IF(KU.GT.2.AND.KU.LE.(NSTORM-1))THEN
          ACCEL=((DISTANCE(KU+1)-DISTANCE(KU-1))/(UTIME(KU+1)-
     1     UTIME(KU-1))-(DISTANCE(KU)-DISTANCE(KU-2))/(UTIME(KU)-
     2     UTIME(KU-2)))*DUTINV
         END IF
        END IF
c
        RHM=MIN(RHM,100.0)
        RHM=MAX(RHM,1.0)
        TM=T600M-273.15
        TMA=T600M
        ESM=6.112*EXP(17.67*TM/(243.5+TM))
        QSM=0.622*ESM/600.0
        SMEX=2.5E6*QSM*(0.01*RHM-1.)-461.0*TMA*QSM*0.01*RHM*
     1    LOG(0.01*RHM)
        XM0=(TS-TO)*SMEX/(CHI*TMA)
        VPBND=MAX(VPOTM,1.0)
        XM0=-(TS-TO)*SMEX/(VPBND*VPBND*TOA)
        XM0=XM0/(1.+XM0)
c        print*, XM0
        XM0=0.01*RHM
        IF(TT.LT.1.5*DT)THEN
	  XM1=1.2*XM0
          XM1=MIN(XM1,1.0)
c          XM1=0.75
          XM2=XM1
        END IF
c
        ACCEL=0.0
        UT0=(DISTANCE(KU)-DISTANCE(KU-1))*DUTINV
        UT=UT0+ACCEL*(TT-0.5*(UTIME(KU)+UTIME(KU-1)))
C
        VMAX2=V2
        VMAX4=MAX(V2,30.0)
        VMAX=VMAX+TRANSFAC*UT
        VMAX3=V2
C
        KH1=NSTORM
        DO 991 J=2,NSTORM
         IF(DISTANCE(J).GT.ADIST)KH1=MIN(KH1,J)
  991   CONTINUE
        AHMLOCAL2=HMIX0(NSTORM)
        SLOCAL=STRAT(NSTORM)
        IF(KH1.LT.NSTORM)THEN
         AHMLOCAL2=(HMIX0(KH1)*(ADIST-DISTANCE(KH1-1))+HMIX0(KH1-1)*
     1    (DISTANCE(KH1)-ADIST))/(DISTANCE(KH1)-DISTANCE(KH1-1))
         SLOCAL=(STRAT(KH1)*(ADIST-DISTANCE(KH1-1))+STRAT(KH1-1)*
     1    (DISTANCE(KH1)-ADIST))/(DISTANCE(KH1)-DISTANCE(KH1-1))
        END IF
        SLOCAL=MAX(SLOCAL,0.1)
c
c  ***  Modify potential intensity using a variant of Schade's ocean feedback formulation. 12/14/2016 ***
c
        ZFAC=0.01*(SLOCAL)**(-0.4)*AHMLOCAL2*ABS(UT)*VPOTM/(0.1+V2)
        ZFAC=MAX(ZFAC,0.0)
        ALPHA=1.-0.87*EXP(-ZFAC)
C
C  ***  Calculate forcings for toy model ***
C
        vpnew2=alpha*1.4*vpotm*vpotm
c        
c  ***  Calculate maximum surface wind speed including background wind
c
        VABS=SQRT(V1*V1+vstar*vstar)
        aff=2.*(vstar/90)*V1/VABS
        dum=1-aff*aff/16.-15.*aff**4/(16.*64.)
        VABS=VABS*dum
c        aff=vstar/(V1+1.0e-6)
c         IF(aff.gt.1.0)then
c           VABS=V1*(1.0+0.47*(aff-1.0)**1.12)
c         END IF
c
        CD1=CDY
c
c   ***  Reduce drag coefficient at very high winds over ocean
c
        CDV1=0.0
        IF(BATHY.LT.0.0)CDV1=CDV
        IF(BATHY.LT.0.0.AND.VABS.GE.50.0)THEN
          CD1=0.7*CDY
        END IF
c
        CVABS=CD1*(1.+MIN(CDV1*VABS,CDMAX))
        CKABS=CD1*(1.+MIN(CDV1*VABS,CKMAX))
        IF(VABS.GT.35.0)CKABS=CKABS*1.20
        CRAT=(1.-BATHFRAC)*CKABS/CVABS
        CDBATH=CDY/8.0E-4
c
c        VFORC=(0.67-0.1*alpha)*CVABS*HATMI*(CRAT*vpnew2*
c     1   XM2*XM2*XM2-CDBATH*V1*V1)
        VFORC=CVABS*HATMI*((0.67-0.1)*CRAT*vpnew2*XM2*XM2*XM2
     1   -(CDBATH-XM2*XM2*XM2*(0.33+0.1*ALPHA))*V1*V1)

        XMFORC=CVABS*CDBATH*HATMI*((1.-XM1)*VABS-2.2*VSHEAR*XM1)
C
C  ***  Print progress report to screen  ***
C
        IF(TT.GT.TIMEPRINT)THEN
c         PRINT*, ' Finished  ',TIMEPRINT,'  hours'
         TIMEPRINT=TIMEPRINT+12.0
        END IF
C
C  Interpolate latitude and longitude to current time
C
        ALATS=LATTE(NSTORM)
        ALOGS=LOGGTE(NSTORM)
        IF(KU.LE.NSTORM)THEN
         ALATS=(LATTE(KU)*(TT-UTIME(KU-1))+LATTE(KU-1)*(UTIME(KU)-
     1    TT))*DUTINV
         ALOGS=(LOGGTE(KU)*(TT-UTIME(KU-1))+LOGGTE(KU-1)*(UTIME(KU)-
     1    TT))*DUTINV
        END IF
c
c   Find maximum wind within city radius
c
        VMTEMP1=V2+0.8*UT
        VMTEMP2=V2*3600./1852.
        VMTEMP1=VMTEMP1*3600./1852.
c
        DELX=AFAC2*COS(ALATS*ADEG)*(ALOGS-CLOG1)
        DELY=AFAC2*(ALATS-CLAT1)
        DISTA=0.001*SQRT(DELX*DELX+DELY*DELY)
c
c    Statment below added 8/1/2009
c
c        IF(VMTEMP1.GE.VMCRIT.AND.VMTEMP2.GE.25.0)IFLAGB=1
c
         IF(((shape.eq.'circ'.and.DISTA.LE.CRAD).or.jint.eq.-1).AND. 
     1    VMTEMP1.GE.VMCRIT.AND.VMTEMP2.GE.25.0)THEN
         IFLAG=1
        END IF
c
        if(shape.eq.'poly'.and.mpoly.eq.0)then
          xtemp(1)=alogs
          xtemp(2)=alogsold
          ytemp(1)=alats
          ytemp(2)=alatsold
          alatsold=alats
          alogsold=alogs
          CALL box(xtemp,ytemp,2,x1vert,y1vert,x2vert,y2vert,
     1     nvert,nmap,xint2,yint2,jint2)
          if(jint2.gt.0.and.VMTEMP1.GE.VMCRIT.AND.VMTEMP2.GE.25.0)then
            iflag=1
          end if
         end if
c
        IF(mpoly.eq.1.and.jint.gt.0)THEN
          n1a=1
          do j=1,nvert
           av1=(alats-y1vert(j))*(x2vert(j)-x1vert(j))
           bv1=(alogs-x1vert(j))*(y2vert(j)-y1vert(j))
           if(av1.gt.bv1)n1a=0
          end do
          if(n1a.eq.1.and.vmtemp1.ge.vmcrit.and.vmtemp2.ge.25.0)then
           iflag=1
          end if
        END IF
C
C          ***  STORE CERTAIN GRAPHICS ARRAYS EVERY DTG TIME UNITS  ***
C
        TT1=TT1+DT
        IF(ABS(TT1-DTG).GT.(0.5*DT).AND.TT.GT.(1.5*DT))GOTO 500
        TT1=0
        NTT=NTT+1
C
C   ***   Create time series    ***
C
         TIMEG(NTT)=TT
         TIMEG(NTT)=TIMEG(NTT)/(24.*3600.)
         LATTIME(NTT)=ALATS
         LOGTIME(NTT)=ALOGS
         IF(LOGTIME(NTT).LT.0.0)LOGTIME(NTT)=LOGTIME(NTT)+360.0
         IF(LOGTIME(NTT).GT.360.0)LOGTIME(NTT)=LOGTIME(NTT)-360.0
c
         UTRANS(NTT)=AFAC2*COS(ADEG*LATTE(KU))*(LOGGTE(KU)-
     1     LOGGTE(KU-1))*DUTINV
         VTRANS(NTT)=AFAC2*(LATTE(KU)-LATTE(KU-1))*DUTINV
C
         VMG(NTT)=V2
         VMG(NTT)=3600.*VMG(NTT)/1852.
         VMGP=V2
C
         VSHG(NTT)=VSHEAR
         U850G(NTT)=U850T*3600./1852.
         V850G(NTT)=V850T*3600./1852.
         VPOTG(NTT)=VPOTM*3600./1852.
         T600G(NTT)=T600M
         RHG(NTT)=RHM
c
c  *** Algorithm to deduce radius of maximum winds from maximum wind speed ***
c
         FC1=1.45E-4*ABS(SIN(ALATS*0.0175)) ! Coriolis parameter
         ro=400.0        ! outer radius in km
         wc=3.0          ! Radiative subsidence speed in mm/s
         cdouter=1.2E-3  ! Outer region drag coefficient
         drouter=0.002   ! Nondimensional radial increment
c
         nouter=1+int(ro*1000.*FC1/(V2*drouter))
c
         allocate(router(nouter),vouter(nouter))
c
         CALL VOUTERNEW(V2,FC1,ro,wc,cdouter,nouter,vouter,router,imin)
         epa=(vouter(imin)-V2)/(vouter(imin)-vouter(imin+1))
         RMG(NTT)=router(imin)+epa*(router(imin+1)-router(imin))
c
         deallocate(router(nouter),vouter(nouter))
c
c
c  *** Formula for central pressure from quadratic curve fit to full CHIPS output  ***
c
         PCG(NTT,1)=1012.0-0.0095*V2*V2-0.7*V2
c
c
  500 CONTINUE
C
C            ***  ADVANCE VARIABLES WITH TIME SMOOTHER ***
C
        V3=V1+2.*DT*VFORC
        V1=V2+DAMP*(V1+V3-2.*V2)
        V2=V3
C
        XM3=XM1+2.*DT*XMFORC
        XM1=XM2+DAMP*(XM1+XM3-2.*XM2)
        XM1=MIN(XM1,1.0)
	XM1=MAX(XM1,0.0)
        XM2=XM3 
        XM2=MIN(XM2,1.0)
        XM2=MAX(XM2,0.0)
C
      GOTO 77
  705 CONTINUE
c
        IF(IFLAG.EQ.0)THEN
          call random_number(gdum)
c
c  Statements added 8/1/2009
c
c         IF(IFLAGB.EQ.0.AND.ICOUNT.LE.20)THEN
c           BACKSPACE(18)
c           ICOUNT=ICOUNT+1
c         ELSE
c           ICOUNT=0
c         END IF
         GOTO 7
        END IF
C
C   ***     Write to output files     ***
C
        IMON=IMONI
        IDATE=IDATEI
        ITIME=ITIMEI
c
        NUMSTORMS=NUMSTORMS+1
c
        if(shape.eq.'poly')then
          write(15,*)NUMSTORMS,clog1,clat1,jint
        end if
c
        if(numstorms.lt.10)then
         write(ch1,1201)numstorms
         DIRECTORY=TRIM(CITYNAME)//'/hurr'//ch1//'.out'
c
        elseif(numstorms.ge.10.and.numstorms.lt.100)then
         write(ch2,1202)numstorms
         DIRECTORY=TRIM(CITYNAME)//'/hurr'//ch2//'.out'
c
        elseif(numstorms.ge.100.and.numstorms.lt.1000)then
         write(ch3,1203)numstorms
         DIRECTORY=TRIM(CITYNAME)//'/hurr'//ch3//'.out'
c
        elseif(numstorms.ge.1000.and.numstorms.lt.10000)then
         write(ch4,1204)numstorms
         DIRECTORY=TRIM(CITYNAME)//'/hurr'//ch4//'.out'
c
        elseif(numstorms.ge.10000.and.numstorms.lt.100000)then
         write(ch5,1205)numstorms
         DIRECTORY=TRIM(CITYNAME)//'/hurr'//ch5//'.out'
c
        end if
c
 1201   format(I1.1)
 1202   format(I2.2)
 1203   format(I3.3)
 1204   format(I4.4)
 1205   format(I5.5)
c
        OPEN(UNIT=12,FILE=TRIM(DIRECTORY),STATUS='UNKNOWN')
C
        DO I=1,NTT
C
          IF(I.GT.1)THEN
           IDATE=IDATEI+TIMEG(I)+FLOAT(ITIMEI)/24.
           ITIME=ITIMEI+24*(TIMEG(I)-INT(TIMEG(I)))
           IF(ITIME.GE.24)ITIME=ITIME-24
           IMON=IMONI
           IMOB=IMON
           IDATB=IDATE
           IF((IMOB.EQ.1.OR.IMOB.EQ.3.OR.IMOB.EQ.5.OR.IMOB.EQ.7.OR.
     1     IMOB.EQ.8.OR.IMOB.EQ.10.OR.IMOB.EQ.12).AND.IDATB.GT.31)THEN
                  IDATE=IDATE-31
              IMON=IMON+1
            END IF
            IF((IMOB.EQ.4.OR.IMOB.EQ.6.OR.IMOB.EQ.9.OR.IMOB.EQ.11)
     1                .AND.IDATB.GT.30)THEN
                  IDATE=IDATE-30
              IMON=IMON+1
            END IF
            IF(IMOB.EQ.2.AND.IDATB.GT.28)THEN
              IDATE=IDATE-28
              IMON=IMON+1
            END IF
            IF(IMON.GT.12)IMON=1
         END IF
C
          WRITE(12,810)IMON,IDATE,ITIME,LATTIME(I),LOGTIME(I),
     1    VMG(I),PCG(I,1),RMG(I),VSHG(I),VPOTG(I),U850G(I),V850G(I),
     2    T600G(I),RHG(I),VMGSE(I),RMGSE(I)
c
          END DO
c
        CLOSE(12)
C
  810 FORMAT(3(1X,I2),2(1X,F8.3),11(1X,F7.2))
c
        if(gmeth.ne.'clim')then
c
         SFREQ=basinfac(ibas)*globalfac*areafac*
     1     float(numstorms)/float(ntotal)
c
        else
c
         SFREQ=1.1*climfac*float(numstorms)/float(ntotal)
c
        end if
c
         OPEN(UNIT=12,FILE=TRIM(CITYNAME)//'/stats'//filenum//'.out',
     1    status='unknown')
c         write(12,1222)BAS,' ',gmeth,' ',shape,' ',TRIM(CITYNAME),CLAT,
c     1    CLOG,CRAD,numstorms,ntotal,SFREQ
         write(12,1222)BAS,gmeth,shape,TRIM(CITYNAME),CLAT,
     1    CLOG,CRAD,numstorms,ntotal,SFREQ
         CLOSE(12)
 1222    FORMAT(1X,A2,1X,A4,1X,A4,1X,A40,1X,F9.4,1X,F9.4,1X,F9.1,1X,
     1    I5,1X,I9,1X,F12.6)
C
        IF(NUMSTORMS.GE.NUMTRACKS)THEN
c
        if(shape.eq.'circ')then
         write(15,*)'1  ',clog1,clat1,jint
        end if
c
         OPEN(UNIT=12,FILE=TRIM(CITYNAME)//'/stats'//filenum//'.txt',
     1    status='unknown')
         write(12,'(1X,A35,T45,A)')'Basin:',BAS
         write(12,'(1X,A35,T45,A)')'Filename:',TRIM(CITYNAME)
         write(12,'(1X,A35,T45,A)')'Seeding method (random or climo):'
     1     , gmeth
         write(12,'(1X,A35,T45,A)')'Filter type:', shape
         write(12,'(1X,a35,T41,F9.4)')'POI Latitude:',CLAT
         write(12,'(1X,a35,T41,F9.4)')'POI Longitude:',CLOG
         write(12,'(1X,A35,T42,F8.1)')'Search radius (km):',CRAD
         write(12,'(1X,A35,T45,F5.1)')'Critical wind speed (kts):',
     1     VMCRIT
         write(12,'(1X,A35,T44,I6)')'Number of tracks:',numstorms
         write(12,'(1X,A35,T44,I8)')'Number of unfiltered tracks:',
     1     ntotal
         write(12,'(1X,A35,T44,F8.2)')'Annual frequency:',sfreq
         CLOSE(12)
c
         STOP
c
        ELSE
         GOTO 5
        END IF
C
        close(15)
c
        END
C
C--------------------------------------------------------------------
C
       SUBROUTINE THEORY(CKCD,H,TS,TO,PA,VM,PC,IFAIL)
C
C       This subroutine calculates the theoretical maximum wind
C         speed and minimum central pressure.
C
       REAL H, LV, PA, PC, PM
       DELTAT=0.0    
       LV=2.5E6
       RD=287.0
       RV=461.0 
       IFAIL=0
C
C      AN is the assumed power dependence of v on radius inside the
C        radius of maximum winds (i.e. v~r**an) used to calculate PC
C
       AN=1.5
C       
        ES=6.112*EXP(17.67*TS/(243.5+TS))
        EP=(TS-TO)/(TO+273.15)
        COEF1=EP*LV*ES/(RV*(TS+273.15))
        COEF2=0.5*CKCD*(1.-H)*(EP+1.)
        COEF3=0.5*CKCD*EP*1000.*DELTAT/(RD*(TS+273.15)*(1.-EP*H))
        PM=PA   
        N=0
   10   CONTINUE
        N=N+1
        PG=PA*EXP(-COEF1*(COEF2/PM+H*(1./PM-1./PA))-COEF3)
        IF(ABS(PG-PM).LT.0.1)THEN
         PM=0.5*(PM+PG)
         VM=SQRT(EP*CKCD*(LV*0.622*ES*(1.-H)/PM+1000.*DELTAT))
         GOTO 20
        END IF 
        IF(N.GT.1000.OR.PG.LE.1.0)THEN
          IFAIL=1
          GOTO 20
        END IF
        PM=0.5*(PM+PG)  
        GOTO 10
   20   CONTINUE
        IF(IFAIL.EQ.1)THEN
         PC=PA
         VM=0.0
        ELSE
         PC=PM*EXP(-VM*VM/(2.*AN*RD*(TS+273.15)))
        END IF
C
      RETURN
      END         
C      
c
        SUBROUTINE choldc(a,n,np,p)
        INTEGER n,np
        REAL*8 a(np,np),p(n)
c
c        Given a positive-definite symmetric matrix a(1:n,1:n), with physical dimension np, this
c        routine constructs its Cholesky decomposition, A = L·LT . On input, only the upper triangle
c        of a need be given; it is not modified. The Cholesky factor L is returned in the lower triangle
c        of a, except for its diagonal elements which are returned in p(1:n).
c
        INTEGER i,j,k
        REAL*8 sum
        do  i=1,n
         do  j=i,n
          sum=a(i,j)
          do  k=i-1,1,-1
           sum=sum-a(i,k)*a(j,k)
          enddo 
          if(i.eq.j)then
           if(sum.le.0.)then
c            print*, 'Cholesky failed with sum = ', sum, 'at', i 
            p(i)=0.0
            return
           end if
           p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
         enddo 
        enddo 
        return
        END
c
        Subroutine Box(xtrack,ytrack,
     1    nd,xa,ya,xb,yb,m,nint,xint,yint,jint)
c
c    This subroutine determines whether a track of dimension nd (given by coordinates xtrack(nd),
c     ytrack(nd)) intersects a series of m-1 line segements whose end points are given by xa(m), ya(m),
c     xb(m), yb(m). Nint=0 means no intersection; nint=1 means intersection. The last point 
c     may equal the first, in which case one has a closed polygon. In this case, the routine also checks
c     for tracks that lie entirely within the polygon. For this reason, the segments must be specified
c     going CLOCKWISE around the polygon.
c
c     The quantities xint and yint constitute the first intersection point, along the 
c     line segment numbered jint. In the case of a closed polygon, if the track is entirely
c     within the polygon, jint is set equal to -1 and xint=yint=0.0.
c
c               K. Emanuel  1/6/2008 revised 1/18/2008
c
        real xtrack(nd), ytrack(nd), xa(m), ya(m), xb(m), yb(m)
        integer nint,im
c
        nint=0
        jint=0
        xint=0.0
        yint=0.0
c
        do i=1,nd-1
         x1=xtrack(i)
         x2=xtrack(i+1)
         y1=ytrack(i)
         y2=ytrack(i+1)
c
         do j=1,m
          x3=xa(j)
          y3=ya(j)
          x4=xb(j)
          y4=yb(j)
c
          denom=(y4-y3)*(x2-x1)-(x4-x3)*(y2-y1)
          if(abs(denom).gt.1e-8)then
           ua=((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))/denom
           ub=((x2-x1)*(y1-y3)-(y2-y1)*(x1-x3))/denom
           if(ua.ge.0.0.and.ua.le.1.0.and.ub.ge.0.0.and.ub.le.1.0)then
            nint=1
            xint=x1+ua*(x2-x1)
            yint=y1+ua*(y2-y1)
            jint=j
            return
           end if
          end if
         end do
        end do
c
        if(nint.eq.0.and.abs(xb(m)-xa(1)).lt.0.01.and.abs(yb(m)-
     1   ya(1)).lt.0.01)then
          call pnpoly(xtrack(1),ytrack(1),xa,ya,m,INOUT)
          if(INOUT.eq.1)then
            nint=1
            jint=-1
          end if
        end if
c
        return
        end
C
       SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)                            
C                                                                       PNP1  20
C     ..................................................................PNP1  30
C                                                                       PNP1  40
C        SUBROUTINE PNPOLY                                              PNP1  50
C                                                                       PNP1  60
C        PURPOSE                                                        PNP1  70
C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            PNP1  80
C                                                                       PNP1  90
C        USAGE                                                          PNP1 100
C           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     PNP1 110
C                                                                       PNP1 120
C        DESCRIPTION OF THE PARAMETERS                                  PNP1 130
C           PX      - X-COORDINATE OF POINT IN QUESTION.                PNP1 140
C           PY      - Y-COORDINATE OF POINT IN QUESTION.                PNP1 150
C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         PNP1 160
C                     VERTICES OF POLYGON.                              PNP1 170
C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           PNP1 180
C                     VERTICES OF POLYGON.                              PNP1 190
C           N       - NUMBER OF VERTICES IN THE POLYGON.                PNP1 200
C           INOUT   - THE SIGNAL RETURNED:                              PNP1 210
C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        PNP1 220
C                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     PNP1 230
C                      1 IF THE POINT IS INSIDE OF THE POLYGON.         PNP1 240
C                                                                       PNP1 250
C        REMARKS                                                        PNP1 260
C           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      PNP1 270
C           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           PNP1 280
C           OPTIONALLY BE INCREASED BY 1.                               PNP1 290
C           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      PNP1 300
C           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    PNP1 310
C           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   PNP1 320
C           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              PNP1 330
C           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         PNP1 340
C           THE SIZE OF THE ARRAYS X AND Y MUST BE INCREASED IF N > 20. PNP1 350
C           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   PNP1 360
C                                                                       PNP1 370
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PNP1 380
C           NONE                                                        PNP1 390
C                                                                       PNP1 400
C        METHOD                                                         PNP1 410
C           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  PNP1 420
C           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        PNP1 430
C           POINT IS INSIDE OF THE POLYGON.                             PNP1 440
C                                                                       PNP1 450
C     ..................................................................PNP1 460
C                                             
       REAL X(N),Y(N),XX(N),YY(N)                                    
       LOGICAL MX,MY,NX,NY                                               
       DO I=1,N                                                        
         X(I)=XX(I)-PX                                                     
         Y(I)=YY(I)-PY  
       End do                                                     
       INOUT=-1                                                          
       DO 2 I=1,N                                                        
       J=1+MOD(I,N)                                                      
       MX=X(I).GE.0.0                                                    
       NX=X(J).GE.0.0                                                    
       MY=Y(I).GE.0.0                                                    
       NY=Y(J).GE.0.0                                                    
       IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2       
       IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
       INOUT=-INOUT                                                      
       GOTO 2                                                           
    3  IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5                       
    4       INOUT=0                                                     
       RETURN                                                            
    5  INOUT=-INOUT                                                      
    2  CONTINUE                                                          
       RETURN                                                            
       END   
c-------------------------------------------------------------------------
c
	SUBROUTINE VOUTERNEW(vm,fc,ro,wc,CD,q,v,r,imin)
c
c Numerical integration of the outer wind profile from simple ODE
c Inputs: vm (maximum wind in m/s); f (Coriolis in s^-1), ro (outer radius
c in km), wc (radiative subsidence rate in mm/s), drag coefficient, and
c number of radial points, q.
c
c Outputs:  V (in m/s) versus radius r (im km); imin is minimum index for which V and r are defined
c
c Created January, 2017
c
        integer q,imin
        real v(q),m(q),r(q)
c
c----------------------------------------------------------------------
	assl=0.2     ! Asselin filter coefficient
c-----------------------------------------------------------------------
	ro=ro*1000.0 ! Convert to meters
	wc=wc*0.001  ! Convert to m/s
	chi=CD*fc*ro/wc  ! definition of chi
	rdim=vm/fc  
	rond=ro/rdim
        dr=rond/(q-1)
c
c Integrate outer wind ODE
c
        rnd=rond-dr
        m(q)=0.0
        v(q)=0.0
        r(q)=rond
	m(q-1)=0.25*(rond**2-rnd**2)
	v(q-1)=m(q-1)/(rond-dr)
        r(q-1)=rond-dr
c
	do i=q-2,1,-1
          r(i)=r(i+1)-dr
          m(i)=m(i+2)-2.*dr*(chi*m(i+1)**2/(rond**2-
     1     r(i+1)**2)-r(i+1))
          m(i+1)=m(i+1)+assl*(m(i)+m(i+2)-2.*m(i+1))
          v(i)=m(i)/r(i)
          if(v(i).gt.1.0)then  ! Stop integration when v exceeds potential intensity
            imin=i
            goto 10
          end if
        end do
c
   10   continue 
c
        do i=1,imin-1  ! Fill in values inside radius where v = potential intensity
          v(i)=0.0
          r(i)=dr*(i-1)
        end do
c
        do i=1,q  !  Re-dimensionalize
          v(i)=vm*v(i)               ! v in m/s
          r(i)=0.001*vm*r(i)/fc      ! r in km
        end do
c
        return
        end
