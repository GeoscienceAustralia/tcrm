      SUBROUTINE PCMIN(SST,PSL,P,T,R,IDISS,SIG,b,CKCD,N,PMIN,VMAX,IFL)
C
C   Revised on 9/24/2005 to fix convergence problems at high pressure
C
C This subroutine calculates the maximum wind speed
C and mimimum central pressure achievable in tropical cyclones, 
C given a sounding and a sea surface temperature.  
C
C  INPUT:   SST: Sea surface temperature in C
C
C           PSL: Sea level pressure (mb)
C
C           P,T,R: One-dimensional arrays of dimension NA
C             containing pressure (mb), temperature (C),
C             and mixing ratio (g/kg). The arrays MUST be
C             arranged so that the lowest index corresponds
C             to the lowest model level, with increasing index
C             corresponding to decreasing pressure. The temperature
C             sounding should extend to at least the tropopause and
C             preferably to the lower stratosphere, however the
C             mixing ratios are not important above the boundary
C             layer. Missing mixing ratios can be replaced by zeros.
C
C
C  OUTPUT:  PMIN is the minimum central pressure, in mb
C
C           VMAX is the maximum surface wind speed, in m/s
C                  (reduced to reflect surface drag)
C
C           IFL is a flag: A value of 1 means OK; a value of 0
C              indicates no convergence (hypercane); a value of 2
C              means that the CAPE routine failed.
C
C--------------------------------------------------------------------
C
C To build using f2py:
C f2py -c -m pcmin pcmin.f
C
C For more information on f2py, see http://www.scipy.org/F2py
C
C--------------------------------------------------------------------

       REAL T(N), P(N), R(N)
       INTEGER IDISS=1
       REAL SIG=1.0
       REAL b=2.0
       REAL CKCD=0.9

Cf2py intent(in) T, P, R
Cf2py intent(in) SST, PSL
Cf2py integer optional,intent(in) IDISS=1
Cf2py real optional,intent(in) SIG=1.0
Cf2py real optional,intent(in) b=2.0
Cf2py real optional,intent(in) CKCD=0.9
Cf2py intent(out) PMIN, VMAX, IFL

C   CKCD: Adjustable constant: Ratio of C_k to C_D:

C   Adjustable constant for buoyancy of displaced parcels:
C   SIG = 0 : Reversible ascent;
C   SIG = 1 : Pseudo-adiabatic ascent

C   Dissipative heating switch:
C   IDISS = 0 : No dissipative heating
C   IDISS = 1 : Dissipative heating allowed
C
C   Exponent, b, in assumed profile of azimuthal velocity in eye:
C   V=V_m(r/r_m)^b. Used only in calculation of central pressure
C
C Set level from which parcels lifted 
C
	NK=1
C
C Factor to reduce gradient wind to 10 m wind
C
	VREDUC=0.8
C
C Normalize certain quantities 
C
	SSTK=SST+273.15
	ES0=6.112*EXP(17.67*SST/(243.5+SST))
	DO 40 I=1,N
	 R(I)=R(I)*0.001
	 T(I)=T(I)+273.15
   40	CONTINUE
C
C Default values 
C
      VMAX=0.0
	PMIN=PSL
	IFL=1
C
	NP=0
	PM=950.0
C
C Find environmental CAPE
C
      TP=T(NK)
      RP=R(NK)
      PP=P(NK)
      CALL CAPE(TP,RP,PP,T,R,P,N,SIG,CAPEA,TOA,IFLAG)
      IF(IFLAG.NE.1)IFL=2
C
C Begin iteration to find mimimum pressure
C
  100 CONTINUE
C
C Find CAPE at radius of maximum winds
C
      TP=T(NK)
      PP=MIN(PM,1000.0)
      RP=0.622*R(NK)*PSL/(PP*(0.622+R(NK))-R(NK)*PSL)
      CALL CAPE(TP,RP,PP,T,R,P,N,SIG,CAPEM,TOM,IFLAG)
      IF(IFLAG.NE.1)IFL=2
      RAT=SSTK/TOM
      IF(IDISS.EQ.1)RAT=1.0
C
C Find saturation CAPE at radius of maximum winds 
C
      TP=SSTK
      PP=MIN(PM,1000.0)
      RP=0.622*ES0/(PP-ES0)
      CALL CAPE(TP,RP,PP,T,R,P,N,SIG,CAPEMS,TOMS,IFLAG)
      IF(IFLAG.NE.1)IFL=2
C
C Initial estimate of minimum pressure
C
      RS0=RP
      TV1=T(1)*(1.+R(1)/0.622)/(1.+R(1))
	TVAV=0.5*(TV1+SSTK*(1.+RS0/0.622)/(1.+RS0))
C	CAT=0.5*CKCD*RAT*(CAPEMS-CAPEM)
	CAT=CAPEM-CAPEA+0.5*CKCD*RAT*(CAPEMS-CAPEM)
	CAT=MAX(CAT,0.0)
	PNEW=PSL*EXP(-CAT/(287.04*TVAV))
C
C   ***  Test for convergence   ***
C
	IF(ABS(PNEW-PM).GT.0.2)THEN
	 PM=PNEW
	 NP=NP+1
	 IF(NP.GT.1000.OR.PM.LT.400.0)THEN
	  PMIN=PSL
	  IFL=0
	  GOTO 900
	 END IF
	 GOTO 100
	ELSE
	 CATFAC=0.5*(1.+1./b)
C	 CAT=CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
	 CAT=CAPEM-CAPEA+CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
	 CAT=MAX(CAT,0.0)
	 PMIN=PSL*EXP(-CAT/(287.04*TVAV))
	END IF
  900	CONTINUE
	FAC=MAX(0.0,(CAPEMS-CAPEM))
	VMAX=VREDUC*SQRT(CKCD*RAT*FAC)
C
C   ***  Renormalize sounding arrays   ***
C
	DO 910 I=1,N
	 R(I)=R(I)*1000.0
	 T(I)=T(I)-273.15
  910	CONTINUE
C
	RETURN
	END
C
      SUBROUTINE CAPE(TP,RP,PP,T,R,P,N,SIG,CAPED,TOB,IFLAG)
C
C This subroutine calculates the CAPE of a parcel with pressure 
C PP (mb), temperature TP (K) and mixing ratio RP (gm/gm) given 
C a sounding of temperature (T in K) and mixing ratio (R in gm/gm)
C as a function of pressure (P in mb). ND is the dimension of the
C arrays T,R and P, while N is the actual number of points in the 
C sounding. CAPED is the calculated value of CAPE and TOB is the 
C temperature at the level of neutral buoyancy.  IFLAG is a flag
C integer. If IFLAG = 1, routine is successful; if it is 0, routine did
C not run owing to improper sounding (e.g.no water vapor at parcel 
C level). IFLAG=2 indicates that routine did not converge.
C
       REAL T(N),R(N),P(N),TVRDIF(100), NA
       REAL SIG=1.0

Cf2py intent(in) T, R, P, TP, RP, PP
Cf2py optional,intent(in) SIG=1.0
Cf2py intent(out) CAPED,TOB,IFLAG

C
C Default values
C
      CAPED=0.0
      TOB=T(1)
      IFLAG=1
C
C Check that sounding is suitable
C
      IF(RP.LT.1.0E-6.OR.TP.LT.200.0)THEN
       IFLAG=0
       RETURN
      END IF
C
C Assign values of thermodynamic constants
C
      CPD=1005.7
      CPV=1870.0
C      CL=4190.0
      CL=2500.0
      CPVMCL=CPV-CL
      RV=461.5
      RD=287.04
      EPS=RD/RV
      ALV0=2.501E6
C
C Define various parcel quantities, including reversible entropy, S.
C
      TPC=TP-273.15
      ESP=6.112*EXP(17.67*TPC/(243.5+TPC))
      EVP=RP*PP/(EPS+RP)
      RH=EVP/ESP
      RH=MIN(RH,1.0)
      ALV=ALV0+CPVMCL*TPC
      S=(CPD+RP*CL)*LOG(TP)-RD*LOG(PP-EVP)+
     1   ALV*RP/TP-RP*RV*LOG(RH)
C
C Find lifted condensation pressure, PLCL 
C
	CHI=TP/(1669.0-122.0*RH-TP)
	PLCL=PP*(RH**CHI)
C
C Begin updraft loop
C
	NCMAX=0
	DO J=1,N
	    TVRDIF(J)=0.0
	END DO
C
	JMIN=1E6
	DO 200 J=1,N
C
C Don't bother lifting parcel above 60 mb and skip sections of 
C sounding below parcel level 
C
          IF(P(J).LT.59.0.OR.P(J).GE.PP)GOTO 200
C
	  JMIN=MIN(JMIN,J)
C
C Parcel quantities below lifted condensation level 
C
	  IF(P(J).GE.PLCL)THEN
	      TG=TP*(P(J)/PP)**(RD/CPD)
	      RG=RP
C
C Calculate buoyancy
C
	      TLVR=TG*(1.+RG/EPS)/(1.+RG)
	      TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
	   ELSE
C
C Parcel quantities above lifted condensation level 
C
	      TG=T(J)
	      TJC=T(J)-273.15
              ES=6.112*EXP(17.67*TJC/(243.5+TJC))
              RG=EPS*ES/(P(J)-ES)
C
C Iteratively calculate lifted parcel temperature and mixing 
C ratio for reversible ascent
C
	      NC=0
  120	  CONTINUE
	  NC=NC+1
C
C Calculate estimates of the rates of change of the entropy 
C with temperature at constant pressure
C
	  ALV=ALV0+CPVMCL*(TG-273.15)
	  SL=(CPD+RP*CL+ALV*ALV*RG/(RV*TG*TG))/TG
	  EM=RG*P(J)/(EPS+RG)
	  SG=(CPD+RP*CL)*LOG(TG)-RD*LOG(P(J)-EM)+
     1      ALV*RG/TG
	  IF(NC.LT.3)THEN
	   AP=0.3
	  ELSE
	   AP=1.0
	  END IF
	  TGNEW=TG+AP*(S-SG)/SL
C
C Test for convergence 
C
	  IF(ABS(TGNEW-TG).GT.0.001)THEN
	   TG=TGNEW
	   TC=TG-273.15
	   ENEW=6.112*EXP(17.67*TC/(243.5+TC))
C
C Bail out if things get out of hand
C
	   IF(NC.GT.500.OR.ENEW.GT.(P(J)-1.0))THEN
            IFLAG=2
            RETURN
	   END IF
	   RG=EPS*ENEW/(P(J)-ENEW)
	   GOTO 120
	  END IF
	  NCMAX=MAX(NC,NCMAX)
C
C Calculate buoyancy 
C
        RMEAN=SIG*RG+(1.-SIG)*RP
	  TLVR=TG*(1.+RG/EPS)/(1.+RMEAN)
	  TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
	 END IF
  200	CONTINUE
C
C Begin loop to find NA, PA, and CAPE from reversible ascent
C
	NA=0.0
	PA=0.0
C
C Find maximum level of positive buoyancy, INB  
C
	INB=1
	DO 550 J=N,JMIN,-1
	 IF(TVRDIF(J).GT.0.0)INB=MAX(INB,J)
  550	CONTINUE
	IF(INB.EQ.1)RETURN
C
C Find positive and negative areas and CAPE 
C
	IF(INB.GT.1)THEN
	 DO 600 J=JMIN+1,INB
	  PFAC=RD*(TVRDIF(J)+TVRDIF(J-1))*(P(J-1)-P(J))/(P(J)+P(J-1))
	  PA=PA+MAX(PFAC,0.0)
	  NA=NA-MIN(PFAC,0.0)
  600	 CONTINUE
C
C Find area between parcel pressure and first level above it
C
	PMA=(PP+P(JMIN))
	PFAC=RD*(PP-P(JMIN))/PMA
	PA=PA+PFAC*MAX(TVRDIF(JMIN),0.0)
	NA=NA-PFAC*MIN(TVRDIF(JMIN),0.0)
C
C Find residual positive area above INB and TO
C
       PAT=0.0
       TOB=T(INB)
       IF(INB.LT.N)THEN
        PINB=(P(INB+1)*TVRDIF(INB)-P(INB)*TVRDIF(INB+1))/
     1   (TVRDIF(INB)-TVRDIF(INB+1))
        PAT=RD*TVRDIF(INB)*(P(INB)-PINB)/(P(INB)+PINB)
	  TOB=(T(INB)*(PINB-P(INB+1))+T(INB+1)*(P(INB)-PINB))/
     1    (P(INB)-P(INB+1))
       END IF
C
C Find CAPE 
C
	 CAPED=PA+PAT-NA
	 CAPED=MAX(CAPED,0.0)
	END IF
C
	RETURN
	END
