FUNCTION IMF(mass,pset,mass_weighted)

  !define IMFs (dn/dM)

  !the following is just a fast-and-loose way to pass things around:
  !if the imf_type var is +10 then we calculate dn/dm*m
  !if the imf_type var is <10 then we calculate dn/dm
  
  USE SPS_VARS_MODULE_NAME
  IMPLICIT NONE

  REAL(SP), DIMENSION(:), INTENT(in) :: mass
  TYPE(PARAMS), INTENT(in) :: pset
  LOGICAL, INTENT(in) :: mass_weighted
  REAL(SP), DIMENSION(SIZE(mass)) :: imf
  INTEGER :: i,n,imf_type_local
  REAL(SP) :: imfcu

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  imf = 0.0
  imf_type_local = pset%imf_type
  IF (mass_weighted) imf_type_local = imf_type_local + 10

  !Salpeter (1955) IMF
  IF (MOD(imf_type_local,10).EQ.0) THEN
     imf = mass**(-salp_ind) 
     IF (imf_type_local.EQ.10) imf = mass*imf
  ENDIF
  
  !Chabrier (2003) IMF
  IF (MOD(imf_type_local,10).EQ.1) THEN
     DO i=1,size(mass)
        IF (mass(i).LT.1) THEN
           imf(i) = exp(-(log10(mass(i))-log10(chab_mc))**2&
                /2/chab_sigma2)
        ELSE
           imf(i) = exp(-log10(chab_mc)**2/2./chab_sigma2)*&
                mass(i)**(-chab_ind)
        ENDIF
     ENDDO
     !convert from dn/dlnM to dn/dM
     imf = imf/mass
     IF (imf_type_local.EQ.11) imf = mass*imf
  ENDIF
  
  !Kroupa (2001) IMF
  IF (MOD(imf_type_local,10).EQ.2) THEN
     DO i=1,size(mass)
        IF (mass(i).GE.0.08.AND.mass(i).LT.0.5) &
             imf(i) = mass(i)**(-pset%imf1)
        IF (mass(i).GE.0.5.AND.mass(i).LT.1.0) &
             imf(i) = 0.5**(-pset%imf1+pset%imf2)*&
             mass(i)**(-pset%imf2)
        IF (mass(i).GE.1.0) &
             imf(i) = 0.5**(-pset%imf1+pset%imf2)*&
             mass(i)**(-pset%imf3)
     ENDDO
     IF (imf_type_local.EQ.12) imf = mass*imf
  ENDIF
  
  !van Dokkum (2008) IMF
  IF (MOD(imf_type_local,10).EQ.3) THEN
     DO i=1,size(mass)
        IF (mass(i).LE.vd_nc*pset%vdmc) THEN
           imf(i) = vd_al*(0.5*vd_nc*pset%vdmc)**(-vd_ind)*&
                exp(-(log10(mass(i))-log10(pset%vdmc))*&
                (log10(mass(i))-log10(pset%vdmc))/2./vd_sigma2)
        ELSE
           imf(i) = vd_ah*mass(i)**(-vd_ind)
        ENDIF
     ENDDO
     !convert from dn/dlnM to dn/dM
     imf = imf/mass
     IF (imf_type_local.EQ.13) imf = mass*imf
  ENDIF
  
  !Dave (2008) IMF
  IF (MOD(imf_type_local,10).EQ.4) THEN
     DO i=1,size(mass)
        IF (mass(i).GE.0.08.AND.mass(i).LT.pset%mdave) &
             imf(i) = mass(i)**(-pset%imf1)
        IF (mass(i).GE.pset%mdave) &
             imf(i) = pset%mdave**(-pset%imf1+pset%imf2)*&
             mass(i)**(-pset%imf2)
     ENDDO
     IF (imf_type_local.EQ.14) imf = mass*imf
  ENDIF
 
  !user-defined IMF
  IF (MOD(imf_type_local,10).EQ.5) THEN
     DO i=1,size(mass)
        IF (mass(i).GE.pset%imf_user_alpha(1,1).AND.&
                mass(i).LT.pset%imf_user_alpha(2,1)) &
                imf(i) = mass(i)**(-pset%imf_user_alpha(3,1))
        imfcu = 1.0
        DO n=2,pset%n_user_imf
           IF (mass(i).GE.pset%imf_user_alpha(1,n).AND.&
                mass(i).LT.pset%imf_user_alpha(2,n)) &
                imf(i) = mass(i)**(-pset%imf_user_alpha(3,n))*&
                pset%imf_user_alpha(1,n)**(-pset%imf_user_alpha(3,n-1)+&
                pset%imf_user_alpha(3,n))*imfcu
           imfcu = imfcu*pset%imf_user_alpha(1,n)**(-pset%imf_user_alpha(3,n-1)+&
                pset%imf_user_alpha(3,n))
        ENDDO
     ENDDO
     IF (imf_type_local.EQ.15) imf = mass*imf
  ENDIF

END FUNCTION IMF

!---------------------------------------------------------------!
!---------------------------------------------------------------!

SUBROUTINE PREPARE_IMF(pset)

  USE SPS_VARS_MODULE_NAME
  IMPLICIT NONE

  ! Populate per-call custom IMF data without mutating shared module state.
  TYPE(PARAMS), INTENT(inout) :: pset
  INTEGER :: i, stat, imf_unit

  pset%n_user_imf = 0
  pset%imf_user_alpha = 0.
  pset%imf_lower_limit = 0.08
  pset%imf_upper_limit = 120.0

  IF (pset%imf_type.LT.0.OR.pset%imf_type.GT.5) THEN
     WRITE(*,*) 'SSP_GEN ERROR: IMF type outside of range',pset%imf_type
     STOP
  ENDIF

  IF (pset%imf_type.EQ.5) THEN
     IF (TRIM(pset%imf_filename).EQ.'') THEN
        OPEN(NEWUNIT=imf_unit,FILE=TRIM(SPS_HOME)//'/data/imf.dat',&
             ACTION='READ',STATUS='OLD')
     ELSE
        OPEN(NEWUNIT=imf_unit,FILE=TRIM(SPS_HOME)//'/data/'//TRIM(pset%imf_filename),&
             ACTION='READ',STATUS='OLD')
     ENDIF
     DO i=1,100
        READ(imf_unit,*,IOSTAT=stat) pset%imf_user_alpha(1,i),&
             pset%imf_user_alpha(2,i),pset%imf_user_alpha(3,i)
        IF (stat.NE.0) GOTO 29
     ENDDO
     WRITE(*,*) 'SSP_GEN ERROR: didnt finish reading in the imf file'
     STOP
29   CONTINUE
     CLOSE(imf_unit)
     pset%n_user_imf = i-1
     pset%imf_lower_limit = pset%imf_user_alpha(1,1)
     pset%imf_upper_limit = pset%imf_user_alpha(2,pset%n_user_imf)
  ENDIF

END SUBROUTINE PREPARE_IMF
