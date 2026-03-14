!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routine to calculate the evolution of a single stellar          !
!  population from a set of input theoretical isochrones and a     !
!  heterogeneous library of stellar spectra.  The code also allows !
!  for variation in the horizontal branch morphology, TP-AGB       !
!  phase, and the blue straggler population.  The output is a      !
!  time-dependent spectrum from the far-UV to the far-IR.          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   PARAMETER RANGES:
!   1.  if (fbhb,sbs)<1E-3 then (fbhb,sbs)=0.0
!   2.  if abs(delt)>0.5 then abs(delt)=0.5
!
!-----------------------------------------------------------!
!-----------------------------------------------------------!

SUBROUTINE SSP_GEN(pset,mass_ssp,lbol_ssp,spec_ssp)

  USE SPS_VARS_MODULE_NAME
  IMPLICIT NONE

  INTEGER :: i, j, ii,klo,khi !,tlo,thi
  !weight given to the entire horizontal branch
  REAL(SP) :: hb_wght,dt,tco
  !array of IMF weights
  REAL(SP), DIMENSION(nm) :: wght
  !SSP spectrum
  REAL(SP), INTENT(inout), DIMENSION(nspec,ntfull) :: spec_ssp
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: tspec_ssp
  !Mass and Lbol info
  REAL(SP), INTENT(inout), DIMENSION(ntfull) :: mass_ssp, lbol_ssp

  !temp arrays for the isochrone data
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: mini,mact,logl,logt,logg,&
       ffco,phase,lmdot
  !arrays holding the number of mass elements for each
  !isochrone and the age of each isochrone
  INTEGER, DIMENSION(nt)     :: nmass
  REAL(SP), DIMENSION(nt)    :: time
  REAL(SP), DIMENSION(nspec) :: tspec
  !structure containing all necessary parameters
  !(TYPE objects defined in sps_vars.f90)
  TYPE(PARAMS), INTENT(in) :: pset
  TYPE(IMF_RUNTIME) :: imf_state
  !CHARACTER(2) :: istr,istr2

  !-----------------------------------------------------------!
  !--------------------------Setup----------------------------!
  !-----------------------------------------------------------!

  IF (check_sps_setup.EQ.0) THEN
     WRITE(*,*) 'SSP_GEN ERROR0: '//&
          'SPS_SETUP must be run once before calling SSP_GEN. '
     STOP
  ENDIF

  !reset arrays
  spec_ssp = 0.
  mass_ssp = 0.
  lbol_ssp = 0.

  !test metallicity range
  IF (pset%zmet.LT.1.OR.pset%zmet.GT.nz) THEN
     WRITE(*,*) 'SSP_GEN ERROR: metallicity outside of range',pset%zmet
     STOP
  ENDIF

  IF (isoc_type.EQ.'bpss') THEN

     !the BPASS SSPs are stored in a master array so simply
     !pull the relevant metallicity model into spec_ssp
     spec_ssp = bpass_spec_ssp(:,:,pset%zmet)
     mass_ssp = bpass_mass_ssp(:,pset%zmet)

  ELSE

     ! PREPARE_IMF fills the per-call custom IMF table and bounds.
     CALL PREPARE_IMF(pset,imf_state)

     ALLOCATE(mini(nt,nm),mact(nt,nm),logl(nt,nm),logt(nt,nm),logg(nt,nm),&
          ffco(nt,nm),phase(nt,nm),lmdot(nt,nm))

     !transfer isochrones into temporary arrays
     mini  = mini_isoc(pset%zmet,:,:)  !initial mass
     mact  = mact_isoc(pset%zmet,:,:)  !actual (present) mass
     logl  = logl_isoc(pset%zmet,:,:)  !log(Lbol)
     logt  = logt_isoc(pset%zmet,:,:)  !log(Teff)
     logg  = logg_isoc(pset%zmet,:,:)  !log(g)
     ffco  = ffco_isoc(pset%zmet,:,:)  !is the TP-AGB star C-rich or O-rich?
     phase = phase_isoc(pset%zmet,:,:) !flag indicating phase of evolution
     lmdot = lmdot_isoc(pset%zmet,:,:) !log Mdot
     nmass = nmass_isoc(pset%zmet,:)   !number of elements per isochrone
     time  = timestep_isoc(pset%zmet,:)!age of each isochrone in log(yr)

     !write for control
     IF (verbose.EQ.1) THEN
        WRITE(*,*)
        WRITE(*,'("   Log(Z/Zsol): ",F6.3)') LOG10(zlegend(pset%zmet)/0.019)
        WRITE(*,'("   Fraction of blue HB stars: ",F6.3)') pset%fbhb
        WRITE(*,'("   Ratio of BS to HB stars  : ",F6.3)') pset%sbss
        WRITE(*,'("   Shift to TP-AGB [log(Teff),log(Lbol)]: ",F5.2,1x,F5.2)') &
             pset%delt, pset%dell
        IF (pset%imf_type.EQ.2) THEN
           WRITE(*,'("   IMF: ",I1,", slopes= ",3F4.1)') &
                pset%imf_type,pset%imf1,pset%imf2,pset%imf3
        ELSE IF (pset%imf_type.EQ.3) THEN
           WRITE(*,'("   IMF: ",I1,", cut-off= ",F4.2)') &
                pset%imf_type,pset%vdmc
        ELSE
           WRITE(*,'("   IMF: ",I1)') pset%imf_type
        ENDIF
     ENDIF

     !-----------------------------------------------------------!
     !---------------------Generate SSPs-------------------------!
     !-----------------------------------------------------------!

     !loop over each isochrone
     DO i=1,nt

        !flag that allows us to compute only a subset of models
        IF (pset%ssp_gen_age(i).EQ.0) CYCLE

        !reset arrays
        hb_wght  = 0.
        wght     = 0.

        IF (verbose.EQ.1) &
             WRITE(*,'("age=",F5.2)') time(i)

        !compute IMF-based weights
        CALL IMF_WEIGHT(mini(i,:),wght,nmass(i),pset,imf_state)

        !modify the horizontal branch
        !need the hb weight for the blue stragglers too
        IF (pset%fbhb.GT.0.0.OR.pset%sbss.GT.1E-3) &
             CALL MOD_HB(pset%fbhb,i,mini,mact,logl,logt,logg,phase,&
             wght,hb_wght,nmass,time(i))

        !add in blue stragglers
        IF (time(i).GE.bhb_sbs_time.AND.pset%sbss.GT.1E-3) &
             CALL ADD_BS(pset%sbss,i,mini,mact,logl,logt,logg,phase,&
             wght,hb_wght,nmass)

        !modify the TP-AGB stars and Post-AGB stars
        CALL MOD_GB(pset%zmet,i,time,pset%delt,pset%dell,pset%pagb,&
             pset%redgb,pset%agb,nmass(i),logl,logt,phase,wght)

        ii = 1 + (i-1)*time_res_incr

        !compute IMF-weighted mass of the SSP
        mass_ssp(ii) = SUM(wght(1:nmass(i))*mact(i,1:nmass(i)))

        !add in remant masses
        IF (add_stellar_remnants.EQ.1) THEN
           CALL ADD_REMNANTS(mass_ssp(ii),MAXVAL(mini(i,:)),pset,imf_state)
        ENDIF

        !compute IMF-weighted bolometric luminosity (actually log(Lbol))
        lbol_ssp(ii) = LOG10(SUM(wght(1:nmass(i))*10**logl(i,1:nmass(i))))

        !compute SSP spectrum
        spec_ssp(:,ii) = 0.
        DO j=1,nmass(i)

           tco = ffco(i,j)
           IF (phase(i,j).EQ.5.AND.tco.GT.1.0) THEN
              !dilute the C star fraction
              !IF (1.0.GE.pset%fcstar) tco = 1.0
           ENDIF

           CALL GETSPEC(pset,mact(i,j),logt(i,j),&
                10**logl(i,j),logg(i,j),phase(i,j),tco,lmdot(i,j),&
                wght(j)/MAXVAL(wght(1:nmass(i))*10**logl(i,1:nmass(i))),tspec)

           !only construct SSPs for particular evolutionary
           !phases if evtype NE -1
           IF ((pset%evtype.EQ.-1.OR.pset%evtype.EQ.phase(i,j))&
                .AND.mini(i,j).LT.pset%masscut) &
                spec_ssp(:,ii) = wght(j)*tspec + spec_ssp(:,ii)

        ENDDO

     ENDDO

  ENDIF

  !-------------------------------------------------------------!
  !-now interpolate the SSPs to fill out the expanded time grid-!
  !-------------------------------------------------------------!

  IF (time_res_incr.GT.1) THEN
     DO j=1,ntfull
        IF (MOD(j-1,time_res_incr).EQ.0) CYCLE
        klo = MAX(MIN(locate(time,time_full(j)),nt-1),1)
        dt  = (time_full(j)-time(klo))/(time(klo+1)-time(klo))
        klo = 1+(klo-1)*time_res_incr
        khi = klo+time_res_incr
        spec_ssp(:,j) = 10**( (1-dt)*LOG10(spec_ssp(:,klo)) + &
             dt*LOG10(spec_ssp(:,khi)))
        lbol_ssp(j)   = (1-dt)*lbol_ssp(klo)   + dt*lbol_ssp(khi)
        mass_ssp(j)   = (1-dt)*mass_ssp(klo)   + dt*mass_ssp(khi)
     ENDDO
  ENDIF

  !-------------------------------------------------------------!
  !-------add the nebular emission model at the SSP level-------!
  !-------------------------------------------------------------!

  IF (add_neb_emission.EQ.2) THEN
     IF (.NOT.ALLOCATED(tspec_ssp)) ALLOCATE(tspec_ssp(nspec,ntfull))
     CALL ADD_NEBULAR(pset,spec_ssp,tspec_ssp)
     spec_ssp = tspec_ssp
  ENDIF

  !-------------------------------------------------------------!
  !---------------add X-ray binaries the SSP level--------------!
  !-------------------------------------------------------------!

  IF (add_xrb_emission.EQ.1) THEN
     IF (.NOT.ALLOCATED(tspec_ssp)) ALLOCATE(tspec_ssp(nspec,ntfull))
     CALL ADD_XRB(pset,spec_ssp,tspec_ssp)
     spec_ssp = tspec_ssp
  ENDIF

  !-------------------------------------------------------------!
  !--------now smooth by an instrumental LSF if provided--------!
  !-------------------------------------------------------------!

  IF (smooth_lsf.EQ.1) THEN
     DO j=1,ntfull
        CALL SMOOTHSPEC(spec_lambda,spec_ssp(:,j),99.d0,lsfinfo%minlam,&
             lsfinfo%maxlam,lsfinfo%lsf)
     ENDDO
  ENDIF


END SUBROUTINE SSP_GEN
