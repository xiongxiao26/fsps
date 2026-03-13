SUBROUTINE ADD_REMNANTS(mass,maxmass,pset)

  !add remnant (WD, NS, BH) masses back into the total mass
  !of the SSP.  These initial-mass-dependent remnant
  !formulae are taken from Renzini & Ciotti 1993.

  USE SPS_VARS_MODULE_NAME
  IMPLICIT NONE

  REAL(SP), INTENT(inout) :: mass
  !maximum mass still alive
  REAL(SP), INTENT(in) :: maxmass
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP) :: minmass, imfnorm

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !normalize the weights
  imfnorm  = FUNCINT(pset%imf_lower_limit,pset%imf_upper_limit,&
       pset,.TRUE.)

  !BH remnants if any
  !40<M_max<imf_up leave behind a 0.5*M BH
  !if imf_upper_limit < 40 or < maxmass, add no mass since no stars could make BH and mlo=mhi=imf_upper_limit
  minmass = MIN(MAXVAL((/mlim_bh,maxmass/)), pset%imf_upper_limit)
  mass = mass + 0.5*FUNCINT(minmass,pset%imf_upper_limit,&
       pset,.TRUE.)/imfnorm

  !Add NS remnants
  !8.5<M_max<40 also eave behind 1.4 Msun NS
  !if imf_upper_limit < 8.5 , add no mass since no stars could make NS, and mlo=mhi=imf_upper_limit
  IF (maxmass.LE.mlim_bh) THEN
     minmass = MIN(MAXVAL((/mlim_ns,maxmass/)), pset%imf_upper_limit)
     mass = mass + 1.4*FUNCINT(minmass,MIN(mlim_bh,pset%imf_upper_limit),&
          pset,.FALSE.)/imfnorm
  ENDIF

  !Add WD remnants
  !M_max<8.5 also leave behind 0.077*M+0.48 WD
  !if imf_upper_limit < 8.5, only add mass if maxmass < imf_upper_limit
  !since otherwise no stars have evolved and mlo=mhi=imf_upper_limit
  IF (maxmass.LE.8.5) THEN
     minmass = MIN(maxmass, pset%imf_upper_limit)
     mass = mass + 0.48*FUNCINT(minmass,MIN(mlim_ns,pset%imf_upper_limit),&
          pset,.FALSE.)/imfnorm
     mass = mass + 0.077*FUNCINT(minmass,MIN(mlim_ns,pset%imf_upper_limit),&
          pset,.TRUE.)/imfnorm

  ENDIF

  RETURN

END SUBROUTINE ADD_REMNANTS
