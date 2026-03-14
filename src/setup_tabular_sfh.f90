! Tabulated SFH defined in a file called sfh.dat that must reside in the data
! directory (or a file specified via the parameter sfh filename; see
! below). The file must contain three columns. The first column is time since
! the Big Bang in Gyr, the second is the SFR in units of solar masses per year,
! the third is the absolute metallicity. An example is provided in the data
! directory. The time grid in this file can be arbitrary (so long as the units
! are correct), but it is up to the user to ensure that the tabulated sfh is
! well-sampled so that the outputs are stable. Obviously, highly oscillatory
! data require dense sampling.
!
! SFRs are clipped to a minimum value of 1e-30 to avoid divide by zero errors.


subroutine setup_tabular_sfh(pset, nzin)

  use, intrinsic :: iso_c_binding, only: c_associated, c_f_pointer, c_loc
  use SPS_VARS_MODULE_NAME, only: nz, tiny_number, tiny30, PARAMS, SPS_HOME, &
       SP, sfh_tab_file_buffer
  implicit none
  type(PARAMS), intent(inout) :: pset
  integer, intent(in) :: nzin
  integer :: stat, n, nrows
  real(SP), pointer :: sfh_tab(:,:)
  real(SP) :: age_value, sfr_value, z_value

  IF (pset%sfh.EQ.2) THEN

     if (pset%sf_start.gt.tiny_number) then
        WRITE(*,*) 'COMPSP ERROR: Tabular sfh, but sf_start > 0'
        STOP
     endif

     ! Read the sfh file
     IF (TRIM(pset%sfh_filename).EQ.'') THEN
        OPEN(3,FILE=TRIM(SPS_HOME)//'/data/sfh.dat',ACTION='READ',STATUS='OLD')
     ELSE
        OPEN(3,FILE=TRIM(SPS_HOME)//'/data/'//TRIM(pset%sfh_filename),&
             ACTION='READ',STATUS='OLD')
     ENDIF
     nrows = 0
     DO
        IF (nzin.EQ.nz) THEN
           READ(3,*,IOSTAT=stat) age_value, sfr_value, z_value
        ELSE
           READ(3,*,IOSTAT=stat) age_value, sfr_value
        ENDIF
        IF (stat.NE.0) EXIT
        nrows = nrows + 1
     ENDDO
     IF (stat.GT.0) THEN
        WRITE(*,*) 'COMPSP ERROR: failed while counting the tabular SFH file'
        STOP
     ENDIF
     IF (nrows.EQ.0) THEN
        WRITE(*,*) 'COMPSP ERROR: tabular SFH file is empty'
        STOP
     ENDIF
     REWIND(3)
     IF (ALLOCATED(sfh_tab_file_buffer)) DEALLOCATE(sfh_tab_file_buffer)
     ALLOCATE(sfh_tab_file_buffer(3,nrows))
     DO n=1,nrows
        IF (nzin.EQ.nz) THEN
           READ(3,*,IOSTAT=stat) sfh_tab_file_buffer(1,n),sfh_tab_file_buffer(2,n),sfh_tab_file_buffer(3,n)
        ELSE
           READ(3,*,IOSTAT=stat) sfh_tab_file_buffer(1,n),sfh_tab_file_buffer(2,n)
           sfh_tab_file_buffer(3,n)=0.0
        ENDIF
        IF (stat.NE.0) THEN
           WRITE(*,*) 'COMPSP ERROR: failed while reading the tabular SFH file'
           STOP
        ENDIF
     ENDDO
     CLOSE(3)

     pset%ntabsfh = nrows
     pset%sfh_tab = c_loc(sfh_tab_file_buffer)
     CALL c_f_pointer(pset%sfh_tab, sfh_tab, [3, pset%ntabsfh])
     sfh_tab(1,1:pset%ntabsfh) = sfh_tab(1,1:pset%ntabsfh)*1E9 !convert to yrs

     !special switch to compute only the last time output
     !in the tabulated file
     ! IF (pset%tage.EQ.-99.) imin=imax

  ELSE IF (pset%sfh.EQ.3) THEN

     !pset%sfh_tab array is supposed to already be filled in, check that it is
     IF ((pset%ntabsfh.EQ.0).OR.(.NOT.c_associated(pset%sfh_tab))) THEN
        WRITE(*,*) 'COMPSP ERROR: sfh=3 but pset%sfh_tab array not initialized!'
        STOP
     ENDIF
     CALL c_f_pointer(pset%sfh_tab, sfh_tab, [3, pset%ntabsfh])

  ENDIF

  ! clip SFR to a minimum of 1e-30
  do n=1, pset%ntabsfh
     sfh_tab(2, n) = max(sfh_tab(2, n), tiny30)
  enddo

end subroutine setup_tabular_sfh
