module datm_datamode_atmw_mod

  use ESMF             , only : ESMF_State, ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO
  use NUOPC            , only : NUOPC_Advertise
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_precip_mod   , only : shr_precip_partition_rain_snow_ramp
  use shr_const_mod    , only : shr_const_tkfrz, shr_const_rhofw, shr_const_rdair
  use dshr_methods_mod , only : dshr_state_getfldptr, chkerr
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_get_stream_pointer
  use dshr_strdata_mod , only : shr_strdata_type
  use dshr_fldlist_mod , only : fldlist_type, dshr_fldlist_add

  implicit none
  private ! except

  public  :: datm_datamode_atmw_advertise
  public  :: datm_datamode_atmw_init_pointers
  public  :: datm_datamode_atmw_advance

  ! export state data
   real(r8), pointer :: Sa_tbot(:)           => null()
   real(r8), pointer :: Sa_u10m(:)           => null()
   real(r8), pointer :: Sa_v10m(:)           => null()

  ! stream data
  real(r8), pointer :: strm_mask(:)         => null()

  real(r8) :: tbotmax ! units detector
  real(r8) :: maskmax ! units detector

  real(r8) , parameter :: tKFrz    = SHR_CONST_TKFRZ
  real(r8) , parameter :: rdair    = SHR_CONST_RDAIR ! dry air gas constant ~ J/K/kg
  real(r8) , parameter :: rhofw    = SHR_CONST_RHOFW ! density of fresh water ~ kg/m^3

  character(*), parameter :: nullstr = 'undefined'
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine datm_datamode_atmw_advertise(exportState, fldsexport, &
       flds_scalar_name, rc)

    ! input/output variables
    type(esmf_State)   , intent(inout) :: exportState
    type(fldlist_type) , pointer       :: fldsexport
    character(len=*)   , intent(in)    :: flds_scalar_name
    integer            , intent(out)   :: rc

    ! local variables
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport, 'Sa_tbot'    )
    call dshr_fldList_add(fldsExport, 'Sa_u10m'    )
    call dshr_fldList_add(fldsExport, 'Sa_v10m'    )

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(datm_comp_advertise): Fr_atm '//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

  end subroutine datm_datamode_atmw_advertise

  !===============================================================================
  subroutine datm_datamode_atmw_init_pointers(exportState, sdat, rc)

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: exportState
    type(shr_strdata_type) , intent(in)    :: sdat
    integer                , intent(out)   :: rc

    ! local variables
    character(len=*), parameter :: subname='(datm_init_pointers): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! initialize pointers for module level stream arrays
    call shr_strdata_get_stream_pointer( sdat, 'Sa_mask'   , strm_mask , rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get export state pointers
     call dshr_state_getfldptr(exportState, 'Sa_tbot'    , fldptr1=Sa_tbot    , rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     call dshr_state_getfldptr(exportState, 'Sa_u10m'    , fldptr1=Sa_u10m    , rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     call dshr_state_getfldptr(exportState, 'Sa_v10m'    , fldptr1=Sa_v10m    , rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine datm_datamode_atmw_init_pointers

  !===============================================================================
  subroutine datm_datamode_atmw_advance(exportstate, mainproc, logunit, mpicom, target_ymd, target_tod, model_calendar, rc)
    use ESMF, only: ESMF_VMGetCurrent, ESMF_VMAllReduce, ESMF_REDUCE_MAX, ESMF_VM
    ! input/output variables
    type(ESMF_State)       , intent(inout) :: exportState
    logical                , intent(in)    :: mainproc
    integer                , intent(in)    :: logunit
    integer                , intent(in)    :: mpicom
    integer                , intent(in)    :: target_ymd
    integer                , intent(in)    :: target_tod
    character(len=*)       , intent(in)    :: model_calendar
    integer                , intent(out)   :: rc

    ! local variables
    logical  :: first_time = .true.
    integer  :: n                   ! indices
    integer  :: lsize               ! size of attr vect
    real(r8) :: rtmp(2)
    type(ESMF_VM) :: vm
    character(len=*), parameter :: subname='(datm_datamode_atmw_advance): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lsize = size(strm_mask)

    if (first_time) then
       call ESMF_VMGetCurrent(vm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! determine tbotmax (see below for use)
       rtmp(1) = maxval(Sa_tbot(:))
       call ESMF_VMAllReduce(vm, rtmp, rtmp(2:), 1, ESMF_REDUCE_MAX, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       tbotmax = rtmp(2)

       if (mainproc) write(logunit,*) trim(subname),' tbotmax = ',tbotmax

       ! determine maskmax (see below for use)
       rtmp(1) = maxval(strm_mask(:))
       call ESMF_VMAllReduce(vm, rtmp, rtmp(2:), 1, ESMF_REDUCE_MAX, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       maskmax = rtmp(2)
       if (mainproc) write(logunit,*) trim(subname),' maskmax = ',maskmax

       ! reset first_time
       first_time = .false.
    end if

    do n = 1, lsize
       !--- temperature ---
       if (tbotmax < 50.0_r8) Sa_tbot(n) = Sa_tbot(n) + tkFrz
       ! Limit very cold forcing to 180K
       Sa_tbot(n) = max(180._r8, Sa_tbot(n))
    end do

  end subroutine datm_datamode_atmw_advance

end module datm_datamode_atmw_mod
