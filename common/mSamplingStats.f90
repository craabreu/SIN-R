module mSamplingStats

use mGlobal

implicit none

integer, parameter :: extra = 1000

type SamplingStats

  logical :: initialized = .false.

  integer :: NS           ! Number of sampled states
  integer :: NRT          ! Number of round trips
  logical :: Downhill     ! Flag indicating a downhill walk
  logical :: Collecting   ! Flag indicating that stats are being collected

  integer,  allocatable :: HRT(:)    ! Histogram collected during a round trip
  integer,  allocatable :: HRTdn(:)  ! Histogram collected during a downhill walk
  real(rb), allocatable :: H(:)      ! Accumulated histogram during round trips
  real(rb), allocatable :: Hdn(:)    ! Accumulated histogram during downhill walks

  integer,  allocatable :: Time(:,:) ! Storage of roundtrip times and downhill walk times

  contains

    procedure :: initialize => SamplingStats_initialize
    procedure :: sample => SamplingStats_sample
    procedure :: save => SamplingStats_save

end type SamplingStats

contains

!===================================================================================================

  subroutine SamplingStats_initialize( me, states )
    class(SamplingStats), intent(inout) :: me
    integer,              intent(in)    :: states

    me%NS = states
    allocate( me%HRT(states), me%HRTdn(states), me%H(states), me%Hdn(states) )
    me%Collecting = .false.
    me%NRT = 0
    allocate( me%Time(0,0) )
    me%initialized = .true.

  end subroutine SamplingStats_initialize

!===================================================================================================

  subroutine SamplingStats_sample( me, state )
    class(SamplingStats), intent(inout) :: me
    integer,              intent(in)    :: state

    if (.not.me%Collecting) then
      if (state == me%NS) then
        me%Collecting = .true.
        me%Downhill = .true.
        me%HRT = 0
        me%HRTdn = 0
        me%H = 0.0_rb
        me%Hdn = 0.0_rb
        me%NRT = 0
      end if
    end if

    if (me%Collecting) then
      if (me%Downhill) then
        me%Downhill = (state /= 1)
      else if (state == me%NS) then
        me%Downhill = .true.
        me%NRT = me%NRT + 1
        me%H = me%H + real(me%HRT,rb)
        me%Hdn = me%Hdn + real(me%HRTdn,rb)
        call Update_Buffer( sum(me%HRT), sum(me%HRTdn) )
        me%HRT = 0
        me%HRTdn = 0
      end if
      me%HRT(state) = me%HRT(state) + 1
      if (me%Downhill) me%HRTdn(state) = me%HRTdn(state) + 1
    end if

    contains

      subroutine Update_Buffer( RoundTrip, DownhillWalk )
        integer, intent(in) :: RoundTrip, DownhillWalk
        integer, allocatable :: Aux(:,:)
        if (me%NRT > size(me%Time,2)) then
          Aux = me%Time
          deallocate( me%Time )
          allocate( me%Time(2,me%NRT+extra) )
          me%Time(:,1:size(Aux,2)) = Aux
        end if
        me%Time(:,me%NRT) = [ RoundTrip, DownhillWalk ]
      end subroutine Update_Buffer

  end subroutine SamplingStats_sample

!===================================================================================================

  subroutine SamplingStats_save( me, file, lambda, eta )
    class(SamplingStats), intent(inout)        :: me
    character(*),         intent(in)           :: file
    real(rb),             intent(in), optional :: lambda(me%NS), eta(me%NS)

    integer  :: unit, i
    character(sl) :: process(2) = ["RoundTrip   ","DownhillWalk"]
    real(rb) :: f, avg, stdev

    open( newunit = unit, file = file, status = "replace" )
    if (present(lambda).and.present(eta)) then
      write(unit,'("state,lambda,eta,H,Hdown,f")')
      do i = 1, me%NS
        f = me%Hdn(i)/me%H(i)
        write(unit,'(A)') trim(join([int2str(i), &
                                     real2str([lambda(i),eta(i),me%H(i),me%Hdn(i),f])],","))
      end do
    else
      write(unit,'("state,H,Hdown,f")')
      do i = 1, me%NS
        f = me%Hdn(i)/me%H(i)
        write(unit,'(A)') trim(join([int2str(i),real2str([me%H(i),me%Hdn(i),f])],","))
      end do
    end if

    write(unit,'(/,"Process,Average,StDev")')
    do i = 1, 2
      avg = sum(me%Time(i,1:me%NRT))/me%NRT
      stdev = sqrt(sum((me%Time(i,1:me%NRT) - avg)**2)/(me%NRT - 1))
      write(unit,'(A)') trim(join([process(i), real2str([avg,stdev])],","))
    end do

    write(unit,'(/,"RoundTrip,TotalTime,DownhillTime")')
    do i = 1, me%NRT
      write(unit,'(A)') trim(join(int2str([i,me%Time(:,i)]),","))
    end do
    close( unit )

  end subroutine SamplingStats_save

!===================================================================================================

end module mSamplingStats
