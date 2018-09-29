module mCorrelate

use mGlobal

implicit none

type, abstract :: tCorrelate

  real(rb), allocatable :: comtab(:,:,:,:), disptab(:,:,:)   ! unwrap positions at each block
  real(rb), allocatable :: msd(:,:)
  real(rb) ::  msd0(3)
  integer,  allocatable :: counter(:,:), samples(:)
  integer,  private     :: nlevels, step, bsize, N

  contains

    procedure :: tCorrelate_setup
    generic :: setup  => tCorrelate_setup
    procedure :: tCorrelate_sample
    generic :: sample => tCorrelate_sample
    procedure :: tCorrelate_save_to_unit, tCorrelate_save_to_file
    generic :: save => tCorrelate_save_to_unit, tCorrelate_save_to_file
    procedure(tCorrelate_operation), nopass, deferred :: operation 

end type tCorrelate

abstract interface

 elemental function tCorrelate_operation(x,y) result(z)
 import :: rb
 real(rb), intent(in) :: x,y
 real(rb) :: z
 end function tCorrelate_operation

end interface

type, extends(tCorrelate) :: tMSD
   
  contains
    
    procedure, nopass :: operation => tMSD_operation

end type

type, extends(tCorrelate) :: tACF

  contains
    
    procedure, nopass :: operation => tACF_operation

end type

contains

!===================================================================================================

  subroutine tCorrelate_setup( me, nevery, blocksize, Nsteps, R)

    class(tCorrelate), intent(inout) :: me
    integer,         intent(in)    :: nevery, Nsteps, blocksize
    real(rb),        intent(in)    :: R(:,:)

    
    integer :: level, jump, m
    
    me%bsize = blocksize
    me%nlevels = 1
    level = Nsteps/me%bsize
    do while (level /= 0)
      me%nlevels = me%nlevels + 1
      level = level/me%bsize
    end do

    me%N = size(R,2)

    allocate( me%comtab(3,me%N,me%bsize,me%nlevels), source = 0.0_rb )
    allocate( me%disptab(3,me%bsize-1,me%nlevels), source = 0.0_rb ) 
    allocate( me%counter(me%bsize-1,me%nlevels), source = 0 ) 
    allocate( me%msd(5,me%nlevels*(me%bsize-1)+1), source = 0.0_rb )
    me%msd0 = 0.0_rb 
    allocate( me%samples(me%nlevels), source = 0 )

    m = 1
    do level = 1,  me%nlevels
      do jump = 1, me%bsize-1
        m = m + 1
        me%msd(1,m) = nevery*jump*me%bsize**(level-1)
      end do
      me%comtab(:,:,1,level) = R
    end do

    me%nlevels = 0
    me%step = 0

  end subroutine tCorrelate_setup

!===================================================================================================

  subroutine tCorrelate_sample( me, R )
    class(tCorrelate), intent(inout) :: me
    real(rb),        intent(out)   :: R(3,me%N)

    integer  :: level, current, origin, jump
    real(rb) :: delta(3,me%N)

    me%step = me%step + 1

    me%nlevels = 1
    level = me%step/me%bsize
    do while (level /= 0)     
      me%nlevels = me%nlevels + 1
      level = level/me%bsize
    end do

    me%msd0 = me%msd0 + sum(me%operation(R,R),2)
 
    do level = 1, me%nlevels
      if (mod(me%step,me%bsize**(level-1)) == 0) then
        current = me%samples(level) + 1
        me%comtab(:,:,mod(current,me%bsize)+1,level) = R
        do origin = max(0, current - me%bsize + 1), current-1
          jump = current - origin
          me%counter(jump,level) = me%counter(jump,level) + 1
          delta = me%operation(R,me%comtab(:,:,mod(origin,me%bsize)+1,level))
          me%disptab(:,jump,level) = me%disptab(:,jump,level) + sum(delta,2)
        end do
        me%samples(level) = current
      end if
    end do

  end subroutine tCorrelate_sample

!===================================================================================================

  subroutine tCorrelate_save_to_file( me, file, append )
    class(tCorrelate), intent(inout)        :: me
    character(*),         intent(in)           :: file
    logical,        intent(in), optional :: append
    integer, parameter :: out = 65
    logical :: app
    app = present(append)
    if (app) app = append
    if (app) then
      open(unit=out,file=file,status="old",position="append")
    else
      open(unit=out,file=file,status="replace")
    end if
    call tCorrelate_save_to_unit( me, out )
    close(out)

  end subroutine tCorrelate_save_to_file

!===================================================================================================

 subroutine tCorrelate_save_to_unit(me, unit)

    class(tCorrelate), intent(inout) :: me
    integer,        intent(in)    :: unit
    integer :: m, level, levelsize, jump
    me%msd(2:4,1) = me%msd0/(me%N*me%counter(1,1))
    me%msd(5,1) =  me%msd(2,1) + me%msd(3,1) + me%msd(4,1)
    write(unit,'("lag x y z total")')
    write(unit,*) me%msd(:,1)
    m = 1 
    do level = 1, me%nlevels
      levelsize = min(me%samples(level) + 1, me%bsize)
      do jump = 1, levelsize-1
        m = m + 1
        me%msd(2:4,m) = me%disptab(:,jump,level)/(me%N*me%counter(jump,level))
        me%msd(5,m) = me%msd(2,m) + me%msd(3,m) + me%msd(4,m)
        write(unit,*) me%msd(:,m)
      end do
    end do
    
 end subroutine tCorrelate_save_to_unit
!===================================================================================================

 elemental function tMSD_operation(x,y) result(z)
 real(rb), intent(in) :: x,y
 real(rb) :: z
 z = (x-y)**2 
 end function tMSD_operation
!===================================================================================================

 elemental function tACF_operation(x,y) result(z)
 real(rb), intent(in) :: x,y
 real(rb) :: z
 z = x*y 
 end function tACF_operation
!===================================================================================================
end module mCorrelate
