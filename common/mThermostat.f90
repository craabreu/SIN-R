module mThermostat

use mGlobal

implicit none

type, abstract :: nhc
  integer,  private :: M                     !> Number of thermostats in the chain
  integer,  private :: nloops                !> Number of RESPA loops for integration
  real(rb), private :: kT                    !> Target temperature in energy units
  real(rb), private :: LKT                   !> kT times the number of degrees of freedom
  real(rb) :: damping
  real(rb) :: meanFactor
  real(rb), allocatable :: InvQ(:)  !> Inverse of thermostat inertial parameters
  real(rb), allocatable :: eta(:)   !> Thermostat "coordinates"
  real(rb), allocatable :: p(:)     !> Thermostat "momenta"
  contains
    procedure :: setup => nhc_setup
    procedure :: energy => nhc_energy
    procedure(nhc_integrate), deferred :: integrate
end type

abstract interface
  elemental subroutine nhc_integrate( me, timestep, TwoKE )
    import :: rb, nhc
    class(nhc), intent(inout) :: me
    real(rb),   intent(in)    :: timestep, TwoKE
  end subroutine nhc_integrate
end interface

type, extends(nhc) :: nhc_pscaling
  contains
    procedure :: integrate => nhc_pscaling_integrate
end type nhc_pscaling

type, extends(nhc) :: nhc_boosting
  contains
    procedure :: integrate => nhc_boosting_integrate
end type nhc_boosting

type, extends(nhc) :: nhc_kamberaj
  contains
    procedure :: integrate => nhc_kamberaj_integrate
end type nhc_kamberaj

contains

  !=================================================================================================

  elemental subroutine nhc_setup( me, nchain, kT, tdamp, dof, nloops )
    class(nhc), intent(inout) :: me
    integer,    intent(in)    :: nchain
    real(rb),   intent(in)    :: kT, tdamp
    integer,    intent(in)    :: dof
    integer,    intent(in)    :: nloops
    me%M = nchain
    me%kT = kT
    me%LKT = dof*kT
    me%nloops = nloops
    allocate( me%InvQ(nchain), me%eta(nchain), me%p(nchain) )
    me%InvQ(1) = one/(me%LKT*tdamp**2)
    me%InvQ(2:nchain) = one/(kT*tdamp**2)
    me%eta = zero
    me%p = zero
  end subroutine nhc_setup

  !=================================================================================================

  elemental function nhc_energy( me ) result( energy )
    class(nhc), intent(in) :: me
    real(rb)               :: energy
    if (me%M /= 0) then
      energy = me%LkT*me%eta(1) + me%kT*sum(me%eta(2:me%M)) + half*sum(me%p**2*me%InvQ)
    else
      energy = zero
    end if
  end function nhc_energy

  !=================================================================================================

  elemental real(rb) function phi( x )
    real(rb), intent(in) :: x
    if (abs(x) > 1E-4_rb ) then
      phi = (one - exp(-x))/x
    else
      phi = one + half*x*(third*x*(one - fourth*x) - one)
    end if
  end function phi

  !=================================================================================================

  elemental subroutine nhc_pscaling_integrate( me, timestep, TwoKE )
    class(nhc_pscaling), intent(inout) :: me
    real(rb),            intent(in)    :: timestep, TwoKE

    integer :: i, j
    real(rb) :: dt, dt_2, twodt, alpha, alphaSum, factor, sumFactor

    dt = timestep/me%nloops
    dt_2 = half*dt
    twodt = two*dt
    alphaSum = zero
    factor = one
    sumFactor = zero
    do i = 1, me%nloops
      me%p(me%M) = me%p(me%M) + (me%p(me%M-1)**2*me%InvQ(me%M-1) - me%kT)*dt_2
      do j = me%M-1, 2, -1
        call integrate( me, j, me%p(j+1)*me%InvQ(j+1), me%p(j-1)**2*me%InvQ(j-1) - me%kT, dt_2 )
      end do
      call integrate( me, 1, me%p(2)*me%InvQ(2), factor**2*twoKE - me%LkT, dt_2 )
      alpha = me%p(1)*me%InvQ(1)
      alphaSum = alphaSum + alpha
      factor = exp(-alphaSum*dt)
      sumFactor = sumFactor + factor
      call integrate( me, 1, me%p(2)*me%InvQ(2), factor**2*twoKE - me%LkT, dt_2 )
      do j = 2, me%M-1
        call integrate( me, j, me%p(j+1)*me%InvQ(j+1), me%p(j-1)**2*me%InvQ(j-1) - me%kT, dt_2 )
      end do
      me%p(me%M) = me%p(me%M) + (me%p(me%M-1)**2*me%InvQ(me%M-1) - me%kT)*dt_2
    end do
    me%eta(1) = me%eta(1) + alphaSum*dt
    me%damping = alphaSum/me%nloops
    me%meanFactor = sumFactor/me%nloops

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine integrate( me, j, alpha, G, dt )
        class(nhc_pscaling), intent(inout) :: me
        integer,             intent(in) :: j
        real(rb),            intent(in) :: alpha, G, dt
        me%p(j) = me%p(j) + (G - alpha*me%p(j))*phi(alpha*dt)*dt
        me%eta(j+1) = me%eta(j+1) + alpha*dt
      end subroutine integrate
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine nhc_pscaling_integrate

  !=================================================================================================

  elemental subroutine nhc_boosting_integrate( me, timestep, TwoKE )
    class(nhc_boosting), intent(inout) :: me
    real(rb),            intent(in)    :: timestep, TwoKE

    integer :: i, j
    real(rb) :: dt, dt_2

    dt = timestep/me%nloops
    dt_2 = half*dt
    do i = 1, me%nloops
      me%p(me%M) = me%p(me%M) + (me%p(me%M-1)**2*me%InvQ(me%M-1) - me%kT)*dt_2
      do j = me%M-1, 2, -1
        call integrate( me, j, me%p(j+1)*me%InvQ(j+1), me%p(j-1)**2*me%InvQ(j-1) - me%kT, dt_2 )
      end do
      call integrate( me, 1, me%p(2)*me%InvQ(2), twoKE - me%LkT, dt )
      do j = 2, me%M-1
        call integrate( me, j, me%p(j+1)*me%InvQ(j+1), me%p(j-1)**2*me%InvQ(j-1) - me%kT, dt_2 )
      end do
      me%p(me%M) = me%p(me%M) + (me%p(me%M-1)**2*me%InvQ(me%M-1) - me%kT)*dt_2
    end do

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine integrate( me, j, alpha, G, dt )
        class(nhc_boosting), intent(inout) :: me
        integer,             intent(in) :: j
        real(rb),            intent(in) :: alpha, G, dt
        me%p(j) = me%p(j) + (G - alpha*me%p(j))*phi(alpha*dt)*dt
        me%eta(j+1) = me%eta(j+1) + alpha*dt
      end subroutine integrate
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine nhc_boosting_integrate

  !=================================================================================================

  elemental subroutine nhc_kamberaj_integrate( me, timestep, TwoKE )
    class(nhc_kamberaj), intent(inout) :: me
    real(rb),            intent(in)    :: timestep, TwoKE

    integer :: i, j
    real(rb) :: dt, dt_2

    dt = timestep/me%nloops
    dt_2 = half*dt
    do i = 1, me%nloops
      me%p(me%M) = me%p(me%M) + (me%p(me%M-1)**2*me%InvQ(me%M-1) - me%kT)*dt_2
      do j = me%M-1, 2, -1
        call integrate( me, j, me%p(j+1)*me%InvQ(j+1), me%p(j-1)**2*me%InvQ(j-1) - me%kT, dt_2 )
      end do
      call integrate( me, 1, me%p(2)*me%InvQ(2), twoKE - me%LkT, dt_2 )
      me%eta = me%eta + me%p*me%InvQ*dt
      call integrate( me, 1, me%p(2)*me%InvQ(2), twoKE - me%LkT, dt_2 )
      do j = 2, me%M-1
        call integrate( me, j, me%p(j+1)*me%InvQ(j+1), me%p(j-1)**2*me%InvQ(j-1) - me%kT, dt_2 )
      end do
      me%p(me%M) = me%p(me%M) + (me%p(me%M-1)**2*me%InvQ(me%M-1) - me%kT)*dt_2
    end do

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine integrate( me, j, alpha, G, dt )
        class(nhc_kamberaj), intent(inout) :: me
        integer,             intent(in) :: j
        real(rb),            intent(in) :: alpha, G, dt
        me%p(j) = me%p(j) + (G - alpha*me%p(j))*phi(alpha*dt)*dt
      end subroutine integrate
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine nhc_kamberaj_integrate

  !=================================================================================================

end module mThermostat
