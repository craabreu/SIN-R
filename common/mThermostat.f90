module mThermostat

use mGlobal
use mRandom

implicit none

type, abstract :: tThermostat
  real(rb) :: damping
  real(rb) :: meanFactor
  contains
    procedure(tThermostat_integrate), deferred :: integrate
    procedure(tThermostat_energy),    deferred :: energy
end type tThermostat

abstract interface
  subroutine tThermostat_integrate( me, timestep, TwoKE )
    import :: rb, tThermostat
    class(tThermostat), intent(inout) :: me
    real(rb),           intent(in)    :: timestep, TwoKE
  end subroutine tThermostat_integrate

  function tThermostat_energy( me ) result( energy )
    import :: rb, tThermostat
    class(tThermostat), intent(in) :: me
    real(rb)                       :: energy
  end function tThermostat_energy
end interface

!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: nve
contains
  procedure :: energy => nve_energy
  procedure :: integrate => nve_integrate
end type nve

!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: nhc
  integer,  private :: M                     !> Number of thermostats in the chain
  integer,  private :: nloops                !> Number of RESPA loops for integration
  real(rb), private :: kT                    !> Target temperature in energy units
  real(rb), private :: LKT                   !> kT times the number of degrees of freedom
  real(rb), allocatable :: InvQ(:)  !> Inverse of thermostat inertial parameters
  real(rb), allocatable :: eta(:)   !> Thermostat "coordinates"
  real(rb), allocatable :: p(:)     !> Thermostat "momenta"
  contains
    procedure :: setup => nhc_setup
    procedure :: energy => nhc_energy
    procedure :: integrate => nhc_integrate
end type


!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: csvr
  integer  :: dof
  real(rb) :: TwoKE
  real(rb) :: tau
  real(rb) :: Hthermo
  type(mt19937) :: random
contains
  procedure :: setup => csvr_setup
  procedure :: energy => csvr_energy
  procedure :: integrate => csvr_integrate
end type csvr

!---------------------------------------------------------------------------------------------------

contains

  !=================================================================================================

  function nve_energy( me ) result( energy )
    class(nve), intent(in) :: me
    real(rb)               :: energy
    energy = zero
  end function nve_energy

  !-------------------------------------------------------------------------------------------------

  subroutine nve_integrate( me, timestep, TwoKE )
    class(nve), intent(inout) :: me
    real(rb),   intent(in)    :: timestep, TwoKE
    me%damping = zero
  end subroutine nve_integrate

  !=================================================================================================

  subroutine nhc_setup( me, nchain, kT, tdamp, dof, nloops )
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

  !-------------------------------------------------------------------------------------------------

  function nhc_energy( me ) result( energy )
    class(nhc), intent(in) :: me
    real(rb)               :: energy
    if (me%M /= 0) then
      energy = me%LkT*me%eta(1) + me%kT*sum(me%eta(2:me%M)) + half*sum(me%p**2*me%InvQ)
    else
      energy = zero
    end if
  end function nhc_energy

  !-------------------------------------------------------------------------------------------------

  subroutine nhc_integrate( me, timestep, TwoKE )
    class(nhc), intent(inout) :: me
    real(rb),   intent(in)    :: timestep, TwoKE

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
        class(nhc), intent(inout) :: me
        integer,    intent(in) :: j
        real(rb),   intent(in) :: alpha, G, dt
        me%p(j) = me%p(j) + (G - alpha*me%p(j))*phi(alpha*dt)*dt
        me%eta(j+1) = me%eta(j+1) + alpha*dt
      end subroutine integrate
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure real(rb) function phi( x )
        real(rb), intent(in) :: x
        if (abs(x) > 1E-4_rb ) then
          phi = (one - exp(-x))/x
        else
          phi = one + half*x*(third*x*(one - fourth*x) - one)
        end if
      end function phi
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine nhc_integrate

  !=================================================================================================

  subroutine csvr_setup( me, kT, tau, dof, seed )
    class(csvr), intent(inout) :: me
    real(rb),    intent(in)    :: kT, tau
    integer,     intent(in)    :: dof, seed
    me%TwoKE = dof*kT
    me%tau = tau
    me%dof = dof
    me%Hthermo = zero
    call me % random % setup( seed )
  end subroutine csvr_setup

  !-------------------------------------------------------------------------------------------------

  function csvr_energy( me ) result( energy )
    class(csvr), intent(in) :: me
    real(rb)                :: energy
    energy = me%Hthermo
  end function csvr_energy

  !-------------------------------------------------------------------------------------------------

  subroutine csvr_integrate( me, timestep, TwoKE )
    class(csvr), intent(inout) :: me
    real(rb),    intent(in)    :: timestep, TwoKE

    real(rb) :: alphaSq, R1, x, sumRs, A, B, C

    R1 = me%random%normal()
    if (mod(me%dof, 2) == 1) then
      x = (me%dof - 1)/2
      sumRs = 2.0*me%random%gamma(x)
    else
      x = (me%dof - 2)/2
      sumRs = 2.0*me%random%gamma(x) + me%random%normal()**2
    end if
    A = exp(-timestep/me%tau)
    B = 1.0 - A
    C = me%TwoKE/(me%dof*TwoKE)
    alphaSq = A + C*B*(R1**2 + sumRs) + 2.0*sqrt(C*B*A)*R1
    me%damping = -half*log(alphaSq)/timestep
    me%Hthermo = me%Hthermo + (one - alphaSq)*half*TwoKE

  end subroutine csvr_integrate

  !=================================================================================================

end module mThermostat
