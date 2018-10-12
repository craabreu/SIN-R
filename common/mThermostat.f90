module mThermostat

use mGlobal
use mRandom
use iso_c_binding
use omp_lib

implicit none

type, abstract :: tThermostat
  integer  :: dof                         ! Number of degrees of freedom
  integer  :: seed                        ! Seed for random numbers
  integer  :: nthreads                    ! Number of parallel threads
  real(rb) :: kT                          ! Temperature in energy units
  real(rb) :: tau                         ! Characteristic time scale
  real(rb),      allocatable :: m(:)      ! Mass associated to each degree of freedom
  integer,       allocatable :: first(:)  ! The first degree of freedom of each thread
  integer,       allocatable :: last(:)   ! The last degree of freedom of each thread
  type(mt19937), allocatable :: random(:) ! A random number generator for each thread
  contains
    procedure :: setup => tThermostat_setup
    procedure :: finish_setup => tThermostat_finish_setup
    procedure :: boost => tThermostat_boost
    procedure :: energy => tThermostat_energy
    procedure(tThermostat_integrate), deferred :: integrate
end type tThermostat

abstract interface
  subroutine tThermostat_integrate( me, timestep, v )
    import :: rb, tThermostat
    class(tThermostat), intent(inout) :: me
    real(rb),           intent(in)    :: timestep
    real(rb),           intent(inout) :: v(*)
  end subroutine tThermostat_integrate
end interface

!---------------------------------------------------------------------------------------------------
!                                  F A K E    T H E R M O S T A T
!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: nve
contains
  procedure :: integrate => nve_integrate
end type nve

!---------------------------------------------------------------------------------------------------
!                                N O S É - H O O V E R     C H A I N
!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: nhc
  integer,  private :: NC           !> Number of thermostats in the chain
  integer,  private :: nloops       !> Number of RESPA loops for integration
  real(rb), private :: LKT          !> kT times the number of degrees of freedom
  real(rb), allocatable :: InvQ(:)  !> Inverse of thermostat inertial parameters
  real(rb), allocatable :: eta(:)   !> Thermostat "coordinates"
  real(rb), allocatable :: p(:)     !> Thermostat "momenta"
  contains
    procedure :: initialize => nhc_initialize
    procedure :: energy => nhc_energy
    procedure :: integrate => nhc_integrate
end type

!---------------------------------------------------------------------------------------------------
!                    S T O C H A S T I C    V E L O C I T Y    R E S C A L I N G
!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: csvr
  real(rb) :: TwoKE
  real(rb) :: Hthermo
contains
  procedure :: finish_setup => csvr_finish_setup
  procedure :: energy => csvr_energy
  procedure :: integrate => csvr_integrate
end type csvr

!---------------------------------------------------------------------------------------------------
!                                          L A N G E V I N
!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: langevin
  real(rb) :: gamma
  real(rb), allocatable :: invSqrtM(:)
contains
  procedure :: finish_setup => langevin_finish_setup
  procedure :: initialize => langevin_initialize
  procedure :: integrate => langevin_integrate
  procedure :: boost => langevin_boost
end type langevin

!---------------------------------------------------------------------------------------------------
!        S T O C H A S T I C    I S O K I N E T I C    N O S E - H O O V E R    R E S P A
!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: sinr
  real(rb) :: Q1       ! Inertial parameter
  real(rb) :: Q2       ! Inertial parameter
  real(rb) :: halfQ1   ! Half inertial parameter
  real(rb) :: gamma
  real(rb), allocatable :: v1(:)
  real(rb), allocatable :: v2(:)
contains
  procedure :: finish_setup => sinr_finish_setup
  procedure :: initialize => sinr_initialize
  procedure :: integrate => sinr_integrate
  procedure :: boost => sinr_boost
end type sinr

!---------------------------------------------------------------------------------------------------

contains

  !=================================================================================================
  !                                A L L    T H E R M O S T A T S
  !=================================================================================================

  subroutine tThermostat_setup( me, dof, mass, kT, tau, seed, threads )
    class(tThermostat), intent(inout) :: me
    integer,            intent(in)    :: dof
    real(rb),           intent(in)    :: mass(*), kT, tau
    integer,            intent(in)    :: seed, threads

    integer :: i, perThread, newseed

    me%dof = dof
    me%m = mass(1:dof)
    me%kT = kT
    me%tau = tau

    me%nthreads =  threads
    perThread = (dof + threads - 1)/threads
    me%first = [((i - 1)*perThread + 1, i=1, threads)]
    me%last = [(min(i*perThread, dof), i=1, threads)]
    allocate( me%random(threads) )
    newseed = seed
    do i = 1, threads
      call me % random(i) % setup( newseed )
      newseed = me%random(i)%i32()
    end do

    call me % finish_setup()

  end subroutine tThermostat_setup

  !-------------------------------------------------------------------------------------------------

  subroutine tThermostat_finish_setup( me )
    class(tThermostat), intent(inout) :: me
  end subroutine tThermostat_finish_setup

  !-------------------------------------------------------------------------------------------------

  subroutine tThermostat_boost( me, dt, a, v )
    class(tThermostat), intent(in)    :: me
    real(rb),           intent(in)    :: dt, a(*)
    real(rb),           intent(inout) :: v(*)
    v(1:me%dof) = v(1:me%dof) + dt*a(1:me%dof)
  end subroutine tThermostat_boost

  !-------------------------------------------------------------------------------------------------

  function tThermostat_energy( me ) result( energy )
    class(tThermostat), intent(in) :: me
    real(rb)                       :: energy
    energy = zero
  end function tThermostat_energy

  !=================================================================================================
  !                                F A K E    T H E R M O S T A T
  !=================================================================================================

  subroutine nve_integrate( me, timestep, v )
    class(nve), intent(inout) :: me
    real(rb),   intent(in)    :: timestep
    real(rb),   intent(inout) :: v(*)
  end subroutine nve_integrate

  !=================================================================================================
  !                              N O S É - H O O V E R     C H A I N
  !=================================================================================================

  subroutine nhc_initialize( me, nchain, nloops )
    class(nhc), intent(inout) :: me
    integer,    intent(in)    :: nchain
    integer,    intent(in)    :: nloops
    me%NC = nchain
    me%LKT = (me%dof - 3)*me%kT
    me%nloops = nloops
    allocate( me%InvQ(nchain), me%eta(nchain), me%p(nchain) )
    me%InvQ(1) = one/(me%LKT*me%tau**2)
    me%InvQ(2:nchain) = one/(me%kT*me%tau**2)
    me%eta = zero
    me%p = zero
  end subroutine nhc_initialize

  !-------------------------------------------------------------------------------------------------

  function nhc_energy( me ) result( energy )
    class(nhc), intent(in) :: me
    real(rb)               :: energy
    if (me%NC /= 0) then
      energy = me%LkT*me%eta(1) + me%kT*sum(me%eta(2:me%NC)) + half*sum(me%p**2*me%InvQ)
    else
      energy = zero
    end if
  end function nhc_energy

  !-------------------------------------------------------------------------------------------------

  subroutine nhc_integrate( me, timestep, v )
    class(nhc), intent(inout) :: me
    real(rb),   intent(in)    :: timestep
    real(rb),   intent(inout) :: v(*)

    integer  :: i, j
    real(rb) :: TwoKE, dt, dt_2, twodt, alpha, alphaSum, factor

    TwoKE = sum(me%m*v(1:me%dof)**2)
    dt = timestep/me%nloops
    dt_2 = half*dt
    twodt = two*dt
    alphaSum = zero
    factor = one
    do i = 1, me%nloops
      me%p(me%NC) = me%p(me%NC) + (me%p(me%NC-1)**2*me%InvQ(me%NC-1) - me%kT)*dt_2
      do j = me%NC-1, 2, -1
        call integrate( me, j, me%p(j+1)*me%InvQ(j+1), me%p(j-1)**2*me%InvQ(j-1) - me%kT, dt_2 )
      end do
      call integrate( me, 1, me%p(2)*me%InvQ(2), factor**2*twoKE - me%LkT, dt_2 )
      alpha = me%p(1)*me%InvQ(1)
      alphaSum = alphaSum + alpha
      factor = exp(-alphaSum*dt)
      call integrate( me, 1, me%p(2)*me%InvQ(2), factor**2*twoKE - me%LkT, dt_2 )
      do j = 2, me%NC-1
        call integrate( me, j, me%p(j+1)*me%InvQ(j+1), me%p(j-1)**2*me%InvQ(j-1) - me%kT, dt_2 )
      end do
      me%p(me%NC) = me%p(me%NC) + (me%p(me%NC-1)**2*me%InvQ(me%NC-1) - me%kT)*dt_2
    end do
    me%eta(1) = me%eta(1) + alphaSum*dt
    v(1:me%dof) = v(1:me%dof)*exp(-dt*alphaSum)

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
  !                  S T O C H A S T I C    V E L O C I T Y    R E S C A L I N G
  !=================================================================================================

  subroutine csvr_finish_setup( me )
    class(csvr), intent(inout) :: me
    me%TwoKE = me%dof*me%kT
    me%Hthermo = zero
  end subroutine csvr_finish_setup

  !-------------------------------------------------------------------------------------------------

  function csvr_energy( me ) result( energy )
    class(csvr), intent(in) :: me
    real(rb)                :: energy
    energy = me%Hthermo
  end function csvr_energy

  !-------------------------------------------------------------------------------------------------

  subroutine csvr_integrate( me, timestep, v )
    class(csvr), intent(inout) :: me
    real(rb),    intent(in)    :: timestep
    real(rb),    intent(inout) :: v(*)

    real(rb) :: TwoKE, alphaSq, R1, x, sumRs, A, B, C

    TwoKE = sum(me%m*v(1:me%dof)**2)
    R1 = me%random(1)%normal()
    if (mod(me%dof, 2) == 1) then
      x = (me%dof - 1)/2
      sumRs = 2.0*me%random(1)%gamma(x)
    else
      x = (me%dof - 2)/2
      sumRs = 2.0*me%random(1)%gamma(x) + me%random(1)%normal()**2
    end if
    A = exp(-timestep/me%tau)
    B = 1.0 - A
    C = me%TwoKE/(me%dof*TwoKE)
    alphaSq = A + C*B*(R1**2 + sumRs) + 2.0*sqrt(C*B*A)*R1
    me%Hthermo = me%Hthermo + (one - alphaSq)*half*TwoKE
    v(1:me%dof) = v(1:me%dof)*sqrt(alphaSq)

  end subroutine csvr_integrate

  !=================================================================================================
  !        S T O C H A S T I C    I S O K I N E T I C    N O S E - H O O V E R    R E S P A
  !=================================================================================================

  subroutine langevin_finish_setup( me )
    class(langevin), intent(inout) :: me

    me%gamma = one/me%tau
    allocate( me%invSqrtM(me%dof) )
    me%invSqrtM = one/sqrt(me%m)

  end subroutine langevin_finish_setup

  !-------------------------------------------------------------------------------------------------

  subroutine langevin_initialize( me, friction )
    class(langevin), intent(inout) :: me
    real(rb),        intent(in)    :: friction
    me%gamma = friction
  end subroutine langevin_initialize

  !-------------------------------------------------------------------------------------------------

  subroutine langevin_integrate( me, timestep, v )
    class(langevin), intent(inout) :: me
    real(rb),    intent(in)    :: timestep
    real(rb),    intent(inout) :: v(*)

    real(rb) :: A, B

    A = exp(-me%gamma*timestep)
    B = sqrt(me%kT*(one - A*A))

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, i
      thread = omp_get_thread_num() + 1
      do i = me%first(thread), me%last(thread)
        v(i) = A*v(i) + B*me%invSqrtM(i)*me%random(thread)%normal()
      end do
    end block
    !$omp end parallel

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine ornstein_uhlenbeck( i, dt, R )
        integer,  intent(in) :: i
        real(rb), intent(in) :: dt, R
        
      end subroutine ornstein_uhlenbeck
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine langevin_integrate

  !-------------------------------------------------------------------------------------------------

  subroutine langevin_boost( me, dt, a, v )
    class(langevin), intent(in)    :: me
    real(rb),    intent(in)    :: dt, a(*)
    real(rb),    intent(inout) :: v(*)
    v(1:me%dof) = v(1:me%dof) + dt*a(1:me%dof)
  end subroutine langevin_boost

  !=================================================================================================
  !        S T O C H A S T I C    I S O K I N E T I C    N O S E - H O O V E R    R E S P A
  !=================================================================================================

  subroutine sinr_finish_setup( me )
    class(sinr), intent(inout) :: me
    me%Q1 = me%kT*me%tau**2
    me%Q2 = me%kT*me%tau**2
    me%halfQ1 = half*me%Q1
    me%gamma = one/me%tau
    allocate( me%v1(me%dof), me%v2(me%dof), source = zero )
  end subroutine sinr_finish_setup

  !-------------------------------------------------------------------------------------------------

  subroutine sinr_initialize( me, friction, v )
    class(sinr), intent(inout) :: me
    real(rb),    intent(in)    :: friction
    real(rb),    intent(inout) :: v(*)

    me%gamma = friction

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, i
      thread = omp_get_thread_num() + 1
      do i = me%first(thread), me%last(thread)
        call initialize( thread, v(i), me%v1(i), me%v2(i), me%m(i) )
      end do
    end block
    !$omp end parallel

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine initialize( thread, v, v1, v2, m )
        integer,  intent(in)  :: thread
        real(rb), intent(out) :: v, v1, v2
        real(rb), intent(in)  :: m
        real(rb) :: factor
        v = sqrt(half*me%kT/m)*me%random(thread)%normal()
        v1 = sqrt(me%kT/me%Q1)*me%random(thread)%normal()
        factor = sqrt(me%kT/(m*v*v + me%halfQ1*v1*v1))
        v = factor*v
        v1 = factor*v1
        v2 = sqrt(me%kT/me%Q2)*me%random(thread)%normal()
      end subroutine initialize
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine sinr_initialize

  !-------------------------------------------------------------------------------------------------

  subroutine sinr_integrate( me, timestep, v )
    class(sinr), intent(inout) :: me
    real(rb),    intent(in)    :: timestep
    real(rb),    intent(inout) :: v(*)

    real(rb) :: dt_2
    real(rb) :: A, B, C

    dt_2 = half*timestep
    A = exp(-me%gamma*timestep)
    B = sqrt(me%kT*(one - A*A)/me%Q2)
    C = (one - A)/(me%gamma*me%Q2)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, i
      thread = omp_get_thread_num() + 1
      do i = me%first(thread), me%last(thread)
        v(i) = v(i)*exp(-me%v1(i)*dt_2)
        me%v1(i) = me%v1(i) + dt_2*(me%m(i)*v(i)**2 - me%kT)/me%Q1
        me%v1(i) = me%v1(i)*exp(-me%v2(i)*dt_2)
        me%v2(i) = A*me%v2(i) + B*me%random(thread)%normal() + C*(me%Q1*me%v1(i)**2 - me%kT)
        me%v1(i) = me%v1(i)*exp(-me%v2(i)*dt_2)
        me%v1(i) = me%v1(i) + dt_2*(me%m(i)*v(i)**2 - me%kT)/me%Q1
        v(i) = v(i)*exp(-me%v1(i)*dt_2)
      end do
    end block
    !$omp end parallel

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine isokinetic( i, dt )
        integer,  intent(in) :: i
        real(rb), intent(in) :: dt
        real(rb) :: expmv2dt, mvv, H
        expmv2dt = exp(-me%v2(i)*dt)
        mvv = me%m(i)*v(i)**2
        H = sqrt(me%kT/(mvv + (me%kT - mvv)*expmv2dt*expmv2dt))
        v(i) = v(i)*H
        me%v1(i) = me%v1(i)*H*expmv2dt
      end subroutine isokinetic
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine sinr_integrate

  !-------------------------------------------------------------------------------------------------

  subroutine sinr_boost( me, dt, a, v )
    class(sinr), intent(in)    :: me
    real(rb),    intent(in)    :: dt, a(*)
    real(rb),    intent(inout) :: v(*)
    v(1:me%dof) = v(1:me%dof) + dt*a(1:me%dof)
  end subroutine sinr_boost

  !=================================================================================================

end module mThermostat
