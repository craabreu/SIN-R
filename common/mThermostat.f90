module mThermostat

use mGlobal
use mRandom
use iso_c_binding
use omp_lib

implicit none

type :: tThermostat
  integer  :: d                            ! Dimensionality of the system
  integer  :: N                            ! Number of particles in the system
  integer  :: dof                          ! Number of degrees of freedom
  integer  :: seed                         ! Seed for random numbers
  integer  :: nthreads                     ! Number of parallel threads
  real(rb) :: kT                           ! Temperature in energy units
  real(rb) :: tau                          ! Characteristic time scale
  logical  :: columns = .true.             ! True if each particle fills a column of a matrix
  integer  :: actualDof                    ! Actual number of degrees of freedom
  real(rb),      allocatable :: invm(:)    ! Inverse mass associated to each degree of freedom
  integer,       allocatable :: first(:)   ! The first degree of freedom of each thread
  integer,       allocatable :: last(:)    ! The last degree of freedom of each thread
  character(sl), allocatable :: keyword(:) ! List of keywords for finishing thermostat setup
  real(rb),      allocatable :: value(:)   ! Values associated to the list of keywords
  type(xrsr128), allocatable :: random(:)  ! A random number generator for each thread
  contains
    procedure :: setup => tThermostat_setup
    procedure :: get_value => tThermostat_get_value
    procedure :: finish_setup => tThermostat_finish_setup
    procedure :: initialize => tThermostat_initialize
    procedure :: zero_total_momenta => tThermostat_zero_total_momenta
    procedure :: compute_kT => tThermostat_compute_kT
    procedure :: external => tThermostat_external
    procedure :: boost => tThermostat_boost
    procedure :: internal => tThermostat_internal
    procedure :: energy => tThermostat_energy
end type tThermostat

!---------------------------------------------------------------------------------------------------
!                                  F A K E    T H E R M O S T A T
!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: nve
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
    procedure :: finish_setup => nhc_finish_setup
    procedure :: energy => nhc_energy
    procedure :: internal => nhc_internal
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
  procedure :: internal => csvr_internal
end type csvr

!---------------------------------------------------------------------------------------------------
!                                          L A N G E V I N
!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: langevin
  real(rb) :: gamma
  real(rb), allocatable :: invSqrtM(:)
contains
  procedure :: finish_setup => langevin_finish_setup
  procedure :: internal => langevin_internal
end type langevin

!---------------------------------------------------------------------------------------------------
!                                        I S O K I N E T I C
!---------------------------------------------------------------------------------------------------

type, extends(tThermostat) :: isokinetic
  real(rb) :: twoKE
contains
  procedure :: finish_setup => isokinetic_finish_setup
  procedure :: initialize => isokinetic_initialize
  procedure :: boost => isokinetic_boost
end type isokinetic

!---------------------------------------------------------------------------------------------------
!        S T O C H A S T I C    I S O K I N E T I C    N O S E - H O O V E R    R E S P A
!---------------------------------------------------------------------------------------------------

type, abstract, extends(tThermostat) :: sinr
  real(rb) :: Q1       ! Inertial parameter
  real(rb) :: Q2       ! Inertial parameter
  real(rb) :: halfQ1   ! Half inertial parameter
  real(rb) :: gamma
  real(rb), allocatable :: v1(:)
  real(rb), allocatable :: v2(:)
contains
  procedure :: finish_setup => sinr_finish_setup
  procedure :: initialize => sinr_initialize
  procedure :: boost => sinr_boost
end type sinr

type, extends(sinr) :: sinr_new
contains
  procedure :: internal => sinr_new_internal
end type sinr_new

type, extends(sinr) :: sinr_xo_respa
contains
  procedure :: internal => sinr_xo_respa_internal
  procedure :: external => sinr_xo_respa_external
end type sinr_xo_respa

!---------------------------------------------------------------------------------------------------

contains

  !=================================================================================================
  !                                A L L    T H E R M O S T A T S
  !=================================================================================================

  subroutine tThermostat_setup( me, d, N, mass, kT, tau, seed, threads, keyword, value, columns )
    class(tThermostat),     intent(inout) :: me
    integer,                intent(in)    :: d, N
    real(rb),               intent(in)    :: mass(*), kT, tau
    integer,                intent(in)    :: seed, threads
    character(*), optional, intent(in)    :: keyword(:)
    real(rb),     optional, intent(in)    :: value(:)
    logical,      optional, intent(in)    :: columns

    integer :: i, perThread, newseed
    real(rb), allocatable :: m(:,:)

    me%d = d
    me%N = N
    me%dof = d*N
    me%kT = kT
    me%tau = tau
    if (present(columns)) me%columns = columns
    if (me%columns) then
      allocate( m(d,N) )
      forall (i=1:N) m(:,i) = mass(i)
    else
      allocate( m(N,d) )
      forall (i=1:N) m(i,:) = mass(i)
    end if
    me%invm = reshape(one/m, [me%dof])
    me%actualDof = me%dof

    me%nthreads =  threads
    perThread = (me%dof + threads - 1)/threads
    me%first = [((i - 1)*perThread + 1, i=1, threads)]
    me%last = [(min(i*perThread, me%dof), i=1, threads)]
    allocate( me%random(threads) )
    newseed = seed
    do i = 1, threads
      call me % random(i) % setup( newseed )
      newseed = me%random(i)%i32()
    end do

    if (present(keyword).and.present(value)) then
      if (size(keyword) /= size(value)) then
        call error( "keywords and values in thermostat setup do not match" )
      end if
      me%keyword = keyword
      me%value = value
    else
      allocate( me%keyword(0), me%value(0) )
    end if

    call me % finish_setup()
  end subroutine tThermostat_setup

  !-------------------------------------------------------------------------------------------------

  subroutine tThermostat_finish_setup( me )
    class(tThermostat), intent(inout) :: me
  end subroutine tThermostat_finish_setup

  !-------------------------------------------------------------------------------------------------

  function tThermostat_get_value( me, keyword, default ) result( answer )
    class(tThermostat), intent(in) :: me
    character(*),       intent(in) :: keyword
    real(rb),           intent(in) :: default
    real(rb)                       :: answer

    integer :: i
    logical :: found

    i = 0
    found = .false.
    do while ((i < size(me%keyword)).and.(.not.found))
      i = i + 1
      found = keyword == me%keyword(i)
    end do
    answer = merge(me%value(i), default, found)
  end function tThermostat_get_value

  !-------------------------------------------------------------------------------------------------

  subroutine tThermostat_zero_total_momenta( me, p )
    class(tThermostat), intent(inout) :: me
    real(rb), target,   intent(inout) :: p(*)

    integer :: i
    real(rb) :: vcm(me%d)
    real(rb), pointer :: pp(:,:)

    if (me%columns) then
      call c_f_pointer(c_loc(p(1)), pp, [me%d, me%N])
      vcm = sum(pp,dim=2)/sum(one/me%invm)
      p(1:me%dof) = p(1:me%dof)*me%invm
      forall (i=1:me%N) pp(:,i) = pp(:,i) - vcm
    else
      call c_f_pointer(c_loc(p(1)), pp, [me%N, me%d])
      vcm = sum(pp,dim=1)/sum(one/me%invm)
      p(1:me%dof) = p(1:me%dof)*me%invm
      forall (i=1:me%N) pp(i,:) = pp(i,:) - vcm
    end if
    p(1:me%dof) = p(1:me%dof)/me%invm
  end subroutine tThermostat_zero_total_momenta

  !-------------------------------------------------------------------------------------------------

  subroutine tThermostat_initialize( me, p )
    class(tThermostat), intent(inout) :: me
    real(rb),           intent(inout) :: p(*)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, i
      thread = omp_get_thread_num() + 1
      do i = me%first(thread), me%last(thread)
        p(i) = sqrt(me%kT/me%invm(i))*me%random(thread)%normal()
      end do
    end block
    !$omp end parallel
    call me % zero_total_momenta( p )
    p(1:me%dof) = p(1:me%dof)*sqrt(me%actualDof*me%kT/sum(me%invm*p(1:me%dof)**2))

  end subroutine tThermostat_initialize

  !-------------------------------------------------------------------------------------------------

  function tThermostat_compute_kT( me, p ) result( kT )
    class(tThermostat), intent(in) :: me
    real(rb),           intent(in) :: p(*)
    real(rb)                       :: kT
    kT = sum(me%invm*p(1:me%dof)**2)/me%actualDof
  end function tThermostat_compute_kT

  !-------------------------------------------------------------------------------------------------

  subroutine tThermostat_external( me, timestep, p )
    class(tThermostat), intent(inout) :: me
    real(rb),           intent(in)    :: timestep
    real(rb),           intent(inout) :: p(*)
  end subroutine tThermostat_external

  !-------------------------------------------------------------------------------------------------

  subroutine tThermostat_boost( me, dt, F, p )
    class(tThermostat), intent(inout) :: me
    real(rb),           intent(in)    :: dt, F(*)
    real(rb),           intent(inout) :: p(*)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, i
      thread = omp_get_thread_num() + 1
      forall (i = me%first(thread):me%last(thread))
        p(i) = p(i) + dt*F(i)
      end forall
    end block
    !$omp end parallel
  end subroutine tThermostat_boost

  !-------------------------------------------------------------------------------------------------

  subroutine tThermostat_internal( me, timestep, p )
    class(tThermostat), intent(inout) :: me
    real(rb),           intent(in)    :: timestep
    real(rb),           intent(inout) :: p(*)
  end subroutine tThermostat_internal

  !-------------------------------------------------------------------------------------------------

  function tThermostat_energy( me ) result( energy )
    class(tThermostat), intent(in) :: me
    real(rb)                       :: energy
    energy = zero
  end function tThermostat_energy

  !=================================================================================================
  !                              N O S É - H O O V E R     C H A I N
  !=================================================================================================

  subroutine nhc_finish_setup( me )
    class(nhc),   intent(inout) :: me

    integer :: nchains, nloops

    nchains = nint(me%get_value("nchains", one))
    nloops = nint(me%get_value("nloops", one))
    me%NC = nchains
    me%actualDof = me%dof - 3
    me%LKT = me%actualDof*me%kT
    me%nloops = nloops
    allocate( me%InvQ(nchains), me%eta(nchains), me%p(nchains) )
    me%InvQ(1) = one/(me%LKT*me%tau**2)
    me%InvQ(2:nchains) = one/(me%kT*me%tau**2)
    me%eta = zero
    me%p = zero

  end subroutine nhc_finish_setup

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

  subroutine nhc_internal( me, timestep, p )
    class(nhc), intent(inout) :: me
    real(rb),   intent(in)    :: timestep
    real(rb),   intent(inout) :: p(*)

    integer  :: i, j
    real(rb) :: TwoKE, dt, dt_2, twodt, alpha, alphaSum, factor

    TwoKE = sum(me%invm*p(1:me%dof)**2)
    dt = timestep/me%nloops
    dt_2 = half*dt
    twodt = two*dt
    alphaSum = zero
    factor = one
    do i = 1, me%nloops
      me%p(me%NC) = me%p(me%NC) + (me%p(me%NC-1)**2*me%InvQ(me%NC-1) - me%kT)*dt_2
      do j = me%NC-1, 2, -1
        call internal( me, j, me%p(j+1)*me%InvQ(j+1), me%p(j-1)**2*me%InvQ(j-1) - me%kT, dt_2 )
      end do
      call internal( me, 1, me%p(2)*me%InvQ(2), factor**2*twoKE - me%LkT, dt_2 )
      alpha = me%p(1)*me%InvQ(1)
      alphaSum = alphaSum + alpha
      factor = exp(-alphaSum*dt)
      call internal( me, 1, me%p(2)*me%InvQ(2), factor**2*twoKE - me%LkT, dt_2 )
      do j = 2, me%NC-1
        call internal( me, j, me%p(j+1)*me%InvQ(j+1), me%p(j-1)**2*me%InvQ(j-1) - me%kT, dt_2 )
      end do
      me%p(me%NC) = me%p(me%NC) + (me%p(me%NC-1)**2*me%InvQ(me%NC-1) - me%kT)*dt_2
    end do
    me%eta(1) = me%eta(1) + alphaSum*dt
    p(1:me%dof) = p(1:me%dof)*exp(-dt*alphaSum)

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine internal( me, j, alpha, G, dt )
        class(nhc), intent(inout) :: me
        integer,    intent(in) :: j
        real(rb),   intent(in) :: alpha, G, dt
        me%p(j) = me%p(j) + (G - alpha*me%p(j))*phi(alpha*dt)*dt
        me%eta(j+1) = me%eta(j+1) + alpha*dt
      end subroutine internal
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
  end subroutine nhc_internal

  !=================================================================================================
  !                  S T O C H A S T I C    V E L O C I T Y    R E S C A L I N G
  !=================================================================================================

  subroutine csvr_finish_setup( me )
    class(csvr),  intent(inout) :: me
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

  subroutine csvr_internal( me, timestep, p )
    class(csvr), intent(inout) :: me
    real(rb),    intent(in)    :: timestep
    real(rb),    intent(inout) :: p(*)

    real(rb) :: TwoKE, alphaSq, R1, x, sumRs, A, B, C

    TwoKE = sum(me%invm*p(1:me%dof)**2)
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
    p(1:me%dof) = p(1:me%dof)*sqrt(alphaSq)

  end subroutine csvr_internal

  !=================================================================================================
  !                                         L A N G E V I N
  !=================================================================================================

  subroutine langevin_finish_setup( me )
    class(langevin), intent(inout) :: me

    me%gamma = me%get_value("friction", one/me%tau)
    allocate( me%invSqrtM(me%dof) )
    me%invSqrtM = one/sqrt(me%invm)

  end subroutine langevin_finish_setup

  !-------------------------------------------------------------------------------------------------

  subroutine langevin_internal( me, timestep, p )
    class(langevin), intent(inout) :: me
    real(rb),    intent(in)    :: timestep
    real(rb),    intent(inout) :: p(*)

    real(rb) :: A, B

    A = exp(-me%gamma*timestep)
    B = sqrt(me%kT*(one - A*A))

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, i
      thread = omp_get_thread_num() + 1
      do i = me%first(thread), me%last(thread)
        p(i) = A*p(i) + B*me%invSqrtM(i)*me%random(thread)%normal()
      end do
    end block
    !$omp end parallel
  end subroutine langevin_internal

  !=================================================================================================
  !                                       I S O K I N E T I C
  !=================================================================================================

  subroutine isokinetic_finish_setup( me )
    class(isokinetic), intent(inout) :: me
    me%twoKE = me%dof*me%kT
  end subroutine isokinetic_finish_setup

  !-------------------------------------------------------------------------------------------------

  subroutine isokinetic_initialize( me, p )
    class(isokinetic), intent(inout) :: me
    real(rb),          intent(inout) :: p(*)

    real(rb) :: twokin(me%nthreads)

    call me % tThermostat % initialize( p )

    !$omp parallel num_threads(me%nthreads)
    block
      integer  :: thread, first, last
      real(rb) :: twoKE
      thread = omp_get_thread_num() + 1
      first = me%first(thread)
      last = me%last(thread)
      twokin(thread) = sum(me%invm(first:last)*p(first:last)**2)
      !$omp barrier
      twoKE = sum(twokin)
      p(first:last) = p(first:last)*sqrt(me%twoKE/twoKE)
    end block
    !$omp end parallel
  end subroutine isokinetic_initialize

  !-------------------------------------------------------------------------------------------------

  subroutine isokinetic_boost( me, dt, F, p )
    class(isokinetic), intent(inout) :: me
    real(rb),    intent(in)    :: dt, F(*)
    real(rb),    intent(inout) :: p(*)

    real(rb) :: Fv(me%nthreads), Fa(me%nthreads)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, first, last
      real(rb) :: lb0, b2, b, coshbt, sinhbt, s, spinv

      thread = omp_get_thread_num() + 1
      first = me%first(thread)
      last = me%last(thread)
      Fv(thread) = sum(me%invm(first:last)*F(first:last)*p(first:last))
      Fa(thread) = sum(me%invm(first:last)*F(first:last)**2)
      !$omp barrier
      lb0 = sum(Fv)/me%twoKE
      b2 = sum(Fa)/me%twoKE
      b = sqrt(b2)
      coshbt = cosh(b*dt)
      sinhbt = sinh(b*dt)
      s = (lb0*(coshbt - one) + b*sinhbt)/b2
      spinv = b/(lb0*sinhbt + b*coshbt)
      p(first:last) = (p(first:last) + F(first:last)*s)*spinv
    end block
    !$omp end parallel

  end subroutine isokinetic_boost

  !=================================================================================================
  !        S T O C H A S T I C    I S O K I N E T I C    N O S E - H O O V E R    R E S P A
  !=================================================================================================

  subroutine sinr_finish_setup( me )
    class(sinr),  intent(inout) :: me

    me%gamma = me%get_value("friction", one/me%tau)
    me%Q1 = me%kT*me%tau**2
    me%Q2 = me%kT*me%tau**2
    me%halfQ1 = half*me%Q1
    allocate( me%v1(me%dof), me%v2(me%dof), source = zero )
  end subroutine sinr_finish_setup

  !-------------------------------------------------------------------------------------------------

  subroutine sinr_initialize( me, p )
    class(sinr), intent(inout) :: me
    real(rb),    intent(inout) :: p(*)

    real(rb) :: twokin(me%nthreads)

    call me % tThermostat % initialize( p )

    !$omp parallel num_threads(me%nthreads)
    block
      integer  :: thread, first, last, i
      real(rb) :: twoKE, sigma1, sigma2, factor
      thread = omp_get_thread_num() + 1
      first = me%first(thread)
      last = me%last(thread)
      twokin(thread) = sum(me%invm(first:last)*p(first:last)**2)
      !$omp barrier
      twoKE = sum(twokin)
      p(first:last) = 0.25_rb*p(first:last)
      sigma1 = sqrt(me%kT/me%Q1)
      sigma2 = sqrt(me%kT/me%Q2)
      do i = first, last
        me%v1(i) = sigma1*me%random(thread)%normal()
        me%v2(i) = sigma2*me%random(thread)%normal()
        factor = sqrt(me%kT/(me%invm(i)*p(i)**2 + me%halfQ1*me%v1(i)**2))
        p(i) = factor*p(i)
        me%v1(i) = factor*me%v1(i)
      end do
    end block
    !$omp end parallel
  end subroutine sinr_initialize

  !-------------------------------------------------------------------------------------------------

  subroutine sinr_boost( me, dt, F, p )
    class(sinr), intent(inout) :: me
    real(rb),    intent(in)    :: dt, F(*)
    real(rb),    intent(inout) :: p(*)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, i
      real(rb) :: lb0, bsq, b, bdt, x, y, s, spinv

      thread = omp_get_thread_num() + 1
      do i = me%first(thread), me%last(thread)
        lb0 = me%invm(i)*F(i)*p(i)/me%kT
        bsq = me%invm(i)*F(i)**2/me%kT
        b = sqrt(bsq)
        bdt = b*dt
        x = coshxm1_x2(bdt)
        y = sinhx_x(bdt)
        s = (lb0*x*dt + y)*dt
        spinv = one/(lb0*y*dt + bdt*bdt*x + one)
        p(i) = (p(i) + F(i)*s)*spinv
        me%v1(i) = me%v1(i)*spinv
      end do
    end block
    !$omp end parallel

  end subroutine sinr_boost

  !-------------------------------------------------------------------------------------------------

  subroutine sinr_new_internal( me, timestep, p )
    class(sinr_new), intent(inout) :: me
    real(rb),        intent(in)    :: timestep
    real(rb),        intent(inout) :: p(*)

    real(rb) :: dt_2, A, B, C

    dt_2 = half*timestep
    A = exp(-me%gamma*timestep)
    B = sqrt(me%kT*(one - A*A)/me%Q2)
    C = (one - A)/(me%gamma*me%Q2)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, i
      thread = omp_get_thread_num() + 1
      do i = me%first(thread), me%last(thread)
        call isokinetic_part( i, dt_2 )
        me%v2(i) = A*me%v2(i) + B*me%random(thread)%normal() + C*(me%Q1*me%v1(i)**2 - me%kT)
        call isokinetic_part( i, dt_2 )
      end do
    end block
    !$omp end parallel

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine isokinetic_part( i, dt )
        integer,  intent(in) :: i
        real(rb), intent(in) :: dt
        real(rb) :: x, v1x, H
        x = exp(-me%v2(i)*dt)
        v1x = me%v1(i)*x
        H = sqrt(me%kT/(me%invm(i)*p(i)**2 + me%halfQ1*v1x*v1x))
        p(i) = p(i)*H
        me%v1(i) = v1x*H
      end subroutine isokinetic_part
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine sinr_new_internal

  !-------------------------------------------------------------------------------------------------

  subroutine sinr_xo_respa_internal( me, timestep, p )
    class(sinr_xo_respa), intent(inout) :: me
    real(rb),             intent(in)    :: timestep
    real(rb),             intent(inout) :: p(*)

    real(rb) :: A, B

    A = exp(-me%gamma*timestep)
    B = sqrt(me%kT*(one - A*A)/me%Q2)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, i
      thread = omp_get_thread_num() + 1
      do i = me%first(thread), me%last(thread)
        me%v2(i) = A*me%v2(i) + B*me%random(thread)%normal()
      end do
    end block
    !$omp end parallel
  end subroutine sinr_xo_respa_internal

  !-------------------------------------------------------------------------------------------------

  subroutine sinr_xo_respa_external( me, timestep, p )
    class(sinr_xo_respa), intent(inout) :: me
    real(rb),             intent(in)    :: timestep
    real(rb),             intent(inout) :: p(*)

    real(rb) :: C

    C = half*timestep/me%Q2

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, i
      thread = omp_get_thread_num() + 1
      do i = me%first(thread), me%last(thread)
        me%v2(i) = me%v2(i) + C*(me%Q1*me%v1(i)**2 - me%kT)
        call isokinetic_part( i, timestep )
        me%v2(i) = me%v2(i) + C*(me%Q1*me%v1(i)**2 - me%kT)
      end do
    end block
    !$omp end parallel

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine isokinetic_part( i, dt )
        integer,  intent(in) :: i
        real(rb), intent(in) :: dt
        real(rb) :: x, v1x, H
        x = exp(-me%v2(i)*dt)
        v1x = me%v1(i)*x
        H = sqrt(me%kT/(me%invm(i)*p(i)**2 + me%halfQ1*v1x*v1x))
        p(i) = p(i)*H
        me%v1(i) = v1x*H
      end subroutine isokinetic_part
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine sinr_xo_respa_external

  !=================================================================================================

end module mThermostat
