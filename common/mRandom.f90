!> This module defines classes for managing random number generation (RNG)
!! algorithms.
!!
!! @author Charlles R. A. Abreu (abreu@eq.ufrj.br)
!! @date Sept 05, 2013
!! @date Sept 13, 2013
module mRandom
implicit none

!> Parameters for the mt19937 (Mersenne Twister) routines:
integer, parameter, private :: allbit_mask =    transfer(Z'ffffffff',1), &
                               upper_mask  =    transfer(Z'80000000',1), &
                               lower_mask  =    transfer(Z'7fffffff',1), &
                               t1_mask     =    transfer(Z'9d2c5680',1), &
                               t2_mask     =    transfer(Z'efc60000',1), &
                               mag01(0:1)  = [0,transfer(Z'9908b0df',1)  ]

!> An abstract class for general random number generation tasks based on
!! algorithms that generate pseudo-random sequences of 32-bit integers.
type, abstract :: i32rng

  integer, private :: kn(0:127)
  real(8), private :: wn(0:127)
  real(8), private :: fn(0:127)

  contains

    !> Performs initialization of the random number generator.
    !! @param[in] seed (integer) a 32-bit integer seed.
    procedure :: setup => i32rng_setup

    !> Seeds the specific random number generator.
    !! @param[in] seed (integer) a 32-bit integer seed.
    procedure(i32rng_init), deferred, private :: init

    !> Generates a random 32-bit integer.
    procedure(i32rng_i32), deferred :: i32

    !> Generates a random 64-bit real number uniformly distributed between
    !! 0 (included) and 1 (excluded).
    procedure :: uniform => i32rng_uniform

    !> Generates a random 64-bit real number normally distributed with
    !! mean 0 and standard deviation 1.
    procedure :: normal => i32rng_normal

    !> Generates a random 64-bit real number distributed according to the
    !! standard gamma distribution with integer order
    procedure :: gamma => i32rng_gamma

    !> Generates a random 32-bit integer number geometrically distributed with
    !! mean n.
    !! @param[in] n (integer) the mean of the geometric distribution used to
    !!            generate the random variate.
    procedure :: geometric => i32rng_geometric

    !> Perfoms a timing test of the random number generator.
    !! @param[in] n (integer) the amount of random numbers to be generated.
    !! @return the time (in seconds) spent to generate the sequence of
    !!         random numbers.
    procedure :: timing => i32rng_timing

    !> Generates a randomly shuffled sequence of integers from 1 to n.
    !! @param[in] n (integer) the size of the sequence to be generated.
    procedure :: sequence => i32rng_sequence

    !> Picks a random label according to a given array of probabilities.
    !! @param[in] P (real array) Probability distribution of labels from 1 to N, where N is the
    !!            number of entries in P. Note: the distribution does not have to be normalized.
    procedure :: label => i32rng_label

end type i32rng

private :: i32rng_init, i32rng_i32

abstract interface

  !> Deferred binding procedure for random number sequence seeding.
  subroutine i32rng_init( a, seed )
    import :: i32rng
    class(i32rng), intent(inout) :: a
    integer,       intent(in)    :: seed
  end subroutine i32rng_init

  !> Deferred binding procedure for 32-bit integer generation.
  function i32rng_i32( a ) result( i32 )
    import :: i32rng
    class(i32rng), intent(inout) :: a
    integer                      :: i32
  end function i32rng_i32

end interface

!> A class for handling pseudo-random number sequences using the SHR3
!! algorithm.
!!
!! @remark Reference: Marsaglia and Tsang, Journal of Statistical
!!         Software 5(8), 2000.
type, extends(i32rng) :: shr3
  integer, private :: jsr
  contains
    ! Deferred bindings:
    procedure :: init => shr3_init
    procedure :: i32 => shr3_i32
end type shr3

!> A class for handling pseudo-random number sequences using the KISS
!! algorithm.
!!
!! @remark Reference: Marsaglia and Zaman, The KISS generator,
!!         Technical report, Departmentof Statistics,
!!         Florida State University, Tallahassee, FL, USA, 1993.
type, extends(i32rng) :: kiss
  integer, private :: x, y, z, w
  contains
    ! Deferred bindings:
    procedure :: init => kiss_init
    procedure :: i32 => kiss_i32
end type kiss

!> A class for handling pseudo-random number sequences using the MT19937
!! (Mersenne Twister) algorithm of Matsumoto and Nishimura (1998).
!!
!! @remark Reference: Matsumoto and Nishimura,
!!         ACM Trans. Model. Comput. Simul. 8, 3-30 (1998).
type, extends(i32rng) :: mt19937
  integer, private :: mti, mt(0:623)
  contains
    ! Deferred bindings:
    procedure :: i32 => mt19937_i32
    procedure :: init  => mt19937_init
end type mt19937

!> A class for handling pseudo-random number sequences using the xoroshiro128+ algorithm.
type, extends(i32rng) :: xoroshiro128plus
  integer(8), private :: s(2) = [123456789_8, 987654321_8]
  logical :: saved = .false.
  real(8) :: rnor
  contains
    procedure :: init => xoroshiro128plus_init
    procedure :: i32 => xoroshiro128plus_i32
    procedure :: i64 => xoroshiro128plus_i64
    procedure :: uniform => xoroshiro128plus_uniform
    procedure :: normal => xoroshiro128plus_normal
end type xoroshiro128plus

private :: i32rng_uniform, i32rng_normal, i32rng_timing
private :: shr3_init, shr3_i32
private :: kiss_init, kiss_i32
private :: mt19937_init, mt19937_i32
private :: xoroshiro128plus_init, xoroshiro128plus_i32

contains
  !-----------------------------------------------------------------------------
  !                                GENERAL i32rng
  !-----------------------------------------------------------------------------
  subroutine i32rng_setup( a, seed )
    class(i32rng), intent(inout) :: a
    integer,       intent(in)    :: seed
    integer :: i
    real(8), parameter :: m1 = 2147483648.0_8
    real(8) :: q, dn, tn, vn
    ! Seed the random number generator:
    call a%init( seed )
    ! Generate tables for the ziggurat algorithm:
    dn = 3.442619855899_8
    tn = 3.442619855899_8
    vn = 0.00991256303526217_8
    q = vn*exp(0.5_8*dn*dn)
    a%kn(0) = int((dn/q)*m1)
    a%kn(1) = 0
    a%wn(0) = q/m1
    a%wn(127) = dn/m1
    a%fn(0) = 1.0_8
    a%fn(127) = exp( -0.5_8*dn*dn )
    do i = 126, 1, -1
      dn = sqrt( -2.0_8 * log( vn/dn + exp( -0.5_8*dn*dn ) ) )
      a%kn(i+1) = int((dn/tn)*m1)
      tn = dn
      a%fn(i) = exp(-0.5_8*dn*dn)
      a%wn(i) = dn/m1
    end do
  end subroutine i32rng_setup
  !-----------------------------------------------------------------------------
  function i32rng_uniform( a ) result( uni )
    class(i32rng), intent(inout) :: a
    real(8)                   :: uni
    uni = 0.5_8 + 0.2328306e-9_8 * a%i32()
  end function i32rng_uniform
  !-----------------------------------------------------------------------------
  function i32rng_normal( a ) result( rnor )
    class(i32rng), intent(inout) :: a
    real(8)                      :: rnor
    real(8), parameter :: r = 3.442620_8
    integer :: hz, iz
    real(8) :: x, y
    hz = a%i32()
    iz = iand( hz, 127 )
    if ( abs( hz ) < a%kn(iz) ) then
      rnor = hz * a%wn(iz)
    else
      do
        if ( iz == 0 ) then
          do
            x = -0.2904764_8* log( a%uniform() )
            y = -log( a%uniform() )
            if ( y+y >= x*x ) exit
          end do
          rnor = r + x
          if ( hz <= 0 ) rnor = -rnor
          return
        end if
        x = hz * a%wn(iz)
        if ( a%fn(iz) + a%uniform()*(a%fn(iz-1) - a%fn(iz)) < exp(-0.5_8*x*x) ) then
          rnor = x
          return
        end if
        hz = a%i32()
        iz = iand( hz, 127 )
        if ( abs( hz ) < a%kn(iz) ) then
          rnor = hz * a%wn(iz)
          return
        end if
      end do
   end if
  end function i32rng_normal
  !-----------------------------------------------------------------------------
  function i32rng_gamma( a, alpha ) result( ans )
    class(i32rng), intent(inout) :: a
    real(8),       intent(in)    :: alpha
    real(8)                      :: ans
    ! Marsaglia and Tsang, "A Simple Method for Generating Gamma Variables"
    ! ACM Transactions on Mathematical Software, Vol 26, No 3, 2000, pp. 363-372.
    real(8) :: u, d, c, x, xx, v
    if (alpha >= 1.0_8) then
      d = alpha - 1.0_8/3.0_8
      c = 1.0_8/sqrt(9.0_8*d)
      do
        x = a%normal()
        v = 1.0_8 + c*x
        if (v <= 0.0_8) cycle
        v = v*v*v
        ans = d*v
        u = a%uniform()
        xx = x*x
        if (u < 1.0_8 - 0.0331_8*xx*xx) return
        if (log(u) < 0.5_8*xx + d - ans + d*log(v)) return
      end do
    else
      ans = a%gamma(alpha + 1.0_8)*a%uniform()**(1.0_8/alpha)
    end if
  end function i32rng_gamma
  !-----------------------------------------------------------------------------
  function i32rng_geometric( a, n ) result( igeo )
    class(i32rng), intent(inout) :: a
    integer,       intent(in)    :: n
    integer                      :: igeo
    igeo = ceiling( log(0.5_8 + 0.2328306e-9_8*a%i32())/log(1.0_8 - 1.0_8/n) )
  end function i32rng_geometric
  !-----------------------------------------------------------------------------
  function i32rng_timing( a, n ) result( time )
    class(i32rng), intent(inout) :: a
    integer,       intent(in)    :: n
    real(8)                      :: time
    integer :: i
    real(8) :: start, i32
    call cpu_time( start )
    do i = 1, n
      i32 = a%uniform()
    end do
    call cpu_time( time )
    time = time - start
  end function i32rng_timing
  !-----------------------------------------------------------------------------
  function i32rng_sequence( a, n ) result( seq )
    class(i32rng), intent(inout) :: a
    integer,       intent(in)    :: n
    integer                      :: seq(n)
    integer :: i, m, iaux
    real(8) :: r(n), raux
    logical :: replaced
    do i = 1, n
      seq(i) = i
      r(i) = a%uniform()
    end do
    replaced = .true.
    m = n - 1
    do while (replaced)
      replaced = .false.
      do i = 1, m
        if ( r(i) > r(i+1) ) then
          raux = r(i)
          r(i) = r(i+1)
          r(i+1) = raux
          iaux = seq(i)
          seq(i) = seq(i+1)
          seq(i+1) = iaux
          replaced = .true.
        end if
      end do
      m = m - 1
    end do
  end function i32rng_sequence
  !-----------------------------------------------------------------------------
  function i32rng_label( a, P ) result( k )
    class(i32rng), intent(inout) :: a
    real(8),       intent(in)    :: P(:)
    integer                      :: k
    real(8) :: r, acc, pi(size(P))
    pi = P/sum(P)
    r = a%uniform()
    k = 1
    acc = pi(1)
    do while (r > acc)
      k = k + 1
      acc = acc + pi(k)
    end do
  end function i32rng_label
  !-----------------------------------------------------------------------------
  !                                  SHR3
  !-----------------------------------------------------------------------------
  subroutine shr3_init( a, seed )
    class(shr3), intent(inout) :: a
    integer,     intent(in)    :: seed
    a%jsr = seed
  end subroutine shr3_init
  !-----------------------------------------------------------------------------
  function shr3_i32( a ) result( i32 )
    class(shr3), intent(inout) :: a
    integer                    :: i32
    integer :: jz
    jz = a%jsr
    a%jsr = ieor(a%jsr,ishft(a%jsr, 13))
    a%jsr = ieor(a%jsr,ishft(a%jsr,-17))
    a%jsr = ieor(a%jsr,ishft(a%jsr,  5))
    i32 = jz + a%jsr
  end function shr3_i32
  !-----------------------------------------------------------------------------
  !                                  KISS
  !-----------------------------------------------------------------------------
  subroutine kiss_init( a, seed )
    class(kiss), intent(inout) :: a
    integer,     intent(in)    :: seed
    a%x = m(m(m(seed,13),-17),5)
    a%y = m(m(m(a%x,13),-17),5)
    a%z = m(m(m(a%y,13),-17),5)
    a%w = m(m(m(a%z,13),-17),5)
    contains
      integer function m( k, n )
        integer, intent(in) :: k, n
        m = ieor(k,ishft(k,n))
      end function m
  end subroutine kiss_init
  !-----------------------------------------------------------------------------
  function kiss_i32( a ) result( i32 )
    class(kiss), intent(inout) :: a
    integer                    :: i32
    a%x = 69069*a%x + 1327217885
    a%y = ieor(a%y,ishft(a%y,13))
    a%y = ieor(a%y,ishft(a%y,-17))
    a%y = ieor(a%y,ishft(a%y,5))
    a%z = 18000*iand(a%z,65535) + ishft(a%z,-16)
    a%w = 30903*iand(a%w,65535) + ishft(a%w,-16)
    i32 = a%x + a%y + ishft(a%z,16) + a%w
  end function kiss_i32
  !-----------------------------------------------------------------------------
  !                          MT19937 (MERSENNE TWISTER)
  !-----------------------------------------------------------------------------
  subroutine mt19937_init( a, seed )
    class(mt19937), intent(inout) :: a
    integer,        intent(in)    :: seed
    integer :: i
    a%mt(0) = iand(seed,allbit_mask)
    do i = 1, 623
      a%mt(i) = 1812433253*ieor(a%mt(i-1),ishft(a%mt(i-1),-30)) + i
      a%mt(i) = iand(a%mt(i),allbit_mask)
    end do
    a%mti = 624
  end subroutine mt19937_init
  !-----------------------------------------------------------------------------
  function mt19937_i32( a ) result( i32 )
    class(mt19937), intent(inout) :: a
    integer                       :: i32
    integer :: y, i
    if (a%mti >= 624) then
      do i = 0, 226
        y = ior(iand(a%mt(i),upper_mask),iand(a%mt(i+1),lower_mask))
        a%mt(i) = ieor(ieor(a%mt(i+397),ishft(y,-1)),mag01(iand(y,1)))
      end do
      do i = 227, 622
        y = ior(iand(a%mt(i),upper_mask),iand(a%mt(i+1),lower_mask))
        a%mt(i) = ieor(ieor(a%mt(i+227),ishft(y,-1)),mag01(iand(y,1)))
      end do
      y = ior(iand(a%mt(623),upper_mask),iand(a%mt(0),lower_mask))
      a%mt(i) = ieor(ieor(a%mt(396),ishft(y,-1)),mag01(iand(y,1)))
      a%mti = 0
    endif
    y = a%mt(a%mti)
    a%mti = a%mti + 1
    y = ieor(y,ishft(y,-11))
    y = ieor(y,iand(ishft(y,7),t1_mask))
    y = ieor(y,iand(ishft(y,15),t2_mask))
    i32 = ieor(y,ishft(y,-18))
  end function mt19937_i32
  !-----------------------------------------------------------------------------
  !                                  Xoroshiro128+
  !-----------------------------------------------------------------------------
  subroutine xoroshiro128plus_init( a, seed )
    class(xoroshiro128plus), intent(inout) :: a
    integer,                 intent(in)    :: seed
    integer(4) :: b(2), c(2)
    integer(8) :: mold
    type(shr3) :: random
    call random % setup(seed)
    b = [random%i32(), random%i32()]
    c = [random%i32(), random%i32()]
    a%s = [transfer(b, mold), transfer(c, mold)]
  end subroutine xoroshiro128plus_init
  !-----------------------------------------------------------------------------
  function xoroshiro128plus_i32( a ) result( i32 )
    class(xoroshiro128plus), intent(inout) :: a
    integer(4)                             :: i32
    i32 = int(a%i64(), 4)
  end function xoroshiro128plus_i32
  !-----------------------------------------------------------------------------
  function xoroshiro128plus_i64( a ) result( i64 )
    class(xoroshiro128plus), intent(inout) :: a
    integer(8)                             :: i64
    integer(8) :: t(2)
    t = a%s
    i64 = t(1) + t(2)
    t(2) = ieor(t(1), t(2))
    a%s(1) = ieor(ieor(ior(shiftl(t(1), 55), shiftr(t(1), 9)), t(2)), shiftl(t(2), 14))
    a%s(2) = ior(shiftl(t(2), 36), shiftr(t(2), 28))
  end function xoroshiro128plus_i64
  !-----------------------------------------------------------------------------
  function xoroshiro128plus_uniform( a ) result( uni )
    class(xoroshiro128plus), intent(inout) :: a
    real(8)                                :: uni
    integer(8), parameter :: z = shiftl(1023_8, 52)
    uni = transfer(ior(z, shiftr(a%i64(), 12)), uni) - 1.0_8
  end function xoroshiro128plus_uniform
  !-----------------------------------------------------------------------------
  function xoroshiro128plus_normal( a ) result( rnor )
    class(xoroshiro128plus), intent(inout) :: a
    real(8)                                :: rnor
    real(8) :: x, y, sum_sq, factor
    if (a%saved) then
      rnor = a%rnor
      a%saved = .false.
    else
      call draw_values
      do while ((sum_sq >= 1.0_8) .or. (sum_sq == 0.0_8))
        call draw_values
      end do
      factor = sqrt(-2.0_8*log(sum_sq)/sum_sq)
      rnor = x*factor
      a%rnor = y*factor
      a%saved = .true.
    end if
    contains
      subroutine draw_values
        x = 2.0_8 * a%uniform() - 1.0_8
        y = 2.0_8 * a%uniform() - 1.0_8
        sum_sq = x*x + y*y
      end subroutine draw_values
  end function xoroshiro128plus_normal
  !-----------------------------------------------------------------------------
end module mRandom
