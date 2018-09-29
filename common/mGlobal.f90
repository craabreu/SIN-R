module mGlobal

integer,      parameter :: rb = 8      !< Default number of bytes for real numbers
integer,      parameter :: sl = 512    !< Default character string length
character(3), parameter :: csl = "512" !< String with default character string length

real(rb), parameter :: zero  = 0.0_rb, &
                       one   = 1.0_rb, &
                       two   = 2.0_rb, &
                       three = 3.0_rb, &
                       six   = 6.0_rb, &
                       half  = 0.5_rb, &
                       third = 0.3333333333333333_rb, &
                       fourth = 0.25_rb

integer :: stdout = 6                  !< Standard output unit
integer :: logunit = 0                 !< Output unit for logging

character(2), parameter, private :: delimiters = achar(32)//achar(9)
character,    parameter, private :: comment_mark = "#"

!> A simple random number generator:
type rng
  integer, private :: jsr
  contains
    procedure :: init => rng_init
    procedure :: i32 => rng_i32
    procedure :: letters => rng_letters
end type rng

interface join
  module procedure :: join_with_space
  module procedure :: join_with_sep
end interface

contains

  !=================================================================================================

  subroutine init_log( file )
    character(*), intent(in) :: file
    logunit = 83
    open( unit = logunit, file = file, status = "replace" )
  end subroutine init_log

  !=================================================================================================

  subroutine stop_log
    close( logunit )
    logunit = 0
  end subroutine stop_log

  !=================================================================================================

  subroutine write_msg( prefix, msg )
    character(*), intent(in) :: prefix, msg
    write(stdout,'(A,A)',advance='no') prefix, msg
    if (logunit /= 0) write(logunit,'(A,A)',advance='no') prefix, msg
  end subroutine write_msg

  !=================================================================================================

  subroutine end_line()
    write(stdout,'()')
    if (logunit /= 0) write(logunit,'()')
  end subroutine end_line

  !=================================================================================================

  subroutine writeln( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9, advance )
    character(*), intent(in), optional :: msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    logical,      intent(in), optional :: advance
    if (present(msg))  call write_msg( "", trim(msg) )
    if (present(msg1)) call write_msg( " ", trim(msg1) )
    if (present(msg2)) call write_msg( " ", trim(msg2) )
    if (present(msg3)) call write_msg( " ", trim(msg3) )
    if (present(msg4)) call write_msg( " ", trim(msg4) )
    if (present(msg5)) call write_msg( " ", trim(msg5) )
    if (present(msg6)) call write_msg( " ", trim(msg6) )
    if (present(msg7)) call write_msg( " ", trim(msg7) )
    if (present(msg8)) call write_msg( " ", trim(msg8) )
    if (present(msg9)) call write_msg( " ", trim(msg9) )
    if (present(advance)) then
      if (advance) call end_line
    else
      call end_line
    end if
  end subroutine writeln

  !=================================================================================================

  subroutine warning( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9 )
    character(*), intent(in)           :: msg
    character(*), intent(in), optional :: msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    call write_msg( "WARNING: ", trim(msg) )
    if (present(msg1)) call write_msg( " ", trim(msg1) )
    if (present(msg2)) call write_msg( " ", trim(msg2) )
    if (present(msg3)) call write_msg( " ", trim(msg3) )
    if (present(msg4)) call write_msg( " ", trim(msg4) )
    if (present(msg5)) call write_msg( " ", trim(msg5) )
    if (present(msg6)) call write_msg( " ", trim(msg6) )
    if (present(msg7)) call write_msg( " ", trim(msg7) )
    if (present(msg8)) call write_msg( " ", trim(msg8) )
    if (present(msg9)) call write_msg( " ", trim(msg9) )
    call end_line
  end subroutine warning

  !=================================================================================================

  subroutine error( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9 )
    character(*), intent(in)           :: msg
    character(*), intent(in), optional :: msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    call write_msg( "ERROR: ", trim(msg) )
    if (present(msg1)) call write_msg( " ", trim(msg1) )
    if (present(msg2)) call write_msg( " ", trim(msg2) )
    if (present(msg3)) call write_msg( " ", trim(msg3) )
    if (present(msg4)) call write_msg( " ", trim(msg4) )
    if (present(msg5)) call write_msg( " ", trim(msg5) )
    if (present(msg6)) call write_msg( " ", trim(msg6) )
    if (present(msg7)) call write_msg( " ", trim(msg7) )
    if (present(msg8)) call write_msg( " ", trim(msg8) )
    if (present(msg9)) call write_msg( " ", trim(msg9) )
    call end_line
    stop
  end subroutine error

  !=================================================================================================

  subroutine delete_files( filename )
    character(sl),   intent(in) :: filename(:)
    integer :: ifile, unit, ierr
    do ifile = 1, size(filename)
      open( newunit = unit, file = filename(ifile), status = "old", iostat = ierr )
      if (ierr == 0) close(unit, status = "delete")
    end do
  end subroutine delete_files

  !=================================================================================================

  subroutine rng_init( a, seed )
    class(rng), intent(inout) :: a
    integer,    intent(in)    :: seed
    a%jsr = seed
  end subroutine rng_init

  !=================================================================================================

  function rng_i32( a ) result( i32 )
    class(rng), intent(inout) :: a
    integer                   :: i32
    integer :: jz
    jz = a%jsr
    a%jsr = ieor(a%jsr,ishft(a%jsr, 13))
    a%jsr = ieor(a%jsr,ishft(a%jsr,-17))
    a%jsr = ieor(a%jsr,ishft(a%jsr,  5))
    i32 = jz + a%jsr
  end function rng_i32

  !=================================================================================================

  function rng_letters( a, n ) result( word )
    class(rng), intent(inout) :: a
    integer,    intent(in)    :: n
    character(sl)             :: word
    integer :: i
    word = ""
    do i = 1, n
      word = trim(word)//achar(97+mod(abs(a % i32()),26))
    end do
  end function rng_letters

  !=================================================================================================

  subroutine clean( str )
    character(*), intent(inout) :: str
    integer :: m, n
    m = scan(trim(str),comment_mark)
    if (m > 0) str = str(1:m-1)
    m = verify(str,delimiters)
    if (m == 0) then
      str = ""
    else
      n = verify(str,delimiters,back=.true.)
      str = str(m:n)
    end if
  end subroutine clean

  !=================================================================================================

  subroutine split( str, narg, arg )
    character(*), intent(in)  :: str
    integer,      intent(out) :: narg
    character(*), intent(out) :: arg(:)
    logical :: letter, word
    integer :: i, wlen
    narg = 0
    wlen = 0
    word = .false.
    do i = 1, len_trim(str)
      letter = scan(str(i:i),delimiters) == 0
      if (word) then
        if (letter) then
          wlen = wlen + 1
          arg(narg)(wlen:wlen) = str(i:i)
        else
          arg(narg) = arg(narg)(1:wlen)
          if (narg == size(arg)) return
          word = .false.
        end if
      else
        if (letter) then
          narg = narg + 1
          wlen = 1
          arg(narg)(wlen:wlen) = str(i:i)
          word = .true.
        end if
      end if
    end do
    if (word) arg(narg) = arg(narg)(1:wlen)
  end subroutine split

  !=================================================================================================

  function join_with_space( arg ) result( str )
    character(*), intent(in) :: arg(:)
    character(sl)            :: str
    integer :: i, narg
    narg = size(arg)
    if (narg == 0) then
      str = ""
    else
      str = arg(1)
      do i = 2, narg
        str = trim(str)//" "//arg(i)
      end do
    end if
  end function join_with_space

  !=================================================================================================

  function join_with_sep( arg, sep ) result( str )
    character(*), intent(in) :: arg(:)
    character,    intent(in) :: sep
    character(sl)            :: str
    integer :: i, narg
    narg = size(arg)
    if (narg == 0) then
      str = ""
    else
      str = arg(1)
      do i = 2, narg
        str = trim(str)//sep//arg(i)
      end do
    end if
  end function join_with_sep

  !=================================================================================================

  subroutine str_swap( a, b )
    character(sl), intent(inout) :: a, b
    character(sl) :: aux
    aux = a
    a = b
    b = aux
  end subroutine str_swap

  !=================================================================================================

  subroutine next_command( unit, narg, arg )
    integer,      intent(in)  :: unit
    integer,      intent(out) :: narg
    character(*), intent(out) :: arg(:)
    narg = 0
    call add_items( narg, arg(1:) )
    contains
      recursive subroutine add_items( narg, arg )
        integer,      intent(inout) :: narg
        character(*), intent(inout) :: arg(:)
        integer       :: ioerr, extra
        character(sl) :: line
        read(unit,'(A'//csl//')',iostat=ioerr) line
        call clean( line )
        do while ((ioerr == 0).and.(line == ""))
          read(unit,'(A'//csl//')',iostat=ioerr) line
          call clean( line )
        end do
        if (ioerr == 0) then
          call split( line, extra, arg )
          if ((arg(extra) == "...").or.(arg(extra) == "&")) then
            narg = narg + extra - 1
            call add_items( narg, arg(extra:) )
          else
            narg = narg + extra
          end if
        end if
      end subroutine add_items
  end subroutine next_command

  !=================================================================================================

  function str2int( str ) result( i )
    character(*), intent(in) :: str
    integer                  :: i
    integer :: ioerr
    read(str,*,iostat=ioerr) i
    if (ioerr /= 0) call error( "bad integer" )
  end function str2int

  !=================================================================================================

  function str2real( str ) result( r )
    character(*), intent(in) :: str
    real(rb)                 :: r
    integer :: ioerr
    read(str,*,iostat=ioerr) r
    if (ioerr /= 0) call error( "bad real number" )
  end function str2real

  !=================================================================================================

  elemental function int2str( i ) result( str )
    integer, intent(in) :: i
    character(sl)       :: str
    write(str,*) i
    str = adjustl(str)
  end function int2str

  !=================================================================================================

  elemental function real2str( a ) result( str )
    real(rb), intent(in) :: a
    character(sl)        :: str
    real(4) :: b
    b = a
    write(str,*) b
    str = adjustl(str)
  end function real2str

  !=================================================================================================

  elemental function is_int( arg ) result( ok )
    character(sl), intent(in) :: arg
    logical                   :: ok
    integer :: ioerr, i
    read(arg,*,iostat=ioerr) i
    ok = ioerr == 0
  end function is_int

  !=================================================================================================

  elemental function is_real( arg ) result( ok )
    character(sl), intent(in) :: arg
    logical                   :: ok
    integer  :: ioerr
    real(rb) :: r
    read(arg,*,iostat=ioerr) r
    ok = ioerr == 0
  end function is_real

  !=================================================================================================

end module mGlobal
