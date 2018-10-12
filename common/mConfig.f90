module mConfig

use mGlobal

implicit none

character, parameter :: element(5) = ["H","C","N","O","S"]
integer,   parameter :: atomic_mass(5) = [1,12,14,16,32]

type tHarmonicModel
  real(rb) :: k, x0
end type tHarmonicModel

type tStructute
  integer :: type, atom1, atom2, atom3
end type tStructute

type tConfig

  real(rb) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(rb), pointer :: Lx, Ly, Lz

  integer  :: ntypes
  real(rb), allocatable :: epsilon(:), sigma(:)

  integer  :: natoms
  integer,  pointer :: Mol(:), Type(:), ntype(:) !Added by Ana
  real(rb), pointer :: Charge(:)

  real(rb), pointer :: R(:,:), F(:,:), P(:,:)

  real(rb), pointer :: mass(:)
  real(rb), pointer :: Rx(:), Ry(:), Rz(:) ! Positions
  real(rb), pointer :: Fx(:), Fy(:), Fz(:) ! Forces
  real(rb), pointer :: Px(:), Py(:), Pz(:) ! Linear momenta
  real(rb), allocatable :: InvMass(:,:)

  logical :: velocity_input = .false.

  integer :: nBondTypes, nAngleTypes
  type(tHarmonicModel), allocatable :: BondModel(:), AngleModel(:)

  integer :: nbonds, nangles
  type(tStructute), allocatable :: Bond(:), Angle(:)

  contains

    procedure :: tConfig_Read, tConfig_Read_from_File
    generic :: Read => tConfig_Read, tConfig_Read_from_File

    procedure :: tConfig_Write, tConfig_Write_to_File
    generic :: Write => tConfig_Write, tConfig_Write_to_File

    procedure :: tConfig_Save_XYZ_to_unit, tConfig_Save_XYZ_to_file
    generic :: Save_XYZ => tConfig_Save_XYZ_to_unit, tConfig_Save_XYZ_to_file

    procedure :: bring_to_central_box => tConfig_bring_to_central_box
end type tConfig

type(tConfig) :: Config

contains

  !=================================================================================================

  subroutine tConfig_Read( me, unit )
    class(tConfig), intent(inout) :: me
    integer,        intent(in)    :: unit
    integer       :: narg, i, k
    integer       :: t !Added by Ana
    character(sl) :: arg(10)
    real(rb) :: mass
    call next_command( unit, narg, arg )
    do while (narg > 0)
      call next_command( unit, narg, arg )

      if ((narg == 3).and.(join(arg(2:3)) == "atom types")) then
        me % ntypes = str2int( arg(1) )

      else if ((narg == 3).and.(join(arg(2:3)) == "bond types")) then
        me % nBondTypes = str2int( arg(1) )

      else if ((narg == 3).and.(join(arg(2:3)) == "angle types")) then
        me % nAngleTypes = str2int( arg(1) )

      else if ((narg == 2).and.(arg(2) == "atoms")) then
        me % natoms = str2int( arg(1) )

      else if ((narg == 2).and.(arg(2) == "bonds")) then
        me % nbonds = str2int( arg(1) )

      else if ((narg == 2).and.(arg(2) == "angles")) then
        me % nangles = str2int( arg(1) )

      else if ((narg == 4).and.(join(arg(3:4)) == "xlo xhi")) then
        me % xmin = str2real(arg(1))
        me % xmax = str2real(arg(2))
        allocate( me % Lx )
        me % Lx = me % xmax - me % xmin

      else if ((narg == 4).and.(join(arg(3:4)) == "ylo yhi")) then
        me % ymin = str2real(arg(1))
        me % ymax = str2real(arg(2))
        allocate( me % Ly )
        me % Ly = me % ymax - me % ymin

      else if ((narg == 4).and.(join(arg(3:4)) == "zlo zhi")) then
        me % zmin = str2real(arg(1))
        me % zmax = str2real(arg(2))
        allocate( me % Lz )
        me % Lz = me % zmax - me % zmin

      else if ((narg == 1).and.(arg(1) == "Masses")) then
        allocate( me % mass(me % ntypes) )
        do i = 1, me % ntypes
          call next_command( unit, narg, arg )
          me % mass(i) = str2real(arg(2))
        end do

      else if ((narg == 2).and.(join(arg(1:2)) == "Pair Coeffs")) then
        allocate( me % epsilon(me % ntypes), me % sigma(me % ntypes) )
        do k = 1, me % ntypes
          call next_command( unit, narg, arg )
          i = str2int(arg(1))
          me % epsilon(i) = str2real(arg(2))
          me % sigma(i) = str2real(arg(3))
        end do

      else if ((narg == 2).and.(join(arg(1:2)) == "Bond Coeffs")) then

        allocate( me % bondModel(me % nBondTypes) )
        do k = 1, me % nBondTypes
          call next_command( unit, narg, arg )
          i = str2int(arg(1))
          me % bondModel(i) % k = str2real(arg(2))
          me % bondModel(i) % x0 = str2real(arg(3))
        end do

      else if ((narg == 2).and.(join(arg(1:2)) == "Angle Coeffs")) then

        allocate( me % angleModel(me % nAngleTypes) )
        do k = 1, me % nAngleTypes
          call next_command( unit, narg, arg )
          i = str2int(arg(1))
          me % angleModel(i) % k = str2real(arg(2))
          me % angleModel(i) % x0 = str2real(arg(3))
        end do

      else if ((narg == 1).and.(arg(1) == "Atoms")) then
        associate( N => me % natoms )
          allocate( me%Mol(N), me%Type(N), me%Charge(N) )
          allocate( me%ntype(me%ntypes), source = 0 ) !Added by Ana
          allocate( me%R(3,N), me%F(3,N), me%P(3,N) )
          me%Rx => me%R(1,:)
          me%Ry => me%R(2,:)
          me%Rz => me%R(3,:)
          me%Fx => me%F(1,:)
          me%Fy => me%F(2,:)
          me%Fz => me%F(3,:)
          me%Px => me%P(1,:)
          me%Py => me%P(2,:)
          me%Pz => me%P(3,:)
          allocate( me%InvMass(3,N) )
        end associate
        do k = 1, me % natoms
          call next_command( unit, narg, arg )
          i = str2int(arg(1))
          me % Mol(i) = str2int(arg(2))
          me % Type(i) = str2int(arg(3))
          me % Charge(i) = str2real(arg(4))
          me % Rx(i) = str2real(arg(5))
          me % Ry(i) = str2real(arg(6))
          me % Rz(i) = str2real(arg(7))
          me % Px(i) = zero
          me % Py(i) = zero
          me % Pz(i) = zero
          t = str2int(arg(3)) !Added by Ana
          me%ntype(t) = me%ntype(t) + 1 !Added by Ana
          if (narg == 10) then
            me % Rx(i) = me % Rx(i) + str2real(arg( 8)) * me % Lx
            me % Ry(i) = me % Ry(i) + str2real(arg( 9)) * me % Ly
            me % Rz(i) = me % Rz(i) + str2real(arg(10)) * me % Lz
          end if
        end do
        forall (k=1:me % natoms) me % InvMass(:,k) = one/me%Mass(me%Type(k))

      else if ((narg == 1).and.(arg(1) == "Velocities")) then

        me % velocity_input = .true.
        do k = 1, me % natoms
          call next_command( unit, narg, arg )
          i = str2int(arg(1))
          mass = me%Mass(me%Type(i))
          me % Px(i) = mass*str2real(arg(2))
          me % Py(i) = mass*str2real(arg(3))
          me % Pz(i) = mass*str2real(arg(4))
        end do

      else if ((narg == 1).and.(arg(1) == "Bonds")) then

        allocate( me%Bond(me%nbonds) )
        do k = 1, me % nbonds
          call next_command( unit, narg, arg )
          i = str2int(arg(1))
          me%Bond(i)%type  = str2int(arg(2))
          me%Bond(i)%atom1 = str2int(arg(3))
          me%Bond(i)%atom2 = str2int(arg(4))
        end do

      else if ((narg == 1).and.(arg(1) == "Angles")) then

        allocate( me%Angle(me%nangles) )
        do k = 1, me % nangles
          call next_command( unit, narg, arg )
          i = str2int(arg(1))
          me%Angle(i)%type  = str2int(arg(2))
          me%Angle(i)%atom1 = str2int(arg(3))
          me%Angle(i)%atom2 = str2int(arg(4))
          me%Angle(i)%atom3 = str2int(arg(5))
        end do

      end if
    end do
    call me % bring_to_central_box()
  end subroutine tConfig_Read

  !=================================================================================================

  subroutine tConfig_Read_from_File( me, file )
    class(tConfig), intent(inout) :: me
    character(*),   intent(in)    :: file
    integer :: inp, stat
    inp = 69
    open( unit = inp, file = file, status = "old", iostat = stat )
    if (stat /= 0) call error( "Configuration file", trim(file), "was not found." )
    call me % Read( inp )
    close(inp)
  end subroutine tConfig_Read_from_File

  !=================================================================================================

  subroutine tConfig_Write( me, unit, velocities )
    class(tConfig), intent(inout) :: me
    integer,        intent(in)    :: unit
    logical,        intent(in), optional :: velocities
    integer :: i
    real(rb) :: xi, yi, zi
    write(unit,'("LAMMPS data file",/)')
    write(unit,'(A," atom types",/)') trim(int2str(me % ntypes))
    write(unit,'(A," atoms",/)') trim(int2str(me % natoms))
    write(unit,'(A," xlo xhi")') trim(join(real2str([me%xmin,me%xmax])))
    write(unit,'(A," ylo yhi")') trim(join(real2str([me%ymin,me%ymax])))
    write(unit,'(A," zlo zhi")') trim(join(real2str([me%zmin,me%zmax])))
    write(unit,'(/,"Masses",/)')
    do i = 1, me % ntypes
      write(unit,'(A)') trim(join([int2str(i),real2str(me%mass(i))]))
    end do
    write(unit,'(/,"Pair Coeffs",/)')
    do i = 1, me % ntypes
      write(unit,'(A)') trim(join([int2str(i),real2str(me%epsilon(i)),real2str(me%sigma(i))]))
    end do
    write(unit,'(/,"Atoms",/)')
    do i = 1, me % natoms
      xi = me%R(1,i)
      yi = me%R(2,i)
      zi = me%R(3,i)
      write(unit,'(A,X,A)') trim(join(int2str([i,me%Mol(i),me%Type(i)]))),  &
                            trim(join(real2str([me%Charge(i),xi,yi,zi])))
    end do
    if (present(velocities)) then
      if (velocities) then
        write(unit,'(/,"Velocities",/)')
        do i = 1, me % natoms
          associate (mass => me%mass(me%Type(i)))
            xi = me%Px(i)/mass
            yi = me%Py(i)/mass
            zi = me%Pz(i)/mass
          end associate
          write(unit,'(A,X,A)') trim(int2str(i)), trim(join(real2str([xi,yi,zi])))
        end do
      end if
    end if
  end subroutine tConfig_Write

  !=================================================================================================

  subroutine tConfig_Write_to_File( me, file, velocities )
    class(tConfig), intent(inout) :: me
    character(*),   intent(in)    :: file
    logical,        intent(in), optional :: velocities
    integer :: out, stat
    out = 69
    open( unit = out, file = file, status = "replace", iostat = stat )
    if (stat /= 0) call error( "Cannot open file", trim(file), "for writing." )
    if (present(velocities)) then
      call me % Write( out, velocities )
    else
      call me % Write( out )
    end if
    close(out)
  end subroutine tConfig_Write_to_File

  !=================================================================================================

  subroutine tConfig_Save_XYZ_to_file( me, file, append )
    class(tConfig), intent(inout)        :: me
    character(*),   intent(in)           :: file
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
    call tConfig_Save_XYZ_to_unit( me, out )
    close(out)
  end subroutine tConfig_Save_XYZ_to_file


  !=================================================================================================

  subroutine tConfig_Save_XYZ_to_unit( me, unit )
    class(tConfig), intent(inout) :: me
    integer,        intent(in)    :: unit
    integer :: i, m
    character(sl) :: itype
    call me % bring_to_central_box()
    write(unit,*) me % natoms
    write(unit,*)
    do i = 1, me % natoms
      m = nint(me%Mass(me%Type(i)))
      if (any(atomic_mass == m)) then
        itype = element(maxloc(transfer(atomic_mass == m, atomic_mass),dim=1))
      else
        itype = int2str(me%Type(i))
      end if
      write(unit,*) trim(itype), me%R(:,i)
    end do
  end subroutine tConfig_Save_XYZ_to_unit

  !=================================================================================================

  subroutine tConfig_bring_to_central_box( me )
    class(tConfig), intent(inout) :: me
    integer  :: i, nmol, imol
    real(rb) :: imass, ri(3), L(3), Lmin(3)
    real(rb), allocatable :: rcm(:,:), molMass(:)
    nmol = maxval(me%mol)
    allocate( rcm(3,nmol), molMass(nmol) )
    rcm = 0.0_rb
    molMass = 0.0_rb
    do i = 1, me%natoms
      imol = me%mol(i)
      imass = me%mass(me%Type(i))
      molMass(imol) = molMass(imol) + imass
      rcm(:,imol) = rcm(:,imol) + imass*me%R(:,i)
    end do
    forall (imol=1:nmol) rcm(:,imol) = rcm(:,imol)/molMass(imol)
    L = [me%Lx, me%Ly, me%Lz]
    Lmin = [me%xmin, me%ymin, me%zmin]
    do imol = 1, nmol
      ri = rcm(:,imol)/L
      rcm(:,imol) = L*(ri - floor(ri))
    end do
    do i = 1, me % natoms
      imol = me%mol(i)
      me%R(:,i) = me%R(:,i) - L*anint((me%R(:,i) - rcm(:,imol))/L)
    end do
  end subroutine tConfig_bring_to_central_box

  !=================================================================================================

end module mConfig
