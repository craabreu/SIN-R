#include "emdee.f03"

program lj_nvt

use mGlobal
use mConfig
use mRandom
use mThermostat
use iso_c_binding
use EmDee

implicit none

real(rb), parameter :: mvv2e = 2390.057364_rb         ! Da*A²/fs² to kcal/mol
real(rb), parameter :: kB = 8.31451E-7_rb             ! Boltzmann constant in Da*A²/(fs²*K)
real(rb), parameter :: Pconv = 1.6388244954E+8_rb     ! Da/(A*fs²) to atm
real(rb), parameter :: kCoul = 0.13893545755135628_rb ! Coulomb constant in Da*A³/(fs*e)²

integer, parameter :: velocity_verlet = 0, &
                      nose_hoover_chain = 1, &
                      stochastic_velocity_rescaling = 2

! Simulation specifications:
character(sl) :: Base
integer       :: N, seed, Nconf, thermo, Nequil, Nprod, rotationMode
real(rb)      :: T, Rc, Rm, dt, skin, kspacePrecision

! System properties:
integer :: dof
real(rb) :: Volume

! Thermostat variables:
integer :: method, M, ndamp, nloops
real(rb) :: tdamp, Hthermo
class(nhc), pointer :: thermostat

! Radial distribution function variables:
logical :: computeRDF
integer :: nevery, bins, npairs, counter
real(rb), allocatable :: gr(:,:), rdf(:,:)
integer , allocatable :: itype(:), jtype(:)

! Other variables:
integer :: step
real(rb) :: dt_2, dt_4, KE_sp, kT
character(256) :: filename, configFile

integer :: threads
type(tEmDee) :: md
type(mt19937) :: random
type(c_ptr), allocatable :: model(:)
type(c_ptr), allocatable :: bond_model(:)
type(c_ptr), allocatable :: angle_model(:)

character(*), parameter :: titles = "Step Temp Press KinEng DispEng CoulEng "// &
                                    "BondEng AngleEng PotEng TotEng Virial H_nhc"

! Executable code:
call writeln( "md/lj/coul ("//__DATE__//")" )

call Get_Command_Line_Args( threads, filename )
call Read_Specifications( filename )

call Config % Read( configFile )
call Setup_Simulation

allocate( model(Config%ntypes), bond_model(Config%nBondTypes), angle_model(Config%nAngleTypes) )
call Configure_System( md, Rm, Rc )
call Config % Save_XYZ( trim(Base)//".xyz" )

call writeln( titles )
step = 0
call writeln( properties() )
do step = 1, NEquil
  call execute_step
  if (mod(step,thermo) == 0) call writeln( properties() )
end do
call writeln( "Loop time of", real2str(md%Time%Total), "s." )
call Report( md )

call writeln( )
call writeln( "Memory usage" )
call writeln( titles )
step = NEquil
call writeln( properties() )

if (computeRDF) then
  allocate(gr(bins,npairs), source = 0.0_rb)
  allocate(rdf(bins,npairs), source = 0.0_rb)
  counter = 1
  call EmDee_rdf(md, bins, npairs, itype, jtype, gr)
end if

do step = NEquil+1, NEquil+NProd
  call execute_step
  if (mod(step,Nconf)==0) then
    call EmDee_download( md, "coordinates"//c_null_char, c_loc(Config%R(1,1)) )
    call Config % Save_XYZ( trim(Base)//".xyz", append = .true. )
  end if
  if (mod(step,thermo) == 0) call writeln( properties() )
  if (computeRDF .and. (mod(step, nevery) == 0)) then
    call EmDee_rdf(md, bins, npairs, itype, jtype, rdf)
    gr = gr + rdf
    counter = counter + 1
  end if
end do
if (computeRDF) call rdf_save( trim(Base)//".rdf" )
call writeln( "Loop time of", real2str(md%Time%Total), "s." )
call Report( md )
call EmDee_download( md, "coordinates"//c_null_char, c_loc(Config%R(1,1)) )
call Config % Write( trim(Base)//"_out.lmp", velocities = .true. )
call stop_log

contains
  !-------------------------------------------------------------------------------------------------
  subroutine execute_step
    select case (method)
      case (velocity_verlet); call Verlet_Step
      case (nose_hoover_chain); call NHC_Step
      case (stochastic_velocity_rescaling); call Bussi_Step
    end select
  end subroutine execute_step
  !-------------------------------------------------------------------------------------------------
  subroutine Configure_System( md, Rm, Rc )
    type(tEmDee), intent(inout) :: md
    real(rb),     intent(in)    :: Rm, Rc

    integer :: i
    real(rb) :: theta0

    md = EmDee_system( threads, 1, Rc, skin, Config%natoms, &
                       c_loc(Config%Type(1)), c_loc(Config%mass(1)), c_null_ptr )

    md%Options%rotationMode = rotationMode

    do i = 1, Config%nBondTypes
      bond_model(i) = EmDee_bond_harmonic(Config%BondModel(i)%k/mvv2e, Config%BondModel(i)%x0)
    end do
    do i = 1, Config%nbonds
      associate (bond => Config%Bond(i))
        call EmDee_add_bond( md, bond%atom1, bond%atom2, bond_model(bond%type) )
      end associate
    end do

    do i = 1, Config%nAngleTypes
      theta0 = Config%AngleModel(i)%x0*pi/180.0_rb
      angle_model(i) = EmDee_angle_harmonic(Config%AngleModel(i)%k/mvv2e, theta0)
    end do
    do i = 1, Config%nangles
      associate (angle => Config%Angle(i))
        call EmDee_add_angle( md, angle%atom1, angle%atom2, angle%atom3, angle_model(angle%type) )
      end associate
    end do

    do i = 1, Config%ntypes
      if (abs(Config%epsilon(i)) < epsilon(1.0_rb)) then
        model(i) = EmDee_pair_none()
      else
        model(i) = EmDee_smoothed( &
                     EmDee_pair_lj_cut( Config%epsilon(i)/mvv2e, Config%sigma(i) ), Rc-Rm )
      end if
      call EmDee_set_pair_model( md, i, i, model(i), kCoul )
    end do

    call EmDee_set_coul_model( md, EmDee_coul_long() )
    call EmDee_set_kspace_model( md, EmDee_kspace_ewald( kspacePrecision ) )
    call EmDee_upload( md, "charges"//c_null_char, c_loc(Config%Charge(1)) )

    call EmDee_upload( md, "box"//c_null_char, c_loc(Config%Lx) )
    call EmDee_upload( md, "coordinates"//c_null_char, c_loc(Config%R(1,1)) )
    call EmDee_random_momenta( md, kT, .true._1, seed )

  end subroutine Configure_System
  !-------------------------------------------------------------------------------------------------
  subroutine Report( md )
    type(tEmDee), intent(in) :: md
    real(rb) :: other
    call writeln( repeat("-",40) )
    call writeln( "Pair time      =", real2str(md%Time%Pair), "s." )
    call writeln( "Neighbor time  =", real2str(md%Time%Neighbor), "s." )
    other = md%Time%Total - (md%Time%Pair + md%Time%Neighbor)
    call writeln( "Neighbor list builds =", int2str(md%builds) )
    call writeln( repeat("-",40) )
  end subroutine Report
  !-------------------------------------------------------------------------------------------------
  character(sl) function properties()
    real(rb) :: Temp, H
    Temp = (md%Energy%Kinetic/KE_sp)*T
    H = md%Energy%Potential + md%Energy%Kinetic
    if (method == nose_hoover_chain) Hthermo = thermostat%energy()
    properties = trim(adjustl(int2str(step))) // " " // &
                 join(real2str([ Temp, &
                                 Pconv*((dof-3)*kB*Temp/3 + md%Virial/3.0_rb)/Volume, &
                                 mvv2e*[md%Energy%Kinetic, &
                                        md%Energy%Dispersion, &
                                        md%Energy%Coulomb, &
                                        md%Energy%Bond, &
                                        md%Energy%Angle, &
                                        md%Energy%Potential, &
                                        H, &
                                        md%Virial, &
                                        H + Hthermo]]))
  end function properties
  !-------------------------------------------------------------------------------------------------
  subroutine Get_Command_Line_Args( threads, filename )
    integer,        intent(out) :: threads
    character(256), intent(out) :: filename
    integer :: argcount
    character(256) :: line
    argcount = command_argument_count()
    if (argcount == 1) then
      threads = 1
      call get_command_argument( 1, line )
    else if (argcount == 2) then
      call get_command_argument( 1, line )
      read(line,*) threads
      call get_command_argument( 2, line )
    else
      write(0,'("Usage: md_lj/md_lj_coul [number-of-threads] input-file")')
      stop
    end if
    filename = line
  end subroutine Get_Command_Line_Args
  !-------------------------------------------------------------------------------------------------
  subroutine Read_Specifications( file )
    character(*), intent(in) :: file
    integer :: inp, i
    open( newunit=inp, file = file, status = "old" )
    read(inp,*); read(inp,*) Base
    read(inp,*); read(inp,*) configFile
    read(inp,*); read(inp,*) T
    read(inp,*); read(inp,*) Rc, Rm
    read(inp,*); read(inp,*) kspacePrecision
    read(inp,*); read(inp,*) seed
    read(inp,*); read(inp,*) dt
    read(inp,*); read(inp,*) skin
    read(inp,*); read(inp,*) Nconf
    read(inp,*); read(inp,*) thermo
    read(inp,*); read(inp,*) Nequil, Nprod
    read(inp,*); read(inp,*) method
    read(inp,*); read(inp,*) ndamp, M, nloops
    read(inp,*); read(inp,*) rotationMode
    read(inp,*); read(inp,*) npairs
    computeRDF = npairs /= 0
    allocate( itype(npairs), jtype(npairs) )
    read(inp,*); read(inp,*) nevery, bins, (itype(i), jtype(i), i=1, npairs)
    close(inp)
    call init_log( trim(Base)//".log" )
    call writeln()
    call writeln( "Base for file names:", Base )
    call writeln( "Name of configuration file:", configFile )
    call writeln( "Temperature:", real2str(T), "K" )
    call writeln( "Cutoff distance:", real2str(Rc), "Å" )
    call writeln( "Seed for random numbers:", int2str(seed) )
    call writeln( "Time step size:", real2str(dt), "fs" )
    call writeln( "Skin size for neighbor lists:", real2str(skin), "Å" )
    call writeln( "Interval for saving configurations:", int2str(Nconf) )
    call writeln( "Interval for printing properties:", int2str(thermo) )
    call writeln( "Number of equilibration steps:", int2str(Nequil) )
    call writeln( "Number of production steps:", int2str(Nprod) )
    call writeln( "Thermostat method:", int2str(method) )
    call writeln( "Thermostat parameters:", int2str(ndamp), int2str(M), int2str(nloops) )
    if (rotationMode == 0) then
      call writeln( "Rotation mode: exact solution" )
    else
      call writeln( "Rotation mode: Miller with", int2str(rotationMode), "respa steps" )
    end if
    call writeln()
  end subroutine Read_Specifications
  !-------------------------------------------------------------------------------------------------
  subroutine Setup_Simulation
    real(rb) :: Lx, Ly, Lz
    N = Config % natoms
    Lx = Config % Lx
    Ly = Config % Ly
    Lz = Config % Lz
    if (Rc+skin >= half*min(Lx,Ly,Lz)) call error( "minimum image convention failed!" )
    dt_2 = half*dt
    dt_4 = 0.25_rb*dt
    dof = 3*N - 3
    kT = kB*T
    KE_sp = half*dof*kT
    Volume = Lx*Ly*Lz
    call random % setup( seed )
    if (all(method /= [velocity_verlet, nose_hoover_chain, stochastic_velocity_rescaling])) then
        call error("Specified method not implemented")
    end if
    if (method == nose_hoover_chain) then
        allocate( nhc_pscaling :: thermostat )
        call thermostat % setup( M, kT, ndamp*dt, dof, nloops )
    end if
    tdamp = ndamp*dt
    Hthermo = 0.0_rb
  end subroutine Setup_Simulation
  !-------------------------------------------------------------------------------------------------
  subroutine Verlet_Step
      call EmDee_boost( md, one, zero, dt_2 )
      call EmDee_displace( md, one, zero, dt )
      call EmDee_boost( md, one, zero, dt_2 )
  end subroutine Verlet_Step
  !-------------------------------------------------------------------------------------------------
  function BussiScale(KE, KE_sp, dof, tau, dt) result( alphaSq )
    real(rb), intent(in) :: KE, KE_sp, tau, dt
    integer,  intent(in) :: dof
    real(rb)             :: alphaSq
    real(rb) :: R1, x, sumRs, A, B, C
    R1 = random%normal()
    if (mod(dof, 2) == 1) then
      x = (dof - 1)/2
      sumRs = 2.0*random%gamma(x)
    else
      x = (dof - 2)/2
      sumRs = 2.0*random%gamma(x) + random%normal()**2
    end if
    A = exp(-dt/tau)
    B = 1.0 - A
    C = KE_sp/(dof*KE)
    alphaSq = A + C*B*(R1**2 + sumRs) + 2.0*sqrt(C*B*A)*R1
  end function BussiScale
  !-------------------------------------------------------------------------------------------------
  subroutine Bussi_Step
    real(rb) :: factor
    factor = BussiScale( md%Energy%Kinetic, KE_sp, dof, tdamp, dt_2 )
    Hthermo = Hthermo + (1.0 - factor)*md%Energy%Kinetic
    call EmDee_boost( md, zero, -log(factor)/dt, dt_2 )
    call Verlet_Step
    factor = BussiScale( md%Energy%Kinetic, KE_sp, dof, tdamp, dt_2 )
    Hthermo = Hthermo + (1.0 - factor)*md%Energy%Kinetic
    call EmDee_boost( md, zero, -log(factor)/dt, dt_2 )
  end subroutine Bussi_Step
  !-------------------------------------------------------------------------------------------------
  subroutine NHC_Step
      call thermostat % integrate( dt_2, two*md%Energy%Kinetic )
      call EmDee_boost( md, zero, thermostat%damping, dt_2 )
      call Verlet_Step
      call thermostat % integrate( dt_2, two*md%Energy%Kinetic )
      call EmDee_boost( md, zero, thermostat%damping, dt_2 )
  end subroutine NHC_Step
!---------------------------------------------------------------------------------------------------
  subroutine rdf_save( file )
    character(*), intent(in) :: file
    integer :: i, unit
    character(3) :: Ci, Cj
    open(newunit=unit, file=file, status="replace")
    write(unit,'(A)',advance='no') "r"
    do i = 1, npairs
      write(Ci,'(I3)') itype(i)
      write(Cj,'(I3)') jtype(i)
      write(unit,'(A)',advance='no') " g("//trim(adjustl(Ci))//","//trim(adjustl(Cj))//")"
    end do
    write(unit,*)
    do i = 1, bins
      write(unit,*) (i - 0.5)*Rc/bins, gr(i,:)/real(counter,rb)
    end do
  end subroutine rdf_save
!===================================================================================================
end program lj_nvt
