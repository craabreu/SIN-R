#include "emdee.f03"

program lj_nvt

use mGlobal
use mConfig
use mRandom
use mThermostat
use iso_c_binding
use EmDee

implicit none

real(rb), parameter :: mvv2e = 2390.057364_rb         ! Da*A^2/fs^2 to kcal/mol
real(rb), parameter :: kB = 8.31451E-7_rb             ! Boltzmann constant in Da*A^2/(fs^2*K)
real(rb), parameter :: Pconv = 1.6388244954E+8_rb     ! Da/(A*fs^2) to atm
real(rb), parameter :: kCoul = 0.138935456_rb         ! Coulomb constant in Da*A^3/(fs*e)^2

integer, parameter :: velocity_verlet = 0, &
                      nose_hoover_chain = 1, &
                      stochastic_velocity_rescaling = 2, &
                      stochastic_isokinetic_nh_respa = 3

! Simulation specifications:
integer  :: seed, Nconf, thermo, Nequil, Nprod
real(rb) :: T, RcIn, Rc, SwitchWidth, dt, skin, kspacePrecision
logical  :: shift
character(sl) :: Base

! RESPA-related variables:
integer  :: RespaN(3)

! System properties:
integer  :: dof
real(rb) :: Volume

! SIN(R) variables:
real(rb), pointer :: F(:,:,:)

! Thermostat variables:
integer  :: method, M, ndamp, nloops
real(rb) :: tdamp
class(tThermostat), pointer :: thermostat

! Radial distribution function variables:
logical :: computeRDF
integer :: nevery, bins, npairs, counter
real(rb), allocatable :: gr(:,:), rdf(:,:)
integer , allocatable :: itype(:), jtype(:)

! Other variables:
integer  :: step
real(rb) :: dt_2, dt_4, KE_sp, kT
character(256) :: filename, configFile

integer :: threads
type(tEmDee)  :: md
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
call Configure_System( md )
call Config % Save_XYZ( trim(Base)//".xyz" )

call writeln( titles )
call writeln( properties(0) )
do step = 1, NEquil
  md%Options%Compute = mod(step,thermo) == 0
  call execute_step
  if (md%Options%Compute) call writeln( properties(step) )
end do
call writeln( "Loop time of", real2str(md%Time%Total), "s." )
call Report( md )

call writeln( )
call writeln( "Memory usage" )
call writeln( titles )
call writeln( properties(NEquil) )

if (computeRDF) then
  allocate(gr(bins,npairs), source = 0.0_rb)
  allocate(rdf(bins,npairs), source = 0.0_rb)
  counter = 1
  call EmDee_rdf(md, bins, npairs, itype, jtype, gr)
end if

do step = NEquil+1, NEquil+NProd
  md%Options%Compute = mod(step,thermo) == 0
  call execute_step
  if (md%Options%Compute) call writeln( properties(step) )
  if (mod(step,Nconf)==0) call Config % Save_XYZ( trim(Base)//".xyz", append = .true. )
  if (computeRDF .and. (mod(step, nevery) == 0)) then
    call EmDee_rdf(md, bins, npairs, itype, jtype, rdf)
    gr = gr + rdf
    counter = counter + 1
  end if
end do
if (computeRDF) call rdf_save( trim(Base)//".rdf" )
call writeln( "Loop time of", real2str(md%Time%Total), "s." )
call Report( md )
call Config % Write( trim(Base)//"_out.lmp", velocities = .true. )
call stop_log

contains
  !-------------------------------------------------------------------------------------------------
  subroutine execute_step
    call RESPA_Step( dt, 3 )
  end subroutine execute_step
  !-------------------------------------------------------------------------------------------------
  subroutine Configure_System( md )
    type(tEmDee), intent(inout) :: md

    integer, parameter :: nlayers = 3
    integer, parameter :: ones(nlayers) = 1
    integer :: i
    real(rb) :: theta0
    type(c_ptr) :: LJ, models(nlayers), address

    ! Initialize a system with 3 model layers:
    md = EmDee_system( threads, nlayers, Rc, skin, Config%natoms, &
                       c_loc(Config%Type(1)), c_loc(Config%mass(1)), c_null_ptr )

    ! Define layer-based parameters:
    ! Layer 1 - Bond and angle interactions only
    ! Layer 2 - Smoothed LJ and Coulomb potentials with small cutoff
    ! Layer 3 - Bond + Angle + Smoothed LJ + Ewald with large cutoff
    call EmDee_layer_based_parameters( md, RcIn, Apply = [0, 1, 0], Bonded = [1, 0, 1] )

    ! Add bond models and bonds:
    do i = 1, Config%nBondTypes
      bond_model(i) = EmDee_bond_harmonic(two*Config%BondModel(i)%k/mvv2e, Config%BondModel(i)%x0)
    end do
    do i = 1, Config%nbonds
      associate (bond => Config%Bond(i))
        call EmDee_add_bond( md, bond%atom1, bond%atom2, bond_model(bond%type) )
      end associate
    end do

    ! Add angle models and angles:
    do i = 1, Config%nAngleTypes
      theta0 = Config%AngleModel(i)%x0*pi/180.0_rb
      angle_model(i) = EmDee_angle_harmonic(two*Config%AngleModel(i)%k/mvv2e, theta0)
    end do
    do i = 1, Config%nangles
      associate (angle => Config%Angle(i))
        call EmDee_add_angle( md, angle%atom1, angle%atom2, angle%atom3, angle_model(angle%type) )
      end associate
    end do

    ! Add van der Waals interaction pairs:
    do i = 1, Config%ntypes
      if (abs(Config%epsilon(i)) < epsilon(1.0_rb)) then
        call EmDee_set_pair_model( md, i, i, EmDee_pair_none(), kCoul )
      else
        LJ = EmDee_pair_lj_cut( Config%epsilon(i)/mvv2e, Config%sigma(i) )
        models = [ EmDee_pair_none(),                                &
                   merge( EmDee_shifted_smoothed( LJ, SwitchWidth ), &
                          EmDee_smoothed( LJ, SwitchWidth ),         &
                          shift ),                                   &
                   EmDee_smoothed( LJ, SwitchWidth )                 ]
        call EmDee_set_pair_multimodel( md, i, i, models, kCoul*ones )
      end if
    end do

    models = [ EmDee_coul_none(),                                               &
               merge( EmDee_shifted_smoothed( EmDee_coul_cut(), SwitchWidth ),  &
                      EmDee_smoothed( EmDee_coul_cut(), SwitchWidth ),          &
                      shift ),                                                  &
               EmDee_coul_long()                                                ]
    call EmDee_set_coul_multimodel( md, models )

    call EmDee_set_kspace_model( md, EmDee_kspace_ewald( kspacePrecision ) )

    call EmDee_upload( md, "charges"//c_null_char, c_loc(Config%Charge(1)) )
    call EmDee_upload( md, "box"//c_null_char, c_loc(Config%Lx) )
    call EmDee_upload( md, "coordinates"//c_null_char, c_loc(Config%R(1,1)) )
    if (Config%velocity_input) then
      call EmDee_upload( md, "momenta"//c_null_char, c_loc(Config%P(1,1)) )
    else
      call EmDee_random_momenta( md, kT, .true._1, seed )
    end if

    address = EmDee_memory_address( md, "coordinates"//c_null_char )
    nullify(Config%R)
    call c_f_pointer( address, Config%R, [3, Config%natoms] )

    address = EmDee_memory_address( md, "momenta"//c_null_char )
    nullify(Config%P)
    call c_f_pointer( address, Config%P, [3, Config%natoms] )

    address = EmDee_memory_address( md, "layerForces"//c_null_char )
    call c_f_pointer( address, F, [3, Config%natoms, nlayers] )

    call EmDee_switch_model_layer( md, nlayers )
    call EmDee_compute_forces( md )
    call discount_forces()

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
  character(sl) function properties(step)
    integer, intent(in) :: step
    real(rb) :: Temp, kineticEnergy, TotalEnergy
    kineticEnergy = half*sum(Config%invMass*Config%P**2)
    Temp = (kineticEnergy/KE_sp)*T
    TotalEnergy = md%Energy%Potential + kineticEnergy
    properties = trim(adjustl(int2str(step))) // " " // &
                 join(real2str([ Temp, &
                                 Pconv*((dof-3)*kB*Temp/3 + md%Virial%Total/3.0_rb)/Volume, &
                                 mvv2e*[kineticEnergy, &
                                        md%Energy%Dispersion, &
                                        md%Energy%Coulomb, &
                                        md%Energy%Bond, &
                                        md%Energy%Angle, &
                                        md%Energy%Potential, &
                                        TotalEnergy, &
                                        md%Virial%Total, &
                                        TotalEnergy + thermostat%energy()]]))
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
    read(inp,*); read(inp,*) RcIn, Rc, SwitchWidth
    read(inp,*); read(inp,*) shift
    read(inp,*); read(inp,*) kspacePrecision
    read(inp,*); read(inp,*) seed
    read(inp,*); read(inp,*) dt
    read(inp,*); read(inp,*) RespaN
    read(inp,*); read(inp,*) skin
    read(inp,*); read(inp,*) Nconf
    read(inp,*); read(inp,*) thermo
    read(inp,*); read(inp,*) Nequil, Nprod
    read(inp,*); read(inp,*) method
    read(inp,*); read(inp,*) ndamp, M, nloops
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
    call writeln( "Internal cutoff distance:", real2str(RcIn), "A" )
    call writeln( "External cutoff distance:", real2str(Rc), "A" )
    call writeln( "Switching region width:", real2str(SwitchWidth), "A" )
    call writeln( "Seed for random numbers:", int2str(seed) )
    call writeln( "Time step size:", real2str(dt), "fs" )
    call writeln( "Numbers of RESPA loops: ", join(int2str(RespaN)) )
    call writeln( "Skin for neighbor lists:", real2str(skin), "A" )
    call writeln( "Interval for saving configurations:", int2str(Nconf) )
    call writeln( "Interval for printing properties:", int2str(thermo) )
    call writeln( "Number of equilibration steps:", int2str(Nequil) )
    call writeln( "Number of production steps:", int2str(Nprod) )
    call writeln( "Thermostat method:", int2str(method) )
    call writeln( "Thermostat parameters:", int2str(ndamp), int2str(M), int2str(nloops) )
    call writeln()
  end subroutine Read_Specifications
  !-------------------------------------------------------------------------------------------------
  subroutine Setup_Simulation
    real(rb) :: L(3)
    L = [Config % Lx, Config % Ly, Config % Lz]
    if (Rc+skin >= half*minval(L)) call error( "minimum image convention failed!" )
    dt_2 = half*dt
    dt_4 = 0.25_rb*dt
    dof = 3*Config%natoms - 3
    kT = kB*T
    KE_sp = half*dof*kT
    Volume = product(L)
    tdamp = ndamp*dt

    call random % setup( seed )

    select case (method)
      case (velocity_verlet)
        allocate( nve :: thermostat )
      case (nose_hoover_chain)
        allocate( nhc :: thermostat )
      case (stochastic_velocity_rescaling)
        allocate( csvr :: thermostat )
      case (stochastic_isokinetic_nh_respa)
        allocate( sinr :: thermostat )
      case default
        call error("Specified thermostat method not implemented")
    end select

    select type(thermo => thermostat)
      class is (nhc)
        call thermo % setup( M, kT, tdamp, dof, nloops )
      class is (csvr)
        call thermo % setup( kT, tdamp, dof, seed )
      class is (sinr)
        call thermo % setup( kT, tdamp, 3, Config%natoms, Config%mass(Config%Type), seed )
    end select

  end subroutine Setup_Simulation
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
  !-------------------------------------------------------------------------------------------------
  subroutine discount_forces()
    F(:,:,3) = F(:,:,3) - (F(:,:,1) + F(:,:,2))
  end subroutine discount_forces
  !-------------------------------------------------------------------------------------------------
  subroutine RESPA_Step( timestep, layer )
    real(rb), intent(in) :: timestep
    integer,  intent(in) :: layer

    integer  :: step
    real(rb) :: dt, dt_2

    dt = timestep/RespaN(layer)
    dt_2 = half*dt
    do step = 1, RespaN(layer)
      call EmDee_switch_model_layer( md, layer )
      Config%P = Config%P + F(:,:,layer)*dt_2
      if (layer == 1) then
        call EmDee_displace( md, one, zero, dt_2 )
        call thermostat % integrate( dt, sum(Config%invMass*Config%P**2) )
        Config%P = Config%P*exp(-thermostat%damping*dt)
        call EmDee_displace( md, one, zero, dt_2 )
      else
        call RESPA_Step( dt, layer-1 )
      end if
      call EmDee_switch_model_layer( md, layer )
      call EmDee_compute_forces( md )
      if (layer == 3) call discount_forces()
      Config%P = Config%P + F(:,:,layer)*dt_2
    end do

  end subroutine
!===================================================================================================
end program lj_nvt
