!#################################################
!#################################################          

! the module defines the global variables for a 
! variational MC for hard spheres in an Harmonic
! potential in 2D.

! The additional subroutine, 'input', read the 
! input gile and allocate the variables accordingly

!#################################################
!#################################################          

module VMC_parameters
    implicit none
    
    !define pi
    real*8, parameter :: PI=dacos(-1.d0)

    !space dimensions 
    integer           :: DIM = 2

    !MCparams
    !{
    integer           :: NMCsteps     !number of MC steps with default value 
    integer           :: NThermSteps  !number of thermalization steps       
    integer           :: NStabSteps   !number of stability steps    
    real*8            :: dt   
    real*8            :: D = 1        !unity of energy 
    real*8            :: sigma        !diffusion variance    
    !}

    !system parameters
    !{
    integer           :: Nwalkers = 1 !number of walkers(in VMC = 1)
    integer           :: Natoms       !number of atoms/particles
    !}

    !other parameters
    !{
    integer           :: NdensProfileSteps
    real*8            :: densProfileStep
    !}


    !dimensional parameter
    !{
    real*8, parameter :: h2over2m = 1 !set to one
    !}

    !verbosity parameters
    !{
    logical           :: PRINT_INITIAL_CONFIGURATION = .FALSE.
    logical           :: PRINT_FINAL_CONFIGURATION   = .FALSE.
    logical           :: PRINT_DENSITY_PROFILE       = .FALSE.
    logical           :: PRINT_ENERGY_EVOLUTION      = .FALSE.
    logical           :: PRINT_INITIAL_PARAMETERS    = .FALSE.
    logical           :: PRINT_TRIAL_WAVEFUNCTION    = .FALSE.

    integer,parameter  :: MAX_FILENAME_LENGHT = 50
    character(MAX_FILENAME_LENGHT) :: outfile_path = "./data/"
    character(MAX_FILENAME_LENGHT) :: init_conf_filename        = "init_conf.dat"
    character(MAX_FILENAME_LENGHT) :: fin_conf_filename         = "fin_conf.dat"
    character(MAX_FILENAME_LENGHT) :: dens_profile_filename     = "density_profile.dat"
    character(MAX_FILENAME_LENGHT) :: energy_evolution_filename = "energy_evolution.dat"
    character(MAX_FILENAME_LENGHT) :: twf_filename              = "twf.dat"
    !}
    
end module 
!######################################################

!######################################################
!               Input Subroutine
!######################################################
subroutine input(filename)
    use VMC_parameters
    use iofile,only: io_open,read_data,is_present

    implicit none

    
    character(*), intent(in) :: filename
    call io_open(input_filename = filename)

    !read MCparams
    call read_data("NMCsteps",NMCsteps)
    call read_data("NTsteps",NThermSteps)
    call read_data("NStabSteps",NStabSteps)
    call read_data("dt",dt)

    !read system parameters
    !call read_data("Nwalkers",Nwalkers) !VMC supposed to be one 
    call read_data("Natoms",Natoms)
    call read_data("NdensProfileSteps",NdensProfileSteps)

    !read verbosity flag
    call read_data("INITIAL_CONFIGURATION",PRINT_INITIAL_CONFIGURATION)
    call read_data("FINAL_CONFIGURATION"  ,PRINT_FINAL_CONFIGURATION  )
    call read_data("DENSITY_PROFILE"      ,PRINT_DENSITY_PROFILE      )
    call read_data("ENERGY_EVOLUTION"     ,PRINT_ENERGY_EVOLUTION     )
    call read_data("INITIAL_PARAMETERS"   ,PRINT_INITIAL_PARAMETERS   )
    call read_data("TRIAL_WF"             ,PRINT_TRIAL_WAVEFUNCTION   )

    if (is_present("INIT_CONF_FILENAME")) then 
        call read_data("INIT_CONF_FILENAME", init_conf_filename)
    end if 
    if (is_present("FIN_CONF_FILENAME")) then 
        call read_data("FIN_CONF_FILENAME", fin_conf_filename)
    end if 
    if (is_present("DENS_PROFILE_FILENAME")) then 
        call read_data("DENS_PROFILE_FILENAME", dens_profile_filename)
    end if 
    if (is_present("ENERGY_EVOLUTION_FILENAME")) then 
        call read_data("ENERGY_EVOLUTION_FILENAME", energy_evolution_filename)
    end if 
    if (is_present("TWF_FILENAME")) then 
        call read_data("TWF_FILENAME", twf_filename)
    end if

    if (PRINT_INITIAL_PARAMETERS) then 
        call print_parameters()
    endif    
end subroutine

subroutine print_parameters()
    use VMC_parameters
    implicit none
    
    print *, "#####################################################"
    print *, " Montecarlo Simulation Using the following parameters"
    print *, "   - Number MC Steps       : ", NMCsteps
    print *, "   - Number Therm Steps    : ", NThermSteps
    print *, "   - Number Stability steps: ", NStabSteps
    print *, "   - dt                    : ", dt
    print *, "   - Atoms Number          : ", Natoms
    print *, "   - Number step dens prof : ", NdensProfileSteps
    print *, " FLAGS:                                              "
    print *, "   - INITIAL CONFIGURATION : ", PRINT_INITIAL_CONFIGURATION
    print *, "   - FINAL CONFIGUTATION   : ", PRINT_FINAL_CONFIGURATION
    print *, "   - DENSITY PROFILE       : ", PRINT_DENSITY_PROFILE
    print *, "   - ENERGY EVOLUTION      : ", PRINT_ENERGY_EVOLUTION
    print *, "   - INITIAL PARAMETERS    : ", PRINT_INITIAL_PARAMETERS 
    print *, "   - TRIAL WAVEFUNCTION    : ", PRINT_TRIAL_WAVEFUNCTION
    
    print *, "######################################################"

end subroutine