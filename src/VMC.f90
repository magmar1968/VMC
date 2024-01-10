!author: Lorenzo Magnoni
!date: 23 oct 2023
!mail: lorenzo.magnoni@studenti.unimi.it or lolemagno99@gmail.com
!github:



module VMC
    implicit none
    private

    integer                             :: MC_step,step
    integer                             :: NEW,OLD,Nacceptances
    real*8                              :: MAX_RADIUS
    real*8                              :: E_avg,E2_avg
    real*8,allocatable,dimension(:)     :: density_profile
    real*8,allocatable,dimension(:,:,:) :: walker
    real*8,dimension(2)                 :: TWF                !to store Trial WF values

    type,public :: VMC_results
        real*8 :: E
        real*8 :: error
        real*8 :: acceptrate
    end type
    

    type,public :: VMC_varParameters
        real*8 Rv
        real*8 alpha
        real*8 k0
        real*8 a_osc
    end type 

    public :: run_VMC
    contains

    !####################################################
    !#           Run Variational MC Simulation          #
    !####################################################    
    subroutine run_VMC(params,results)
        use VMC_print
        use VMC_parameters
        use HS_puregas
        use random,Only:init_random_seed
        implicit none

        type(VMC_varParameters),intent(in)  :: params
        type(VMC_results),intent(out)       :: results
        real*8  :: E
        real*8  :: c1
        real*8  :: difference
        integer :: COUNTER = 0
        
        allocate(walker(Natoms,DIM,2))
        allocate(density_profile(NdensProfileSteps))
        call init_parameters(params) !initialize the simulation parameters
        call gen_TWF_tables(TWFNPartitions)
        call init_random_seed()

        call gen_initial_configuration(walker(:,:,OLD))
        if (PRINT_INITIAL_CONFIGURATION) then 
            call print_inital_conf_toFile(walker(:,:,OLD))    
        end if 
        ! main cycle 
        do MC_step = - NStabSteps, NMCsteps
            do step = 1, NThermSteps!thermalization steps to avoid correlated results
                COUNTER = COUNTER + 1 !total cycles
                
                call diffuse(walker(:,:,OLD),walker(:,:,NEW),sigma)
                if (check_hcore_crosses(walker(:,:,NEW)) .eqv. .FALSE.) then
                    TWF(OLD) = trial_WF(walker(:,:,OLD)) 
                    TWF(NEW) = trial_WF(walker(:,:,NEW))

                    !metropolis question
                    difference = 2*(TWF(NEW) - TWF(OLD))
                    if (difference > 0) then !new wavefunction greater than old 
                        OLD = 3 - OLD; NEW = 3 - NEW !swap NEW <--> OLD
                        Nacceptances = Nacceptances + 1
                    else
                        call random_number(c1)
                        if( exp(difference) > c1 ) then 
                            OLD = 3 - OLD; NEW = 3 - NEW !swap NEW <--> OLD
                            Nacceptances = Nacceptances + 1
                        end if 
                    end if
                else 
                    continue
                end if 

                 
            end do 
            !update accumulators
            if(MC_step > 0) then 
                E      = Elocal(walker(:,:,OLD))
                E_avg  = E_avg  * real(MC_step-1,kind=8)/real(MC_step,8) + (E)/real(MC_step,8)
                E2_avg = E2_avg * real(MC_step-1,kind=8)/real(MC_step,8) + (E**2)/real(MC_step,8)

                call print_states(walker(:,:,OLD)) !energy acc/ density profile
            end if 
        end do 

        results%E          = E_avg                            !FIRST RETURNED VALUE
        results%error      = sqrt( (E2_avg - E_avg**2))       !SECOND RETURNED VALUE  
        results%acceptRate = (100. * Nacceptances )/ COUNTER  !THIRD RETURNED VALUE  
        
        if(PRINT_FINAL_CONFIGURATION) then 
            call print_final_conf_toFile(walker(:,:,OLD))
        end if 

        if(PRINT_DENSITY_PROFILE) then 
            call print_density_profile_toFile(density_profile)
        end if 

        deallocate(walker);deallocate(density_profile)
        call free_TWF_tables()
    end subroutine 
    !################################################################


    !####################################################
    !#                  Init                            #
    !#################################################### 
    subroutine init_parameters(params)
        use VMC_parameters
        use HS_puregas,Only:gen_initial_configuration
        use HS_puregas,Only:set_HS_parameters
        use VMC_print,Only:print_inital_conf_toFile
        use VMC_print,Only:print_TWF_tofile
        implicit none 
        type(VMC_varParameters), intent(in) :: params
        
        !initialize everything
        call set_HS_parameters(params%alpha,params%Rv,params%k0,params%a_osc)
       
       
        !setting the max distance from the center to a reasonable ammount
        MAX_RADIUS      = 3*params%a_osc 
        densProfileStep = MAX_RADIUS/NdensProfileSteps
        density_profile = 0
        sigma           = sqrt(2*h2over2m*dt)
        Nacceptances    = 0;
        E_avg           = 0.
        E2_avg          = 0.
        NEW = 2; OLD = 1  !init swappers

        if (PRINT_TRIAL_WAVEFUNCTION) then 
            call print_TWF_tofile()
        end if 
    end subroutine init_parameters 


    !####################################################
    !#               Update States                      #
    !####################################################
    subroutine print_states(R)
        use VMC_parameters
        use VMC_print
        use HS_puregas,Only: get_energies
        real*8,dimension(:,:) :: R
        real*8 :: epot,ekin,ekinfor,E

        if(PRINT_DENSITY_PROFILE) then 
            call update_density_profile(R)
        end if

        if(PRINT_ENERGY_EVOLUTION) then
            call get_energies(R,epot,ekin,ekinfor) 
            E = epot + ekin
            call print_E_evolution_toFile(MC_step,E_avg,E,epot,ekin,ekinfor)
        end if 

        if(PRINT_ATOMS_PATH) then 
            call print_atoms_path_toFile(R,MC_step)
        end if 
        
    end subroutine print_states






    !####################################################
    !#               Update Density Profile             #
    !#################################################### 
    subroutine update_density_profile(R)
        use VMC_parameters
        use HS_puregas
        use array_utility,Only:increase_size
        implicit none
        real*8,intent(in),dimension(Natoms,DIM) :: R
        real*8     :: radius 
        integer    :: i_atom,i_step
        integer    :: Nenlarging,new_size
        
        do i_atom = 1, Natoms
            radius = norm2(R(i_atom,:)) 

            if (radius > MAX_RADIUS) then 
                Nenlarging = floor( (radius - MAX_RADIUS)/densProfileStep)
                Nenlarging = Nenlarging + 2 !adding some extra space 
                new_size   = size(density_profile) + Nenlarging
                
                call increase_size(array=density_profile, new_size=new_size)
                
                NdensProfileSteps = size(density_profile)
                MAX_RADIUS        = densProfileStep*new_size
            end if 
            i_step = int(radius/densProfileStep) + 1
            density_profile(i_step) = density_profile(i_step) + 1
        end do
    end subroutine

end module
    

