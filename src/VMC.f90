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
    real*8                              :: E_acc,E2_acc
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
        use VMC_print,only:print_inital_conf_toFile,print_final_conf_toFile
        use VMC_parameters
        use HS_puregas
        implicit none

        type(VMC_varParameters),intent(in)  :: params
        type(VMC_results),intent(out)    :: results
        real*8  :: E,E2
        integer :: COUNTER = 0
        
        allocate(walker(Natoms,DIM,2))
        allocate(density_profile(NdensProfileSteps))
        call init(params)
        
        ! main cycle 
        do MC_step = - NStabSteps, NMCsteps
            do step = 1, NThermSteps!thermalization steps to avoid correlated results
                call diffuse(walker(:,:,OLD),walker(:,:,NEW),sigma)
                TWF(OLD) = trial_WF(walker(:,:,OLD)) !compute the initial value of the TWF
                TWF(NEW) = trial_WF(walker(:,:,NEW))

                !metropolis question
                if( (TWF(NEW)/TWF(OLD))**2 > rand()  ) then 
                    OLD = 3 - OLD; NEW = 3 - NEW !swap NEW <--> OLD
                    Nacceptances = Nacceptances + 1
                end if 
                COUNTER = COUNTER + 1
            end do 
            !update accumulators
            if(MC_step > 0) then 
                call update_energy_accumulator(walker(:,:,OLD))
                call update_density_profile(walker(:,:,OLD))
            end if 
        end do 
        E  = E_acc/NMCsteps                 
        E2 = E2_acc/NMCsteps

        results%E          = E                           !FIRST RETURNED VALUE
        results%error      = sqrt( (E2 - E**2)/NMCsteps) !SECOND RETURNED VALUE  
        results%acceptRate = (100. * Nacceptances )/ COUNTER  !THIRD RETURNED VALUE  
          
        
        if(PRINT_FINAL_CONFIGURATION) then 
            call print_final_conf_toFile(walker(:,:,OLD))
        end if 

        deallocate(walker);deallocate(density_profile)
    end subroutine 
    !################################################################


    !####################################################
    !#                  Init                            #
    !#################################################### 
    subroutine init(params)
        use VMC_parameters
        use HS_puregas,Only:gen_initial_configuration
        use HS_puregas,Only:set_HS_parameters
        use VMC_print,Only:print_inital_conf_toFile
        implicit none 
        type(VMC_varParameters), intent(in) :: params
        
        !initialize everything
        call set_HS_parameters(params%alpha,params%Rv,params%k0,params%a_osc)
       
       
        MAX_RADIUS      = 4*params%a_osc !setting the max distance from the center
        densProfileStep = MAX_RADIUS/NdensProfileSteps
        sigma           = sqrt(2*h2over2m*dt)
        Nacceptances    = 0;
        E_acc           = 0.
        E2_acc          = 0.
        NEW=2;OLD=1            !init swappers

        call gen_initial_configuration(walker(:,:,OLD))
        if (PRINT_INITIAL_CONFIGURATION) then 
            call print_inital_conf_toFile(walker(:,:,OLD))    
        end if 
    end subroutine init 



    !####################################################
    !#           Update Energy Accumulators             #
    !####################################################
    subroutine update_energy_accumulator(R)
        use VMC_parameters
        use VMC_print
        use HS_puregas,Only:get_energies,Elocal
        implicit none

        real*8 :: ekin, ekinfor,epot,E
        real*8,intent(in),dimension(Natoms,DIM) :: R

        E      = Elocal(R)
        E_acc  = E_acc + E
        E2_acc = E2_acc + E**2
        
        if(PRINT_ENERGY_EVOLUTION) then
            call get_energies(R,epot,ekin,ekinfor) 
            call print_E_evolution_toFile(MC_step,E,epot,ekin,ekinfor)
        end if 

    end subroutine


    !####################################################
    !#               Update Density Profile             #
    !#################################################### 
    subroutine update_density_profile(R)
        use VMC_parameters
        use HS_puregas
        implicit none
        real*8,intent(in),dimension(Natoms,DIM) :: R
        real*8 :: delta_r,radius
        integer :: i_atom,i_step
        
        do i_atom = 1, Natoms
            radius = norm2(R(i_atom,:)) 
            i_step = int(radius/delta_r)
            !possible source of segsev because I'm too lazy 
            !to implement a better code that adapt the array 
            ! size to the furthest particle 
            density_profile(i_step) = density_profile(i_step) + 1
        end do
    end subroutine
end module
    

