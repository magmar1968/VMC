module VMC_print
    use VMC_parameters
    implicit none

    integer , private :: FID = 11

    contains

    !####################################################
    !#               Print TWF to File                  #
    !####################################################
    subroutine print_TWF_tofile()
        use HS_puregas,Only:get_TWF_terms
        implicit none
        integer,parameter :: NFUNC = 3 !func,funcprime,funcdoubleprime
        real*8, dimension(NFUNC,TWFNPartitions) :: TWF_corr,TWF_harm
        real*8  :: rstep,r
        integer :: i_step
        character(MAX_FILENAME_LENGHT) :: filename
        filename = trim(outfile_path)//trim(twf_filename)
        
        call get_TWF_terms(TWF_corr=TWF_corr,TWF_harm=TWF_harm,rstep=rstep,&
                           Npartitions=TWFNPartitions)

        open(FID,file=filename,status="replace")
        write(FID,"(A)") "step,r,corr,corrprime,corrdoubleprime,harm,harmprime,harmdoubleprime"
        do i_step=1,TWFNPartitions
            r = (i_step-1)*rstep
            write(FID,*) i_step,r, &
                         TWF_corr(1,i_step),TWF_corr(2,i_step),TWF_corr(3,i_step), &
                         TWF_harm(1,i_step),TWF_harm(2,i_step),TWF_harm(3,i_step)
        end do 
        
        close(FID)
        
    end subroutine

    

    !####################################################
    !#           Print Configuration to File            #
    !####################################################
    subroutine print_configuration_toFile(R,filename)
        use VMC_parameters
        implicit none
        real*8,dimension(Natoms,DIM),intent(in) :: R 
        character(*) :: filename
        integer :: i_atom

        open(FID,file=filename,status="replace")
        write(FID,"(A)") "atom,x,y"
        do i_atom=1, Natoms
            write(FID,*) i_atom,",",R(i_atom,1),",",R(i_atom,2)
        end do
        close(FID)
    end subroutine


    subroutine print_inital_conf_toFile(R)
        use VMC_parameters
        implicit none 
        real*8,intent(in),dimension(Natoms,DIM) :: R
        character(MAX_FILENAME_LENGHT) :: filename 
        filename = trim(outfile_path)//trim(init_conf_filename)
        call print_configuration_toFile(R,filename)
    end subroutine print_inital_conf_toFIle

    subroutine print_final_conf_toFile(R)
        use VMC_parameters
        implicit none
        real*8,intent(in),dimension(Natoms,DIM) :: R
        character(MAX_FILENAME_LENGHT) :: filename
        filename = trim(outfile_path)//trim(fin_conf_filename)
        call print_configuration_toFile(R,filename)
    end subroutine
    !####################################################

    !####################################################
    !#         Print Energies Evolution to file         #
    !####################################################       
    subroutine print_E_evolution_toFile(MC_step,E,Epot,Ekin,Ekinfor)
        use VMC_parameters
        implicit none
        integer,intent(in) :: MC_step
        real*8, intent(in) :: E,Epot,Ekin,Ekinfor
        logical,save ::FirstTime = .TRUE. 
        character(MAX_FILENAME_LENGHT),save :: filename
        
        if (FirstTime) then 
            filename = trim(outfile_path)//trim(energy_evolution_filename)
            open(FID,file=filename,status="replace")
            write(FID,"(A)") "step,E,E^2,epot,ekin,ekinfor"
            FirstTime = .FALSE. 
            close(FID)
        end if

        open(FID,file=filename,access="append") 
        write(FID,*) MC_step,E,E**2,Epot,Ekin,Ekinfor
        close(FID)
    end subroutine

    !##################################################
    !#           Print Density Profile to File        #
    !##################################################
    subroutine print_density_profile_toFile(density_profile)
        use VMC_parameters
        implicit none
        real*8,dimension(NdensProfileSteps) :: density_profile
        integer :: i_step 
        real*8 :: rmin,rmax,area,dens
        character(MAX_FILENAME_LENGHT) :: filename 
        
        filename= trim(outfile_path)//trim(dens_profile_filename)

        open(FID,file=filename,status="replace")
        write(FID,"(A)")"rmin,rmax,dens"

        !normalize for number of evaluations
        density_profile = density_profile/(NMCsteps)
        !normalize for number of particles
        density_profile = density_profile/Natoms
        !normalize for area and printing 
        do i_step = 1, NdensProfileSteps
            rmin = (i_step-1)*densProfileStep
            rmax = (i_step)  *densProfileStep
            area = PI*(rmax**2-rmin**2)
            dens = density_profile(i_step)/area

            write(FID,*) rmin,rmax,dens
        end do 
        close(FID)
    end subroutine
end module 