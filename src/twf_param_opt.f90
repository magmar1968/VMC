!author: Lorenzo Magnoni
!data: 23 october 2023
!mail: lorenzo.magnoni@studenti.unimi.it or lolemagno99@gmail.com
!code: The aim of this code is to optimize the parameters of the 
!      the trial function via a VMC algorithm. The trial function is 
!      composed by a bessel function linked with a negative exponential
!      function of the form 1-c*exp{-r/alpha}. The linking point Rv and 
!      alpha are the variational parameters evaluated by this program.
!      The code make use of a VMC module for the variatonal algorithm part. 

program parameters_optimization
    call main()
     
contains 
subroutine main()    
    use VMC,Only:run_VMC
    use VMC,Only:VMC_varParameters,VMC_results
    use iofile, Only: io_open
    implicit none

    real    :: T1,T2
    integer :: secs,min, delta
    type(VMC_varParameters) ::  params
    type(VMC_results) :: results
    real*8,allocatable,dimension(:) :: alpha_array
    real*8,allocatable,dimension(:) :: Rv_array
    real*8,allocatable,dimension(:) :: k0_array
    integer :: count,N_lines
    character(len= 60) :: filename =  "./data/variational_parameters.dat"
    
    !read variational parameters
    call io_open(filename,N_lines)
    allocate(alpha_array(N_lines-1))
    allocate(Rv_array(N_lines-1   ))
    allocate(k0_array(N_lines-1   ))
    
    
    call read_variational_parameters(filename,alphas=alpha_array,&
                                             Rvs=Rv_array,&
                                             k0s=k0_array)
    call input("./data/inputfile.dat")
    do count = 1,size(alpha_array)
         call cpu_time(T1)
       
         params = VMC_varParameters(Rv=Rv_array(count),&
                                 alpha=alpha_array(count),&
                                 k0   =k0_array(count),&
                                 a_osc=50)
       
         print *,params
         ! call run_VMC(params=params,results=results)
         ! call cpu_time(T2)
         ! delta = int(T2-T1)
         ! min = delta/60
         ! secs = int(delta - min*60)
         ! print *, "######################################################"
         ! print *, "RESULTS:"
         ! print *, "  E:       ", results%E
         ! print *, "  error:   ", results%error
         ! print *, "  accrate: ", results%acceptrate 
         ! print *, "  comp time ", delta
         ! print *, "          -> mins: ", min
         ! print *, "          -> secs: ", secs
         ! print *, "######################################################"
     end do
    
    deallocate(alpha_array)
    deallocate(k0_array)
    deallocate(Rv_array)
end subroutine



subroutine read_variational_parameters(filename,alphas,Rvs,k0s)
   
    Use, intrinsic :: iso_fortran_env, Only : iostat_end
    use strings, only: parse,value
    implicit none
    
    character(*),intent(in) :: filename
    real*8,dimension(:),intent(out) :: alphas,RVs,k0s 
    integer :: FID = 33,IERROR,nargs
    character(60) :: ioerrmsg
    character(60) :: line
    character(50),dimension(10) :: args
    integer :: count, ios
    
    open(FID,file=filename,status='old',action='read')
    IERROR = 0
    read(FID,"(A)",iostat= IERROR,iomsg=ioerrmsg) !read first line
    count = 1
    do while(IERROR == 0)
        read(FID,"(A)",iostat= IERROR,iomsg=ioerrmsg) line

        Select Case(IERROR)
        Case(0)
            call parse(line,delims=" ",args= args,nargs=nargs )
            call value(args(2),RVs(count),ios=ios)
            call value(args(3),alphas(count),ios=ios)
            call value(args(4),k0s(count),ios=ios)
            count = count + 1
        Case(iostat_end)
            exit
        Case Default 
            print *, "IERROR: ", IERROR 
            print *, ioerrmsg
        End Select
    end do 
    close(FID)
end subroutine


end program parameters_optimization

