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
    use VMC,Only:run_VMC
    use VMC,Only:VMC_varParameters,VMC_results
    implicit none

    real    :: T1,T2
    integer :: secs,min, delta
    type(VMC_varParameters) ::  params
    type(VMC_results) :: results

    call cpu_time(T1)
    call input("./data/inputfile.dat")
    params = VMC_varParameters(Rv=20,&
                            alpha=36.0,&
                            k0   =0.010459459459459461,&
                            a_osc=50)

    call run_VMC(params=params,results=results)
    call cpu_time(T2)
    
    delta = int(T2-T1)
    min = delta/60
    secs = int(delta - min*60)

    print *, "######################################################"
    print *, "RESULTS:"
    print *, "  E:       ", results%E
    print *, "  error:   ", results%error
    print *, "  accrate: ", results%acceptrate 
    print *, "  comp time ", delta
    print *, "          -> mins: ", min
    print *, "          -> secs: ", secs
    print *, "######################################################"

end program parameters_optimization

