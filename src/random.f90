!author: Lorenzo Magnoni
!date:   28/09/23
!code: This module is intended to contain usefull functions for the 
!      generation of random number for all possible distribution
!
!#############################################################
! File contents:
!
! Subroutine                           Label               Line    
! ----------                           ---                 ----
! Initial Random Seed                  init_random_seed
!
! Uniform Random Number Generator      uniform
!
! Gaussian Random Number Generator     gauss         
!#######################################################################
!Last Update 07/11/2023 L.Magnoni lorenzo.magnoni@studenti.unimi.it

module random
    contains

    !###################################
    !       Initialize Random Seed
    !###################################
    subroutine init_random_seed()

        integer              :: i, n, clock
        integer, allocatable :: seed(:)
        
        call random_seed(size = n)
        allocate(seed(n))
        
        call system_clock(count=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /) + node
        
        call random_seed(put = seed)
        
        deallocate(seed)
        
    end subroutine init_random_seed
        

    !##############################################
    !     Uniform random number generator
    !##############################################
    !usage: generate a uniform distributed random number in [min, max] interval
    !input: - min (default 0) left limit of the interval
    !       - max (default 1) right limit of the inteval
    function uniform(min,max) result(w)
        implicit none

        real :: w                    !output 
        real, optional :: min, max   ! input (optional)
        real :: a_min,a_max ! a is short for "actual..."

        if (present(min)) then
            a_min = min
        else 
            a_min = 0
        end if 

        if (present(max)) then
            a_max = max
        else 
            a_max = 1
        end if

        call random_number(w)
        w = w*(a_max-a_min) + a_min

    end function uniform


    !##############################################
    !     Gaussian random number generator
    !##############################################
    !usage: generate a gaussian distributed random number using the BOX-MULLER method
    !input: - mu (default 0): the mu of the normal distribution
    !       - variance (default 1): the standar deviation of the normal distribution  
    double precision function gauss(sigma,mu)

        implicit none
        real*8, intent(in),optional :: sigma,mu
        real*8                  :: m_sigma,m_mu
        real*8                  :: w, v1, v2, l
        real*8                  :: s1, s2
        w = 2.d0

        if(present(sigma)) then 
            m_sigma = sigma
        else 
            m_sigma = 1. 
        end if 

        if(present(mu)) then
            m_mu = mu
        else 
            m_mu = 0.
        end if 

        do
            call random_number(s1)
            call random_number(s2)

            v1 = 2.d0*s1 - 1.d0
            v2 = 2.d0*s2 - 1.d0
            w = v1*v1 + v2*v2

            if (w.lt.1.d0) exit
        end do

        l     = v1*sqrt(-2.d0*log(w)/(w))
        l     = m_sigma*l + m_mu 
        gauss = l
    end function gauss
    !##############################################
end module