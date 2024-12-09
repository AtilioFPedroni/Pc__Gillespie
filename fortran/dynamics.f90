! ## File: dynamics.f90
! ## - main program: running the SIS dynamics, based on OGA (Optimized Gillespie Algorithm).
! ## See README.md for more information and use
!-----------------------------------------------------------------------------
! SIS epidemic model algorithm based on the article
!           Computer Physics Communications 219C (2017) pp. 303-312
!           "Optimized Gillespie algorithms for the simulation of 
!            Markovian epidemic processes on large and heterogeneous networks"
! Copyright (C) 2017 Wesley Cota, Silvio C. Ferreira
! 
! Please cite the above cited paper (available at <http://dx.doi.org/10.1016/j.cpc.2017.06.007> ) 
! as reference to our code.
! 
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------
! Author    : Wesley Cota
! Email     : wesley@wcota.me
! Date      : 27 Mar 2017
! Version   : 1.0
!-----------------------------------------------------------------------------
! See README.md for more details
! This code is available at <https://github.com/wcota/dynSIS>
! For pure Python, see <https://github.com/wcota/dynSIS-py>
! For NetworkX library, see <https://github.com/wcota/dynSIS-networkx> (NetworkX implementation)

module mod_SIS_OGA
use mod_netdata
use mod_random
use mod_read_tools
implicit none
    
    ! File read/write arguments and parameter
    character*1024                :: f_output, f_temp
    integer, parameter            :: und_output = 10
    
    
    ! Dynamics samples variables
    integer                       :: dynp_i, dynp_sam       ! samples vars
    real*8                        :: dynp_lb                ! lambda infection rate. mu is defined as = 1
    integer                       :: dynp_tmax              ! Maximum time steps
    real*8                        :: dynp_pINI              ! fraction of network first infected (random)
    
    ! Dynamics    
    real*8                        :: dyn_t, dyn_dt          ! Times and time step variables
    
    ! SIS-OGA - Dynamics Variables
    real*8                        :: dyn_m                  ! m = M/R. 1 - m = w = W/R
    real*8                        :: dyn_R                  ! Total rate
    integer, allocatable          :: dyn_VI(:), dyn_sig(:)  ! Lists V^I and sigma
    integer                       :: dyn_NI         ! # of infected vertex N_I and # of infected edges N_k
    
    ! SIS-OGA - Network structure variables
    integer                       :: net_kmax               ! Used in the rejection probability
    
    integer                       :: aux_expoente
    real*8                        :: aux_pg_do_tempo_atual
    real*8, parameter             :: aux_escala = 1.05d0
    
contains
    
    subroutine read_dyn_parameters()
        call read_f(dynp_lb,"Value of infection rate lambda (mu is defined as equal to 1): ")
        call read_i(dynp_tmax,"Maximum time steps (it stops if the absorbing state is reached): ")
        call read_f(dynp_pINI,"Fraction of infected vertices on the network as initial condition (is random for &
 & each sample): ")
        
        ! Allocate the SIS-OGA lists V^I
        allocate(dyn_VI(net_N),dyn_sig(net_N))
    end subroutine
    
    subroutine random_initial_condition()
        integer :: ver, vti
        
        dyn_sig = 0 ! sigma
        !dyn_VI = 0 ! list V^I (not needed)
        dyn_NI = 0 ! N_I
        
        ! Sort vertices and apply the initial condition
        do vti = 1, int(net_N*dynp_pINI)
            vti_ver : do
                ver = random_int(1,net_N)
                if (dyn_sig(ver) == 0) then
                    dyn_NI = dyn_NI + 1
                    dyn_VI(dyn_NI) = ver
                    dyn_sig(ver) = 1
                    exit vti_ver
                endif
            enddo vti_ver
        enddo
        
        dyn_t = 0d0
        aux_expoente = 0
        aux_pg_do_tempo_atual = aux_escala**aux_expoente
    end subroutine    

    subroutine dyn_run()
    
        integer :: pos_inf
        integer :: ver, pos_nei
        real*8  :: rnd
    
        call random_initial_condition()
        
        dyn_time_loop : do while (dyn_t <= dynp_tmax)
            
            ! Calculate the total rate
            dyn_R = (dyn_NI + 1d0*dynp_lb * dyn_NI)
            
            ! Select the time step
            rnd = max(random_d(), 1e-12) ! Avoid rnd = 0
            dyn_dt = -log(rnd) / dyn_R
            
            ! Update the time
            dyn_t = dyn_t + dyn_dt
            
            ! Probability m to heal
            dyn_m = 1d0*dyn_NI/ dyn_R
            
            ! Try to heal
            rnd = random_d()
            if (rnd < dyn_m) then
                ! Select one infected vertex from V^I
                pos_inf = random_int(1,dyn_NI) 
                ver = dyn_VI(pos_inf)
                
                ! Then, heal it
                dyn_sig(ver) = 0
                dyn_VI(pos_inf) = dyn_VI(dyn_NI) ! Swap positions
                dyn_NI = dyn_NI - 1             ! Then, short the list
                
            ! If not, try to infect: w = 1 - m
            else
                ! Select the infected vertex i with prob. proportional to k_i
                select_infec : do
                    pos_inf = random_int(1,dyn_NI)
                    ver = dyn_VI(pos_inf)
                    if (random_d() < 1d0*net_k(ver)/(1d0*net_kmax)) exit select_infec
                enddo select_infec
                
                ! select one of its neighbors
                pos_nei = random_int(net_ini(ver) , net_ini(ver) + net_k(ver) - 1)
                ver = net_con(pos_nei)
                
                ! if not a phantom process, infect
                if (dyn_sig(ver) == 0) then
                    dyn_sig(ver) = 1
                    dyn_NI = dyn_NI + 1     ! Increase by 1 the list
                    dyn_VI(dyn_NI) = ver    ! Add one element to list
                endif
            endif
            
            ! if a absorbing state is reached, exit
            if (dyn_NI == 0) exit dyn_time_loop
            
            ! Try to save the dynamics by time unit
            do while (dyn_t >= aux_pg_do_tempo_atual)

                write(und_output, *) dyn_t, aux_pg_do_tempo_atual, 1d0*dyn_NI/net_N

                aux_expoente = aux_expoente + 1
                aux_pg_do_tempo_atual = aux_escala**aux_expoente
            enddo
            
        enddo dyn_time_loop
        
    end subroutine

    subroutine run_prog()
        call print_info('################################################################################')
        call print_info('### Optimized Gillespie algorithms for the simulation of Markovian epidemic  ###')
        call print_info('############ processes on large and heterogeneous networks: SIS-OGA ############')
        call print_info('##============ Copyright (C) 2017 Wesley Cota, Silvio C. Ferreira ============##')
        call print_info('##===== Paper available at <http://dx.doi.org/10.1016/j.cpc.2017.06.007> =====##')
        call print_info('##======= The codes are available at <https://github.com/wcota/dynSIS> =======##')
        call print_info('##======== Please cite the above cited paper as reference to our code ========##')
        call print_info('##=== This code is under GNU General Public License. Please see README.md. ===##')
        call print_info('################################################################################')
    
        ! initial value of input counter to read command arguments, used by mod_read_tools
        inp_pos = 1 
    
        ! Files arguments read
        call read_arg(f_input)
        call read_arg(f_output)
    
        ! Read network to memory
        call readEdges()
    
        ! To be used in the SIS-OGA algorithm. Calculate the k_max of the network
        net_kmax = maxval(net_k)
    
        ! Initate the random generator
        call print_info('')
        call print_progress('Generating random seed')
        call random_ini()
        call print_done()
    
        ! We are ready! All network data is here, now we need to read the dynamical parameters.
        call print_info('')
        call print_info('Now we need to read the dynamical parameters.')
        call read_dyn_parameters()
    
        ! Let's run the dynamics.
        call print_info('')
        call print_info('Running dynamics...')
    
        ! Loop over all the samples
    
        ! Open file and write info
        open(und_output,file=f_output)
    
        ! Run dynamics
        call dyn_run()
            
        ! Close output file
        close(und_output)
        call print_done()
    
        call print_info('')
        call print_info('Everything ok!')
        call print_info('Input file: '//trim(adjustl(f_input)))
        call print_info('Output file: '//trim(adjustl(f_output)))
        call print_info('')
        call print_info('*****Algorithm used: Optimized Gillespie Algorithm for SIS (SIS-OGA, Fortran)*****')
        call print_info('Codes available at <https://github.com/wcota/dynSIS>.')
        
    end subroutine

end module

program dynSIS
use mod_SIS_OGA
implicit none
    
    call run_prog()
    
end program
