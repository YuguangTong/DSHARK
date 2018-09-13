module solver

  use param_mod
  implicit none

contains

  !> Scans through the requested wavenumber interval, computes corresponding frequencies and writes them to output file
  subroutine solve_disp(Nspecies_in, theta_in, delta_in, &
                                ! rf_error_in, eps_error_in, &
                                ! q_in, mu_in, dens_in, &
                                ! beta_para_in, beta_perp_in, &
                                ! kappa_in, drift_in, &
                                ! omega_r_in, omega_i_in, &
                                ! increment_r, increment_i, &
       kstart, kend, nk, &
       solution)

    !--------------------------------------------------------------------------------
    ! input params
    !--------------------------------------------------------------------------------
    integer, intent(in) :: Nspecies_in
    real, intent(in) :: theta_in
    real, intent(in) :: delta_in 
    ! real, optional, intent(in) :: rf_error_in
    ! real, optional, intent(in) :: eps_error_in

    ! real, intent(in), dimension (Nspecies_in) :: q_in
    ! real, intent(in), dimension (Nspecies_in) :: mu_in
    ! real, intent(in), dimension (Nspecies_in) :: dens_in
    ! real, intent(in), dimension (Nspecies_in) :: beta_perp_in
    ! real, intent(in), dimension (Nspecies_in) :: beta_para_in
    ! real, intent(in), dimension (Nspecies_in) :: kappa_in
    ! real, intent(in), dimension (Nspecies_in) :: drift_in

    ! real, intent(in) :: omega_r_in, omega_i_in
    ! real, intent(in) :: increment_r, increment_i
    real, intent(in) :: kstart, kend
    integer, intent(in) :: nk

    !--------------------------------------------------------------------------------
    ! output params
    !--------------------------------------------------------------------------------
    complex, intent(out) :: solution(nk)

    !--------------------------------------------------------------------------------
    ! local variables
    !--------------------------------------------------------------------------------
    real :: omega_r, omega_i
    real :: start, finish
    integer :: ik, n
    real :: dk
    complex :: omega_start, increment
    real, allocatable, dimension (:) :: krange

    !--------------------------------------------------------------------------------
    ! Allocate memory & assign values for variables in param_mod
    !--------------------------------------------------------------------------------
    Nspecies = Nspecies_in
    theta = theta_in*2.0*pi/360.0
    delta = delta_in
    ! omega_r = omega_r_in
    ! omega_i = omega_i_in

    ! rf_error = 1.0e-3
    ! eps_error = 1.0e-4
    ! if (present(rf_error_in)) rf_error = rf_error_in
    ! if (present(eps_error_in)) eps_error = eps_error_in

    write(*,'(A12,F12.8)') 'pi:', pi
    write(*,'(A12,I8)') 'Nspecies:', Nspecies_in
    ! write(*,'(A12,E20.10)') 'rf_error:', rf_error
    ! write(*,'(A12,E20.10)') 'eps_error:', eps_error


    ! allocate(mu(Nspecies),drift(Nspecies),q(Nspecies))
    ! allocate(dens(Nspecies),kappa(Nspecies))
    ! allocate(beta_para(Nspecies),beta_perp(Nspecies),beta_ratio(Nspecies))

    ! do n=1,Nspecies

    !    write(*,'(A12,I8)') 'species:', n
    !    write(*,'(A12,F8.6)') 'q:', q_in(n)
    !    write(*,'(A12,F8.6)') 'mu:', mu_in(n)
    !    write(*,'(A12,F8.6)') 'dens:', dens_in(n)
    !    write(*,'(A12,F8.6)') 'drift:', drift_in(n)
    !    write(*,'(A12,F8.6)') 'beta_para:', beta_para_in(n)
    !    write(*,'(A12,F8.6)') 'beta_per:', beta_perp_in(n)

    !    q(n)=q_in(n)
    !    mu(n)=mu_in(n)
    !    dens(n)=dens_in(n)
    !    drift(n)=drift_in(n)
    !    beta_para(n)=beta_para_in(n)
    !    beta_perp(n)=beta_perp_in(n)
    !    kappa(n)=kappa_in(n)

    ! enddo

    ! beta_ratio=beta_perp/beta_para
    ! !================================================================================

    ! if(omega_r.eq.0.0) omega_r=10.0**(-12) !algorithm cannot handle exact zero
    ! if(omega_i.eq.0.0) omega_i=10.0**(-12) !algorithm cannot handle exact zero
    ! if(theta.eq.0.0) theta=0.0001          !algorithm cannot handle exact zero

    ! omega_start=omega_r+i*omega_i
    ! increment=increment_r+i*increment_i

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! write(*,'(A12,F12.8,A9,F12.8)')  'omega_r:', omega_r, '   omega_i:', omega_i
    write(*,'(A12,F12.8,A10,F12.8)') 'kstart:', kstart, '    kend:', kend
    write(*,'(A12, F12.8)') 'theta:', theta_in
    write(*,'(A12,I8)') 'nk:', nk

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call cpu_time(start)

    allocate(krange(nk))

    dk=(kend-kstart)/(1.0*nk)

    do ik=1,nk
       krange(ik)=kstart+(ik-1)*dk
    enddo

    !scan through wavenumber interval
    ! do ik=1,nk

    !    write(*,*) ' '
    !    write(*,'(A7,I6,A10,F12.8)') '-------',ik,'------- k=', krange(ik)
    !    write(*,'(A12,E20.10,A9,E20.10)')  'guess omega:', real(omega_start), '   gamma:',aimag(omega_start)

    !    !use Muller method to iterate root of dispersion relation
    !    call muller(omega_start,krange(ik),solution(ik))

    !    write(*,'(A9,E20.10,A9,E20.10)')  '   omega:', real(solution(ik)), '   gamma:',aimag(solution(ik))


    !    if ((ik.ge.3).and.(ik.lt.nk))  then

    !       !if three subsequent solutions omega(k) are found, use quadratic polynomial fit 
    !       !to guess next starting frequency for Muller iteration
    !       call polyfit(krange(ik-2:ik+1),solution(ik-2:ik),omega_start)

    !    else

    !       !for the first two solution omega(k) guess next starting frequency for Muller iteration
    !       !by raising the computed omega by an increment which is provided by the user
    !       omega_start=solution(ik)+increment

    !    end if

    ! enddo

    call cpu_time(finish)

    write(*,*) 'Time elapsed: ', finish-start


    ! deallocate(krange)
    ! deallocate(drift)
    ! deallocate(mu,q,kappa,dens)
    ! deallocate(beta_para,beta_perp,beta_ratio)

  end subroutine solve_disp

end module solver
