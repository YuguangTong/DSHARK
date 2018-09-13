module solver

  use param_mod
  implicit none

contains

  subroutine solve_disp(Nspecies_in, theta_in, delta_in, &
       rf_error_in, eps_error_in, &
       q_in, mu_in, dens_in, &
       beta_para_in, beta_perp_in, &
       kappa_in, drift_in, &
       omega_r_in, omega_i_in, &
       increment_r, increment_i, &
       kstart, kend, nk, solution)

    !--------------------------------------------------------------------------------
    ! input params
    !--------------------------------------------------------------------------------
    integer, intent(in) :: Nspecies_in
    real*8, intent(in) :: theta_in
    real*8, intent(in) :: delta_in 
    real*8, intent(in) :: kstart, kend
    integer, intent(in) :: nk
    real*8, optional, intent(in) :: rf_error_in
    real*8, optional, intent(in) :: eps_error_in

    real*8, intent(in), dimension (Nspecies_in) :: q_in
    real*8, intent(in), dimension (Nspecies_in) :: mu_in
    real*8, intent(in), dimension (Nspecies_in) :: dens_in
    real*8, intent(in), dimension (Nspecies_in) :: beta_perp_in
    real*8, intent(in), dimension (Nspecies_in) :: beta_para_in
    real*8, intent(in), dimension (Nspecies_in) :: kappa_in
    real*8, intent(in), dimension (Nspecies_in) :: drift_in
    real*8, intent(in) :: omega_r_in, omega_i_in
    real*8, intent(in) :: increment_r, increment_i
    ! --------------------------------------------------------------------------------
    ! output params
    !--------------------------------------------------------------------------------
    complex, intent(out) :: solution(nk)

    !--------------------------------------------------------------------------------
    ! local variables
    !--------------------------------------------------------------------------------
    real*8 :: omega_r, omega_i
    real*8 :: start, finish
    integer :: ik, n
    real*8 :: dk
    complex*16 :: omega_start, increment
    real*8, allocatable, dimension (:) :: krange

    !--------------------------------------------------------------------------------
    ! Allocate memory & assign values for variables in param_mod
    !--------------------------------------------------------------------------------
    Nspecies = Nspecies_in
    theta = theta_in*2.0*pi/360.0
    delta = delta_in
    rf_error = 1.0e-4
    eps_error = 1.0e-3
    if (present(rf_error_in)) rf_error = rf_error_in
    if (present(eps_error_in)) eps_error = eps_error_in

    write(*,'(A12,F12.8)') 'pi:', pi
    write(*,'(A12,I8)') 'Nspecies:', Nspecies_in
    write(*,'(A12,E20.10)') 'rf_error:', rf_error
    write(*,'(A12,E20.10)') 'eps_error:', eps_error

    write(*,'(A12, F12.8)') 'theta:', theta_in
    write(*,'(A12, F12.8)') 'delta:', delta_in
    write(*,'(A12,F12.8,A10,F12.8)') 'kstart:', kstart, '    kend:', kend
    write(*,'(A12,I8)') 'nk:', nk

    print *, q_in
    print *, mu_in
    print *, dens_in
    print *, beta_para_in
    print *, beta_perp_in

    allocate(mu(Nspecies),dens(Nspecies),q(Nspecies))
    allocate(drift(Nspecies), kappa(Nspecies))
    allocate(beta_para(Nspecies),beta_perp(Nspecies),beta_ratio(Nspecies))

    do n=1,Nspecies

       ! print *, 'species:', n
       ! print *, 'q:', q_in(n)
       ! print *, 'mu:', mu_in(n)
       ! print *, 'dens:', dens_in(n)
       ! print *, 'drift:', drift_in(n)
       ! print *, 'beta_para:', beta_para_in(n)
       ! print *, 'beta_per:', beta_perp_in(n)

       q(n)=q_in(n)
       mu(n)=mu_in(n)
       dens(n)=dens_in(n)
       drift(n)=drift_in(n)
       beta_para(n)=beta_para_in(n)
       beta_perp(n)=beta_perp_in(n)
       kappa(n)=kappa_in(n)

    enddo

    beta_ratio=beta_perp/beta_para
    !================================================================================

    print *, 'omega_r:', omega_r_in, '   omega_i_in:', omega_i_in
    omega_r = omega_r_in
    omega_i = omega_i_in
    
    if(omega_r.eq.0.0) omega_r=10.0**(-12) !algorithm cannot handle exact zero
    if(omega_i.eq.0.0) omega_i=10.0**(-12) !algorithm cannot handle exact zero
    if(theta.eq.0.0) theta=0.0001          !algorithm cannot handle exact zero

    omega_start=omega_r+i*omega_i
    
    print *, 'omega_start:', omega_start
    
    increment=increment_r+i*increment_i
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, 'omega_r:', omega_r, '   omega_i:', omega_i
    print *, 'kstart:', kstart, '    kend:', kend
    print *, 'theta:', theta_in
    print *, 'nk:', nk
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
    !    !call muller(omega_start,krange(ik),solution(ik))

    !    write(*,'(A9,E20.10,A9,E20.10)')  '   omega:', real(solution(ik)), '   gamma:',aimag(solution(ik))


    !    if ((ik.ge.3).and.(ik.lt.nk))  then

    !       !if three subsequent solutions omega(k) are found, use quadratic polynomial fit 
    !       !to guess next starting frequency for Muller iteration
    !       !call polyfit(krange(ik-2:ik+1),solution(ik-2:ik),omega_start)
    !       print *, 'dummy line'
    !    else

    !       !for the first two solution omega(k) guess next starting frequency for Muller iteration
    !       !by raising the computed omega by an increment which is provided by the user
    !       omega_start=solution(ik)+increment

    !    end if

    ! enddo
    
    call cpu_time(finish)

    write(*,*) 'Time elapsed: ', finish-start

    deallocate(krange)
    deallocate(drift)
    deallocate(mu,q,kappa,dens)
    deallocate(beta_para,beta_perp,beta_ratio)
    
  end subroutine solve_disp
  
end module solver
