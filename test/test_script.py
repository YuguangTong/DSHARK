import py_dshark

Nspecies_in = 3
theta_in = 3.0
theta_in = 80.
delta_in = 1e-4
kstart = 0.2
kend = 0.6
nk = 30
rf_error_in = 1e-2
eps_error_in = 1.0

q_in = [1, -1, -1]
mu_in = [1, 1836, 1836]
dens_in = [1, 0.8, 0.2]
beta_para_in = [0.2, 1.6, 0.9]
beta_perp_in = [0.7, 0.3, 1.1]
kappa_in = [4, 5, 7.6]
drift_in = [0, -0, 36.7]

omega_r_in = 0.335
omega_i_in = 0.556
increment_r = 0.0
increment_i = 0.37

sol = py_dshark.solver.solve_disp(
    nspecies_in=Nspecies_in,
    theta_in=theta_in,
    delta_in=delta_in,
    kstart=kstart,
    kend=kend,
    nk=nk,
    rf_error_in=rf_error_in,
    eps_error_in=eps_error_in,
    q_in=q_in,
    mu_in=mu_in,
    dens_in=dens_in,
    beta_para_in=beta_para_in,
    beta_perp_in=beta_perp_in,
    kappa_in=kappa_in,
    drift_in=drift_in,
    omega_r_in=omega_r_in,
    omega_i_in=omega_i_in,
    increment_r=increment_r,
    increment_i=increment_i
    )

