import numpy as np
import py_dshark

kstart = .2
kend = 0.4
nk = 20
theta = 80.
omega_r = 0.035
omega_i = -1e-6

krange, sol = py_dshark.solver.solve_disp(
    beta_para_in=np.array([0.25, 0.45, 0.3 ]),
    beta_perp_in=np.array([0.25, 0.45, 0.3 ]),
    delta_in=0.0001,
    dens_in=np.array([1. , 0.9, 0.1]),
    drift_in=np.array([ 0.,  0., 0.]),
    eps_error_in=1e-6,
    increment_i=0,
    increment_r=0,
    kappa_in=np.array([50, 50, 50]),
    mu_in=np.array([1.000e+00, 1.836e+03, 1.836e+03]),
    nspecies_in=3,
    q_in=np.array([ 1., -1., -1.]),
    rf_error_in=0.0001,
    kstart=kstart,
    kend=kend,
    nk=nk,
    theta_in=theta,
    omega_r_in=omega_r,
    omega_i_in=omega_i
    )
print('sol:')
for k, s in zip(krange, sol):
    print('k = {0:.3g}, w = {1:.3g} + {2:.3g} i'.format(k, s.real, s.imag))
