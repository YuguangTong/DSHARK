import pandas as pd

def get_wave_number(kstart=.1, kend=.4, nk=30):
    string = \
f'''
&wavenumber
kstart = {kstart}
kend = {kend}
nk = {nk}
/
'''
    return string

def get_initial_guess(omega_r=0., omega_i=.013, increment_r=1e-5, increment_i=3e-3):
    string = \
f'''
&initial_guess
omega_r = {omega_r}
omega_i = {omega_i}
increment_r = {increment_r}
increment_i = {increment_i}
/
'''
    return string

def get_setup(nspecies=2, theta=45., delta=1e-5):
    string = \
f'''
&setup
Nspecies = {nspecies}
theta = {theta}
delta = {delta}
/
'''
    return string

def get_accuracy(rf_error=1e-3, eps_error=1e-6):
    string = \
f'''
&accuracy
rf_error = {rf_error}
eps_error = {eps_error}
/
'''
    return string

def get_species(q_in=1., mu_in=1., dens_in=1., drift_in=0., beta_para_in=4., beta_perp_in=2., kappa_in=7.):
    string = \
f'''
&species
q_in = {q_in}
mu_in = {mu_in}
dens_in = {dens_in}
drift_in = {drift_in}
beta_para_in = {beta_para_in}
beta_perp_in = {beta_perp_in}
kappa_in = {kappa_in}
/
'''
    return string

def write_input_data(input_file_text, input_file_name="input.dat"):
    with open(input_file_name,"w") as f:
        f.write(input_file_text)
    
def read_disp_data(dispersion_file_name="omega.dat"):
    return pd.read_csv(dispersion_file_name, sep='\s+', names=['k', 'wr', 'wi'])