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

def gen_input_file_text(input_params):
    scalar_keys = ['kstart', 'kend', 'nk', 'omega_r', 'omega_i', 'increment_r', 'increment_i', 'nspecies',
            'theta', 'delta'
           ]
    vector_keys = ['q_in', 'mu_in', 'dens_in', 'drift_in', 'beta_para_in', 'beta_perp_in', 'kappa_in']
    for key in scalar_keys:
        assert key in input_params.keys()
    for key in vector_keys:
        assert key in input_params.keys()
        assert input_params['nspecies'] == len(input_params[key])
    
    input_file_text = ''
    input_file_text += get_wave_number(
        input_params['kstart'], 
        input_params['kend'], 
        input_params['nk']
    )
    input_file_text += get_initial_guess(
        omega_r=input_params['omega_r'], 
        omega_i=input_params['omega_i'], 
        increment_r=input_params['increment_r'],
        increment_i=input_params['increment_i']
    )
    input_file_text += get_setup(
        nspecies=input_params['nspecies'], 
        theta=input_params['theta'], 
        delta=input_params['delta']
    )
    input_file_text += get_accuracy(
        rf_error=1e-4, 
        eps_error=1e-6
    )

    for j in range(len(input_params['q_in'])):
        input_file_text += get_species(
            q_in=input_params['q_in'][j],
            mu_in=input_params['mu_in'][j],
            dens_in=input_params['dens_in'][j],
            drift_in=input_params['drift_in'][j],
            beta_para_in=input_params['beta_para_in'][j],
            beta_perp_in=input_params['beta_perp_in'][j],
            kappa_in=input_params['kappa_in'][j],
        )    
    return input_file_text

def write_input_data(input_file_text, input_file_name="input.dat"):
    with open(input_file_name,"w") as f:
        f.write(input_file_text)
    
def read_disp_data(dispersion_file_name="omega.dat"):
    return pd.read_csv(dispersion_file_name, sep='\s+', names=['k', 'wr', 'wi'])