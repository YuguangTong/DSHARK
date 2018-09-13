objects=param_mod.o disp_det.o intgrnd.o Zk_func.o Z_func.o dZ_func.o  integrator.o intdeiini.o intdei.o dBessel_Jn.o muller.o polyfit.o  exp_Bessel_In.o exp_dBessel_In.o cerror.o cont_frac.o F21.o F21_cont_frac.o
f90comp = gfortran
options = -fdefault-real-8 -fPIC
f2py_options = --fcompiler=gfortran

build: $(objects) param_mod.mod solve_disp.f90
	f2py -c $(f2py_options) -I $(objects) -m py_dshark solve_disp.f90 --f90flags="-fdefault-real-8"

test_build: param_mod.o solve_disp.f90
	f2py -c $(f2py_options) -I. param_mod.o -m py_dshark solve_disp.f90

param_mod.mod: param_mod.o param_mod.f90
	$(f90comp) -c $(options) param_mod.f90

param_mod.o: param_mod.f90
	$(f90comp) -c $(options) param_mod.f90

disp_det.o: param_mod.mod disp_det.f90
	$(f90comp) -c $(options) disp_det.f90

intgrnd.o: intgrnd.f90
	$(f90comp) -c $(options)  intgrnd.f90

integrator.o: param_mod.mod integrator.f90
	$(f90comp) -c $(options)  integrator.f90

intdeiini.o: intdeiini.f90
	$(f90comp) -c $(options)  intdeiini.f90

intdei.o: intdei.f90
	$(f90comp) -c $(options)  intdei.f90

dBessel_Jn.o: dBessel_Jn.f90
	$(f90comp) -c $(options)  dBessel_Jn.f90

Zk_func.o: param_mod.mod Zk_func.f90
	$(f90comp) -c $(options)  Zk_func.f90

Z_func.o: param_mod.mod Z_func.f90
	$(f90comp) -c $(options)  Z_func.f90

dZ_func.o: dZ_func.f90
	$(f90comp) -c $(options)  dZ_func.f90

cerror.o: param_mod.mod cerror.f90
	$(f90comp) -c $(options)  cerror.f90

cont_frac.o: cont_frac.f90
	$(f90comp) -c $(options)  cont_frac.f90

muller.o: param_mod.mod muller.f90
	$(f90comp) -c $(options)  muller.f90

polyfit.o: polyfit.f90
	$(f90comp) -c $(options)  polyfit.f90

exp_Bessel_In.o: param_mod.mod  exp_Bessel_In.f90
	$(f90comp) -c $(options)  exp_Bessel_In.f90

exp_dBessel_In.o: exp_dBessel_In.f90
	$(f90comp) -c $(options)  exp_dBessel_In.f90

F21.o: param_mod.mod  F21.f90
	$(f90comp) -c $(options)  F21.f90

F21_cont_frac.o: F21_cont_frac.f90
	$(f90comp) -c $(options)  F21_cont_frac.f90

clean:
	rm param_mod.mod $(objects) *.so
	rm -rf py_dshark.*.so.dSYM

test:
	python test_script.py
