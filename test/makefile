build:
	gfortran -c -fdefault-real-8 -fPIC param_mod.f90
	f2py -c --fcompiler=gfortran -I. param_mod.o -m py_dshark solve_disp.f90 --f90flags="-fdefault-real-8"
clean:
	rm *.o *.mod *.so *~

test:
	python test_script.py
