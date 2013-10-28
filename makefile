# master makefile 

default:
	@echo Specify what to make:
	@echo " test - compile and run a simple test"
	@echo " test-openmp - compile and run a simple OpenMP test"
	@echo " mex - matlab/octave mex files"
	@echo " mex-openmp - matlab/octave mex files (OpenMP)"
	@echo " mwrap - matlab wrapper generator"
	@echo " all - all of the above"
	@echo " clean - remove object files"

SRC = ./src
DOC = ./doc
MATLAB = ./matlab
EXAMPLES = ./examples
MWRAP = ./packages/mwrap-0.33

test: 
	cd $(EXAMPLES); $(MAKE) HOST=linux-gfortran 
	
test-openmp: 
	cd $(EXAMPLES); $(MAKE) HOST=linux-gfortran-openmp 
	
mex: 
	cd $(MATLAB); $(MAKE) -f makefile.mwrap TARGET=octave-linux
	
mex-openmp: 
	cd $(MATLAB); $(MAKE) -f makefile.mwrap TARGET=octave-linux-openmp
	
mwrap: 
	cd $(MWRAP); $(MAKE) 
	cp -f $(MWRAP)/mwrap ./bin

all: test test-openmp mex mex-openmp mwrap


clean:
	cd $(SRC); $(MAKE) clean
	cd $(EXAMPLES); $(MAKE) clean
	cd $(MATLAB); $(MAKE) clean

distclean:
	cd $(SRC); $(MAKE) distclean
	cd $(EXAMPLES); $(MAKE) distclean
	cd $(MATLAB); $(MAKE) distclean



