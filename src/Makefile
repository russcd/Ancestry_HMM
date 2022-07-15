##  A MAKEFILE FOR The Ancestry HMM. TO USE, RUN THE COMMAND "make" VIA COMMAND LINE ##

TCFLAGS = -ltcmalloc 
ARMAFLAGS = -larmadillo
CONDAFLAGS = -fvisibility-inlines-hidden -std=c++17 -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem ${CONDA_PREFIX}/include

all:
	$(CXX) -O3 $(CXXFLAGS) -o ancestry_hmm ancestry_hmm.cpp  $(ARMAFLAGS)
	$(CXX) -O3 $(CXXFLAGS) -o ahmm-s ahmms.cpp  $(ARMAFLAGS)
	
conda:
	$(CXX) -O3 $(CONDAFLAGS) -o ahmm-s ahmms.cpp -L ${CONDA_PREFIX}/lib -I ${CONDA_PREFIX}/include $(ARMAFLAGS)
	$(CXX) -O3 $(CONDAFLAGS) -o ancestry_hmm ancestry_hmm.cpp -L ${CONDA_PREFIX}/lib -I ${CONDA_PREFIX}/include $(ARMAFLAGS)
	
ahmms:
	$(CXX) -O3 $(CXXFLAGS) -o ahmm-s ahmms.cpp $(ARMAFLAGS)
	
ahmm:
	$(CXX) -O3 $(CXXFLAGS) -o ancestry_hmm ancestry_hmm.cpp $(ARMAFLAGS)

	
## if you have a local install of google perftools, please add  a TCFlag link. 
###  $(LINK.cc) -std=c++11 -O3 -o ahmm-s ahmms.cpp $(ARMAFLAGS)  $(TCFLAGS)

## if you have a local armadillo installation, you will need to provide the directory during compile time and possible also link lblas and lapack 
## our recommendation is to use miniconda3 to do the installation
##	$ conda install -c conda-forge armadillo
## then you will have the appropriate lib and include files in your home directory under subdirectory miniconda3/
## so, replace USERNAME with your unix id on the following line and try this

##	$(CXX) -std=c++11 -O3 -o ahmm-s ahmms.cpp -L ${CONDA_PREFIX}/lib -I ${CONDA_PREFIX}/include $(ARMAFLAGS)
      
## if it builds correctly, you may also need to link the library during runtime
## to do this, add the following line to your ~/.bash_profile or ~/.bashrc
##     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib

