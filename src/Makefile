##  A MAKEFILE FOR The Ancestry HMM. TO USE, RUN THE COMMAND "make" VIA COMMAND LINE ##

TCFLAGS = -ltcmalloc 
ARMAFLAGS = -larmadillo

all:
	$(LINK.cc) -O3 -o ancestry_hmm ancestry_hmm.cpp $(TCFLAGS) $(ARMAFLAGS)


## if you do not have a local install of google perftools, please simply remove the TCFlag link. 
###    $(LINK.cc) -O3 -o ancestry_hmm ancestry_hmm_V0.94.cpp $(ARMAFLAGS)

## if you have a local armadillo installation, you will need to provide the directory during compile time and possible also link lblas and lapack 
