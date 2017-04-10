##  A MAKEFILE FOR The Ancestry HMM. TO USE, RUN THE COMMAND "make" VIA COMMAND LINE ##

TCFLAGS = -ltcmalloc 
ARMAFLAGS = -larmadillo

all:
	$(LINK.cc) -O3 -o ancestry_hmm ancestry_hmm_V0.93.cpp $(TCFLAGS) $(ARMAFLAGS)


## if you do not have a local install of google perftools, please simply remove the TCflag. 
###    $(LINK.cc) -O3 -o ancestry_hmm ancestry_hmm_V0.93.cpp $(ARMAFLAGS)