##  A MAKEFILE FOR The Ancestry HMM. TO USE, RUN THE COMMAND "make" VIA COMMAND LINE ##

TCFLAGS = -ltcmalloc 
ARMAFLAGS = -larmadillo

all:
	$(LINK.cc) -O3 -o ancestry_hmm ancestry_hmm.cpp $(ARMAFLAGS)
### if tcmalloc is installed, this can be linked using the command below instead
#	$(LINK.cc) -O3 -o ancestry_hmm ancestry_hmm.cpp $(TCFLAGS) $(ARMAFLAGS)

## if you have a local armadillo installation, you will need to provide the directory during compile time and possible also link lblas and lapack 
## our recommendation is to use miniconda3 to do the installation
##	$ conda install -c conda-forge armadillo
## then you will have the appropriate lib and include files in your home directory under subdirectory miniconda3/
## so, replace USERNAME with your unix id on the following line and try this

       #$(LINK.cc) -O3 -o ancestry_hmm ancestry_hmm.cpp -L /home/USERNAME/miniconda3/lib/ -I /home/USERNAME/miniconda3/include/ $(ARMAFLAGS)
      
## if it builds correctly, you will also need to link the library during runtime
## to do this, add the following line to your ~/.bash_profile or ~/.bashrc
##     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/USERNAME/miniconda3/lib/
