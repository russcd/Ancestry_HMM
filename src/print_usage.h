#ifndef __PRINT_USAGE_H
#define __PRINT_USAGE_H

void print_usage() {
    
    cerr << endl << endl << "ancestry_hmm usage:" << endl << endl ;
    cerr << "\trequired:" << endl ;
    cerr << "\t\t-i [string]\t\tinput file name" << endl ;
    cerr << "\t\t-s [string]\t\tsample id and ploidy file" << endl ;
    cerr << "\t\t-a [int] [float] [float] ..." << endl ;
    cerr << "\t\t\tnumber of ancestral populations and ancestry proportion attributable to each" << endl ;
    cerr << "\t\t-p [int] [int] [float]" << endl ;
    cerr << "\t\t\tancestry pulse with format, ancestral population, time," << endl ;
    cerr << "\t\t\tand proportion of final ancestry from this pulse" << endl ;
    cerr << "\t\t\tnegative time or proportions indicate that parameters are to be estimated" << endl << endl ;
    
    cerr << "\toptional:" << endl ;
    cerr << "\t\t--help\t\t\tprint this help statement" << endl ;
    cerr << "\t\t--ne [int]\t\teffective population size of the admixed population" << endl ;
    cerr << "\t\t-g\t\t\tsamples are specified with genotypes rather than read counts" << endl ;
    cerr << "\t\t--precision [int]\tmodify float and double precision to int" << endl ;
    cerr << "\t\t-v\t\t\tviterbi decoding" << endl ;
    cerr << "\t\t-b [int] [int]\t\tnumber of bootstraps and bootstrap block size in number of SNPs" << endl ;
    cerr << "\t\t--tmax [int]\t\tmaximum time of an admixture pulse" << endl ;
    cerr << "\t\t--tmin [int]\t\tminimum time of an admixture pulse" << endl ;
    cerr << "\t\t--tolerance [float]\tdistance in lnL units to just convergence" << endl ;
    cerr << "\t\t-e [float]\t\terror rates" << endl ;
    cerr << "\t\t-E\t\t\tsite specific error rates are included" << endl ;
    cerr << "\t\t--fix\t\t\tancestral allele frequencies are certain" << endl << endl ;
    
    cerr << "\toptional and relevant only for multiple pulse models:" << endl ;
    cerr << "\t\t--output-ancestry\toutput ancestry posteriors rather than pulses" << endl ;
    cerr << "\t\t-r [int]\t\tnumber of random restarts during nelder-mead optimization" << endl ;
    cerr << "\t\t--pmax [int]\t\tmaximum proportion ancestry in an admixture pulse" << endl ;
    cerr << "\t\t--pmin [int]\t\tminimum proportion ancestry in an admixture pulse" << endl ;
    
}

#endif

