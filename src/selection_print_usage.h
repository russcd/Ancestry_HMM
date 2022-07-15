#ifndef __SELECTION_PRINT_USAGE_H
#define __SELECTION_PRINT_USAGE_H

void print_usage() {
    
    cerr << endl << endl << "ahmm-s usage:" << endl << endl ;
    cerr << "\trequired:" << endl ;
    cerr << "\t\t-i [string]\n\t\t\tinput file name" << endl ;
    cerr << "\t\t-s [string]\n\t\t\tsample id and ploidy file" << endl ;

    cerr << "\t\t-p [int] [int] [float]" << endl ;
    cerr << "\t\t\t ancestry pulse with format, ancestral population, time," << endl ;
    cerr << "\t\t\t and proportion of final ancestry from this pulse" << endl ;
    cerr << "\t\t--ne [int]\n\t\t\teffective population size of the admixed population" << endl ;

    cerr << "\n\tselect one of the following working modes:" << endl ;
    
    cerr << "\t\t--gss [int] [int] [int] [float] [float]" << endl ;
    cerr << "\t\t\t golden section search for optimal selection coeffient at each site." << endl ;
    cerr << "\t\t\t parameters: chromosomal position start, stop, step, selection coefficient start, stop" << endl ;

    cerr << "\t\t--grid [int] [int] [int] [float] [float] [float]" << endl ;
    cerr << "\t\t\t calculate likelihood ratios in a grid." << endl ;
    cerr << "\t\t\t parameters: chromosomal position start, stop, step, selection coefficient start, stop, step." << endl ;

    cerr << "\t\t--site [int] [float]" << endl ;
    cerr << "\t\t\t calculate likelihood ratios for a single value of s at a single site."  << endl ;
    cerr << "\t\t\t parameters: chromosomal position, selective coeffient" << endl ;


    cerr << "\n\toptional:" << endl ;
    cerr << "\t\t--help\n\t\t\tprint this help statement" << endl ;
    cerr << "\t\t-g\n\t\t\tsamples are specified with genotypes rather than read counts" << endl ;
    
    cerr << "\t\t--chr [string]" << endl ;
    cerr << "\t\t\t specify chromosome that will be analyzed" << endl ;
    cerr << "\t\t\t (only necessary when there are multiple chromosomes in input file)" << endl ;
    cerr << "\t\t--chr_win [int] [int]" << endl ;
    cerr << "\t\t\t limit region on chromosome that will be analyzed" << endl ;
    
    cerr << "\t\t--gss_precision [float]" << endl ;
    cerr << "\t\t\t specify precision in finding optimal value of s using golden section search. default: 1e-5" << endl ;
    cerr << "\t\t--unit_coords" << endl ;
    cerr << "\t\t\t unit for start and stop position in grid and gss search can be defined as chromosome" << endl ;
    cerr << "\t\t\t coordinates rather than as line number in input file. default off" << endl ;
    cerr << "\t\t--window [string] [float]" << endl ;
    cerr << "\t\t\t specify size of Markov chain in percent or Morgans." << endl ;
    cerr << "\t\t\t \"p 10\" extends the markov chain 10% of chromosome length on each side of selected site." << endl ;
    cerr << "\t\t\t \"m 0.1\" extends the windows 0.1 Morgan on each side of the selected site." << endl ;
    cerr << "\t\t\t default: \"p 100\"" << endl ;
    cerr << "\t\t--traj [int]" << endl ;
    cerr << "\t\t\t change algorithm for generating selection trajectories." << endl ;
    cerr << "\t\t\t 4: 4-point approximation, 3: 3-point approximation (legacy option, not recommended)." << endl ;
    cerr << "\t\t\t default: forward iteration." << endl ;
    cerr << "\t\t--stochastic" << endl ;
    cerr << "\t\t\t enables the stochastic method for generation selection trajectory." << endl ;
    cerr << "\t\t\t (Experimental. Slow. Use for small values of s.)" << endl ;
    cerr << "\t\t--stochastic_reps [int]" << endl ;
    cerr << "\t\t\t specifies number of simulations for the stochastic trajectory algorithm." << endl ;
    cerr << "\t\t\t default: 10000" << endl ;
    cerr << "\t\t--full_selection_space" << endl ;
    cerr << "\t\t\t turns off optimization of the selection coeffient search space. (Experimental)" << endl ;
}

#endif

