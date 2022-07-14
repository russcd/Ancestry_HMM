#ifndef __SELECTION_CMD_LINE_H
#define __SELECTION_CMD_LINE_H

/// command line information and global parameters
class cmd_line {
public:
    
    /// terms to bound the optimization
    double t_max ;
    double t_min ;
    
    /// to bound proportion search
    double p_max ;
    double p_min ;
    
    /// to create intial simplex points
    double t_length ;
    double p_length ;
    
    /// number of restarts
    int n_restarts ;
    
    /// proportion ancestry for 0-n must sum to 1
    /// these are therefore the final ancestry proportion, not necessary the proportion that fluxed if additional pulses occured closer to the present
    vector<double> ancestry_proportion ;
    
    /// store relevant ancestry information
    vector<pulse> ancestry_pulses ;
    
    /// diploid effective population size ( i.e. 2n )
    double ne ;
    
    /// tolerance for parameter search
    double tolerance ;
    
    /// minimum recombinational distance between markers
    double minimum_distance ;
    
    /// error rates for reads (if read based) or genotypes (if genotype based)
    double error_rate ;
    
    /// bool sample is expressed as genotypes, not read counts
    bool genotype ;
    
    /// ancestral genotype frequencies are fixed
    bool ancestral_fixed ;
    
    /// viterbi output
    /// caution: not recommended for more samples of ploidy > 1
    bool viterbi ;
    
    /// number of digits of precision to include
    int precision ;
    
    /// output actual pulses rather than ancestry states
    bool output_pulses ;
    
    /// error rates specifed
    bool error_rates ; 
    
    /// input file name
    string input_file ;
    
    /// sample file
    string sample_file ;
    
    /// bootstrap
    int n_bootstraps ;
    int block_size ; 

    /// +=+=+=+=+=+=+ selection +=+=+=+=+=+=+


    bool is_limit; // limits to only one chromosome
    string limit_chr ; // specifies which chromosome
    int limit_win_start ; // window to analyse start
    int limit_win_end ; // window to analyze end

    // grid search
    bool calc_grid; // set grid search
    int grid_pstart; // position start
    int grid_pstop;
    int grid_pstep;
    double grid_sstart; // selection coeffient start
    double grid_sstop;
    double grid_sstep;

    // single site test
    bool test_point;
    int test_pos;
    double test_sel;

    // golden section search
    bool run_gss; // set gss search
    int gs_pstart;
    int gs_pstop;
    int gs_pstep;
    double gs_sstart;
    double gs_sstop;
    double gs_sstep;
    int gs_max_iterations;
    double gs_precision;
    
    bool is_coord; // use chromosome coordinates

    bool limit_sel_space; // use full search space for s

    int traj_function = 0; // set trajectory function

    // HMM chain window size
    string win_unit;
    double win_morgan;
    double win_percent;
    // int win_bp;

    // stochastic trajectory
    bool use_stochastic; // use stochastic trajectory function
    int stochastic_reps; // how many repeats

    /// read relevant information
    void read_cmd_line ( int argc, char *argv[] ) ;

} ;

#endif

