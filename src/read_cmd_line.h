#ifndef __READ_CMD_LINE_H
#define __READ_CMD_LINE_H

void cmd_line::read_cmd_line ( int argc, char *argv[] ) {
    
    ///defaults
    ancestral_fixed = false ;                  /// set to true for qtl or experimental evolution application if ancestral genotypes are known and at fixed frequencies.
    
    /// ideally we recommend pruning LD in advance
    minimum_distance = 0 ;                          /// minimum distance in morgans between sites to consider
    ne = 2e4 ;                                      /// actually 2ne
    
    // time params to bound our search
    t_max = 10000 ;
    t_min = 1 ;
    p_max = 0.99999 ;
    p_min = 0.00001 ;
    t_length = 0.8 ;
    p_length = 0.8 ;
    
    /// if set, we clear once
    bool clear = false ;
    
    /// error rates
    error_rates = false ;
    
    // the default behavior is a single pulse of ancestry 1 into ancestry 0
    ancestry_pulses.resize( 2 ) ;
    ancestry_pulses[0].type = 0 ;
    ancestry_pulses[1].type = 1 ;
    
    /// default is 50:50 with single pulse of 1 into 0
    ancestry_pulses[0].proportion = 0.5 ;
    ancestry_pulses[1].proportion = 0.5 ;
    ancestry_pulses[0].proportion_fixed = true ;
    ancestry_pulses[1].proportion_fixed = true ;
    
    /// also the ancestry proportions are known
    ancestry_proportion.assign(2,0.5) ;
    
    /// time is not fixed by default, pulse of 1 into 0
    /// does not matter, really since 0>1 would be identical in formulation
    ancestry_pulses[0].time = 3000 ;
    ancestry_pulses[0].time_fixed = true ;
    ancestry_pulses[1].time = 10 ;
    ancestry_pulses[1].time_fixed = false ;

    /// end parameter this will be in lnl units
    /// i.e. must obtain <= this amount of improvement between all vertices to quit
    tolerance = 1e-5 ;
    
    /// restart number
    n_restarts = -1 ;

    /// per site per read error rate
    error_rate = 0.01 ;
    
    /// genotype data rather than read data?
    genotype = false ;

    // viterbi
    viterbi = false ;
    
    /// output pulses rather than ancestry counts
    output_pulses = true ;
    
    /// set output precision
    precision = 10 ;
    
    /// sample file
    sample_file = "null" ;
    
    // intput file
    input_file = "null" ;
    
    /// bootstraps
    n_bootstraps = 0 ;
    block_size = 0 ;
    
	/// accept command line parameters
	for (int i=1; i<argc; i++) {
        
        /// set this first, unless defaults are used
        if ( strcmp(argv[i],"-p") == 0 ) {
            
            /// if we deviate from default, clear the vector once
            if ( clear == false ) {
                ancestry_pulses.clear() ;
                clear = true ;
            }

            /// create new pulse, and record the parameters
            pulse new_ancestry_pulse ;
            new_ancestry_pulse.type = atoi(argv[++i]) ;
            new_ancestry_pulse.time = atof(argv[++i]) ;
            new_ancestry_pulse.proportion = atof(argv[++i]) ;
            
            // if time is set, we are not estimating it
            ///// set time with a negative number to provide the starting guess for this parameter
            if ( new_ancestry_pulse.time > 0 ) {
                new_ancestry_pulse.time_fixed = true ;
            }
            else {
                new_ancestry_pulse.time = new_ancestry_pulse.time * -1 ;
                new_ancestry_pulse.time_fixed = false ;
            }
            
            // if proporion is set, we are not estimating it
            ////// set proporiton with a negative number to provide the starting guess for this parameter
            if ( new_ancestry_pulse.proportion > 0 ) {
                new_ancestry_pulse.proportion_fixed = true ;
            }
            else {
                new_ancestry_pulse.proportion_fixed = false ;
                new_ancestry_pulse.proportion = -1 * new_ancestry_pulse.proportion ;
            }
            ancestry_pulses.push_back( new_ancestry_pulse ) ;
            ancestry_pulses.back().entry_order = ancestry_pulses.size() - 1 ;
        }
    
        //// for each ancestry type, set the total ancestry fraction
        //// this must be set and equal to all the ancestry types listed above
        if ( strcmp(argv[i],"-a") == 0 ) {
            ancestry_proportion.clear() ;
            int stop = atoi(argv[++i]) ;
            float sum = 0 ;
            for ( int l = 0 ; l < stop ; l ++ ) {
                ancestry_proportion.push_back( atof(argv[++i]) ) ;
                sum += ancestry_proportion.back() ;
            }
            
            //// check that ancestry proportions sum to one
            if ( sum < 0.9999 || sum > 1.0001 ) {
                cerr << "\n\n\t\t ERROR: ancestry proportions must sum to one\n\n" ;
                print_usage() ;
                exit(1) ;
            }
        }
        
        if ( strcmp(argv[i],"--help") == 0 ) {
            print_usage() ;
            exit(1) ; 
        }
        
        if ( strcmp(argv[i],"-g") == 0 ) {
            genotype = true ;
        }
        
        if ( strcmp(argv[i],"--output-ancestry") == 0 ) {
            output_pulses = false ;
        }
        if ( strcmp(argv[i],"--precision") == 0 ) {
            precision = atoi(argv[++i]) ;
            cout.precision(precision) ;
            cerr.precision(precision) ;
        }
        
        if ( strcmp(argv[i],"-v") == 0 ) {
            viterbi = true ;
        }
        
        if ( strcmp(argv[i],"-r") == 0 ) {
            n_restarts = atoi(argv[++i]) ;
        }
        
        /// bootstraps supplied as '-b <int, number> <int, block size>
        if ( strcmp(argv[i],"-b") == 0 ) {
            n_bootstraps = atoi(argv[++i]) ;
            block_size = atoi(argv[++i]) ;
        }
        
        /// to bound possible pulse times
        if ( strcmp(argv[i],"--tmax") == 0 ) {
            t_max = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--tmin") == 0 ) {
            t_min = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--pmin") == 0 ) {
            p_min = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--pmax") == 0 ) {
            p_max = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--tlength") == 0 ) {
            t_length = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--plength") == 0 ) {
            p_length = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--tolerance") == 0 ) {
            tolerance = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i], "-e" ) == 0 ) {
            error_rate = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i], "-E" ) == 0 ) {
            error_rates = true ; 
        }
        if ( strcmp(argv[i],"--ne") == 0 ) {
            ne = 2 * atof(argv[++i]) ;
        }

        /// this version will allow inputting all samples in a single file with separate posterior output files
        if ( strcmp(argv[i],"-i") == 0 ) {
            input_file = string(argv[++i]) ;
        }
        
        /// sample file 
        if ( strcmp(argv[i],"-s") == 0 ) {
            sample_file = string(argv[++i]) ;
        }
        
        if ( strcmp(argv[i],"-d") == 0 ) {
            minimum_distance = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--fix") == 0 ) {
            ancestral_fixed = true ;
        }
    }
    
    if ( input_file == "null" ) {
        cerr << "\n\n\t\tERROR: must provide input file\n\n\t\t\t-i [path/to/input_file]\n\n" ;
        print_usage() ;
        exit(1) ;
    }
    if ( sample_file == "null" ) {
        cerr << "\n\n\t\tERROR: must provide sample file\n\n\t\t\t-s [path/to/sample_file]\n\n" ;
        print_usage() ;
        exit(1) ;
    }
        
    return ;
}

#endif

