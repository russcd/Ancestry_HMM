#ifndef __GOLDEN_SEARCH_H
#define __GOLDEN_SEARCH_H

/// golden section search for single parameter optimization
/// not currently included, but may be useful to legacy versions 
vector<pulse> golden_search ( cmd_line &options, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, map<int,vector<vector<int> > > &state_changes ) {
    
    /// now do golden search until we reach tolerance threshhold and stop
    double phi = ( sqrt(5) - 1 ) / 2 ;
    
    /// parameter values to hold during search
    double low_bracket ;
    double high_bracket ;
    double param_low ;
    double param_high ;
    
    /// now figure out the parameters that vary
    for ( int p = 0 ; p < options.ancestry_pulses.size() ; p ++ ) {
        if ( options.ancestry_pulses[p].time_fixed == false ) {
            high_bracket = options.t_max ;
            low_bracket = options.t_min ;
            param_low = options.t_max + phi * ( options.t_min - options.t_max ) ;
            param_high = options.t_min + phi * ( options.t_max - options.t_min ) ;

        }
        else if ( options.ancestry_pulses[p].proportion_fixed == false ) {
            high_bracket = options.p_max ;
            low_bracket = options.p_min ;
            param_low = options.p_max + phi * ( options.p_min - options.p_max ) ;
            param_high = options.p_min + phi * ( options.p_max - options.p_min ) ;
        }
    }

    /// likelihood information for test points
    double lnl_low = 0 ;
    double lnl_high = 0 ;
    double lnl_diff = 1000 ;
    int iteration = 0 ;
    
    cerr << "\titeration\tlow bound\tlow test\thigh test\thigh bound\tlnl low\tlnl high\n" ;
    while ( options.tolerance < lnl_diff ) {
        
        /// compute probabilty of low point
        if ( lnl_low == 0 ) {
            vector<pulse> v_low = options.ancestry_pulses ;
            for ( int p = 0 ; p < v_low.size() ; p ++ ) {
                if ( options.ancestry_pulses[p].time_fixed == false ) {
                    v_low[p].time = param_low ;
                }
                else if ( options.ancestry_pulses[p].proportion_fixed == false ) {
                    v_low[p].fraction_of_remainder = param_low ;
                }
            }
            
            lnl_low += evaluate_vertex( v_low, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes ) ;
        }
                
        /// compute probability of high point
        if ( lnl_high == 0 ) {
            vector<pulse> v_high = options.ancestry_pulses ;
            for ( int p = 0 ; p < v_high.size() ; p ++ ) {
                if ( options.ancestry_pulses[p].time_fixed == false ) {
                    v_high[p].time = param_high ;
                }
                else if ( options.ancestry_pulses[p].proportion_fixed == false ) {
                    v_high[p].fraction_of_remainder = param_high ;
                }
            }
            lnl_high = evaluate_vertex( v_high, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes ) ;
        }
        
        /// print update
        cerr << "\t" << iteration << "\t" << low_bracket << "\t" << param_low << "\t" << param_high << "\t" << high_bracket << "\t" << lnl_low << "\t" << lnl_high << endl ;
        iteration ++ ;
        
        /// record dfiference
        lnl_diff = abs( lnl_low - lnl_high ) ;
        
        /// if true, we know that the maximum is between param_low and high_bracket
        if ( lnl_high >= lnl_low ) {
            low_bracket = param_low ;
            param_low = param_high ;
            param_high = low_bracket + ( high_bracket - low_bracket ) * phi ;
            lnl_low = lnl_high ;
            lnl_high = 0 ;
        }
        
        /// otherwise, the maximum is between low_bracket and param_high
        else {
            high_bracket = param_high ;
            param_high = param_low ;
            param_low = high_bracket + ( low_bracket - high_bracket ) * phi ;
            lnl_high = lnl_low ;
            lnl_low = 0 ;
        }
    }
    
    vector<pulse> optimum = options.ancestry_pulses ;
    for ( int p = 0 ; p < optimum.size() ; p ++ ) {
        if ( options.ancestry_pulses[p].time_fixed == false ) {
            optimum[p].time = ( param_high + param_low ) / 2 ;
        }
        else if ( options.ancestry_pulses[p].proportion_fixed == false ) {
            optimum[p].fraction_of_remainder = ( param_high + param_low ) / 2 ;
        }
    }
    
    return optimum ;
}


#endif

