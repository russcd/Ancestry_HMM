#ifndef __INBRED_H
#define __INBRED_H

void fix_ibd_transitions( vector<vector< map< vector<transition_information>, double > > > &transition_matrix_information, vector<vector<int> > &states, vector<pulse> &ancestry_pulses, double inbreeding_transition_rate ) {
    
    /// number of ancestry pulses
    int number_pulses = ancestry_pulses.size() ;
    
    /// this is the first pulse
    for ( int i = 0 ; i < transition_matrix_information.size() ; i ++ ) {
        /// this is the pulse we're transitioning into
        for ( int j = 0 ; j < transition_matrix_information[i].size() ; j ++ ) {
            
            /// skip situations where we are comparing two outbred states
            if ( i < states.size() - number_pulses && j < states.size() - number_pulses ) {
                continue ;
            }
            
            /// skip situtions where we compare two inbred states
            if ( i >= states.size() - number_pulses && j >= states.size() - number_pulses ) {
                continue ;
            }
            
            /// okay then we need to fix it for ibd comparisons
            /// clear the existing transisions, they are wrong
            transition_matrix_information[i][j].clear() ;
            
            ////
            // ask how many are identical in state in inbred versus outbred? if 0,
            // then there are no ways to transitions between the states, we can record the
            // proportion of possible donor chromsomes and multiply by the rate of IBD
            vector<int> sum ( number_pulses, 0 ) ;
            for ( int p = 0 ; p < states[i].size() ; p ++ ) {
                sum[0] += states[i][p] ;
                sum[1] += states[j][p] ;
            }
            
            /// find and store valid transitions
            /// must have one chromsome of same type as ibd state
            /// then prob is related to the proportion of chromosomes of the same type
            for ( int p = 0 ; p < states[i].size() ; p ++ ) {
                if ( states[i][p] > 0 && states[j][p] > 0 ) {
                    vector<transition_information> new_transition (1) ;
                    new_transition[0].start_state = p ;
                    new_transition[0].end_state = p ;
                    new_transition[0].transition_count = 1 ;
                    new_transition[0].ibd_transition = true ; 
                    double prob_chrom = (double)states[j][p]/(double)sum[1]* (double)states[i][p]/(double)sum[0] ;
                    transition_matrix_information[i][j][new_transition] = inbreeding_transition_rate * prob_chrom ;
                }
            }
        }
    }
}

#endif




