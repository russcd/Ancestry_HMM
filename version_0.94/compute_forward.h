#ifndef __COMPUTE_FORWARD_H
#define __COMPUTE_FORWARD_H

double markov_chain::compute_forward_probabilities( map<int, vector<mat> > &transition_probabilites, vector<mat> &interploidy_transitions ) {
        
    /// return log likelihood which is sum of cts
    double lnl = 0 ;
    
    /// clear the fw probs matrix
    alphas.resize( transition_probabilites[ploidy_switch[0]].size() ) ;

    /// ploidy index to tract where in path we are
    int ploidy_index = 0 ;
    
    //// set all values to zero, but mostly just reize
    alphas[0].resize( transition_probabilites[ploidy_switch[0]][1].n_cols ) ;
    
    /// get initial state set
    alphas[0] = emission_probabilities[0] * start_prob ;
    lnl += normalize( alphas[0] ) ;
    
    /// do all other sites
    for ( int i = 1 ; i < emission_probabilities.size() ; i ++ ) {
        
        /// if we're at or past the next switch position
        bool ploidy_change = false ;
        if ( i >= ploidy_switch_position[ploidy_index+1] ) {
            ploidy_index ++ ;
            if ( ploidy_switch[ploidy_index] != ploidy_switch[ploidy_index-1] ) {
                ploidy_change = true ;
            }
        }
        
        /// resize matrix
        alphas[i].resize( transition_probabilites[ploidy_switch[ploidy_index]][1].n_cols ) ;
        
        /// requires slightly different math if we are transitioning in ploidy between two adjacent sites
        if ( ploidy_change == true ) {
            
            /// transitions across a chromosome boundary will have low self-self rates
            if ( transition_probabilites[ploidy_switch[ploidy_index]][i](0,0) < 0.75 ) {
                alphas[i].fill( 1 ) ;
            }
            
            //// otherwise, this is a transition across ploidy types on the same chromosome use the interploidy transition rates
            else {
                alphas[i] = interploidy_transitions[ploidy_switch[ploidy_index-1]-1] * alphas[i-1] % emission_probabilities[i] ;
            }
        }
        
        /// otehrwise business as ususal
        else {
            alphas[i] = transition_probabilites[ploidy_switch[ploidy_index]][i] * alphas[i-1] % emission_probabilities[i] ;
        }

        /// normalize and updated likelihood
        lnl += normalize( alphas[i] ) ;
        
    }

    return lnl ;
}

#endif
