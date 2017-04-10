#ifndef __COMPUTE_BACKWARD_H
#define __COMPUTE_BACKWARD_H

/// backward probabilities are same except transpose transition matrix
void markov_chain::compute_backward_probabilities( map<int, vector<mat> > &transition_probabilites, vector<mat> &interploidy_transitions ) {

    // resize betas matrix
    betas.resize(transition_probabilites[number_chromosomes].size()) ;

    /// index to tract position in ploidy path will iterate through backwards for backward probs
    int ploidy_index = ploidy_switch_position.size() - 2 ;
        
    /// set last states to one
    betas.back().resize(alphas.back().size()) ;
    
    /// start with last position and multiply
    betas.back() = emission_probabilities.back() * end_prob ;
    normalize(betas.back()) ;
    
    /// now go from t-1 to 0 to iterate through backwards
    for ( int i = emission_probabilities.size()-2 ; i > -1 ; i -- ) {
                        
        /// if we're at or before the last switch position
        bool ploidy_change = false ;
        if ( i < ploidy_switch_position[ploidy_index] ) {
            ploidy_index -- ;
            if ( ploidy_switch[ploidy_index] != ploidy_switch[ploidy_index+1] ) {
                ploidy_change = true ;
            }
        }
        
        /// resize vector
        betas[i].zeros( transition_probabilites[ploidy_switch[ploidy_index]][1].n_cols ) ;
        
        //// if there was a transition in ploidy
        if ( ploidy_change == true ) {
            
            /// if self transition is this low, we switched across a chromosome boundary or marker densities are so low that LAI is pointless, add a check for this
            if ( transition_probabilites[ploidy_switch[ploidy_index]][i+1](0,0) < 0.75 ) {
                /// new chromosome just starts with flat probs, will normalize to one
                betas[i].fill( 1 ) ;
            }
            //// otherwise, this is a transition across ploidy types on the same chromosome use the interploidy transition rates
            else {
                betas[i] = interploidy_transitions[ploidy_switch[ploidy_index]-1].t() * betas[i+1] % emission_probabilities[i] ;
            }
        }
        
        /// multiply matrices to produce latest betas
        else {
            betas[i] = transition_probabilites[ploidy_switch[ploidy_index]][i+1].t() * betas[i+1] % emission_probabilities[i] ;
        }
        
        /// normalize vector to prevent underflow issues
        normalize(betas[i]) ;
    }
}

#endif

