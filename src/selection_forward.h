#ifndef __SELECTION_FORWARD_H
#define __SELECTION_FORWARD_H



// Forward algoritm modified for selection inferrence.
double markov_chain::selection_forward_probabilities_genotypes( map<int, vector<mat> > &transition_probabilites, vector<mat> &interploidy_transitions, selection &point, bool go_downstream, vector<double> &genofreq, vector<int> &position ) {
    //cerr << "cp2_1 " << genofreq[0] << endl;

    /// return log likelihood which is sum of cts
    double lnl = 0 ;
    
    /// clear the fw probs matrix
    alphas.resize( transition_probabilites[ploidy_switch[0]].size() ) ;

    /// ploidy index to tract where in path we are
    int ploidy_index = 0 ;
    
    //// set all values to zero, but mostly just reize
    alphas[0].resize( transition_probabilites[ploidy_switch[0]][1].n_cols ) ;
    
    /// get initial state set
    //alphas[0] = emission_probabilities[point.pos] * start_prob ;

    // Populate starting conditions

    // Check how to specify nn. The current way is a bit of a hack.
    double nn = transition_probabilites[ploidy_switch[0]][1].n_cols - 1;
    for (int k = nn; k >= 0; k--) {
        alphas[0][nn-k] = binomial(nn, k, genofreq[0]);
    }

    lnl += normalize( alphas[0] ) ;
    
    /// do all other sites
    /// Checks if going upstream or downstream from the selected site.
    if (go_downstream == true) {
        selection_forward_loop_reverse(transition_probabilites, interploidy_transitions, point, lnl, ploidy_index, position) ;
    }
    else {
        selection_forward_loop(transition_probabilites, interploidy_transitions, point, lnl, ploidy_index, position) ;
    }
    
    return lnl ;
}

// Loop in forward algorithm going downstream from the selected site
void markov_chain::selection_forward_loop( map<int, vector<mat> > &transition_probabilites, vector<mat> &interploidy_transitions, selection &point, double &lnl, int ploidy_index, vector<int> &position) {
    // WARNING: Check what +1 index does. May be unnecessary.
    //cerr << "Emission probabilities: " << emission_probabilities.size() << " " << point.pos << endl;

    int j = 1;
    int k;
    double normalpha;

    for ( int i = 1 ; i < transition_probabilites[ploidy_switch[ploidy_index]].size() ; i ++ ) {
        k = point.pos + i;
        /// if we're at or past the next switch position
        bool ploidy_change = false ;
        if ( i >= ploidy_switch_position[ploidy_index+1] ) {
            ploidy_index ++ ;
            if ( ploidy_switch[ploidy_index] != ploidy_switch[ploidy_index-1] ) {
                ploidy_change = true ;
            }
        }
        /// resize matrix
        alphas[j].resize( transition_probabilites[ploidy_switch[ploidy_index]][1].n_cols ) ;
        
        /// requires slightly different math if we are transitioning in ploidy between two adjacent sites
        if ( ploidy_change == true ) {
            /// transitions across a chromosome boundary will have low self-self rates
            if ( transition_probabilites[ploidy_switch[ploidy_index]][i](0,0) < 0.75 ) {
                alphas[j].fill( 1 ) ;
            }
            
            //// otherwise, this is a transition across ploidy types on the same chromosome use the interploidy transition rates
            else {
                alphas[j] = interploidy_transitions[ploidy_switch[ploidy_index-1]-1] * alphas[j-1] % emission_probabilities[k] ;
            }
        }
        
        /// otehrwise business as ususal

        else {
            alphas[j] = transition_probabilites[ploidy_switch[ploidy_index]][j] * alphas[j-1] % emission_probabilities[k] ;
            normalpha = normalize( alphas[j] ) ;
        }

        /// normalize and updated likelihood
        lnl += normalpha ;
        j++;
    }
}

// Loop in forward algorithm going upstream from the selected site
void markov_chain::selection_forward_loop_reverse( map<int, vector<mat> > &transition_probabilites, vector<mat> &interploidy_transitions, selection &point, double &lnl, int ploidy_index, vector<int> &position) {
    // WARNING: Check what +1 index does. May be unnecessary.

    int j = 1;
    int k;
    double normalpha;

    for ( int i = 0 ; i < transition_probabilites[ploidy_switch[ploidy_index]].size()-1 ; i ++ ) {
        k = point.pos - i;
        /// if we're at or past the next switch position
        bool ploidy_change = false ;
        if ( i >= ploidy_switch_position[ploidy_index+1] ) {
            ploidy_index ++ ;
            if ( ploidy_switch[ploidy_index] != ploidy_switch[ploidy_index-1] ) {
                ploidy_change = true ;
            }
        }
        /// resize matrix
        alphas[j].resize( transition_probabilites[ploidy_switch[ploidy_index]][1].n_cols ) ;
        
        /// requires slightly different math if we are transitioning in ploidy between two adjacent sites
        if ( ploidy_change == true ) {
            /// transitions across a chromosome boundary will have low self-self rates
            if ( transition_probabilites[ploidy_switch[ploidy_index]][i](0,0) < 0.75 ) {
                alphas[j].fill( 1 ) ;
            }
            
            //// otherwise, this is a transition across ploidy types on the same chromosome use the interploidy transition rates
            else {
                alphas[j] = interploidy_transitions[ploidy_switch[ploidy_index-1]-1] * alphas[j-1] % emission_probabilities[k] ;
            }
        }
        
        /// otehrwise business as ususal
        else {
            alphas[j] = transition_probabilites[ploidy_switch[ploidy_index]][j] * alphas[j-1] % emission_probabilities[k] ;
            normalpha = normalize( alphas[j] ) ;
        }

        /// normalize and updated likelihood
        lnl += normalpha ;
        j++;
    }
}

#endif
