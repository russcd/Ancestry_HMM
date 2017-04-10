#ifndef __VITERBI_H
#define __VITERBI_H

void markov_chain::viterbi( vector<int> &position, vector<double> &recombination_rate, map<int, vector<vector<int> > > &states, vector<string> &chrom, map<int, vector<mat> > &transition_probabilites, vector<mat> &interploidy_transitions, bool output_pulses, vector<pulse> &pulses ) {
    
    //// find largest states group that's possible for this sample
    int max_n_states = 2 ;
    for ( std::map<int,vector<vector<int> > >::iterator p = states.begin() ; p != states.end() ; ++ p ) {
        if ( p->second.size() > max_n_states ) {
            max_n_states = p->second.size() ;
        }
    }
    
    //// probs and state paths as we transition through
    vector<vector<int> > viterbi_path(max_n_states) ;
    vector<double> path_probs(max_n_states,0) ;
    
    //// starting probs
    for ( int s = 0 ; s < emission_probabilities[0].size() ; s ++ ) {

        path_probs[s] += log( emission_probabilities[0](s) ) + log( start_prob );
        viterbi_path[s].push_back(s) ;
    }
    
    /// ploidy index to tract where in path we are
    int ploidy_index = 0 ;
    
    /// now go across all incoming paths for each state, and select highest prob path
    for ( int i = 1 ; i < emission_probabilities.size() ; i ++ ) {
    
        /// if we're at or past the next switch position
        bool ploidy_change = false ;
        if ( i >= ploidy_switch_position[ploidy_index+1] ) {
            ploidy_index ++ ;
            if ( ploidy_switch[ploidy_index] != ploidy_switch[ploidy_index-1] ) {
                ploidy_change = true ;
            }
        }
        
        /// to record best paths and their probabilities for each state
        vector<int> max_paths(emission_probabilities[i].size()) ;
        vector<double> max_path_probs(emission_probabilities[i].size(), -1.7976931348623157E+308) ;
        
        /// if the ploidy between adjacent states is identical, we can do a fairly normal transition
        if ( ploidy_change == false ) {
            for ( int s = 0 ; s < emission_probabilities[i].size() ; s ++ ) {
                for ( int p = 0 ; p < emission_probabilities[i].size() ; p ++ ) {
                    /// find all ways of getting to the current state
                    double path_prob_state = path_probs[p] + log( transition_probabilites[ploidy_switch[ploidy_index]][i](s,p) ) + log( emission_probabilities[i](s) ) ;
                    /// check if this path prob is greater and update to that
                    if ( path_prob_state > max_path_probs[s] ) {
                        max_paths[s] = p ;
                        max_path_probs[s] = path_prob_state ;
                    }
                }
            }
        }
        
        /// if we have transitioned between adjancent ploidy tracts
        else {
            for ( int s = 0 ; s < emission_probabilities[i].size() ; s ++ ) {
                for ( int p = 0 ; p < emission_probabilities[i-1].size() ; p ++ ) {
                    
                    /// do normal interploidy transition
                    double path_prob_state = path_probs[p] + log( interploidy_transitions[ploidy_switch[ploidy_index-1]-1](s,p) ) + log( emission_probabilities[i](s) ) ;
                    
                    /// if this is a between chromosomes transition all states are equally likely
                    if ( transition_probabilites[ploidy_switch[ploidy_index]][i](0,0) < 0.75 ) {
                        path_prob_state = path_probs[p] + log( emission_probabilities[i](s) ) ;
                    }
                    
                    if ( path_prob_state > max_path_probs[s] ) {
                        max_paths[s] = p ;
                        max_path_probs[s] = path_prob_state ;
                    }
                }
            }
        }
        
        /// update recorded paths and probabilities
        vector<vector<int> > new_viterbi_paths ( emission_probabilities[i].size() ) ;
        vector<double> new_path_probs ( emission_probabilities[i].size() ) ;
        for ( int s = 0 ; s < emission_probabilities[i].size() ; s ++ ) {
            new_viterbi_paths[s] = viterbi_path[max_paths[s]] ;
            new_viterbi_paths[s].push_back(s) ;
            new_path_probs[s] = max_path_probs[s] ;
        }
        
        /// now replace viterbi paths/probs
        swap( viterbi_path, new_viterbi_paths ) ;
        swap( path_probs, new_path_probs ) ;
    }

    /// find optimal viterbi path
    int optimal_path = 0 ;
    for ( int p = 0 ; p < path_probs.size() ; p ++ ) {
        if ( path_probs[p] > path_probs[optimal_path] ) {
            optimal_path = p ;
        }
    }

    /// print optimal viterbi path for all states
    ofstream out ( output_file.c_str() ) ;
    string current_chrom = chrom[0] ;
    int start = 0 ;
    ploidy_index = 0 ;
    int current_state = viterbi_path[optimal_path][0] ;
    
    /// to keep tract of distance in morgans
    double recombination_distance = 0 ;
    double recombination_start = 0 ;
    
    /// iterate across paths
    for ( int p = 0 ; p < chrom.size() ; p ++ ) {
            
        /// if we're at or past the next switch position
        bool ploidy_change = false ;
        if ( p >= ploidy_switch_position[ploidy_index+1] ) {
            ploidy_index ++ ;
            if ( ploidy_switch[ploidy_index] != ploidy_switch[ploidy_index-1] ) {
                ploidy_change = true ;
            }
        }
        
        //// no ploidy change but switch in chromosome
        if ( ploidy_change == false && ( current_chrom != chrom[p] || viterbi_path[optimal_path][p] != current_state ) ) {
            if ( current_chrom != chrom[p] ) {
                out << current_chrom << "\t" << start << "\t" << position[p-1] << "\t" << recombination_start << "\t" << recombination_distance << "\t" ;
                start = 0 ;
                recombination_start = 0 ;
                recombination_distance = 0 ;
            }
            else {
                out << current_chrom << "\t" << start << "\t" << position[p]-1 << "\t" << recombination_start << "\t" << recombination_distance << "\t" ;
                start = position[p] ;
                recombination_start = recombination_distance ;
            }
            for ( int c = 0 ; c < states[ploidy_switch[ploidy_index]][current_state].size() - 1 ; c ++ ) {
                out << states[ploidy_switch[ploidy_index]][current_state][c] << "," ;
            }
            out << states[ploidy_switch[ploidy_index]][current_state].back() << endl ;

            current_state = viterbi_path[optimal_path][p] ;
            current_chrom = chrom[p] ;
        }
        
        else if ( ploidy_change == true ) {
            if ( current_chrom != chrom[p] ) {
                out << current_chrom << "\t" << start << "\t" << position[p-1] << "\t" << recombination_start << "\t" << recombination_distance << "\t" ;
                start = 0 ;
                recombination_start = 0 ;
                recombination_distance = 0 ;
            }
            else {
                out << current_chrom << "\t" << start << "\t" << position[p] - 1 << "\t" << recombination_start << "\t" << recombination_distance << "\t" ;
                start = position[p] ;
                recombination_start = recombination_distance ;
            }
            
            for ( int c = 0 ; c < states[ploidy_switch[ploidy_index-1]][current_state].size() - 1 ; c ++ ) {
                out << states[ploidy_switch[ploidy_index-1]][current_state][c] << "," ;
            }
            out << states[ploidy_switch[ploidy_index-1]][current_state].back() << endl ;
            current_chrom = chrom[p] ;
            current_state = viterbi_path[optimal_path][p] ;
        }
        
        /// add additional recombination
        if ( recombination_rate[p] < 0.4 ) {
            if ( p > 0 ) {
                recombination_distance += recombination_rate[p] * ( position[p] - position[p-1] ) ;
            }
            else {
                recombination_distance += recombination_rate[p] * position[p] ;
            }
        }
    }
    out << current_chrom << "\t" << start << "\t" << position.back() << "\t" << recombination_start << "\t" << recombination_distance << "\t" ;
    for ( int c = 0 ; c < states[ploidy_switch[ploidy_index]][current_state].size() - 1 ; c ++ ) {
        out << states[ploidy_switch[ploidy_index]][current_state][c] << "," ;
    }
    out << states[ploidy_switch[ploidy_index]][current_state].back() << endl ;
}

#endif
