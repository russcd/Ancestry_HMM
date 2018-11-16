#ifndef __BOOTSTRAP_H
#define __BOOTSTRAP_H

/// create_bootstraps
vector<vector<pulse> > bootstraps (vector<vector<pulse> > &vertices, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, cmd_line &options, map<int,vector<vector<int> > > &state_list, vector<string> &chromosomes ) {
    
    /// vector to store optimum admixture model for each bootstrap
    vector<vector<pulse> > bootstrap_models ;
    
    /// rng to select block positions
    int max = floor(recombination_rate.size()/ options.block_size) ;

    /// create b total bootstraps
    for ( int b = 0 ; b < options.n_bootstraps ; b ++ ) {
        
        //// bootstrap data objects
        vector<markov_chain> bootstrap_chain ( markov_chain_information.size() ) ;
        vector<double> bootstrap_recombination_rate ;
        vector<int> bootstrap_positions ;
        
        /// populate these data objects
        while ( bootstrap_positions.size() < position.size() ) {
            
            /// find admissable start positions
            int start = ( rand() % max ) * options.block_size ;
            int end = start + options.block_size ;
            if ( end > position.size() ) {
                end = position.size() - 1 ;
            }
            
            /// create positions in bootstrap region
            bootstrap_positions.insert( bootstrap_positions.end(), position.begin() + start , position.begin() + end ) ;
            
            //// create recombination rate in bootstrap regions
            bootstrap_recombination_rate.push_back( 0.5 ) ;
            bootstrap_recombination_rate.insert( bootstrap_recombination_rate.end(), recombination_rate.begin() + start + 1 , recombination_rate.begin() + end ) ;
            
            //// create markov chain information for each sample
            for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
                
                /// add emissions
                bootstrap_chain[m].emission_probabilities.insert( bootstrap_chain[m].emission_probabilities.end() , markov_chain_information[m].emission_probabilities.begin() + start ,markov_chain_information[m].emission_probabilities.begin() + end ) ;

                /// start probs
                bootstrap_chain[m].start_prob = 1 ;
                bootstrap_chain[m].end_prob = 1 ;
                
                /// ploidy transitions
                bootstrap_chain[m].ploidy_switch_position.push_back( bootstrap_chain[m].emission_probabilities.size() - options.block_size ) ;
                bool start_found = false ;
                for ( int p = 0 ; p < markov_chain_information[m].ploidy_switch_position.size() ; p ++ ) {
                    
                    /// find and insert start ploidy
                    if ( markov_chain_information[m].ploidy_switch_position[p] > start && start_found == false ) {
                        bootstrap_chain[m].ploidy_switch.push_back( markov_chain_information[m].ploidy_switch[p-1] ) ;
                        start_found = true ;
                    }
                    
                    /// insert all ploidy changes across that block
                    if ( start_found == true && markov_chain_information[m].ploidy_switch_position[p] <= end && markov_chain_information[m].ploidy_switch_position[p] >= start ) {
                        bootstrap_chain[m].ploidy_switch.push_back( markov_chain_information[m].ploidy_switch[p] ) ;
                        bootstrap_chain[m].ploidy_switch_position.push_back( bootstrap_chain[m].emission_probabilities.size() - options.block_size + markov_chain_information[m].ploidy_switch_position[p] - start ) ;
                    }
                }
            }
        }
        
        //// to avoid look ahead errors
        for ( int m = 0 ; m < bootstrap_chain.size() ; m ++ ) {
            bootstrap_chain[m].ploidy_switch_position.push_back(2147483647) ;
            bootstrap_chain[m].ploidy_switch.push_back( bootstrap_chain[m].ploidy_switch.back() ) ;
        }
        
        /// if there are params to estimate, do amoeba search
        if ( vertices.size() > 2 ) {
            cerr << "starting nelder-mead search\t\n\tbootstrap no:\t" << b << endl ;
            bootstrap_models.push_back( nelder_mead_search( vertices, options, bootstrap_chain, transition_matrix_information, bootstrap_recombination_rate, bootstrap_positions, state_list ) ) ;
        }
        
        /// or do golden section line search for single parameter optimization
        else {
            cerr << "starting golden section search\t\n\tbootstrap no:\t" << b << endl << endl ;
            bootstrap_models.push_back( golden_search( options, bootstrap_chain, transition_matrix_information, bootstrap_recombination_rate, bootstrap_positions, state_list ) ) ;
        }
    }
    
    return bootstrap_models ;
}

#endif
