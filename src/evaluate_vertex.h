#ifndef __EVALUATE_VERTEX_H
#define __EVALUATE_VERTEX_H

/// evaluate likelihood for a single vertex
double evaluate_vertex( vector<pulse> &vertex, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, cmd_line &options, map<int,vector<vector<int> > > &state_changes ) {
    
    /// create single chromosome transition matrix for single chromosomes
    mat transition_rates = create_transition_rates( vertex, options.ne, options.ancestry_proportion ) ;

    /// create transition matrix for all ploidies
    map<int,vector<mat> > transition_matrix ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        create_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information.at(m).number_chromosomes], recombination_rate, position, markov_chain_information.at(m).number_chromosomes, transition_rates ) ;
        for ( int p = 0 ; p < markov_chain_information[m].ploidy_switch.size() ; p ++ ) {
            create_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information[m].ploidy_switch[p]], recombination_rate, position, markov_chain_information[m].ploidy_switch[p], transition_rates ) ;
        }
    }
    
    /// compute transitions within a state
    vector<mat> interploidy_transitions ;
    interploidy_transitions = create_interploidy_transitions( state_changes, vertex, options.ancestry_proportion ) ;
    
    /// now compute the forward probabilities
    double lnl = 0 ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        lnl += markov_chain_information[m].compute_forward_probabilities( transition_matrix, interploidy_transitions ) ;
    }
    return lnl ;
}

#endif

