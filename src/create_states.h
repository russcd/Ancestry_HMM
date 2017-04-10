#ifndef __CREATE_STATES_H
#define __CREATE_STATES_H

/// create the set of states that are permissable for a given ploidy
void create_initial_states ( double &number_chromosomes, vector<pulse> &ancestry_pulses, map<int,vector<vector<int> > > &state_list ) {
    
    /// if we have already computed the state list for this sample ploidy, move on
    if ( state_list.find( number_chromosomes ) != state_list.end() ) {
        return ;
    }
    
    /// list of all possible states for single chromosomes
    vector<int> states ;
    for ( int p = 0 ; p < ancestry_pulses.size() ; p ++ ) {
        states.push_back( p ) ;
    }
    
    /// now get all possible arrangements and store them in our state list
    vector< vector<int> > results = multichoose( number_chromosomes, states ) ;
    for ( int i = 0 ; i < results.size() ; i ++ ) {
        vector<int> state_vector( ancestry_pulses.size(), 0 ) ;
        for ( int j = 0 ; j < results[i].size() ; j ++ ) {
            state_vector.at(results[i][j]) ++ ;
        }
        state_list[number_chromosomes].push_back( state_vector ) ;
    }
}

#endif
