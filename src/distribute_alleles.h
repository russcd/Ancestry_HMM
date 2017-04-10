#ifndef __DISTRIBUTE_ALLELES_H
#define __DISTRIBUTE_ALLELES_H

/// distribute alleles (if genotype space) across all possible states
/// read counts is the number of chromosomes from each pulse
void distribute_alleles ( vector<int> &read_counts, double &A, double &read_total, map <vector<double>, double > &A_counts ) {
    
    /// read origins
    vector<int> read_origins ;
    for ( int r = 0 ; r < read_counts.size() ; r ++ ) {
        for ( int l = 0 ; l < read_counts.at(r) ; l ++ ) {
            read_origins.push_back( r ) ;
        }
    }
    
    /// create vector of all reads 0 and 1
    vector<int> counts (A, 1) ;
    for ( int l = A ; l < read_total ; l ++ ) {
        counts.push_back( 0 ) ;
    }
    
    /// get all unique arrangements
    vector< vector<int> > results = multipermute( counts ) ;
    
    /// now distribute into class bins by alinging all results with the read origins
    for ( int i = 0 ; i < results.size() ; i ++ ) {
        vector<double> A_count( read_counts.size(), 0 ) ;
        for ( int r = 0 ; r < results.at(i).size() ; r ++ ) {
            if ( results[i][r] == 1 ) {
                A_count[read_origins[r]] ++ ;
            }
        }
        A_counts[A_count] ++ ;
    }
}

/// distribute reads (if pileup space) across all possible states
void distribute_reads ( vector<int> &state, double &read_total, vector<vector<int> > &read_counts ) {
    
    /// list of all possible states to be drawn from where we ignore states with zero chromosomes in that state
    vector<int> states ;
    for ( int i = 0 ; i < state.size() ; i ++ ) {
        if ( state[i] != 0 ) {
            states.push_back( i ) ;
        }
    }
    
    /// now get all possible arrangements of reads and store them in our list
    vector< vector<int> > results = multichoose( read_total, states ) ;
    for ( int i = 0 ; i < results.size() ; i ++ ) {
        vector<int> state_vector( state.size(), 0 ) ;
        for ( int j = 0 ; j < results[i].size() ; j ++ ) {
            state_vector.at(results[i][j]) ++ ;
        }
        read_counts.push_back( state_vector ) ;
    }
}

#endif

