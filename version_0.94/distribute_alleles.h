#ifndef __DISTRIBUTE_ALLELES_H
#define __DISTRIBUTE_ALLELES_H

/// solution to computing all ways of distributing A reads among vector of read counts
/// using recursion to generate all possible lists of outcomes stored in nA
void compute_allele_counts( vector<vector<double> > &nA, vector<double> nA_iteration, vector<int> &read_counts, double A, double t, int pop ) {
   
    int min = 0 ;
    if ( t - read_counts[pop] < A ) {
        min = A - ( t - read_counts[pop] ) ;
    }
    int max = read_counts[pop] ;
    if ( A < read_counts[pop] ) {
        max = A ;
    }
    
    for ( int c = min ; c <= max ; c ++ ) {
        nA_iteration.push_back(c) ;
        if ( pop == read_counts.size() - 1 ) {
            nA.push_back( nA_iteration ) ;
        }
        else {
            compute_allele_counts( nA, nA_iteration, read_counts, A - c, t - read_counts[pop], pop + 1 ) ;
        }
        nA_iteration.pop_back() ;
    }
}

/// distribute alleles (if genotype space) across all possible states
/// read counts is the number of chromosomes from each pulse
void distribute_alleles ( vector<int> &read_counts, double &A, double &read_total, map <vector<double>, double > &A_counts ) {
    
    vector<double> nA_iteration ;
    vector<vector<double> > nA ;
    compute_allele_counts( nA, nA_iteration, read_counts, A, read_total, 0 ) ;
    
    /// now compute ways to get to those arrangements
    double total = 0 ;
    for ( int n = 0 ; n < nA.size() ; n ++ ) {
        double sum = 1 ;
        for ( int l = 0 ; l < nA[n].size() ; l ++ ) {
            sum *= nCk[read_counts[l]][nA[n][l]] ;
        }
        A_counts[nA[n]] = sum ;
        total += sum ;
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

