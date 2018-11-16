#ifndef __READ_EMISSIONS_H
#define __READ_EMISSIONS_H

/// this will calculate the emissions probabilities using read based methods
void create_emissions_matrix( double n, input_line &new_line, bool &ancestral_fixed, vector<vector<int> > &states, int sample_index, vector<pulse> &pulses, vec &emission_matrix ) {
    
    /// if no data in sample, return flat probs
    if ( new_line.sample_counts[sample_index][2] == 0 ) {
        emission_matrix.ones(states.size()) ;
        return ;
    }
    
    /// indexes will correspond to states vector that we passed here
    emission_matrix.zeros(states.size()) ;
    
    /// i is the index of the state we're in
    for ( double i = 0 ; i < states.size() ; i ++ ) {
        
        /// ancestry states
        vector<int> ancestry_states = pulses_to_ancestry( states[i], pulses ) ;
        
        /// find all distributions of reads across existing states using multichoose
        /// this data object is the number of A's derived from each ancestry type
        vector<vector<int> > read_counts ;
        distribute_reads( ancestry_states, new_line.sample_counts[sample_index][2], read_counts ) ;
        
        /// now the probability that you sample reads from each ancestry type
        for ( int rc = 0 ; rc < read_counts.size() ; rc ++ ) {
            
            // the probability of this particular sampling of reads given the state i
            double p_reads = multinomial( n, new_line.sample_counts[sample_index][2], read_counts[rc], ancestry_states ) ;
                                    
            /// now we need to distribute our alleles
            /// this data object is the number of A reads from each class, conditional on the number of reads from each class, i.e. A_counts[i][j] < read_counts[r][j]
            map<vector<double>, double > A_counts ;
            distribute_alleles( read_counts[rc], new_line.sample_counts[sample_index][0], new_line.sample_counts[sample_index][2], A_counts ) ;
            
            /// now compute probability of each sampling arrangement
            for ( std::map<vector<double>,double>::iterator a = A_counts.begin() ; a != A_counts.end() ; ++ a ) {
                
                double prob_counts = 1 ;
                for ( int c = 0 ; c < a->first.size() ; c ++ ) {
                    double sum = 0 ;
                    /// whether ancestral frequencies are not fixed [default]
                    if ( ancestral_fixed == false ) {
                        /// number of samples that are allele A is j
                        for ( double j = 0 ; j <= ancestry_states[c] + 0.01 ; j ++ ) {
                            double pA = j/ancestry_states[c]*(1-new_line.error_1)+(1-j/ancestry_states[c])*new_line.error_2 ;
                            
                            /// ancestral pop probs + choose j probs ( all genotype probs )
                            sum += 1/( new_line.reference_counts[c][2] + ancestry_states[c] + 1 )
                            * 1/nCk[ new_line.reference_counts[c][2] + ancestry_states[c] ][ new_line.reference_counts[c][0] + j ]
                            
                            /// remaining sample genotype prob
                            * nCk[ ancestry_states[c] ][ j ]
                            
                            /// read probs now, i.e. prob of number of A's given the sampling probs
                            * pow( pA, a->first[c] )
                            * pow( 1-pA, read_counts[rc][c] - a->first[c] ) ;
                        }
                    }
                    
                    /// if ancestral frequencies are fixed
                    else {
                        double f = new_line.reference_counts[c][0] / new_line.reference_counts[c][2] ;
                        for ( double j = 0 ; j <= ancestry_states[c] + 0.01 ; j ++ ) {
                            /// error probs are same
                            double pA = j/ancestry_states[c]*(1-new_line.error_1)+(1-j/ancestry_states[c])*new_line.error_2 ;
                            
                            /// no ancestral probs required
                            
                            /// prob of sampling these A's
                            sum += nCk[ ancestry_states[c] ][ j ]
                            * pow( f, j ) * pow ( 1 - f, ancestry_states[c] - j )
                            
                            //// prob of smapling these reads
                            * pow( pA, a->first[c] ) * pow( 1-pA, read_counts[rc][c] - a->first[c] ) ;
                        }
                    }
                    prob_counts *= sum ;
                }
                emission_matrix(i) += p_reads * prob_counts * a->second ;
            }
        }
    }
    return ;
}
#endif

