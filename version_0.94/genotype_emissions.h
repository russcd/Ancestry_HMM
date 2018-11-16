#ifndef __GENOTYPE_EMISSIONS_H
#define __GENOTYPE_EMISSIONS_H

/// each inner vector will be length 4 with, A's that are really A's, a's that are really a's, A's that are really a's, and a's that are really A's
double genotype_total_prob( input_line &new_line, const double &A_count, int &chrom_count, double &real_A_count ) {
    
    if ( chrom_count == 0 ) {
        return 1 ;
    }
    
    /// vector of K balls to go into
    vector<int> observed_genotypes ;
    for ( int i = 0 ; i < A_count ; i ++ ) {
        observed_genotypes.push_back(1) ;
    }
    for ( int i = A_count ; i < chrom_count ; i ++ ) {
        observed_genotypes.push_back(0) ;
    }
    
    /// N bins corresponding to the real (unobserved) genotypes, these are not permuted
    vector<int> real_genotypes ;
    for ( int i = 0 ; i < real_A_count ; i ++ ) {
        real_genotypes.push_back(1) ;
    }
    for ( int i = real_A_count ; i < chrom_count ; i ++ ) {
        real_genotypes.push_back(0) ;
    }
    
    /// multipermute those observed genotypes across the chromosomes
    vector<vector<int> > observed_arrangements = multipermute( observed_genotypes ) ;
    
    /// record total prob of all possible genotype observed/real arrangements
    double total_prob = 0 ;
    
    /// count differences and sames
    for ( int o = 0 ; o < observed_arrangements.size() ; o ++ ) {
        /// record prob of this specific arrangement
        double prob = 1 ;
        for ( int s = 0 ; s < observed_arrangements[o].size() ; s ++ ) {
            if ( observed_arrangements[o][s] == 1 ) {
                if ( real_genotypes[s] == 1 ) {
                    prob *= ( 1 - new_line.error_1 ) ;
                }
                else {
                    prob *= new_line.error_1 ;
                }
            }
            else {
                if ( real_genotypes[s] == 1 ) {
                    prob *= new_line.error_2 ;
                }
                else {
                    prob *= ( 1 - new_line.error_2 ) ;
                }
            }
        }
        total_prob += prob ;
    }
    
    return total_prob ;
}

/// this will calculate the emissions probabilities using genotype based methods
void create_emissions_matrix_genotype( double n, input_line &new_line, bool &ancestral_fixed, vector<vector<int> > &states, int &sample_index, vector<pulse> &pulses, vec &emission_matrix ) {
        
    /// if no data in sample at site, return flat probs
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
        
        /// genotype sampling probabilities are 1 by definition, so no p_reads term
        /// instead just distribute A alleles across sampled chromosomes wihtout replacement
        map<vector<double>, double > A_counts ;
        distribute_alleles( ancestry_states, new_line.sample_counts[sample_index][0], new_line.sample_counts[sample_index][2], A_counts ) ;
        
        /// now compute probability of each sampling arrangement
        for ( std::map<vector<double>,double>::iterator a = A_counts.begin() ; a != A_counts.end() ; ++ a ) {
            
            /// update probs as we iterate across pulses
            double prob_counts = 1 ;
            
            //// c is the ancestry pulse from which we have sampled the genotype
            for ( int c = 0 ; c < a->first.size() ; c ++ ) {
                
                double sum = 0 ;
                
                /// ancestral frequencies are not fixed [default]
                if ( ancestral_fixed == false ) {
                    
                    /// sum across all chromosomes for producing a given set of genotypes which allows us to accomodate error in genotyping similarly to how we accomodate reads, but still assumes one chromosome == one genotype, i.e. permutation not binomial
                    /// j is the number of real A alleles on those chromosomes
                    for ( double j = 0 ; j <= ancestry_states[c] ; j ++ ) {
                        
                        //// then compute probability of ancestral alleles and our sample
                        sum += 1/( new_line.reference_counts[c][2] + ancestry_states[c] + 1 )
                            * 1/nCk[ new_line.reference_counts[c][2] + ancestry_states[c] ][ new_line.reference_counts[c][0] + j ]
                        
                        /// and include total prob of the observed genotypes given true genotypes, j
                        * genotype_total_prob( new_line, a->first[c], ancestry_states[c], j ) ;
                    
                    }
                }
                
                /// if ancestral frequencies are fixed
                else {
                    
                    //// allele frequency is just the allele frequency in the reference panel
                    double f = new_line.reference_counts[c][0] / new_line.reference_counts[c][2] ;

                    /// to accomodate error in genotypes, we still sum across all possible genotype combinations
                    for ( double j = 0 ; j <= ancestry_states[c] ; j ++ ) {
                        
                        /// just figure probabilty of getting j A's given ancestry state and ancestry freq, no integral required
                        sum += nCk[ ancestry_states[c] ][ j ]
                        * pow( f, j )
                        * pow ( 1 - f, ancestry_states[c] - j )

                        /// and include total prob of the observed genotypes given true genotypes, j
                        * genotype_total_prob( new_line, a->first[c], ancestry_states[c], j ) ;
                    }
                }
                
                prob_counts *= sum ;
            }
            
            /// recall that a second records the possible permutations of all genotypes into a given ancestry state and is therefore equivalent to the multinomial combinatoric multiplier
            emission_matrix[i] +=  prob_counts * a->second ;
        }
    }
}

#endif

