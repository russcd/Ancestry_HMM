#ifndef __MARKOV_CHAIN_H
#define __MARKOV_CHAIN_H

/// will include all basic information for input data and functions to compute forward, backward, forward-backward, and viterbi
class markov_chain {
public:
    
    /// sample attributes
    string output_file ;
    double number_chromosomes ;
    
    /// file describing ploidy path across the genome
    string path_file ;
    
    /// data object storing ploidy paths to be looked up during emissions computation with chromosome as key
    vector<ploidy_entry> sample_ploidy_path ;
    vector<int> ploidy_switch_position ;
    vector<double> ploidy_switch ;
    
    /// read from input file
    vector<vec> emission_probabilities ;
    
    /// create initial states to be stored
    double start_prob ;
    double end_prob ;
    
    /// forward probs
    vector<vec> alphas ;
    double compute_forward_probabilities( map<int, vector<mat> > &transition_matrix, vector<mat> &interploidy_transitions  ) ;
    
    /// backward probs
    vector<vec> betas ;
    void compute_backward_probabilities( map<int, vector<mat> > &transition_matrix, vector<mat> &interploidy_transitions ) ;
    
    /// combine probs
    void combine_prob( vector<int> &position, map<int, vector<vector<int> > > &states, vector<string> &chrom, bool output_pulses, vector<pulse> &pulses ) ;
    
    /// output viterbi paths
    void viterbi( vector<int> &position, vector<double> &recombination_rate, map<int,vector<vector<int> > > &states, vector<string> &chrom, map<int,vector<mat> > &transition_probabilites, vector<mat> &interploidy_transitions, bool output_pulses, vector<pulse> &pulses ) ;
    
} ;

#endif
