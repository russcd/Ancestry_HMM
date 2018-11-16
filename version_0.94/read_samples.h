#ifndef __READ_SAMPLES_H
#define __READ_SAMPLES_H

void read_samples( vector<markov_chain> &markov_chain_information, string &input_file, bool viterbi ) {

    ifstream in ( input_file.c_str() ) ;
    while ( !in.eof() ) {
        markov_chain new_sample ;
        in >> new_sample.output_file ;
        
        if ( new_sample.output_file == "" ) {
            continue ;
        }
        
        if ( viterbi == false ) {
            new_sample.output_file.append( ".posterior" ) ;
        }
        else {
            new_sample.output_file.append( ".viterbi" ) ;
        }
        
        in >> new_sample.number_chromosomes ;
        
        new_sample.path_file = "null" ;
        
        /// if ploidy path is present, store file and read later
        if ( new_sample.number_chromosomes < 0 ) {
            in >> new_sample.path_file ;
        }

        /// store new samples
        markov_chain_information.push_back( new_sample ) ;
    }
    in.close() ;
    
    /// read or generate ploidy pahts
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        if ( markov_chain_information[m].path_file != "null" ) {
            read_ploidy_file( markov_chain_information[m].path_file, markov_chain_information[m].sample_ploidy_path ) ;
            markov_chain_information[m].ploidy_switch_position.push_back( 0 ) ;
            markov_chain_information[m].ploidy_switch.push_back( markov_chain_information[m].sample_ploidy_path[0].ploidy ) ;
        }
        else {
            markov_chain_information[m].ploidy_switch_position.push_back( 0 ) ;
            markov_chain_information[m].ploidy_switch.push_back( markov_chain_information[m].number_chromosomes ) ;
            ploidy_entry new_entry ;
            new_entry.ploidy = markov_chain_information[m].number_chromosomes ;
            markov_chain_information[m].sample_ploidy_path.push_back( new_entry ) ;
        }
    }
    
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        markov_chain_information[m].end_prob = 1 ;
        markov_chain_information[m].start_prob = 1 ;
    }
    
    return ;
}

#endif
