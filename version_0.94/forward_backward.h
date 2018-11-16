#ifndef __FORWARD_BACKWARD_H
#define __FORWARD_BACKWARD_H

void markov_chain::combine_prob( vector<int> &position, map<int, vector<vector<int> > > &states, vector<string> &chrom, bool output_pulses, vector<pulse> &pulses ) {

    /// find unique ploidy entries in the sample and obtain a sorted list of those entries
    vector<double> unique_ploidy_entries = ploidy_switch ;
    sort( unique_ploidy_entries.begin(), unique_ploidy_entries.end() ) ;
    vector<double>::iterator u = unique ( unique_ploidy_entries.begin(), unique_ploidy_entries.end() ) ;
    unique_ploidy_entries.resize( distance( unique_ploidy_entries.begin(), u ) ) ;
    
    /// output stream
    ofstream out ( output_file.c_str() ) ;

    //// output ancestry states, rather than pulses
    if ( output_pulses == false ) {
        
        /// pulse to ancestry
        /// ploidy >> index of ancestry pulse vector >> ancestry state vector
        map<int,map<int,vector<int> > > ploidy2pulse2ancestry = create_pulse_map( states, pulses, unique_ploidy_entries ) ;
        
        /// get unique and sorted list of vectors of each ancestry type
        map<int,map<vector<int>,int > > ancestry_states ;
        for ( std::map<int,map<int,vector<int> > >::iterator p = ploidy2pulse2ancestry.begin() ; p != ploidy2pulse2ancestry.end() ; ++ p ) {
            for ( std::map<int,vector<int> >::iterator i = p->second.begin() ; i != p->second.end() ; ++ i ) {
                ancestry_states[p->first][i->second] = 0 ;
            }
        }
        
        /// print those unique lists for each ploidy
        for ( map<int,map<vector<int>,int > >::iterator p = ancestry_states.begin() ; p != ancestry_states.end() ; ++ p ) {
            out << "chrom\tposition" ;
            for ( map<vector<int>,int >::reverse_iterator s = p->second.rbegin() ; s != p->second.rend() ; ++ s ) {
                out << "\t" ;
                for ( int c = 0 ; c < s->first.size() - 1 ; c ++ ) {
                    out << s->first[c] << "," ;
                }
                out << s->first.back() ;
            }
            out << endl ;
        }
        
        /// finally create a state count to ploidy map, this way we can just look up the appropriate map from the alpha matrix size
        map<int, int> statecount2ploidy ;
        for ( map<int, vector<vector<int> > >::iterator s = states.begin() ; s != states.end() ; s ++ ) {
            statecount2ploidy[s->second.size()] = s->first ;
        }
        
        //// now iterate through function and print appropriate ancestry states given ploidy at each site
        for ( int i = 0 ; i < alphas.size() ; i ++ ) {
            
            vec smoothed_probs = alphas[i] % betas[i] ;
            normalize( smoothed_probs ) ;
            
            /// current ploidy
            int ploidy = statecount2ploidy[smoothed_probs.size()] ;

            /// need to translate pulse state probs to ancestry state probs
            map<vector<int>,double> ancestry_states ;
            for ( int l = 0 ; l < smoothed_probs.n_rows ; l ++ ) {
                ancestry_states[ploidy2pulse2ancestry[ploidy][l]] += smoothed_probs(l) ;
            }
            
            out << chrom.at(i) << "\t" << position.at(i) ;
            for ( map<vector<int>,double>::reverse_iterator l = ancestry_states.rbegin() ; l != ancestry_states.rend() ; ++ l ) {
                out << "\t" << l->second ;
            }
            out << endl ;
        }
        
        out.close() ;
    }
    
    /// output ancestry pulses rather than ancestry states [default]
    else {
        
        for ( int p = 0 ; p < unique_ploidy_entries.size() ; p ++ ) {
            out << "chrom\tposition" ;
            for ( int s = 0 ; s < states[unique_ploidy_entries[p]].size() ; s ++ ) {
                out << "\t" ;
                for ( int c = 0 ; c < states[unique_ploidy_entries[p]][s].size() - 1 ; c ++ ) {
                    out << states[unique_ploidy_entries[p]][s][c] << "," ;
                }
                out << states[unique_ploidy_entries[p]][s].back() ;
            }
            out << endl ;
        }
        
        for ( int i = 0 ; i < alphas.size() ; i ++ ) {
                        
            vec smoothed_probs = alphas[i] % betas[i] ;
            normalize( smoothed_probs ) ;
            
            out << chrom.at(i) << "\t" << position.at(i) ;
            for ( int l = 0 ; l < smoothed_probs.n_rows ; l ++ ) {
                out << "\t" << smoothed_probs(l) ;
            }
            out << endl ;
        }
    }
    out.close() ;
    return ; 
}

#endif
