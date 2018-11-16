#ifndef __TRANSITION_INFORMATION_H
#define __TRANSITION_INFORMATION_H

/// binomial information for transitions to look up during computation
class transition_information {
public:
    int start_state ;                                     /// will be start state
    int end_state ;
    int transition_count ;
    
    bool ibd_transition ;
    
    friend bool operator < ( transition_information a, transition_information b ) {
        if ( a.start_state > b.start_state ) return a.start_state > b.start_state ;
        else if ( a.end_state > b.end_state ) return a.end_state > b.end_state ;
        else return a.transition_count > b.transition_count ;
    }
} ;

//// find all possible transitions between two adjacent markers and record
void create_transition_information( double &ploidy, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<vector<int> > &states ) {
    
    /// skip if we already have the ploidy figured out
    if ( transition_matrix_information.find(ploidy) != transition_matrix_information.end() ) {
        return ;
    }
    
    //// resize vectors
    transition_matrix_information[ploidy].resize(states.size()) ;
    for ( int i = 0 ; i < states.size() ; i ++ ) {
        transition_matrix_information[ploidy][i].resize(states.size()) ;
    }
    
    /// iterate aross all possible start states
    for ( int i = 0 ; i < states.size() ; i ++ ) {
        
        /// the vector of places we start at
        vector<int> chrom_start ;
        for ( int l = 0 ; l < states[i].size() ; l ++ ) {
            for ( int r = 0 ; r < states[i][l] ; r ++ ) {
                chrom_start.push_back(l) ;
            }
        }
        
        // iterate across all possible destinations
        for ( int j = 0 ; j < states.size() ; j ++ ) {
            
            /// vector of places we end at
            vector<int> chrom_end ;
            for ( int l = 0 ; l < states[j].size() ; l ++ ) {
                for ( int r = 0 ; r < states[j][l] ; r ++ ) {
                    chrom_end.push_back(l) ;
                }
            }
            
            /// acquire all possible end arrangements
            vector<vector<int> > end_arrangements = multipermute( chrom_end ) ;
            
            /// now create binomial transition infromation
            /// iterate across all end sets
            for ( int e = 0 ; e < end_arrangements.size() ; e ++ ) {
                
                /// store data as we go
                map<pair<int,int>, int > transition_counts ;
                
                /// figure out what changes are required
                for ( int l = 0 ; l < end_arrangements[e].size() ; l ++ ) {
                    pair <int, int> trans (chrom_start[l],end_arrangements[e][l]) ;
                    if ( transition_counts.find( trans ) == transition_counts.end() ) {
                        transition_counts[trans] = 0 ;
                    }
                    transition_counts[trans] ++ ;
                }
                
                vector<transition_information> transition_paths ;
                
                for ( std::map<pair<int,int>,int>::iterator m = transition_counts.begin() ; m != transition_counts.end() ; ++m ) {
                    transition_information new_trans ;
                    new_trans.start_state = m->first.first ;
                    new_trans.end_state = m->first.second ;
                    new_trans.transition_count = m->second ;
                    new_trans.ibd_transition = false ;
                    transition_paths.push_back( new_trans ) ;
                }
                
                /// count occurences
                transition_matrix_information[ploidy][i][j][transition_paths] ++ ;
            }
        }
    }
}

#endif
