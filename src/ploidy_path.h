#ifndef __PLOIDY_PATH_H
#define __PLOIDY_PATH_H

class ploidy_entry {
public:
    string chrom ;
    int start ;
    int stop ;
    double ploidy ;
} ;

void read_ploidy_file ( string path_file, vector<ploidy_entry> &path ) {
    
    ifstream in ( path_file.c_str() ) ;
    while ( !in.eof() ) {
        
        ploidy_entry new_ploidy_entry ;
        
        in >> new_ploidy_entry.chrom >> new_ploidy_entry.start >> new_ploidy_entry.stop >> new_ploidy_entry.ploidy ;
        
        if ( new_ploidy_entry.chrom == "" ) {
            continue ;
        }
        
        path.push_back( new_ploidy_entry ) ;
    }
}

/// create the interploidy transition matrix for each model
vector<mat> create_interploidy_transitions ( map<int,vector<vector<int> > > &state_list, vector<pulse> &vertex, vector<double> &ancestry_proportion ) {
    
    // need to compute what proportion of final is this ancestry type
    vector<double> a = ancestry_proportion ;
    for ( int p = 0 ; p < vertex.size() ; p ++ ) {
        vertex[p].proportion = a[vertex[p].type] * vertex[p].fraction_of_remainder ;
        a[vertex[p].type] -= vertex[p].proportion ;
    }
    
    /// create interploidy transitions from ploidy one to ploidy two
    vector<mat> interploidy_transition_rates ;
    interploidy_transition_rates.resize(2) ;
    interploidy_transition_rates[0].zeros( state_list[2].size(), state_list[1].size() ) ;
    interploidy_transition_rates[1].zeros( state_list[1].size(), state_list[2].size() ) ;
    
    //// iterate through all possible transitions into interploidy states
    for ( int i = 0 ; i < state_list[1].size() ; i ++ ) {
        for ( int j = 0 ; j < state_list[2].size() ; j ++ ) {
            for ( int s = 0 ; s < state_list[1][i].size() ; s ++ ) {
                if ( state_list[2][j][s] - state_list[1][i][s] == 1 ) {
                    interploidy_transition_rates[0]( j, i ) = vertex[s].proportion ;
                    if ( state_list[2][j][s] == 2 ) {
                        interploidy_transition_rates[1]( i, j ) = 1 ;
                    }
                    else {
                        interploidy_transition_rates[1]( i, j ) = 0.5 ;
                    }
                }
            }
        }
    }

    return interploidy_transition_rates ;
}

#endif
