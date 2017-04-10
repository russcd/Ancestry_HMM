#ifndef __PULSES_TO_ANCESTRY_H
#define __PULSES_TO_ANCESTRY_H

vector<int> pulses_to_ancestry ( vector<int> counts, vector<pulse> &pulses ) {
    
    /// create ancestry types
    vector<int> ancestry (2,0) ;
    
    /// now add up the counts
    for ( int c = 0 ; c < counts.size() ; c ++ ) {
        if ( pulses[c].type + 1 > ancestry.size() ) {
            ancestry.resize( pulses[c].type + 1 ) ;
        }
        ancestry[pulses[c].type] += counts[c] ;
    }
    
    return ancestry ;
}

/// create a map that links the states to their ancestry type counts
map<int,map<int,vector<int> > > create_pulse_map ( map<int, vector<vector<int> > > states, vector<pulse> pulses, vector<double> ploidy_list ) {
    
    map<int,map<int,vector<int> > > ploidy2pulses2ancestry ;
    for ( int p = 0 ; p < ploidy_list.size() ; p ++ ) {
        for ( int i = 0 ; i < states[ploidy_list[p]].size() ; i ++ ) {
            ploidy2pulses2ancestry[ploidy_list[p]][i] = pulses_to_ancestry( states[ploidy_list[p]][i], pulses ) ;
        }
    }
    return ploidy2pulses2ancestry ; 
}

#endif
