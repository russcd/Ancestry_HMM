#ifndef __SORT_VERTICES_H
#define __SORT_VERTICES_H

void sort_vertices( vector<vector<pulse> > &vertices, vector<double> &lnl ) {
    
    vector<vector<pulse> > new_pulses ;
    vector<double> new_lnl ;
    
    while ( lnl.size() > 0 ) {
        
        double max = lnl[0] ;
        int index = 0 ;
        
        for ( int i = 0 ; i < lnl.size() ; i ++ ) {
            if ( lnl[i] > max ) {
                max = lnl[i] ;
                index = i ;
            }
        }
        new_pulses.push_back( vertices[index] ) ;
        new_lnl.push_back( max ) ;
        lnl.erase( lnl.begin() + index ) ;
        vertices.erase( vertices.begin() + index ) ;
    }
    
    /// swap out for ordered sets
    swap ( lnl, new_lnl ) ;
    swap ( vertices, new_pulses ) ;
    
}

#endif

