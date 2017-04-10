#ifndef __CREATE_PULSES_H
#define __CREATE_PULSES_H

/// create pulses and create data points for optimization function
int create_pulses ( vector<vector<pulse> > &vertices, cmd_line &options ) {
    
    //// count parameters to estimate and figure out which are dependent on each other
    int nparams = 0 ;
    
    /// count dependent pulses
    map<int, int> dependent_pulses ;
    
    /// figure out how many ancestry types pulse more than once
    /// record as pulse and count of pulse
    vector<double> proportion_accounted ( options.ancestry_proportion.size(), 0 ) ;
    for ( int p = options.ancestry_pulses.size() -1 ; p > -1 ; p-- ) {
        proportion_accounted[options.ancestry_pulses[p].type] += options.ancestry_pulses[p].proportion ;
        if ( options.ancestry_pulses[p].proportion_fixed == false ) {
            if ( dependent_pulses.find( options.ancestry_pulses[p].type ) == dependent_pulses.end() ) {
                dependent_pulses[ options.ancestry_pulses[p].type ] = 1 ;
                options.ancestry_pulses[p].fraction_of_remainder = 1 ;
                options.ancestry_pulses[p].proportion_fixed = true ;
            }
            else {
                options.ancestry_pulses[p].fraction_of_remainder = options.ancestry_pulses[p].proportion/proportion_accounted[options.ancestry_pulses[p].type] ;
                dependent_pulses[ options.ancestry_pulses[p].type ] ++ ;
                nparams ++ ;        /// only have to estimate non-fixed second - n of an ancestry type
            }
        }
        else {
            options.ancestry_pulses[p].fraction_of_remainder = options.ancestry_pulses[p].proportion / proportion_accounted[options.ancestry_pulses[p].type] ;
        }
        if ( options.ancestry_pulses[p].time_fixed == false ) {
            nparams ++ ;
        }
    }
    
    /// if there are model parameters to estimate, create vertices for the starting simplex
    if ( nparams > 0 ) {
        
        /// create the starting point x0
        vector<pulse> x0 = options.ancestry_pulses ;
        vertices.push_back( x0 ) ;
        
        /// need nparm + 1 points in simplex
        for ( int n = 0 ; n < nparams ; n ++ ) {
            int param = 0 ;
            vector<pulse> vertex = x0 ;
            for ( int p = 0 ; p < options.ancestry_pulses.size() ; p ++ ) {
                if ( options.ancestry_pulses[p].time_fixed == false ) {
                    if ( param == n ) {
                        if ( options.t_max - vertex[p].time > vertex[p].time - options.t_min ) {
                            vertex[p].time /= ( 1 - options.t_length ) ;
                        }
                        else {
                            vertex[p].time *= ( 1 - options.t_length ) ;
                        }
                    }
                    param ++ ;
                }
                if ( options.ancestry_pulses[p].proportion_fixed == false ) {
                    if ( param == n ) {
                        if ( options.p_max - vertex[p].fraction_of_remainder > vertex[p].fraction_of_remainder - options.p_min ) {
                            vertex[p].fraction_of_remainder *= options.p_length ;
                        }
                        else {
                            vertex[p].fraction_of_remainder *= ( 1 - options.p_length ) ;
                        }
                    }
                    param ++ ;
                }
            }
            check_vertex( vertex, options ) ;
            vertices.push_back( vertex ) ;
        }
    }
    return nparams ;
}

#endif


