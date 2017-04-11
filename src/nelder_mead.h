#ifndef __NELDER_MEAD_H
#define __NELDER_MEAD_H

/// identify centroid for simplex
vector<pulse> create_centroid ( vector<vector<pulse> > &vertices ) {
    vector<pulse> centroid = vertices[0] ;
    for ( int p = 0 ; p < vertices[0].size() ; p ++ ) {
        if ( centroid[p].time_fixed == false ) {
            for ( int v = 1 ; v < vertices.size()-1 ; v ++ ) {
                centroid[p].time += vertices[v][p].time ;
            }
            centroid[p].time /= ( vertices.size() - 1 ) ;
        }
        if ( centroid[p].proportion_fixed == false ) {
            for ( int v = 1 ; v < vertices.size()-1 ; v ++ ) {
                centroid[p].fraction_of_remainder += vertices[v][p].fraction_of_remainder ;
            }
            centroid[p].fraction_of_remainder /= ( vertices.size() - 1 ) ;
        }
    }
    return centroid ;
}

/// create reflection point
vector<pulse> create_reflection ( vector<pulse> &centroid, vector<pulse> &worst, double &alpha ) {
    vector<pulse> reflection = centroid ;
    for ( int p = 0 ; p < centroid.size() ; p ++ ) {
        if ( centroid[p].time_fixed == false ) {
            reflection[p].time = centroid[p].time + alpha * ( centroid[p].time - worst[p].time ) ;
        }
        if ( centroid[p].proportion_fixed == false ) {
            reflection[p].fraction_of_remainder = centroid[p].fraction_of_remainder + alpha * ( centroid[p].fraction_of_remainder - worst[p].fraction_of_remainder ) ;
        }
    }
    return reflection ;
}

/// create expansion, contraction, etc
vector<pulse> create_expansion ( vector<pulse> &centroid, vector<pulse> &worst, double &scale ) {
    
    vector<pulse> expansion = centroid ;
    for ( int p = 0 ; p < centroid.size() ; p ++ ) {
        if ( centroid[p].time_fixed == false ) {
            expansion[p].time = centroid[p].time + scale * ( worst[p].time - centroid[p].time ) ;
        }
        if ( centroid[p].proportion_fixed == false ) {
            expansion[p].fraction_of_remainder = centroid[p].fraction_of_remainder + scale * ( worst[p].fraction_of_remainder - centroid[p].fraction_of_remainder ) ;
        }
    }
    return expansion ;
}

/// retstart nelder_mead with random simplex
void random_restart( vector<vector<pulse> > &vertices, cmd_line &options ) {
    for ( int v = 0 ; v < vertices.size() ; v ++ ) {
        for ( int p = 0 ; p < vertices[v].size() ; p ++ ) {
            if ( vertices[v][p].time_fixed == false ) {
                vertices[v][p].time = options.t_min + rand()/((double)RAND_MAX) * ( options.t_max - options.t_min ) ;
            }
            if ( vertices[v][p].proportion_fixed == false ) {
                vertices[v][p].fraction_of_remainder = options.p_min + rand()/((double)RAND_MAX) * ( options.p_max - options.p_min ) ;
            }
        }
        check_vertex( vertices[v], options ) ;
    }
}

/// optimation via amoeba search
vector<pulse> nelder_mead_search( vector<vector<pulse> > &vertices, cmd_line &options, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, map<int,vector<vector<int> > > &state_changes ) {
    
    /// amoeba search parameters
    double alpha = 1 ;
    double gamma = 2 ;
    double rho = 0.5 ;
    double sigma = 0.5 ;
    
    /// record best optimum so far
    vector<pulse> best_optimum = vertices[0] ;
    double best_lnl = -1.7976931348623157E+308 ; // set initial best to minimum double value, if we underflow that, shiiittttt
    
    for ( int r = 0 ; r < options.n_restarts+1 ; r ++ ) {
    
        /// lnl of each vertex to be stored in vector
        vector<double> lnls (vertices.size(),0) ;
    
        /// evaluate likelihood of initial points
        for ( int v = 0 ; v < vertices.size() ; v ++ ) {
            lnls[v] = evaluate_vertex( vertices[v], markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes ) ;
        }
        
        /// improvement in new position relative to old position both in lnl units
        int iteration = 0 ;
        while ( lnls[0]-lnls.back() > options.tolerance ) {
            
            cerr << r << "\t" << iteration << "\t" << lnls[0]-lnls.back() << "\t" ;
            iteration ++ ;
            
            cerr << endl ;
            for ( int v = 0 ; v < vertices.size() ; v ++ ) {
                cerr << "\tvertex: " << v << "\t" << "lnl: " << lnls[v] << endl ;
                for ( int p = 0 ; p < vertices[v].size() ; p ++ ) {
                    vertices[v][p].print() ;
                }
            }
            cerr << endl ;
            
            /// create centroid and evaluate
            vector<pulse> centroid = create_centroid ( vertices ) ;
            
            /// reflection point
            vector<pulse> reflection = create_reflection( centroid, vertices.back(), alpha ) ;
            check_vertex( reflection, options ) ;
            double lnl_reflection = evaluate_vertex( reflection, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes ) ;
            
            cerr << "\treflection\t" << "lnl: " << lnl_reflection << endl ;
            for ( int p = 0 ; p < reflection.size() ; p ++ ) {
                reflection[p].print() ;
            }
            cerr << endl ;
            
            if ( lnl_reflection < lnls[0] && lnl_reflection > lnls[lnls.size()-2] ) {
                cerr << "\t\t\t\tREFLECTION ACCEPTED\n" ;
                lnls.back() = lnl_reflection ;
                vertices.back() = reflection ;
                sort_vertices( vertices, lnls ) ;
                continue ;
            }
            
            /// expansion point
            else if ( lnl_reflection >= lnls[0] ) {
                vector<pulse> expansion = create_expansion( centroid, reflection, gamma ) ;
                check_vertex( expansion, options ) ;
                double lnl_expansion = evaluate_vertex( expansion, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes ) ;
                
                cerr << "\texpansion\t" << "lnl: " << lnl_expansion << endl ;
                for ( int p = 0 ; p < expansion.size() ; p ++ ) {
                    expansion[p].print() ;
                }
                cerr << endl ;
                
                if ( lnl_expansion > lnl_reflection ) {
                    cerr << "\t\t\t\tEXPANSION ACCEPTED\n" ;
                    lnls.back() = lnl_expansion ;
                    vertices.back() = expansion ;
                    sort_vertices( vertices, lnls ) ;
                    continue ;
                }
                else {
                    cerr << "\t\t\t\tREFLECTION ACCEPTED\n" ;
                    lnls.back() = lnl_reflection ;
                    vertices.back() = reflection ;
                    sort_vertices( vertices, lnls ) ;
                    continue ;
                }
            }
            
            /// contraction
            vector<pulse> contraction = create_expansion( centroid, vertices.back(), rho ) ;
            check_vertex( contraction, options ) ;
            double lnl_contraction = evaluate_vertex( contraction, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes ) ;
            
            cerr << "\tcontraction\t" << "lnl: " << lnl_contraction << endl ;
            for ( int p = 0 ; p < contraction.size() ; p ++ ) {
                contraction[p].print() ;
            }
            cerr << endl ;
            
            if ( lnl_contraction > lnls.back() ) {
                cerr << "\t\t\t\tCONTRACTION ACCEPTED\n" ;
                lnls.back() = lnl_contraction ;
                vertices.back() = contraction ;
                sort_vertices( vertices, lnls ) ;
                continue ;
            }
            
            /// shrink
            else {
                cerr << "\t\t\t\tSHRINK ACCEPTED\n" ;
                for ( int v = 1 ; v < vertices.size() ; v ++ ) {
                    vector<pulse> replacement = create_expansion( vertices[0], vertices[v], sigma ) ;
                    check_vertex( replacement, options ) ;
                    lnls[v] = evaluate_vertex( replacement, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes ) ;
                    vertices[v] = replacement ;
                }
                sort_vertices( vertices, lnls ) ;
                continue ;
            }
        }
        
        if ( best_lnl < lnls[0] ) {
            best_optimum = vertices[0] ;
            best_lnl = lnls[0] ;
        }
        
        cerr << "restart: " << r << "\t" << lnls[0]-lnls.back() << endl ;
        for ( int p = 0 ; p < vertices[0].size() ; p ++ ) {
            vertices[0][p].print() ;
        }
        cerr << endl ;
        
        random_restart( vertices, options ) ;
        for ( int v = 0 ; v < vertices.size() ; v++ ) {
            lnls[v] = evaluate_vertex( vertices[v], markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes ) ;
        }
        sort_vertices( vertices, lnls ) ;
    }
    
    /// sort points by lnl, with point 0 having highest lnl
    /// sort_vertices( vertices, lnls ) ;
    return best_optimum ;
}


#endif

