#ifndef __SELECTION_NELDER_MEAD_H
#define __SELECTION_NELDER_MEAD_H


// The Nelder mead method is no longer used in AHMM-S
// This file is therefor not necessary and is included out of laziness

// sorts vertex. can probably be optimized
void selection_sort_vertex( vector<selection> &v ) {
    sort(v.begin(), v.end());
    reverse(v.begin(), v.end());
}

// generates random double between fmin and fmax
double double_rand(double fmin, double fmax) {
    double f = (double)rand() / RAND_MAX;
    return fmin + f * (fmax - fmin);
}

/*
// calculates 2 vectors with transition rates going away from the site of interest
vector<vector<mat>> selection_transition_rates(selection point, vector<double> &recombination_rate, cmd_line &options) {
    vector<double> vecf ;
    vector<double> vecb ;
    //split_vector(point.pos, recombination_rate, vecb, vecf, options) ;

    double m = options.ancestry_pulses[1].proportion;
    int generations = options.ancestry_pulses[1].time ;
    int n = options.ne ; /// DOUBLE CHECK HAPLOID/DIPLOID!!    
    int tt = 0;

    // generates vector with allele frequencies of selected allele over time
    vector<double> sel_traject ;
    //double halfsel = 0.5 * point.sel; // test. remove
    selection_trajectory(sel_traject, point.sel, tt, m, generations, n) ; // change tt
    
    //cerr << endl << "fwd_iter" << endl;

    //vector<mat> fwd_trans = fwd_iter(vecf, sel_traject, tt, m, generations, n) ;
    //vector<mat> back_trans = fwd_iter(vecb, sel_traject, tt, m, generations, n) ;

    // generates two vectors of transition rates, both going away from the site of interest in different directions

    //cerr << "fwd_vector" << endl;
    vector<mat> fwd_trans = fwd_iter(vecf, sel_traject, m, options.ne) ; //options.ne
    //cerr << endl << "back_vector" << endl;
    vector<mat> back_trans = fwd_iter(vecb, sel_traject, m, options.ne) ;

    vector<vector<mat>> tr_vector;
    tr_vector.push_back(fwd_trans);
    tr_vector.push_back(back_trans);
    return tr_vector;
}


// function for calculating likelihood of parameters
double selection_evaluate_point(selection &point, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, cmd_line &options, map<int,vector<vector<int> > > &state_changes ) {
    cerr << "BP2: Before transition rates." << endl;
    vector<vector<mat>> t_rates = selection_transition_rates(point, recombination_rate, options);
    //vector<vector<mat>> t_rates = selection_transition_rates_genotypes(point, recombination_rate, options, position); // test. remove
    
    //cerr << "BP3: After transition rates." << endl;
    
    double comb_lnl = 0;
    bool go_backwards = false;
    //go_backwards = true;

    for (int i=0 ; i < 2 ; i++) {
        // transition matrix
        map<int,vector<mat> > transition_matrix ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            selection_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information.at(m).number_chromosomes], recombination_rate, position, markov_chain_information.at(m).number_chromosomes, t_rates[i] ) ;
            // Delete maybe
            for ( int p = 0 ; p < markov_chain_information[m].ploidy_switch.size() ; p ++ ) {
                selection_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information[m].ploidy_switch[p]], recombination_rate, position, markov_chain_information[m].ploidy_switch[p], t_rates[i] ) ;
            }
        }
        //cerr << "BP4: After transition matrix." << endl;
        /// compute transitions within a state
        vector<mat> interploidy_transitions ;
        //interploidy_transitions = create_interploidy_transitions( state_changes, vertex, options.ancestry_proportion ) ;
        
        /// now compute the forward probabilities
        double lnl = 0 ;
        cerr << "markov_chain_information.size()  " << markov_chain_information.size() << endl;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            //cerr << "Sample#: " << m << endl;
            lnl += markov_chain_information[m].selection_forward_probabilities( transition_matrix, interploidy_transitions, point, go_backwards ) ;
        }
        //cerr << "BP5: After compute forward. " << i << " " << lnl << endl;
        comb_lnl += lnl;
        go_backwards = true;
    }
    point.lnl = comb_lnl;
    return comb_lnl ;
    // forward probabilities
    // other probabilities ??
} 
 */


// function for determining if point is within bounds
// iterates until it is. (Somewhat dangerous. May cause slow-down or infinite loop if peak is outsite of bounds)
void selection_check_point(selection &point, cmd_line &options) {
    double pos_range_term = 0.1 ;
    double sel_range_term = 0.05 ;
    bool changed;
    int reps = 0;
    do {
        changed = false;
        if (point.pos > options.pos_max) {
            cout << "pos too large: " << point.pos << "\t";
            point.pos = point.pos - ((double) rand() / (RAND_MAX))*pos_range_term * ( options.pos_max - options.pos_min ) ;
            cout << point.pos << endl;
            changed = true;
        }
        else if (point.pos < options.pos_min)
        {
            cout << "pos too small: " << point.pos << "\t";
            point.pos = point.pos + ((double) rand() / (RAND_MAX))*pos_range_term * ( options.pos_max - options.pos_min ) ;
            cout << point.pos << endl;
            changed = true;
        }
        
        if (point.sel > options.sel_max) {
            cout << "sel too large: " << point.sel << "\t";
            point.sel = point.sel - ((double) rand() / (RAND_MAX))*sel_range_term * ( options.sel_max - options.sel_min ) ;
            cout << point.sel << endl;
            changed = true;
        }
        else if (point.sel < options.sel_min) {
            cout << "sel too small: " << point.sel << "\t";
            point.sel = point.sel + ((double) rand() / (RAND_MAX))*sel_range_term * ( options.sel_max - options.sel_min ) ;
            cout << point.sel << endl;
            changed = true;
        }
        if ( reps > 20) {
            cout << "selection_check_point run more than 20 times. Breaking." << endl;
            break;
        }
        reps++;
    } while (changed == true);
}

// function for reflecting, extending or contracting vertex. What it does depends on "mod" argument
selection selection_reflection(vector<selection> &vertex, double mod, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, cmd_line &options, map<int,vector<vector<int> > > &state_changes ) {
    selection centroid;
    selection newpoint;

    // calculate centroid point
    centroid.sel = vertex[0].sel - (vertex[0].sel - vertex[1].sel)/2;
    centroid.pos = vertex[0].pos - (vertex[0].pos - vertex[1].pos)/2;

    // calculate new point to be used
    newpoint.sel = centroid.sel + mod*(centroid.sel - vertex[2].sel);
    newpoint.pos = centroid.pos + mod*(centroid.pos - vertex[2].pos);

    // checks if new point is within bounds
    selection_check_point(newpoint, options);

    // calculates likelihood
    vector <vector<double> > split_vecs;
    map <double,vector<double> > sel_trajectories;
    selection_evaluate_point_genotypes( newpoint, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;

    return newpoint;
}

// function for shinking vertex
vector<selection> selection_shrink(vector<selection> &vertex, double mod, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, cmd_line &options, map<int,vector<vector<int> > > &state_changes ) {
    vector<selection> new_vertex;
    new_vertex.push_back(vertex[0]);
    
    selection point1;
    selection point2;

    point1.sel = vertex[0].sel - (vertex[0].sel - vertex[1].sel)/2;
    point1.pos = vertex[0].pos - (vertex[0].pos - vertex[1].pos)/2;

    point2.sel = vertex[0].sel - (vertex[0].sel - vertex[2].sel)/2;
    point2.pos = vertex[0].pos - (vertex[0].pos - vertex[2].pos)/2;
    
    selection_check_point(point1, options);
    vector <vector<double> > split_vecs;
    map <double,vector<double> > sel_trajectories;
    selection_evaluate_point_genotypes( point1, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;
    new_vertex.push_back(point1);

    selection_check_point(point2, options);
    vector <vector<double> > split_vecs2;
    map <double,vector<double> > sel_trajectories2;
    selection_evaluate_point_genotypes( point2, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs2, sel_trajectories2 ) ;
    new_vertex.push_back(point2);

    return new_vertex;
}

// generates start vertex
vector<selection> selection_start_vertex(cmd_line &options) {
    vector<selection> vertex;

    // Chromosomal position of first point in vertex is within 0.2<x<0.8 of total sequence length.
    int pos_center_seq = (int)((options.pos_max+options.pos_min)*0.5);
    int pos_start_range = (int)((pos_center_seq-options.pos_min)*options.pos_limit); // NOTE: Maybe change 0.8 to option
    
    cout << "pos_center_seq: " << pos_center_seq << " pos_start_range: " << pos_start_range << " options.pos_min: " << options.pos_min << " options.pos_max: " << options.pos_max << " options.pos_limit: " << options.pos_limit << endl;

    double s_min_lim = options.sel_min + (options.sel_max - options.sel_min)*(1-options.sel_limit)*0.5;
    double s_max_lim = options.sel_min + (options.sel_max - options.sel_min)*(options.sel_limit)*0.5;

    selection point1;
    point1.pos = rand()%(2*pos_start_range + 1) + pos_center_seq-pos_start_range;
    point1.sel = double_rand(s_min_lim, s_max_lim);
    vertex.push_back(point1);

    // Chr. position for point 2 and 3 can be within 0.1 of seq.length compared to first point
    int vertex_range = (int)(pos_center_seq*0.2); // As above

    for (int i = 0; i < 2; i++) {
        selection point;
        point.pos = rand()%(2*vertex_range + 1) + point1.pos-vertex_range;
        point.sel = double_rand(s_min_lim, s_max_lim);
        vertex.push_back(point);
    }
    return vertex;
}

/*
vector<selection> old_selection_start_vertex(cmd_line &options) {
    vector<selection> vertex;

    for (int i; i < 3; i++) {
        selection point;
        point.pos = rand()%(options.pos_max-options.pos_min + 1) + options.pos_min;
        point.sel = double_rand(options.sel_min, options.sel_max);
        vertex.push_back(point);
    }
    return vertex;
}
*/


// function for using Nelder-Mead to find optimal values for selection and position
selection selection_nelder_mead(cmd_line &options, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, map<int,vector<vector<int> > > &state_changes) {

    int r = 1;

    double alpha = 1 ; // for reflecting vertex
    double gamma = 2 ; // for extension vertex
    double rho = -0.5 ; // used for contracting vertex. Note negative value
    double sigma = 0.5 ; // for shrinking vertex

    double best_lnl = -1.7976931348623157E+308;
    selection best_optimum;
    vector<selection> vertex;

    // loop to repeat nelder mead search 20 times (should be changed)
    for (int reps=0 ; reps<20 ; reps++) {
        
        // generates random start vertex
        vertex = selection_start_vertex(options);
        best_optimum = vertex[0];

        // calculates likelihood for start vertex
        // Check parameters
        for ( int v = 0 ; v < vertex.size() ; v ++ ) {
                vector <vector<double> > split_vecs;
                map <double,vector<double> > sel_trajectories; // WARNING: should probably only be defined once in the beginning of Nelder-Mead and then sent as an argument to all subfunctions
                vertex[v].lnl = selection_evaluate_point_genotypes( vertex[v], markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;
                cout << vertex[v] << endl;
            }
        
        // sorts start vertex
        selection_sort_vertex(vertex) ;
        
        int iteration = 0 ;

        // loop for Nelder-Mead search. Breaks when the change between each iteration is smaller than a tolerance
        while ( vertex[0].lnl-vertex.back().lnl > options.tolerance ) {
            /*if (vertex[0].pos == vertex[1].pos && vertex[0].pos == vertex[2].pos) {
                cout << "ABORTING. Degenerate position: " << vertex[0].pos << endl;
                break;
            }*/

            cout << "Vertex: " << r << "\t" << iteration << "\t" << vertex[0].lnl-vertex.back().lnl << "\t" ;
            iteration ++ ;

            cout << vertex[0] << "\t" << vertex[1] << "\t" << vertex[2] << "\t" << endl;

            // creates reflected vertex
            selection reflection = selection_reflection(vertex, alpha, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes);
            //cout << "R1" << endl;

            // checks new vertex and if an extension should be tried instead
            if (reflection < vertex[0] && vertex[1] < reflection) {
                vertex.back() = reflection ;
                cout << "Reflection1" << endl;
                selection_sort_vertex(vertex) ;
                continue ;
            }
            else if (vertex[0] < reflection) {
                selection extension = selection_reflection(vertex, gamma, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes);
                if ( reflection < extension ) {
                    vertex.back() = extension ;
                    cout << "Extension" << endl;
                }
                else {
                    vertex.back() = reflection ; 
                    cout << "Reflection2" << endl;
                }
                selection_sort_vertex(vertex) ;
                continue;
            }
        
            // if reflected or extended vertex is not good enough, tries contracting or shrinking vertex
            selection contraction = selection_reflection(vertex, rho, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes);
            
            if (vertex.back() < contraction) {
                vertex.back() = contraction ;
                cout << "Contraction" << endl;
            }
            else {
                vertex = selection_shrink(vertex, sigma, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes);
                cout << "Shrink" << endl;
            }
            selection_sort_vertex(vertex) ;
        }

        selection_sort_vertex(vertex) ;

        // outputs best point. This will be compared between each iteration of the N-M algorithm
        if ( best_lnl < vertex[0].lnl ) {
            best_optimum = vertex[0] ;
            best_lnl = vertex[0].lnl ;
        }
        cout << "Max point: " << vertex[0] <<  " best_lnl: " << best_lnl << endl;
    }
    
    cout << "END" << endl;
    cout << vertex[0] << "\t" << vertex[1] << "\t" << vertex[2] << "\t" << endl;
    return best_optimum;

}


#endif