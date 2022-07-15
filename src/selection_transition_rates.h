#ifndef __SELECTION_TRANSITION_RATES_H
#define __SELECTION_TRANSITION_RATES_H


// generates transition matrix.
void selection_transition_matrix(map<int,vector<mat> > &transition_matrix , vector<vector< map< vector<transition_information>, double > > > &transition_info, vector<double> &recombination_rate, vector<int> &positions, double &number_chromosomes, vector<mat> &transition_rates ) {
    
    /// check if we already computed this for this sample ploidy
    if ( transition_matrix.find( number_chromosomes ) != transition_matrix.end() ) {
        return ;
    }
    
    /// else, have to create entire matrix
    /// first create data object of approporate size
    //transition_matrix[number_chromosomes].resize(recombination_rate.size()) ;
    transition_matrix[number_chromosomes].resize(transition_rates.size()) ;
    
    //// iterate across all positions and compute transition matrixes
    for ( int p = 0 ; p < transition_rates.size() ; p ++ ) {

        /// create actual transition matrix
        transition_matrix[number_chromosomes][p].set_size( transition_info.size(), transition_info.size() ) ;
        transition_matrix[number_chromosomes][p].fill( 0 ) ;
        
        /// population transitions by summing across all routes
        for ( int i = 0 ; i < transition_info.size() ; i ++ ) {
            for ( int j = 0 ; j < transition_info[i].size() ; j ++ ) {
                for ( std::map<vector<transition_information>,double>::iterator t = transition_info[i][j].begin() ; t != transition_info[i][j].end() ; ++ t ) {
                    double prob_t = 1 ;
                    for ( int r = 0 ; r < t->first.size() ; r ++ ) {
                        prob_t *= pow( transition_rates[p](t->first[r].start_state,t->first[r].end_state), t->first[r].transition_count ) ;
                    }
                    transition_matrix[number_chromosomes][p](j,i) += prob_t * t->second ;
                }
            }
        }        
    }
}

// as above, but to use when printing expected genotype frequences across the chromosome
vector<vector<mat> > selection_transition_rates_genotypes(selection point, vector<double> &recombination_rate, cmd_line &options, vector<int> &position, vector < vector<double> > &genofreqs, vector <vector<double> > &split_vecs, map <double,vector<double> > &sel_trajectories) {
    
    // The selection coefficient is halved, since the model uses the haploid coefficient,
    // but the user specifies the diploid coefficient from the command line
    point.sel = 0.5 * point.sel;

    // splits the chromosome into two vectors and trims them based on window size
    if (split_vecs.size() == 0) {
        split_vecs = split_vector(point.pos, recombination_rate, options) ;
    }

    double m = options.ancestry_pulses[1].proportion;
    int generations = options.ancestry_pulses[1].time ;
    int n = options.ne ; /// DOUBLE CHECK HAPLOID/DIPLOID!!    
    int tt = 0;

    // generates vector with allele frequencies of selected allele over time
    vector<double> sel_traject ;
    map <double,vector<double> >::iterator it;
    it = sel_trajectories.find(point.sel);

    if (options.use_stochastic == true) {
        if (it == sel_trajectories.end()) {
            selection_stochastic_trajectory(sel_traject, point.sel, m, generations, n, options.stochastic_reps) ; // change tt
            sel_trajectories[point.sel] = sel_traject;
        }
        else {
            sel_traject = it->second;
        }
    }
    else {
        if (it == sel_trajectories.end()) {
            selection_trajectory(sel_traject, point.sel, tt, m, generations, n) ; // change tt
            sel_trajectories[point.sel] = sel_traject;
        }
        else {
            sel_traject = it->second;
        }
    }
    
    
    vector<double> gf1;
    vector<double> gf2;
    genofreqs.push_back(gf1);
    genofreqs.push_back(gf2);

    vector<mat> fwd_trans;
    vector<mat> back_trans;

    // checks which trajectory function to use (4point, 3point or forward iteration)
    if (options.traj_function == 4) {
        if (point.sel == 0.0) {
            fwd_trans = fwd_iter_genotype_freq(split_vecs[0], sel_traject, m, options.ne, genofreqs[0]) ; //options.ne
            back_trans = fwd_iter_genotype_freq(split_vecs[1], sel_traject, m, options.ne, genofreqs[1]) ;
        }
        else {
            genofreqs[0].push_back(sel_traject.back());
            genofreqs[1].push_back(sel_traject.back());
            fwd_trans = approx_curve(split_vecs[0], sel_traject, m) ; //options.ne
            back_trans = approx_curve(split_vecs[1], sel_traject, m) ;
        }
    }
    else if (options.traj_function == 3) {
        genofreqs[0].push_back(sel_traject.back());
        genofreqs[1].push_back(sel_traject.back());
        fwd_trans = approx_curve_3point(split_vecs[0], sel_traject, m) ; //options.ne
        back_trans = approx_curve_3point(split_vecs[1], sel_traject, m) ;
    }
    else { 
        fwd_trans = fwd_iter_genotype_freq(split_vecs[0], sel_traject, m, options.ne, genofreqs[0]) ; //options.ne
        back_trans = fwd_iter_genotype_freq(split_vecs[1], sel_traject, m, options.ne, genofreqs[1]) ;
    }
    
    vector<vector<mat> > tr_vector;
    tr_vector.push_back(fwd_trans);
    tr_vector.push_back(back_trans);
    return tr_vector;
}

double selection_evaluate_point_genotypes(selection &point, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, cmd_line &options, map<int,vector<vector<int> > > &state_changes, vector <vector<double> > &split_vecs, map <double,vector<double> > &sel_trajectories) {

    vector < vector<double> > genofreqs ;
    
    vector<vector<mat> > t_rates = selection_transition_rates_genotypes(point, recombination_rate, options, position, genofreqs, split_vecs, sel_trajectories); // test. remove
    
    double comb_lnl = 0;
    bool go_backwards = false;

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
        /// compute transitions within a state
        vector<mat> interploidy_transitions ;
        //interploidy_transitions = create_interploidy_transitions( state_changes, vertex, options.ancestry_proportion ) ;
        
        /// now compute the forward probabilities
        double lnl = 0 ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            lnl += markov_chain_information[m].selection_forward_probabilities_genotypes( transition_matrix, interploidy_transitions, point, go_backwards, genofreqs[i], position ) ;
        }
        comb_lnl += lnl;
        go_backwards = true;
    }

    // tries to handle the the situation where the likelihood is undefined = nan
    if ( options.limit_sel_space == false ) {
        if (isnan(comb_lnl) == true ) {
            comb_lnl = -1000000000;
        }
    }
    point.lnl = comb_lnl;
    return comb_lnl ;
}

// function for calculating likelihoods in a grid
// takes start, stop and step values for selection and position
void selection_grid(int p_start, int p_stop, int p_step, double s_start, double s_stop, double s_step, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, cmd_line &options, map<int,vector<vector<int> > > &state_changes ) {

    map <double,vector<double> > sel_trajectories;

    // check if limits are specified in chromosome coordinates or per site (ie from 1000bp to 2000bp or SNP #1 to SNP #10)
    if (options.is_coord ==  true) {
        int p_start = get_position(options.grid_pstart, position);
        int p_stop = get_position(options.grid_pstop, position);

        cerr << "p_start " << p_start << "p_stop " << p_stop << "p_step " << p_step << endl;

        if (p_start == -1) {
            cerr << "ERROR: specified start coordinate for grid not found on chromosome" << endl;
            exit(1);
        }
        if (p_stop == -1) {
            cerr << "ERROR: specified STOP coordinate for Golden section search not found on chromosome" << endl;
            exit(1);
        }
    }
    if ( p_start > p_stop ) {
        cerr << "ERROR: specified stop coordinate for grid is located before start coordinate." << endl;
        exit(1);
    }
    

    // loop over all sites in specified region with p_step steps
    for (int p = p_start; p < p_stop; p+=p_step) {
        
        vector <vector<double> > split_vecs;

        // generate neutral transition rate for normalization / calculating likelihood ratio
        selection point0;
        point0.pos = p;
        point0.sel = 0;
        selection_evaluate_point_genotypes( point0, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;

        // loop over selective coeffients
        for (double s = s_start; s < s_stop; s=s+s_step) {
            selection point;
            point.pos = p;
            point.sel = s;
            selection_evaluate_point_genotypes( point, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;

            cout << position[point.pos] << "\t" << point.sel << "\t" << setprecision(12) << point.lnl-point0.lnl << endl;
        }
    }
}


void selection_golden_section(vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, cmd_line &options, map<int,vector<vector<int> > > &state_changes) {
    map <double,vector<double> > sel_trajectories;
    double GR = (sqrt(5) + 1) / 2;

    int pstart;
    int pstop;

    double gs_start = options.gs_sstart;
    double gs_stop = options.gs_sstop;

    // check if using full selection space, if not use s that reaches 0.99 frequency in specified number of generations
    if ( options.limit_sel_space == true ) {
        gs_stop = selection_get_max_sel(options.gs_sstart, options.gs_sstop, options.gs_sstep, options.ancestry_pulses[1].proportion, options.ancestry_pulses[1].time, options.ne);
    }
    cerr << "Golden section search search. Likelihood calculated for values of selection between " << gs_start << " and " << gs_stop << endl;

    // check which coordinate format is used, and convert bp to snp number
    if (options.is_coord ==  true) {
        pstart = get_position(options.gs_pstart, position);
        pstop = get_position(options.gs_pstop, position);

        if (pstart == -1) {
            cerr << "ERROR: specified START coordinate for Golden section search not found on chromosome" << endl;
            exit(1);
        }
        if (pstop == -1) {
            cerr << "ERROR: specified STOP coordinate for Golden section search not found on chromosome" << endl;
            exit(1);
        }
    }

    else {
        pstart = options.gs_pstart;
        pstop = options.gs_pstop;
    }

    if ( pstart > pstop ) {
        cerr << "ERROR: specified STOP coordinate for Golden section search is located before START coordinate." << endl;
        exit(1);
    }

    // Golden section searchj
    for (int p = pstart; p < pstop; p+=options.gs_pstep) {
        vector <vector<double> > split_vecs;

        // calculate likelihood for neutral case
        selection point0;
        point0.pos = p;
        point0.sel = 0;
        selection_evaluate_point_genotypes( point0, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;

        // do golden section search
        selection point1;
        selection point2;
        selection point3;
        selection point4;

        point1.pos = p;
        point2.pos = p;
        point3.pos = p;
        point4.pos = p;

        point1.sel = gs_start;
        point2.sel = gs_stop;
        point3.sel = point2.sel - (point2.sel - point1.sel) / GR;
        point4.sel = point1.sel + (point2.sel - point1.sel) / GR;

        selection_evaluate_point_genotypes( point1, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;
        selection_evaluate_point_genotypes( point2, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;
        selection_evaluate_point_genotypes( point3, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;
        selection_evaluate_point_genotypes( point4, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;

        int i = 0;

        while (abs(point3.sel - point4.sel) > options.gs_precision) {
            if (point3.lnl > point4.lnl) {
                point2 = point4;
                point4 = point3;
                point3.sel = point2.sel - (point2.sel - point1.sel) / GR;
                selection_evaluate_point_genotypes( point3, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;
            }
            else {
                point1 = point3;
                point3 = point4;
                point4.sel = point1.sel + (point2.sel - point1.sel) / GR;
                selection_evaluate_point_genotypes( point4, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_changes, split_vecs, sel_trajectories ) ;
            }
            i++;
        }

        // verbose output. used for debugging
        //cout << position[point0.pos] << "\t" << (point3.sel+point4.sel)/2 << "\t" << setprecision(12) << ((point3.lnl+point4.lnl)/2)-point0.lnl << "\t" << i << "\t" << point3.lnl << "\t" << point4.lnl << "\t" << (point3.lnl+point4.lnl)/2 << "\t" << point0.lnl << endl;
        
        // cout << position[point0.pos] << "\t" << (point3.sel+point4.sel)/2 << "\t" << setprecision(12) << ((point3.lnl+point4.lnl)/2)-point0.lnl << "\tsel_trajectories.size() " << sel_trajectories.size() << endl;

        cout << position[point0.pos] << "\t" << (point3.sel+point4.sel)/2 << "\t" << setprecision(12) << ((point3.lnl+point4.lnl)/2)-point0.lnl << endl;
    }
}


#endif