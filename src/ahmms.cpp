/*

 copyright: Russ Corbett-Detig
            rucorbet@ucsc.edu

            Jesper Svedberg
            jsvedber@ucsc.edu

 This is software distributed under the gnu public license version 3.
 
 */

/// headers
#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <string>
#include <fstream>
#include <algorithm>

// Includes specific for Ancestry_HMM-S
#include <cmath>  
#include <cstring>
#include <utility>
#include <iomanip>
#include <cstdlib>
#include <random> // ++++ REQUIRES C++11 ++++


/// linear algebra library is armadillo
#define ARMA_NO_DEBUG
#include <armadillo>
using namespace arma ;
using namespace std ;

/// our header files in /src directory
#include "selection_print_usage.h" // JS
#include "factorial.h"
#include "nchoosek.h" 
#include "selection_subsample.h" 
#include "multichoose.h"
#include "multipermute.h"
#include "normalize.h"
#include "ancestry_pulse.h"
#include "ploidy_path.h" 
#include "selection_class.h" // JS
#include "selection_markov_chain.h"
#include "read_samples.h" 
#include "pulses_to_ancestry.h" 
#include "compute_forward.h"
#include "compute_backward.h"
#include "forward_backward.h"
#include "viterbi.h" 
#include "transition_information.h"
#include "exponentiate_matrix.h"
#include "selection_cmd_line.h"
#include "create_transition_rates.h"
#include "selection_read_cmd_line.h"
#include "evaluate_vertex.h"
#include "check_vertex.h"
#include "sort_vertices.h"
#include "create_pulses.h" 
#include "create_states.h"
#include "input_line.h"
#include "distribute_alleles.h" 
#include "binomial.h"
#include "read_emissions.h"
#include "genotype_emissions.h"
#include "selection_read_input.h"
#include "nelder_mead.h"
#include "golden_search.h"
#include "bootstrap.h"

// Includes specific for Ancestry_HMM-S
#include "selection_get_position.h"
#include "selection_optimize_test_func.h" // Function for testing Nelder-Mead. Remove?
#include "selection_fwd_iter.h"
#include "selection_trajectory.h"
#include "selection_split_vector.h"
#include "selection_forward.h"
#include "selection_stochastic_traj.h"
#include "selection_transition_rates.h"




int main ( int argc, char *argv[] ) {
    
    /// time tracking
    clock_t t = clock() ;
    clock_t total = clock() ;
    
    /// seed prng
    srand (t) ;
    
	// read cmd line 
	cmd_line options ;
    cerr << "reading command line" ; t = clock();
	options.read_cmd_line( argc, argv ) ;

    /// chain objects for each sample
    vector<markov_chain> markov_chain_information ;
    
    /// get sample ids and ploidy from input file
    cerr << "\t\t\t\t" << (double) (clock() - t) << " ms\n" << "reading sample ids and ploidy" ; t = clock();
    read_samples( markov_chain_information, options.sample_file, options.viterbi ) ;

    /// create states matrix
    cerr << "\t\t\t" << (double) (clock() - t) << " ms\n" << "creating states matrix" ; t = clock();
    /// store all possible state space arranged by ploidy and then vector of state counts
    map<int,vector<vector<int> > > state_list ;
    /// now create initial state list
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        for ( int p = 0 ; p < markov_chain_information[m].sample_ploidy_path.size() ; p ++ ) {
            create_initial_states( markov_chain_information.at(m).sample_ploidy_path[p].ploidy, options.ancestry_pulses, state_list ) ;
        }
    }
    
	/// read in panels and update matrices
    cerr << "\t\t\t\t" << (double) (clock() - t) << " ms\n" << "reading data and creating emissions matrices\t" ; t = clock() ;
    /// store recombination rates and positions
    vector<int> position ;
    vector<double> recombination_rate ;
    vector<string> chromosomes ;
    int sel_pos ;
    read_file( options, markov_chain_information, state_list, position, recombination_rate, chromosomes, sel_pos ) ;

    

    /// create basic transition information
    cerr << (double) (clock() - t) << " ms" << endl << "computing transition routes\t\t\t" ; t = clock() ;
    /// 3d map to look up by ploidy, start state, end state, and then relevant transition information
    map<int, vector<vector< map< vector<transition_information>, double > > > > transition_matrix_information ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        for ( int p = 0 ; p < markov_chain_information[m].sample_ploidy_path.size() ; p ++ ) {
            create_transition_information( markov_chain_information.at(m).sample_ploidy_path[p].ploidy, transition_matrix_information, state_list[markov_chain_information.at(m).sample_ploidy_path[p].ploidy] ) ;
        }
    }
    cerr << endl;


    // Below are ahmm-s specific options

    // If using grid search with --grid flag
    if (options.calc_grid == true) {
        int p_start = options.grid_pstart;
        int p_stop = options.grid_pstop;
        int p_step = options.grid_pstep;

        double s_start = options.grid_sstart;
        double s_stop = options.grid_sstop;

        if ( options.limit_sel_space == true ) {
            s_stop = selection_get_max_sel(options.grid_sstart, options.grid_sstop, options.grid_sstep, options.ancestry_pulses[1].proportion, options.ancestry_pulses[1].time, options.ne);
        }
        double s_step = options.grid_sstep;

        cerr << "Grid search. Likelihood calculated for values of selection between " << s_start << " and " << s_stop << endl;

        selection_grid(p_start, p_stop, p_step, s_start, s_stop, s_step, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_list);
        return 0;
    }



    // If testing a single point using --site flag.
    if (options.test_point == true) {
        cerr << "Evaluating point: " << options.test_pos << ", " << options.test_sel << endl;

        map <double,vector<double> > sel_trajectories;
        vector <vector<double> > split_vecs;
        int testpos;
        selection point0;

        if (options.is_coord ==  true) {
            testpos = get_position(options.test_pos, position);

            if (testpos == -1) {
                cerr << "ERROR: specified site not found on chromosome" << endl;
                exit(1);
            }

        }
        else {
            testpos = options.test_pos;
        }

        point0.pos = testpos;
        point0.sel = 0;
        selection_evaluate_point_genotypes( point0, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_list, split_vecs, sel_trajectories ) ;

        selection point1;
        point1.pos = testpos;
        point1.sel = options.test_sel;
        selection_evaluate_point_genotypes( point1, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_list, split_vecs, sel_trajectories ) ;

        cout << "lnL for a selected site s=" << options.test_sel << " at position " << position[point0.pos] << " is: " << point1.lnl-point0.lnl << endl;

        return 0;
    }



    // If using Golden Section Search with --gss flag
    if (options.run_gss == true) {
        selection_golden_section(markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_list);
        return 0;
    }


	return 0 ; 
}
	

