/*

 copyright: Russ Corbett-Detig
            rucorbet@ucsc.edu
 
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
using namespace std ;

/// linear algebra library is armadillo
#define ARMA_NO_DEBUG
#include <armadillo>
using namespace arma ;

/// our header files in /src directory
#include "print_usage.h"
#include "factorial.h"
#include "nchoosek.h" 
#include "subsample.h" 
#include "multichoose.h"
#include "multipermute.h"
#include "normalize.h"
#include "ancestry_pulse.h"
#include "ploidy_path.h" 
#include "markov_chain.h"
#include "read_samples.h" 
#include "pulses_to_ancestry.h" 
#include "compute_forward.h"
#include "compute_backward.h"
#include "forward_backward.h"
#include "viterbi.h" 
#include "transition_information.h"
#include "exponentiate_matrix.h"
#include "cmd_line.h"
#include "create_transition_rates.h"
#include "read_cmd_line.h"
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
#include "read_input.h"
#include "nelder_mead.h"
#include "golden_search.h"
#include "bootstrap.h" 

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
    read_file( options, markov_chain_information, state_list, position, recombination_rate, chromosomes ) ;
        
    /// create basic transition information
    cerr << (double) (clock() - t) << " ms" << endl << "computing transition routes\t\t\t" ; t = clock() ;
    /// 3d map to look up by ploidy, start state, end state, and then relevant transition information
    map<int, vector<vector< map< vector<transition_information>, double > > > > transition_matrix_information ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        for ( int p = 0 ; p < markov_chain_information[m].sample_ploidy_path.size() ; p ++ ) {
            create_transition_information( markov_chain_information.at(m).sample_ploidy_path[p].ploidy, transition_matrix_information, state_list[markov_chain_information.at(m).sample_ploidy_path[p].ploidy] ) ;
        }
    }
    
    /// create admixture model(s)
    cerr << (double) (clock() - t) << " ms" << endl << "creating initial admixture model(s)\t\t" ; t = clock();
    vector<vector<pulse> > vertices ;
    int nparams = create_pulses( vertices, options ) ;
    
    /// set number of restarts if unspecified default is factorial * 2
    cerr << (double) (clock() - t) << " ms" << endl << "estimating " << nparams << " parameters\n" ;
    if ( options.n_restarts < 0 ) {
        options.n_restarts = factorial[nparams] * 2 ;
    }
    
    /// vector of models to be evaluated and optimized
    vector<pulse> optimum ;

    /// if there are params to estimate, do amoeba search
    if ( nparams > 1 ) {
        cerr << "starting nelder-mead search\t\t" << endl ;
        optimum = nelder_mead_search( vertices, options, markov_chain_information, transition_matrix_information, recombination_rate, position, state_list ) ;
        cerr << "\n\t\t\t\t\tSEARCH TIME: " << (double) (clock() - t) << " ms" << endl << endl << "optimal model found:\n\n" ;
    }
    
    /// or do golden section line search for single parameter optimization
    else if ( nparams == 1 ) {
        cerr << "starting golden section search\t\t" << endl ;
        optimum = golden_search( options, markov_chain_information, transition_matrix_information, recombination_rate, position, state_list ) ;
        cerr << "\n\t\t\t\t\tSEARCH TIME: " << (double) (clock() - t) << " ms" << endl << endl << "optimal model found:\n\n" ;
    }
    
    /// otherwise just evaluate the supplied model
    else {
        optimum = options.ancestry_pulses ;
        cerr << endl << endl << "evaluating supplied model:\n\n" ;
    }
    
    /// print model
    cerr << "\ttype\ttime\tproportion\n" ;
    cout << "optimum: \n" ;
    cout << "\ttype\ttime\tproportion\n" ;

    vector<double> a = options.ancestry_proportion ;
    for ( int p = 0 ; p < optimum.size() ; p ++ ) {
        optimum[p].proportion = a[optimum[p].type] * optimum[p].fraction_of_remainder ;
        a[optimum[p].type] -= optimum[p].proportion ;
        cerr << "\t" << optimum[p].type << "\t" << optimum[p].time << "\t" << optimum[p].proportion << endl ;
        cout << "\t" << optimum[p].type << "\t" << optimum[p].time << "\t" << optimum[p].proportion << endl ;
    }
    
    /// bootstrap models as necessary
    if ( options.n_bootstraps > 0 ) {
        cerr << "computing " << options.n_bootstraps << " bootstrap models" << endl ;
        
        vector<vector<pulse> > bootstrap = bootstraps( vertices, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_list, chromosomes ) ;
        
        //// print out bootstrapped admixture models 
        for ( int b = 0 ; b < options.n_bootstraps ; b ++ ) {
         
            cout << "bootstrap: " << b << endl ;
            cout << "\ttype\ttime\tproportion\n" ;
            
            cerr << endl << "bootstrap: " << b << endl ;
            cerr << "\ttype\ttime\tproportion\n" ;
            
            vector<double> a = options.ancestry_proportion ;
            for ( int p = 0 ; p < optimum.size() ; p ++ ) {
                bootstrap[b][p].proportion = a[optimum[p].type] * bootstrap[b][p].fraction_of_remainder ;
                a[bootstrap[b][p].type] -= bootstrap[b][p].proportion ;
                cerr << "\t" << bootstrap[b][p].type << "\t" << bootstrap[b][p].time << "\t" << bootstrap[b][p].proportion << endl ;
                cout << "\t" << bootstrap[b][p].type << "\t" << bootstrap[b][p].time << "\t" << bootstrap[b][p].proportion << endl ;
            }
        }
    }
    
    /// create transition rates for the optimal or supplied set of pulses
    cerr << endl << "creating per morgan transition rates\t\t" ; t = clock();
    mat transition_rates = create_transition_rates( optimum, options.ne, options.ancestry_proportion ) ;
    
    /// create transition information
    cerr << (double) (clock() - t) << " ms" << endl << "creating transition matrices\t\t\t" ; t = clock();
    map<int,vector<mat> > transition_matrix ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        create_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information.at(m).number_chromosomes], recombination_rate, position, markov_chain_information.at(m).number_chromosomes, transition_rates ) ;
        for ( int p = 0 ; p < markov_chain_information[m].ploidy_switch.size() ; p ++ ) {
            create_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information[m].ploidy_switch[p]], recombination_rate, position, markov_chain_information[m].ploidy_switch[p], transition_rates ) ;
        }
    }
    
    //// create interploidy transition matrix
    vector<mat> interploidy_transitions = create_interploidy_transitions ( state_list, optimum, options.ancestry_proportion ) ;
    
    /// output viterbi path for optimized model
    if ( options.viterbi == true ) {
        cerr << (double) (clock() - t) << " ms" << endl << "viterbi posterior decoding and printing\t" ; t = clock() ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            markov_chain_information[m].viterbi( position, recombination_rate, state_list, chromosomes, transition_matrix, interploidy_transitions, options.output_pulses, optimum ) ;
        }
    }

    /// output forward-backward full probability distribution by default
    else {
        cerr << (double) (clock() - t) << " ms" << endl << "computing forward probabilities\t" ; t = clock() ;
        double lnl = 0 ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            lnl += markov_chain_information[m].compute_forward_probabilities( transition_matrix, interploidy_transitions ) ;
        }
        
        cerr << "lnl: " << lnl << "\t\t" << (double) (clock() - t) << " ms" << endl << "computing backward probabilities\t\t\t" ; t = clock() ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            markov_chain_information[m].compute_backward_probabilities( transition_matrix, interploidy_transitions ) ;
        }
        
        cerr << (double) (clock() - t) << " ms" << endl << "forward-backward posterior decoding and printing\t\t\t" ; t = clock() ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            markov_chain_information[m].combine_prob( position, state_list, chromosomes, options.output_pulses, optimum ) ;
        }
    }
        
    cerr << (double) (clock() - t) << " ms" << endl ;
    cerr << "total run time:\t\t\t" << (double) (clock() - total) << " ms" << endl ;

	return 0 ; 
}
	

