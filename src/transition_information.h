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

/// object will be sets of sets of outcomes 
void enumerate_rows ( double ploidy, int row, vector<int> &pulse_types, vector<vector<int> > &output ) { 

	// if we have no chromosomes of this type to transition, skip
	if ( row == 0 ) { 
		vector<int> end (1, 0) ; 
		output.push_back( end ) ; 
		return ; 
	}

	/// find all possible outcomes for a given input 
	vector<vector<int> > end_states = multichoose( row, pulse_types ) ; 

	/// end is then stored as the output pulse (index) and count (value)
	for ( int p = 0 ; p < end_states.size() ; p ++ ) { 
		vector<int> end ( pulse_types.size(), 0 ) ; 
		for ( int m = 0 ; m < end_states[p].size() ; m ++ ) { 
			end[end_states[p][m]] ++ ; 
		}
		output.push_back(end) ; 
	}
}

/// for each row, compute the binomial coefficient n!/k1!k2!...kn! 
/// possible_transitions is the sets for each pulse of possible final transitions 
/// binom_coeffs is the final data object to store it
void compute_binomial_coefficients ( vector<vector<int> > &transitions, vector<double> &binomial_coefficients ) {

	for ( int t = 0 ; t < transitions.size() ; t ++ ) { 
		double multiplier = 1 ; 
		double sum = 0 ;
		for ( int pulse = 0 ; pulse < transitions[t].size() ; pulse ++ ) { 
			sum += transitions[t][pulse] ; 
			multiplier *= 1/factorial[transitions[t][pulse]] ; 
		}
		multiplier *= factorial[sum] ; 

		binomial_coefficients.push_back( multiplier ) ; 
	}
}

void enumerate_combinations ( vector<vector<int> > &groups, vector<vector<vector<int> > > &possible_transitions ) { 

	/// to replace groups with updated groups
	vector<vector<int> > groups_so_far ; 

	/// scroll across all pulses
	for ( int p = 0 ; p < possible_transitions.size() ; p ++ ) { 

		/// record the groups that we have so far
		groups_so_far = groups ; 

		/// then empty groups and rebuild the vector 
		groups.clear() ; 

		/// if empty we just make a vector with everything we've got
		if ( groups_so_far.empty() ) {

			/// just create a vector of these
			for ( int row = 0 ; row < possible_transitions[p].size() ; row ++ ) { 
				vector<int> new_group ( 1, row ) ;
				groups.push_back( new_group ) ; 
			}
		}

		/// otherwise, just append to all existing groups and re-record
		else { 
			for ( int g = 0 ; g < groups_so_far.size() ; g ++ ) { 
				for ( int row = 0 ; row < possible_transitions[p].size() ; row ++ ) {
					groups.push_back( groups_so_far[g] ) ; 
					groups.back().push_back( row ) ; 
				}
			}
		}
	}
}

/// multiply across groups bcs to get ways of producing each matrix
void compute_group_bcs ( vector<vector<int> > &groups, vector<vector<double> > &binomial_coefficients, vector<double> &group_bcs ) { 

	for ( int g = 0 ; g < groups.size() ; g ++ ) { 
		group_bcs.push_back( 1 ) ; 
		for ( int pulse = 0 ; pulse < groups[g].size() ; pulse ++ ) {
			group_bcs.back() *= binomial_coefficients[pulse][groups[g][pulse]] ; 
		}
	}
}

/// figure out the end state for each group 
void determine_end_state( vector<vector<vector<int> > > &possible_transitions, vector<transition_information> &new_transition, vector<int> &group, vector<int> &end ) {
	
	for ( int p = 0 ; p < possible_transitions.size() ; p ++ ) { 
		for ( int e = 0 ; e < possible_transitions[p][group[p]].size() ; e ++ ) {
			transition_information t ; 
			t.start_state = p ; 
			t.end_state = e ; 
			t.transition_count = possible_transitions[p][group[p]][e] ; 
			end[e] += possible_transitions[p][group[p]][e] ; 
			t.ibd_transition = false ;  
			new_transition.push_back(t) ; 
		}
	}
}

/// create map of states back onto their position in the state vector 
void create_state_map ( vector<vector<int> > &states, map<vector<int>, int> &state_map ) {
	for ( int s = 0 ; s < states.size() ; s ++ ) { 
		state_map[states[s]] = s ; 
	}
}

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

	/// generate vector of all possible pulse types 
	vector<int> pulse_types ; 
	for ( int p = 0 ; p < states[0].size() ; p ++ ) { 
		pulse_types.push_back( p ) ; 
	}	

	/// create state map for quick lookup at the end
	map<vector<int>,int> state_map ;  
	create_state_map( states, state_map ) ; 

    /// iterate aross all possible start states
    for ( int i = 0 ; i < states.size() ; i ++ ) {

    	/// all combinations 
    	vector<vector<vector<int> > > possible_transitions ;
    	vector<vector<double> > binomial_coefficients ; 

    	/// iterate through all pulses to find possible outputs (rows) 
    	for ( int row = 0 ; row < states[i].size() ; row ++ ) { 

	    	/// find all possible outputs for the pulse (row)  
    		vector<vector<int> > output ;  
	    	enumerate_rows( ploidy, states[i][row], pulse_types, output ) ;
	    	possible_transitions.push_back( output ) ; 

		    /// find binomial coefficients 
		    vector<double> pulse_bcs ; 
	    	compute_binomial_coefficients( output, pulse_bcs ) ;
	    	binomial_coefficients.push_back( pulse_bcs ) ; 
    	}

	    /// find all possible combinations among the possible transitions
	    vector<vector<int> > groups ; 
	    enumerate_combinations( groups, possible_transitions ) ; 

	    /// compute group bcs
	    vector<double> group_bcs ; 
	    compute_group_bcs( groups, binomial_coefficients, group_bcs ) ; 

	    /// determine the end state fo each combination, count transitions 
	    for ( int g = 0 ; g < groups.size() ; g ++ ) { 
	    	vector<int> end ( pulse_types.size(), 0 ) ; 
	    	vector <transition_information> new_transition ;
	    	determine_end_state( possible_transitions, new_transition, groups[g], end ) ; 
	    	transition_matrix_information[ploidy][i][state_map[end]][new_transition] += group_bcs[g] ;
	    }
	}
}

#endif

























