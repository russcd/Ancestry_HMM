#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "factorial.h"

using namespace std;

//// facotorial look up table
const vector<long double> factorial = create_factorial() ;


/// binomial information for transitions
class binomial_info {
public:
    int up ;
    int down ;
    long double bc ;
} ;

/// hmm input lines
class input_line {
public:
    int pos ;
    long double r1_A_count ;
    long double r1_a_count ;
    long double r2_A_count ;
    long double r2_a_count ;
    long double A_count ;
    long double a_count ;
    long double recombination_rate ;
} ;

/// binomial and combinatoric functions
long double nCk ( long double n, long double k ) {
    return factorial[n]/( factorial[k]*factorial[n-k] ) ;
}
/// binomial probs
long double binomial ( long double n, long double k, long double prob ) {
    return nCk(n,k)*pow((1-prob),n-k)*pow(prob,k) ;
}

/// normalize vector
long double normalize( vector<long double> &vec, long double minimum_probability ) {
    long double sum = 0 ;
    long double max = 0 ;
    for ( int i = 0 ; i < vec.size() ; i ++ ) {
        sum += vec.at(i) ;
    }
    for ( int i = 0 ; i < vec.size() ; i ++ ) {
        vec.at(i) /= sum ;
        if ( vec.at(i) < minimum_probability ) {
            vec.at(i) = minimum_probability ;
        }
        if ( vec.at(i) > 1 - minimum_probability ) {
            vec.at(i) = 1 - minimum_probability ;
        }
    }
    return sum ;
}

/// template
class cmd_line ;

/// will include all basic information for input data
/// to run forward backward algorithm on each sample separately
class markov_chain {
public:
    
    /// cmd line
    string file ;
    string output_file ;
    long double number_chromosomes ;
    
    /// read from input file
    vector<int> positions ;
    vector<long double> recombination_rates ;
    vector<vector<long double> > emission_probabilities ;
    void read_file ( cmd_line &options ) ;
    
    /// create initial states vectors
    vector<long double> start_probs ;
    vector<long double> end_probs ;
    void create_initial_states () ;
    
    /// estimate ancestry proportion
    long double ancestry_weight ;
    long double estimate_ancestry_proportion( vector<vector<long double> > &emission_probabilities ) ;
    
    /// binomials for transition matrix
    vector<vector<vector<binomial_info> > > binomial_coefficients ;
    void compute_binomial_coefficients() ;
    
    /// transition information
    vector<vector<vector<long double> > > transition_probabilites ;
    void resize_transition_matrix() ;
    
    /// forward probs
    vector<vector<long double> > alphas ;
    void compute_transitions_probabilities( long double t, cmd_line &options ) ;
    long double compute_forward() ;
    
    /// backward probs
    vector<vector<long double> > betas ;
    void compute_backward() ;
    
    /// combine probs
    void combine_prob() ;
    
} ;

/// command line information and global parameters
class cmd_line {
public:
    
    // whether parental genotype frequencies are fixed and known
    bool ancestral_fixed ;
    
    /// terms to bound the golden section search
    long double t_max ;
    long double t_min ;
    
    /// user may also specify t
    long double t_fixed ;
    
    /// diploid effective population size ( i.e. 2n )
    long double ne ;
    
    /// tolerance for golden section search
    long double tolerance ;
    
    /// proportion ancestry 1
    bool fixed_a ;
    long double ancestry_proportion ;
    bool initial_ancestry ;
    
    /// minimum recombinational distance between markers
    long double minimum_distance ;
    
    /// error rates
    long double error_rate ;
    long double genotype_error_rate_sample ;
    
    /// bool sample is expressed as genotypes, not read counts
    bool genotype ;
    
    /// minimum emissions probabilities
    long double minimum_emissions_probability ;
    
    /// input states
    string initial_states_input ;
    
    /// read relevant information
    void read_cmd_line ( int argc, char *argv[], vector<markov_chain> &markov_chain_information ) ;

} ;

void markov_chain::combine_prob( ) {
    
    ofstream out ( output_file.c_str() ) ;
    
    for ( int i = 0 ; i < alphas.size() ; i ++ ) {
        
        vector<long double> smoothed_probs ;
        for ( int l = 0 ; l < alphas.at(i).size() ; l ++ ) {
            long double sprob = alphas.at(i).at(l) * betas.at(i).at(l) ;
            smoothed_probs.push_back( sprob ) ;
        }
        normalize( smoothed_probs, 0 ) ;
        
        out << positions.at(i) ;
        for ( int l = 0 ; l < smoothed_probs.size() ; l ++ ) {
            out << "\t" << smoothed_probs.at(l) ;
        }
        out << endl ;
    }
    
    out.close() ;
}

/// this will calculate the
vector<long double> create_emissions_matrix_genotype( long double n, input_line &new_line, bool ancestral_fixed, long double genotype_error_rate_sample ) {
    
    vector<long double> emissions_matrix ;
    
    /// numbers to make algebra easiesr later
    long double p0 = new_line.r1_A_count + new_line.r1_a_count ;
    long double p1 = new_line.r2_A_count + new_line.r2_a_count ;
    
    for ( long double i = 0 ; i <= n ; i ++ ) {
        
        /// total prob
        long double total_prob_state = 0 ;
        
        /// starting condidtions
        int g0A_min = 0 ;
        if ( new_line.A_count - i > 0 ) g0A_min = new_line.A_count - i ;
        int g0A_max = n - i ;
        if ( new_line.A_count < g0A_max ) g0A_max = new_line.A_count ;
        
        /// now sum across all possible genotypes within the sample
        if ( ancestral_fixed == false ) {
            for ( long double g0A = g0A_min ; g0A <= g0A_max ; g0A ++ ) {
                long double g1A = new_line.A_count - g0A ;
                total_prob_state +=  nCk( n - i , g0A ) * 1/(p0+n-i+1) * 1/nCk( p0 + n - i, new_line.r1_A_count + g0A ) * nCk( i , g1A ) * 1/(p1+i+1)*1/nCk( p1 + i, new_line.r2_A_count + g1A ) ;
            }
        }
        else {
            long double f0 = new_line.r1_A_count/p0*(1-genotype_error_rate_sample)+new_line.r1_a_count/p0*genotype_error_rate_sample ;
            long double f1 = new_line.r2_A_count/p1*(1-genotype_error_rate_sample)+new_line.r2_a_count/p1*genotype_error_rate_sample  ;
            for ( long double g0A = g0A_min ; g0A <= g0A_max ; g0A ++ ) {
                long double g1A = new_line.A_count - g0A ;
                total_prob_state += binomial( n-i, g0A, f0 ) * binomial( i , g1A, f1 ) ;
            }
        }
        
        emissions_matrix.push_back( total_prob_state ) ;
    }
    
    normalize( emissions_matrix, 0 ) ;
    return emissions_matrix ;
}

/// this will calculate the emissions probabilities
vector<long double> create_emissions_matrix( long double n, input_line &new_line, long double minimum_probability, long double error, bool &ancestral_fixed ) {
    vector<long double> emission_matrix ;

    /// numbers to make algebra easiesr later
    long double p0 = new_line.r1_A_count + new_line.r1_a_count ;
    long double p1 = new_line.r2_A_count + new_line.r2_a_count ;
    
    /// i is the number of anestry 1
    /// n-i is the number of ancestry 0
    for ( long double i = 0 ; i <= n ; i ++ ) {
        
        /// total prob
        long double total_prob_state = 0 ;
        
        /// now the probability that you sample reads from each ancestry type
        for ( long double r0 = 0 ; r0 <= new_line.A_count+new_line.a_count ; r0 ++ ) {
            long double r1 = new_line.A_count+new_line.a_count - r0 ;
            
            // the probability of this particular sampling of reads
            long double p_reads = binomial( new_line.A_count+new_line.a_count, r1, i/n ) ;
            
            /// skip impossible combinations
            if ( p_reads == 0 ) {
                continue ;
            }
            
            /// now iterate through possible genotypes within those read counts
            long double start = 0 ;
            if ( new_line.A_count - r1 > start ) {
                start = new_line.A_count - r1 ;
            }
            long double stop = r0 ;
            if ( new_line.A_count < stop ) {
                stop = new_line.A_count ;
            }
            
            for ( long double r0A = start ; r0A <= stop ; r0A ++ ) {
                long double r1A = new_line.A_count - r0A ;
                
                /// to store sum across all possible sample genotypes
                long double sum0 = 0 ;
                long double sum1 = 0 ;
                
                if ( ancestral_fixed == false ) {
                    for ( int j = 0 ; j <= n - i ; j ++ ) {
                        long double p_A = j/(n-i)*(1-error)+(1-j/(n-i))*error ;
                        sum0 += 1/( p0 + n - i + 1 ) * 1/nCk( p0 + n - i, new_line.r1_A_count + j ) * nCk( n-i, j ) * pow( p_A, r0A ) * pow( 1-p_A, r0 - r0A ) ;
                    }
                    for ( int j = 0 ; j <= i ; j ++ ) {
                        long double p_A = j/i*(1-error) + (1-j/i)*error ;
                        sum1 += 1/( p1 + i + 1 ) * 1/nCk( p1 + i, new_line.r2_A_count + j ) * nCk( i, j ) * pow( p_A, r1A ) * pow( 1-p_A, r1-r1A ) ;
                    }
                }
                else {
                    long double f0 = new_line.r1_A_count/p0 ;
                    long double f1 = new_line.r2_A_count/p1 ;
                    for ( int j = 0 ; j <= n - i ; j ++ ) {
                        long double p_A = j/(n-i)*(1-error)+(1-j/(n-i))*error ;
                        sum0 += nCk( n-i, j ) * pow( f0, j ) * pow( 1-f0, n-i-j ) * pow( p_A, r0A ) * pow( 1-p_A, r0 - r0A ) ;
                    }
                    for ( int j = 0 ; j <= i ; j ++ ) {
                        long double p_A = j/i*(1-error) + (1-j/i)*error ;
                        sum1 += nCk( i, j ) * pow(f1, j) * pow( 1-f1, i-j ) * pow( p_A, r1A ) * pow( 1-p_A, r1-r1A ) ;
                    }
                }
                
                total_prob_state += p_reads * nCk( r0, r0A ) * sum0 * nCk( r1, r1A ) * sum1 ;
            }
        }
        emission_matrix.push_back( total_prob_state ) ;
    }
    
    normalize( emission_matrix, minimum_probability ) ;
    return emission_matrix ;
}

void markov_chain::compute_backward () {
    
    /// get final state
    vector<long double> last_position ;
    for ( int i = 0 ; i < emission_probabilities.back().size() ; i ++ ) {
        last_position.push_back( emission_probabilities.back().at(i) * end_probs.at(i) ) ;
    }
    normalize( last_position, 0 ) ;
    betas.push_back( last_position ) ;
    
    /// now got from t-1 to 0 to iterate through backwards
    for ( int i = emission_probabilities.size()-2 ; i > -1 ; i -- ) {
        vector<long double> beta ;
        for ( int s1 = 0 ; s1 < emission_probabilities.at(i).size() ; s1 ++ ) {
            long double total_prob = 0 ;
            for ( int s2 = 0 ; s2 < transition_probabilites.at(i+1).at(s1).size() ; s2 ++ ) {
                /// this one iterates through s1 >> s2 since it's backwards probs
                total_prob += betas.back().at(s2) * transition_probabilites.at(i+1).at(s2).at(s1) * emission_probabilities.at(i).at(s1) ;
            }
            beta.push_back(total_prob) ;
        }
        normalize( beta, 0 ) ;
        betas.push_back( beta ) ;
    }
    
    /// reverse betas so they are in the same order as alphas for gamma equations
    reverse( betas.begin(), betas.end() ) ;
}

long double markov_chain::compute_forward() {
    
    /// return log norm
    long double log_normalization = 0 ;
    
    /// clear the fw probs matrix
    alphas.clear() ;
    
    /// get initial state set
    vector<long double> first_position ;
    for ( int i = 0 ; i < emission_probabilities.at(0).size() ; i ++ ) {
        first_position.push_back( emission_probabilities.at(0).at(i) * start_probs.at(i) ) ;
    }
    normalize( first_position, 0 ) ;
    alphas.push_back( first_position ) ;
    
    /// do all other sites
    for ( int i = 1 ; i < emission_probabilities.size() ; i ++ ) {
        vector<long double> alpha(emission_probabilities.at(0).size(),0) ;
        for ( int s1 = 0 ; s1 < emission_probabilities.at(i).size() ; s1 ++ ) {
            for ( int s2 = 0 ; s2 < transition_probabilites.at(i).at(s1).size() ; s2 ++ ) {
                /// s2 to s1 for fwd probs since it's the probability of getting to s1
                alpha.at(s1) += alphas.back().at(s2) * transition_probabilites.at(i).at(s1).at(s2) * emission_probabilities.at(i).at(s1) ;
            }
        }
        
        long double sum = normalize( alpha, 0 ) ;
        log_normalization += log( sum ) ;
        alphas.push_back( alpha ) ;
    }
    
    return log_normalization ;
}


/// compute segment transition rates
void compute_segment_transitions ( cmd_line &options, long double recombination_rate, long double &prob_01, long double &prob_10, long double distance, long double t ) {
    
    /// to get this back to per bp recombination rates
    recombination_rate = recombination_rate / distance ;
    
    /// per basepair rates of transition into and out of segments
    /// based in large part on liang and nielsen 2014
    long double p1 = recombination_rate * options.ne * options.ancestry_proportion * ( 1 - exp(-t/options.ne) ) ;
    long double p2 = recombination_rate * options.ne * ( 1 - options.ancestry_proportion ) * ( 1 - exp(-t/options.ne) ) ;
    
    /// now we have to solve the continuous time MCMC approximation for transition across the entire segment
    prob_01 = p1/(p1+p2) - p1/(p1+p2)*exp(-(p1+p2)*distance ) ;
    prob_10 = p2/(p1+p2) - p2/(p1+p2)*exp(-(p1+p2)*distance ) ;
    
}

void markov_chain::compute_transitions_probabilities( long double t, cmd_line &options ) {
    
    for ( int r = 0 ; r < recombination_rates.size() ; r ++ ) {
        
        //// acquire transition probabilities
        long double prob_01 ;
        long double prob_10 ;
        long double distance = 1 ;
        if ( r != 0 ) {
            long double distance = positions.at(r) - positions.at(r-1) ;
        }
        
        compute_segment_transitions( options, recombination_rates.at(r), prob_01, prob_10, distance, t ) ;
        
        // now compute all possible path probabilities and sum for each start and end point
        for ( int to = 0 ; to < transition_probabilites.at(r).size() ; to ++ ) {
            for ( int from = 0 ; from < transition_probabilites.at(r).at(to).size() ; from ++ ) {
                transition_probabilites.at(r).at(to).at(from) = 0 ;
                
                for ( int b = 0 ; b < binomial_coefficients[to][from].size(); b++ ) {
                    
                    long double new_prob = binomial_coefficients[to][from][b].bc         ///product of binomial coefficients
                    * pow( prob_01, binomial_coefficients[to][from][b].up )              /// how many lienages are going 0 > 1
                    * pow( prob_10, binomial_coefficients[to][from][b].down )            /// how many lineages are going 1 > 0
                    * pow( 1 - prob_01, number_chromosomes - from - binomial_coefficients[to][from][b].up ) /// how many stay 0
                    * pow( 1 - prob_10, from - binomial_coefficients[to][from][b].down ) ;  /// how many stay 1
                    
                    transition_probabilites.at(r).at(to).at(from) += new_prob ;
                }
            }
        }
    }
}

void markov_chain::resize_transition_matrix () {
    /// resize transition matrix vector
    transition_probabilites.resize(emission_probabilities.size()) ;
    for ( int i = 0 ; i < transition_probabilites.size() ; i ++ ) {
        transition_probabilites.at(i).resize( number_chromosomes + 1 ) ;
        for ( int n = 0 ; n <= number_chromosomes ; n ++ ) {
            transition_probabilites.at(i).at(n).resize( number_chromosomes + 1 ) ;
        }
    }
}

void markov_chain::compute_binomial_coefficients( ) {
    
    binomial_coefficients.resize( number_chromosomes + 1 ) ;
    
    for ( int i = 0 ; i <= number_chromosomes ; i ++ ) {                        /// i will be to
        binomial_coefficients.at( i ).resize( number_chromosomes + 1 ) ;
        for ( int l = 0 ; l <= number_chromosomes ; l ++ ) {                    /// l will be from
            for ( int a = 0 ; a <= l ; a ++ ) {                 /// number of lineages that changed downwards
                for ( int b = 0 ; b <= number_chromosomes - l ; b ++ ) {        /// number of lineages that changed upwards
                    if ( l - a + b != i ) continue ;  /// if this is not a doable combination, skip
                    binomial_info new_bc ;
                    new_bc.up = b ;
                    new_bc.down = a ;
                    new_bc.bc = nCk( number_chromosomes - l, b ) * nCk( l, a ) ;
                    binomial_coefficients[i][l].push_back( new_bc ) ;
                }
            }
        }
    }
}

long double markov_chain::estimate_ancestry_proportion ( vector<vector<long double> > &emission_probabilities ) {
    
    long double proportion = 0 ;
    for ( int l = 0 ; l < emission_probabilities.size() ; l ++ ) {
        for ( long double i = 0 ; i < emission_probabilities.at(l).size() ; i ++ ) {
            proportion += i * emission_probabilities.at(l).at(i)/(emission_probabilities.at(l).size()-1) ;
        }
    }
    proportion /= emission_probabilities.size() ;
    return ( proportion * ancestry_weight ) ;
}

void markov_chain::create_initial_states() {
    for ( int i = 0 ; i < number_chromosomes + 1 ; i ++ ) {
        start_probs.push_back( 1/(number_chromosomes+1) ) ;
        end_probs.push_back( 1 ) ;
    }
}

void markov_chain::read_file ( cmd_line &options ) {
    
    ifstream in ( file.c_str() ) ;
    long double extra_recombination = 1 ; //// since the first site transition matrxi does not matter, we can print whatever here.
    
    while( !in.eof() ) {
        input_line new_line ;
        in >> new_line.pos >> new_line.r1_A_count >> new_line.r1_a_count >> new_line.r2_A_count >> new_line.r2_a_count >> new_line.A_count >> new_line.a_count >> new_line.recombination_rate ;
        
        /// ignore lines where recombination may not be suffiicent to make sites independent
        extra_recombination += new_line.recombination_rate ;
        if ( extra_recombination < options.minimum_distance || new_line.a_count + new_line.A_count == 0 ) {
            continue ;
        }
        new_line.recombination_rate = extra_recombination ;
        extra_recombination = 0 ;
        
        /// if recombination rate is 0, set to 1e-9
        if ( new_line.recombination_rate == 0 ) {
            new_line.recombination_rate = 1e-9 ;
        }
        
        /// create transition matrix
        recombination_rates.push_back( new_line.recombination_rate ) ;
        
        /// create emissions
        vector<long double> emissions ;
        if ( options.genotype == false ) {
            emissions = create_emissions_matrix( number_chromosomes, new_line, options.minimum_emissions_probability, options.error_rate, options.ancestral_fixed ) ;
        }
        else {
            emissions = create_emissions_matrix_genotype( number_chromosomes, new_line, options.ancestral_fixed, options.genotype_error_rate_sample ) ;
        }
        emission_probabilities.push_back( emissions ) ;
        
        /// record position
        positions.push_back( new_line.pos ) ;
    }
}

void cmd_line::read_cmd_line ( int argc, char *argv[], vector<markov_chain> &markov_chain_information ) {
    
    ///defaults
    ancestral_fixed = false ;                       /// set to true for qtl or experimental evolution application s
    minimum_emissions_probability = 0 ;             /// minimum probabiliity of emission state
    initial_states_input = "NULL" ;                 /// if nothign is specified, this is a uniform prior
    minimum_distance = 0 ;                          /// minimum distance in morgans between sites to consider
    ne = 2e4 ;                                      /// actually 2ne
    t_max = 10000 ;
    t_min = 1 ;
    t_fixed = 0 ;
    ancestry_proportion = -1 ;
    fixed_a = false ;
    initial_ancestry = false ;
    tolerance = 10 ;
    error_rate = 0.01 ;                             /// this is the error rate per read per site
    genotype_error_rate_sample = 1e-5 ;             /// error rate for sample genotypes if genotypes are provided to hmm
    genotype = false ;
    
	/// accept command line parameters
	for (int i=1; i<argc; i++) {
        if ( strcmp(argv[i],"-g") == 0 ) {
            genotype = true ;
        }
        if ( strcmp(argv[i],"-a") == 0 ) {
            ancestry_proportion = atof(argv[++i]) ;
            fixed_a = true ;
        }
        if ( strcmp(argv[i],"--fix") == 0 ) {
            ancestral_fixed = true ;
        }
        if ( strcmp(argv[i],"--ai") == 0 ) {
            ancestry_proportion = atof(argv[++i]) ;
            initial_ancestry = true ;
        }
        if ( strcmp(argv[i],"--tmax") == 0 ) {
            t_max = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--tolerance") == 0 ) {
            tolerance = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i], "-e" ) == 0 ) {
            error_rate = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i], "--ge" ) == 0 ) {
            genotype_error_rate_sample = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i], "-t" ) == 0 ) {
            t_fixed = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--tmin") == 0 ) {
            t_min = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--ne") == 0 ) {
            ne = 2 * atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-m") == 0 ) {
            minimum_emissions_probability = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-i") == 0 ) {
            markov_chain new_sample ;
            new_sample.number_chromosomes = atof(argv[++i]) ;
            new_sample.file = argv[++i] ;
            new_sample.output_file = argv[++i] ;
            markov_chain_information.push_back(new_sample) ;
        }
        if ( strcmp(argv[i],"-s") == 0 ) {
            initial_states_input = argv[++i] ;
        }
        if ( strcmp(argv[i],"-d") == 0 ) { 
            minimum_distance =  atof(argv[++i]) ; 
        }
    }
}

long double golden_search ( cmd_line &options, vector<markov_chain> &markov_chain_information ) {
     
    /// now do golden search until we reach tolerance threshhold and stop
    long double phi = ( sqrt(5) - 1 ) / 2 ;
    long double t_low = options.t_max + phi * ( options.t_min - options.t_max ) ;
    long double t_high = options.t_min + phi * ( options.t_max - options.t_min ) ;
    long double low_bracket = options.t_min ;
    long double high_bracket = options.t_max ;
    long double lnl_low = 0 ;
    long double lnl_high = 0 ;
    long double a_low = 0 ;
    long double a_high = 0 ;
    int iteration = 0 ;
    
    cerr << "\titeration\tlow bound\tlow test\thigh test\thigh bound\tlnl low\tlnl high\n" ;
    while ( options.tolerance < t_high - t_low ) {
        
        /// compute probabilty of low point
        if ( lnl_low == 0 ) {
            for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
                markov_chain_information.at(m).compute_transitions_probabilities( t_low, options );
                lnl_low += markov_chain_information.at(m).compute_forward() ;
            }
            
            if ( options.fixed_a == false ) {
                a_low = 0 ;
                for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
                    a_low += markov_chain_information.at(m).estimate_ancestry_proportion( markov_chain_information.at(m).alphas ) ;
                }
            }
        }
        
        /// compute probability of high point
        if ( lnl_high == 0 ) {
            for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
                markov_chain_information.at(m).compute_transitions_probabilities( t_high, options );
                lnl_high += markov_chain_information.at(m).compute_forward() ;
            }
            
            if ( options.fixed_a == false ) {
                a_high = 0 ;
                for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
                    a_high += markov_chain_information.at(m).estimate_ancestry_proportion( markov_chain_information.at(m).alphas ) ;
                }
            }
        }
        
        cerr << "\t" << iteration << "\t" << low_bracket << "\t" << t_low << "\t" << t_high << "\t" << high_bracket << "\t" << lnl_low << "\t" << lnl_high << endl ;
        
        iteration ++ ;
        
        /// if true, we know that the maximum is between t_low and high_bracket
        if ( lnl_high >= lnl_low ) {
            low_bracket = t_low ;
            t_low = t_high ;
            t_high = low_bracket + ( high_bracket - low_bracket ) * phi ;
            lnl_low = lnl_high ;
            lnl_high = 0 ;
            
            if ( options.fixed_a == 0 ) {
                a_low = a_high ;
                a_high = 0 ;
                options.ancestry_proportion = a_low ;
            }
        }
    
        /// otherwise, the maximum is between low_bracket and t_high
        else {
            high_bracket = t_high ;
            t_high = t_low ;
            t_low = high_bracket + ( low_bracket - high_bracket ) * phi ;
            lnl_high = lnl_low ;
            lnl_low = 0 ;
            
            if ( options.fixed_a == 0 ) {
                a_high = a_low ;
                a_low = 0 ;
                options.ancestry_proportion = a_high ;
            }
        }
    }
    
    return ( ( low_bracket + high_bracket ) / 2 ) ;
}

int main ( int argc, char *argv[] ) {
	
    /// chain objects
    vector<markov_chain> markov_chain_information ;
    
	// read cmd line 
	cmd_line options ; 
	options.read_cmd_line( argc, argv, markov_chain_information ) ;
    
    /// read in initial state_matrix if exists
    cerr << "creating initial states\n" ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        markov_chain_information.at(m).create_initial_states() ;
    }
    
	/// read in panels and update matrices
    cerr << "reading data and creating emissions matrices\n" ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        markov_chain_information.at(m).read_file( options ) ;
    }
    
    /// ancestry estimte weights
    cerr << "computing ancestry estimation weights\n" ;
    long double total_chroms = 0 ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        total_chroms += markov_chain_information.at(m).number_chromosomes ;
    }
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        markov_chain_information.at(m).ancestry_weight = markov_chain_information.at(m).number_chromosomes / total_chroms ;
    }

    /// estimate ancestry proportion from emissions matrix if not specified
    if ( options.fixed_a == false && options.initial_ancestry == false ) {
        options.ancestry_proportion = 0 ;
        cerr << "estimating initialial ancestry proportion\n" ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            options.ancestry_proportion += markov_chain_information.at(m).estimate_ancestry_proportion( markov_chain_information.at(m).emission_probabilities ) ;
        }
        cerr << "\tinitial ancestry estimate:\t" << options.ancestry_proportion << endl ;
    }
    
    /// acquire binomial coefficients for summation of transition probabilities
    cerr << "producing binomial coefficient matrices\n" ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        markov_chain_information.at(m).compute_binomial_coefficients() ;
        markov_chain_information.at(m).resize_transition_matrix() ;
    }
    
    /// estimate time of admixture using golden section search
    if ( options.t_fixed == 0 ) {
        cerr << "estimating admixture time based on golden section search\n" ;
        options.t_fixed = golden_search( options, markov_chain_information ) ;
        cerr << "\testimated_admixture time:\t" << options.t_fixed << endl ;
        cerr << "\testimated ancestry:\t" << options.ancestry_proportion << endl ;
    }
    else {
        cerr << "admixture time assumed to be:\t" << options.t_fixed << endl ;
    }
    
    /// compute foward probabilities
    cerr << "computing forward probabilities\n" ;
    long double log_normalization = 0 ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        markov_chain_information.at(m).compute_transitions_probabilities( options.t_fixed, options );
        log_normalization += markov_chain_information.at(m).compute_forward() ;
    }
    cerr << "\tlnL:\t" << log_normalization << endl ;
    
    /// compute backwards probabilities
    cerr << "computing backwards probabilities\n" ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        markov_chain_information.at(m).compute_backward( ) ;
    }

    // now combine the two to obtain smoothed probabilities
    cerr << "combining probabilities\n" ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        markov_chain_information.at(m).combine_prob( ) ;
    }
	return 0 ; 
}
	

