#ifndef __INPUT_LINE_H
#define __INPUT_LINE_H

/// hmm input lines
class input_line {
public:
    int pos ;
    string chrom ;
    vector<vector<double> > reference_counts ;
    vector<vector<double> > sample_counts ;
    double recombination_rate ;
    double error_1 ;
    double error_2 ; 
} ;

#endif
