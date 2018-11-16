#ifndef __BINOMIAL_H
#define __BINOMIAL_H

/// binomial probs
double binomial ( double n, double k, double prob ) {
    return nCk[n][k]*pow((1-prob),n-k)*pow(prob,k) ;
}

/// multinomial probs
double multinomial ( double &nc, double &n, vector<int> &k, vector<int> &prob ) {
    double mult = factorial[n] ;
    for ( int x = 0 ; x < k.size() ; x ++ ) {
        mult *= pow( (double)prob.at(x)/nc, k.at(x) )/factorial[k.at(x)] ;
    }
    return mult ;
}

#endif
