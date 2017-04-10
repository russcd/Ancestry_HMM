#ifndef __NORMALIZE_H
#define __NORMALIZE_H

/// normalize vector
double normalize( vector<double> &vec ) {
    double sum = 0 ;
    for ( int i = 0 ; i < vec.size() ; i ++ ) {
        sum += vec.at(i) ;
    }
    for ( int i = 0 ; i < vec.size() ; i ++ ) {
        vec.at(i) /= sum ;
    }
    return log(sum) ;
}

/// normalize vector
double normalize( vec &vector ) {
    double sum = accu( vector ) ;
    for ( int i = 0 ; i < vector.size() ; i ++ ) {
        vector(i) /= sum ;
    }
    return log(sum) ;
}

#endif
