#ifndef __EXPONENTIATE_MATRIX_H
#define __EXPONENTIATE_MATRIX_H

/// exponentiation by squaring until we find a better (working) solution
mat exp_matrix( mat &matrix, int exp ) {
    if ( exp == 2 ) {
        return ( matrix * matrix ) ;
    }
    else if ( exp == 1 ) {
        return matrix ;
    }
    else if ( exp%2 == 1 ) {
        mat result = exp_matrix( matrix, exp-1 ) ;
        return ( matrix * result ) ;
    }
    else {
        mat result = exp_matrix( matrix, exp/2 ) ;
        return result * result ;
    }
}

#endif
