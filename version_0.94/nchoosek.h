#ifndef __NCHOOSEK_H
#define __NCHOOSEK_H

/// this will just create a matrix through which we can look up any nck within the range that can be represented using a double
const int max_n = 1001 ;

vector<vector<double> > create_nck_table () {
    vector<vector<double> > nck(max_n,vector<double>(max_n) ) ;
    for ( double n = 0 ; n < max_n ; n ++ ) {
        for ( double k = 0 ; k <= n ; k ++ ) {
         
            if ( n == k || k == 0 ) {
                nck[n][k] = 1 ;
                continue ;
            }
            
            double l = k ;
            if ( n-k < k ) {
                l = n - k ;
            }
            
            double lf = 2 ;
            double nf = n - l + 1 ;
            double result = 1 ;
    
            while ( lf < l + 1 || nf < n + 1 ) {
                
                /// check if >= 1 because we want to stay in the dynamic range of the double as long as possible to retain precision and make computation possible
                if ( ( result >= 1 || nf == n + 1 ) && lf < l + 1 ) {
                    result /= lf ;
                    lf ++ ;
                }
                else {
                    result *= nf ;
                    nf ++ ;
                }
            }
            nck[n][k] = result ;
        }
    }
    
    return nck ;
}

//// function to create the table is called automatically such that this matrix is available as long as the header is included
const vector<vector<double> > nCk = create_nck_table() ;

#endif
