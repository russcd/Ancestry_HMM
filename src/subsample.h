#ifndef __SUBSAMPLE_H
#define __SUBSAMPLE_H

void subsample_reads ( double &c1, double &c2, int sample_max ) {
    while ( c1 + c2 > sample_max ) {
        double r = ((double) rand() / (RAND_MAX)) ;
        if ( r < c1/(c1+c2) ) {
            c1 -- ;
        }
        else {
            c2 -- ;
        }
    }
}

#endif
