#ifndef __FACTORIAL_H
#define __FACTORIAL_H

#include <vector> 

using namespace std ; 

vector<long double> create_factorial () {

	vector<long double> factorial(1755) ;
	factorial[0] = 1 ; 
	long double fact = 1 ; 
	for ( long double base = 1 ; base < 1755 ; base ++ ) { 
		fact *= base ; 
		factorial[base] = fact ; 
	}
	return factorial ; 
}

#endif
