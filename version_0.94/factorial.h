#ifndef __FACTORIAL_H
#define __FACTORIAL_H

vector<double> create_factorial () {

	vector<double> factorial(1755) ;
	factorial[0] = 1 ; 
	double fact = 1 ; 
	for ( double base = 1 ; base < 1755 ; base ++ ) { 
		fact *= base ; 
		factorial[base] = fact ; 
	}
	return factorial ; 
}

const vector<double> factorial = create_factorial() ;

#endif
