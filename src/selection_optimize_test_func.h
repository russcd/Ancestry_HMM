#ifndef __OPTIMIZE_TEST_FUNC_H
#define __OPTIMIZE_TEST_FUNC_H

double optimize_test_func(int x, double y) {
    return 2000-1000*(pow(sin(x/1000), 10) + cos(10 + y * x/1000) * cos(x/1000));
    //2000-1000*(np.sin(x/1000) ** 10 + np.cos(10 + y * x/1000) * np.cos(x/1000))

}
double optimize_test_func2(int x, double y) {
    double xx = x;
    return 1000 - 50*(pow((xx/1000)-3,2) + pow(y-3,2));
}

double selection_evaluate_point2(selection &point, vector<markov_chain> &markov_chain_information, map<int, vector<vector< map< vector<transition_information>, double > > > > &transition_matrix_information, vector<double> &recombination_rate, vector<int> &position, cmd_line &options, map<int,vector<vector<int> > > &state_changes ) {
    //cout << "Evaluate point1: " << point<< endl;
    point.lnl = optimize_test_func2(point.pos, point.sel);
    //cout << "Evaluate point2: " << point<< endl;
    return point.lnl;
}

#endif