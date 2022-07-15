#ifndef __SPLIT_VECTOR_H
#define __SPLIT_VECTOR_H

/// splits vector (of recombination rates) into two vectors at the selected site.
/// the back vector is generated in reverse order

// min function because of namespace collision between std and arma
int int_min(int a, int b){
    int minout;
    if (a < b) {
        minout = a;
    }
    else {
        minout = b;
    }
    return minout;
}

// splits a vector of chromosomal positions into two vectors going away from a focal site.
// one vector is reversed going back to 0
// also trims the hmm window used
vector <vector<double> > split_vector(int sel_site, vector<double> &whole_vec, cmd_line &options) 
{
    vector <vector<double> > split_vecs;
    vector<double> fwd_vec;
    vector<double> back_vec;

    // trim vector if size is specified in morgans
    if (options.win_unit == "m") {
        double sum_morgans_fwd;
        double sum_morgans_back;
        for (int i = sel_site; i < whole_vec.size(); i++) {
            sum_morgans_fwd += whole_vec[i];
            if (sum_morgans_fwd <= options.win_morgan) {
                fwd_vec.push_back(whole_vec[i]) ;
            }
            else {
                break;
            }
        }
        for (int i = sel_site; i > 0; i--) {
            sum_morgans_back += whole_vec[i];
            if (sum_morgans_back <= options.win_morgan) {
                back_vec.push_back(whole_vec[i]) ;
            }
            else {
                break;
            }
        }
    }

    // trim vector if size is specified in percent
    else if (options.win_unit == "p") {
        int percent_size = whole_vec.size() * (options.win_percent/100);
        int trim_size_fwd = int_min((whole_vec.size() - sel_site) , percent_size); // Check for off by 1 error
        int trim_size_back = int_min(sel_site, percent_size);

        for (int i = sel_site; i < (sel_site + trim_size_fwd); i++) {
            fwd_vec.push_back(whole_vec[i]) ;
        }

        for (int i = sel_site; i > (sel_site - trim_size_back); i--) {
            back_vec.push_back(whole_vec[i]) ;
        }

    }

    split_vecs.push_back(fwd_vec);
    split_vecs.push_back(back_vec);
    return split_vecs;
    //cout << "Split vector lengths: " << fwd_vec.size() << ", " << back_vec.size() << endl;
}



#endif