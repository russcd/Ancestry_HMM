#ifndef __SELECTION_STOCHASTIC_TRAJECTORY_H
#define __SELECTION_STOCHASTIC_TRAJECTORY_H


// generates a frequency trajectory of the selected site using the stochastic method
// the stochastic method will generate many random trajectories and calculate and return the average trajectory
void selection_stochastic_trajectory(vector<double> &trajectory, double s, double m, int generations, int ne, int reps) 
{
    vector<double> traj_sum(generations+1,0); //// +1??????
    default_random_engine rand_gen;

    double fixed_freq;
    double n_sel;
    double n_nonsel;
    int fix_gen = 0;
    double ns_fit;
    double nns_fit;
    double ns_freq;
    double new_ns;

    bool not_lost_fixed = true;
    int rcount = 0;

    
    while (rcount < reps) {
        n_sel = ne * m;
        n_nonsel = ne * (1 - m);
    
        fixed_freq = 0;
        vector<double> traj(generations+1,0);
        traj[0] += m;

        for (int g = 1; g <= generations; g++) { //// not sure about start and end
            ns_fit = n_sel * (1 + s);
            nns_fit = n_nonsel;
            ns_freq = ns_fit / (ns_fit + nns_fit);

            if ( ns_freq == 0.0) {
                not_lost_fixed = false;
                break;
            }

            binomial_distribution<> repopulate(ne, ns_freq);
            new_ns = repopulate(rand_gen);
            traj[g]  = new_ns/ne;

            n_sel = new_ns;
            n_nonsel = ne - new_ns;
        }
        
        if (not_lost_fixed = true) {
            for (int f = 0; f < traj.size(); f++) {
                traj_sum[f] += traj[f];
            }
            rcount++;
        }
        not_lost_fixed = true;
    }
    
    for (int i = 0; i < traj_sum.size(); i++) {
        trajectory.push_back(traj_sum[i]/reps);
    }
    
}

#endif