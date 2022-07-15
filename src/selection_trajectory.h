#ifndef __SELECTION_TRAJECTORY_H
#define __SELECTION_TRAJECTORY_H


// generates vector with allele frequency of selected allele over time/generations
void selection_trajectory(vector<double> &freq, double s, int tt, double m, int generations, int n) 
{
    // returns flat vector if selection is 0 (ie, no change in ellele frequency over time)
    if ( s == 0) {
        freq.assign(generations,m);
        return;
    }

    int t0 = 0;
    int t0min = 0;
    bool found = false;
    double f;

    // loops over generations until initial frequency (m) is reached + number of generations (gen) has passed
    while (found == false) {
        f = 1 / (1 + 2 * n * s * exp(-s*t0));
        if (f > m) {
            if (tt == 0) {
                tt = t0;
            }
            else if (t0 == (tt + generations)) {
                found = true;
            }
            freq.push_back(f);
        }
        t0++;
    }

}   

// checks if selective coeffient causes site to go to fixation in the time since introgression
bool selection_reaches_fixation(double s, double m, int generations, int n) 
{
    double max_freq = 0.99;
    int tt = 0;
    int t0 = 0;
    int t0min = 0;
    bool found = false;
    double f;

    s *= 0.5;

    // loops over generations until initial frequency (m) is reached + number of generations (gen) has passed
    while (found == false) {
        f = 1 / (1 + 2 * n * s * exp(-s*t0));
        if (f > max_freq) {
            cerr << "Frequency " << f << " generation " << t0-tt << " selection " << s <<endl;
            found = true;
            break;
        }
        else if (f > m) {
            if (tt == 0) {
                tt = t0;
            }
            else if (t0 == (tt + generations)) {
                break;
            }
        }
        t0++;
    }
    return found;
}

// returns selection coeffient that reaches 0.99 in the time since introgression
// used to speed up calculations and prevent division errors
double selection_get_max_sel(double min_s, double max_s, double step_s, double m, int generations, int n)
{
    cerr << "selection_get_max_sel " << min_s << " " << max_s << " " << step_s << " " << m << " " << generations << " " << n << endl;
    double last_s = 0;
    for (double s = min_s; s <= max_s; s += step_s) {
        if (selection_reaches_fixation(s, m, generations, n) == true) {
            return last_s;
        }
        last_s = s;
    }
    return max_s;
}

#endif
