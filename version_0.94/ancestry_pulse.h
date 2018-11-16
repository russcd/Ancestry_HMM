#ifndef __ANCESTRY_PULSE_H
#define __ANCESTRY_PULSE_H

/// ancestry pulses
class pulse {
    
public:
    
    /// time of pulse
    double time ;
    
    /// is time fixed?
    bool time_fixed ;
    
    /// ancestry type that pulses
    double type ;
    
    /// proportion of the individuals replaced with this pulse
    double proportion ;
    
    /// is the proportion fixed
    bool proportion_fixed ;
    
    /// if there are multiple pulses of the same ancestry type, proportion that remains that this pulse will take
    double fraction_of_remainder ;
    
    /// order pulses were entered in
    int entry_order ;
    
    /// print a pulse
    void print () ;
    
    /// sort pulses by time
    friend bool operator < ( const pulse &a, const pulse &b ) {
        return a.time > b.time ;
    }
} ;

//// sort vector by type
void sort_pulse_vector ( vector<pulse> &ancestry_pulses, int ancestry_types ) {
    vector<pulse> return_pulses ;
    for ( int a = 0 ; a < ancestry_types ; a++ ) {
        for ( int p = 0 ; p < ancestry_pulses.size() ; p ++ ) {
            if ( ancestry_pulses[p].type == a ) {
                return_pulses.push_back( ancestry_pulses[p] ) ;
            }
        }
    }
    ancestry_pulses = return_pulses ;
}

/// print pulse information
void pulse::print () {
    cerr << "\t\t" << time << "\t" << type << "\t" << fraction_of_remainder << endl ;
}

#endif
