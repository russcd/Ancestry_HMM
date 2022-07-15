#ifndef __SELECTION_CLASS_H
#define __SELECTION_CLASS_H

class selection {
public:
    int pos;
    double sel;
    double lnl;

    /// sort pulses by time
    friend bool operator < ( const selection &a, const selection &b ) {
        return a.lnl < b.lnl ;
    }
} ;

ostream& operator<< (ostream &out, selection const& point) {
    out << "Selection point. pos:" << point.pos << "  sel:" << setprecision(15) << point.sel << "  lnL: "  << point.lnl;
    return out;
}

#endif
