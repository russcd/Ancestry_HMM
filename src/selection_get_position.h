#ifndef __GET_POSITION_H
#define __GET_POSITION_H

int get_position(int pos, vector<int> &positions) {
    for (int i = 0; i < positions.size(); i++) {
        if ( positions[i] >= pos ) {
            return i;
        }
    }
    return -1;
}

#endif