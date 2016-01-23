//
// Created by Mingqiu Wang on 1/22/16.
//

#include "../include/DataStructure.h"
using namespace std;

std::ostream& operator<<(std::ostream& os, const ID & idnum) {
    os << idnum.strandID() << "\t" << idnum.baseID();
    return os;
}