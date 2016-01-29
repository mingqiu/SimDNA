//
// Created by Mingqiu Wang on 1/22/16.
//


#include "../include/Origami.h"
using namespace std;

void Origami::testInput() {

    std::ofstream myfile;
    myfile.open ("test.txt");
    for (auto strand = 0; strand < _strandNum; ++strand)
        for (auto base = _resNumInEachStrand[strand]; base >=1 ; --base) {
            auto a = _nucleotide[{strand, base}];
            if (a.isDiscont()) myfile << a.id() << endl;
        }

    myfile.close();
}

void Origami::testhbpAssign() {
    int hb, hbold, hbp;
    ID id, idleft, idright;

    std::ofstream myfile;
    myfile.open ("test2.txt");
    for (int i = 0; i < _strandNum; i++) {
        vector<int> strand = _breaksOnEachStand[i];

        for (int item = strand.size()-1; item >=0; --item) {
            id = {i,strand[item]};
            hbp = _axialDiscons.findIndexFromID(id);
            myfile << id << "\t" << _nucleotide[id].pairID() << "\t" <<
            (_axialDiscons.findTypeFromIndex(hbp).is_typeA() &&
             _axialDiscons.findTypeFromIndex(hbp).is_typeB()) << endl;
        }

    }
    myfile.close();

}
