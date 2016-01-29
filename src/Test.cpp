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
