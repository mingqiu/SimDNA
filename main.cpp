#include <iostream>
#include "./include/Origami.h"

using namespace std;

int main() {
//    Origami aOrigami = {"/Users/mingqiuwang/Workspace/DNAOrigami/TestCases/" +
//                        inputJsonFileName + '.' + "pairs"};
//    Origami aOrigami = {"/Users/mingqiuwang/Workspace/DNAOrigami/TestCases/rectangle_3x3_4x4_5x5_for_simulation/" + inputJsonFileName + '.' + "pairs"};
//    Origami aOrigami = {"/Users/mingqiuwang/Workspace/DNAOrigami/TestCases/results_20160208/" +
//                        inputJsonFileName + '.' + "pairs"};

    Origami aOrigami = {"/Users/mingqiuwang/Workspace/DNAOrigami/TestCases/sim2/" + inputJsonFileName + '.' + "pairs"};


    aOrigami.identifyDiscontinuity();
    aOrigami.processNodes();
    aOrigami.processCrossovers();
    aOrigami.processStackedJuncs();
    aOrigami.connecting();

    aOrigami.toPDB("/Users/mingqiuwang/Workspace/DNAOrigami/SimDNA/results/Debug/" +
                   inputJsonFileName + '.' + "pdb");
    aOrigami.toXML("/Users/mingqiuwang/Workspace/DNAOrigami/SimDNA/results/Debug/" +
                   inputJsonFileName + '.' + "xml");
    aOrigami.test();


    return 0;
}