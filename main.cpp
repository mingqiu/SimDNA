#include <iostream>
#include "./include/Origami.h"

using namespace std;

int main() {
    Origami aOrigami = {"/Users/mingqiuwang/Workspace/DNAOrigami/TestCases/" +
                        inputJsonFileName + '.' + "pairs"};

    vector<pair<ID, ID>> crossovers = aOrigami.identifyDiscontinuity();
    aOrigami.processNodes();
    aOrigami.processCrossovers(crossovers);
    aOrigami.connecting();
    aOrigami.toPDB("/Users/mingqiuwang/Workspace/DNAOrigami/SimDNA/results/Debug/" +
                   inputJsonFileName + '.' + "pdb");
    aOrigami.toXML("/Users/mingqiuwang/Workspace/DNAOrigami/SimDNA/results/Debug/" +
                   inputJsonFileName + '.' + "xml");
//    aOrigami.test();


    return 0;
}