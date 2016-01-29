#include <iostream>
#include "./include/Origami.h"

using namespace std;

int main() {
    Origami aOrigami = {"/Users/mingqiuwang/Workspace/DNAOrigami/TestCases/" +
                        inputJsonFileName + '.' + "pairs"};

    vector<pair<ID, ID>> crossovers = aOrigami.makeDiscontinuity();
    aOrigami.processNodes();
    aOrigami.processCrossovers(crossovers);
    aOrigami.connecting();
    aOrigami.toPDB("/Users/mingqiuwang/Workspace/DNAOrigami/SimDNA/results/Debug/" +
                   inputJsonFileName + '.' + "pdb");
    return 0;
}