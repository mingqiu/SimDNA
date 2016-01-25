#include <iostream>
#include "./include/Origami.h"

using namespace std;

int main() {
    Origami aOrigami = {"/Users/mingqiuwang/Workspace/DNAOrigami/JunctionDetect4DNAOrigami/test/" +
                        inputJsonFileName + '.' + "pairs"};

    vector<pair<ID, ID>> crossovers = aOrigami.makeHBPs();
    aOrigami.processNodes();
    aOrigami.processCrossovers(crossovers);
    aOrigami.connecting();
    aOrigami.toPDB("/Users/mingqiuwang/Workspace/DNAOrigami/SimDNA/results/Debug/" +
                   inputJsonFileName + '.' + "pdb");
    return 0;
}