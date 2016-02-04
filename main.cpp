#include <iostream>
#include "./include/Origami.h"

using namespace std;

int main() {
    Origami aOrigami = {"/Users/mingqiuwang/Workspace/DNAOrigami/TestCases/" +
                        inputJsonFileName + '.' + "pairs"};

    aOrigami.identifyDiscontinuity();
    aOrigami.processNodes();
    aOrigami.processCrossovers();
    aOrigami.connecting();
    aOrigami.toPDB("/Users/mingqiuwang/Workspace/DNAOrigami/SimDNA/results/Debug/" +
                   inputJsonFileName + '.' + "pdb");
    aOrigami.toXML("/Users/mingqiuwang/Workspace/DNAOrigami/SimDNA/results/Debug/" +
                   inputJsonFileName + '.' + "xml");
//    aOrigami.test();


    return 0;
}