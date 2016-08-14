#include <iostream>
#include "./include/Origami.h"

using namespace std;

int main() {
    Origami aOrigami = {"/Users/mingqiuwang/Workspace/DNAOrigami/TestCases/" +
                        inputJsonFileName + '.' + "pairs"};

    aOrigami.identifyDiscontinuity();
    aOrigami.processNodes();
    aOrigami.processCrossovers();
    aOrigami.processStackedJuncs();
    aOrigami.connecting();
    aOrigami.vitualSite();

    aOrigami.toPDB("/Users/mingqiuwang/Workspace/DNAOrigami/SimDNA/results/Debug/" +
                   inputJsonFileName + '.' + "pdb");
    aOrigami.toXML("/Users/mingqiuwang/Workspace/DNAOrigami/SimDNA/results/Debug/" +
                   inputJsonFileName + '.' + "xml");


    return 0;
}