#include <iostream>
#include "./include/Origami.h"

using namespace std;

int main() {
    Origami aOrigami = {"/Users/mingqiuwang/Workspace/DNAOrigami/JunctionDetect4DNAOrigami/test/" +
                        inputJsonFileName + '.' + "pairs"};
    vector<pair<ID, ID>> crossovers = aOrigami.makeHBPs();
//    for (const auto& item : crossovers)
//        if (item.first.baseID()==2501||item.second.baseID()==2501)
//            cout << item.first << "\t" << item.second << endl;
    aOrigami.processNodes();
    aOrigami.processCrossovers(crossovers);
    return 0;
}