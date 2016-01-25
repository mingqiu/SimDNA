//
// Created by Mingqiu Wang on 1/25/16.
//

#include "../include/Origami.h"

void Origami::toPDB(std::string str) {

    FILE *pdb;
    const char *pdbName = str.c_str();
    if ((pdb = fopen(pdbName, "w")) == NULL){ printf("\nerror on open pdb file!"); exit(0); }
    fprintf(pdb, "MODEL     1\n");

    Vector3Dd coordinate;
    std::vector<std::vector<int>> origamiAdjacencyList = _graph.showGraph();

    for (auto aHelicalBreakPair = 1; aHelicalBreakPair < _graph.howManyNodes()+1; ++aHelicalBreakPair) {
        coordinate = _graph.nodeCenter(aHelicalBreakPair);
        fprintf(pdb, "ATOM  %5d %4d ALL     1    %8.2lf%8.2lf%8.2lf\n",
                aHelicalBreakPair, aHelicalBreakPair, coordinate.x(), coordinate.y(), coordinate.z());
    }

    for (auto aHelicalBreakPair = 1; aHelicalBreakPair < origamiAdjacencyList.size(); ++aHelicalBreakPair) {
        fprintf(pdb, "CONECT%5d", aHelicalBreakPair);
        for (auto connectingHelicalPair : origamiAdjacencyList[aHelicalBreakPair]) {
            fprintf(pdb, "%5d", (int) connectingHelicalPair);
        }
        fprintf(pdb, "\n");
    }

    fprintf(pdb, "ENDMDL\n");
    fclose(pdb);
}