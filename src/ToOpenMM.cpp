//
// Created by Mingqiu Wang on 1/25/16.
//

#include "../include/Origami.h"
using namespace std;

void Origami::toPDB(string str) {

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
            fprintf(pdb, "%5d", connectingHelicalPair);
        }
        fprintf(pdb, "\n");
    }

    fprintf(pdb, "ENDMDL\n");
    fclose(pdb);
}

void Origami::toXML(string str) {

    FILE *xml;
    const char *xmlName = str.c_str();
    if ((xml = fopen(xmlName, "w")) == NULL){ printf("\nerror on open xml file!"); exit(0); }
    fprintf(xml, "<ForceField>\n\n");


    fprintf(xml, " <AtomTypes>\n");
    for (const auto & item : _graph.get_nodes().member())

        fprintf(xml, "  <Type name=\"%d\" class=\"C\" mass=\"%lf\"/>\n", item.second.get_num(), item.second.get_mass());
    fprintf(xml, " </AtomTypes>\n\n");

// Residues
    fprintf(xml, " <Residues>\n");
    fprintf(xml, "  <Residue name=\"ALL\">\n");
    for (auto aHelicalBreakPair = 1; aHelicalBreakPair < _graph.howManyNodes()+1; ++aHelicalBreakPair)
        fprintf(xml, "   <Atom name=\"%d\" type=\"%d\"/>\n", aHelicalBreakPair, aHelicalBreakPair);
    for (const auto &item: _graph.get_edges().member()) {
        fprintf(xml, "   <Bond from=\"%d\" to=\"%d\"/>\n", item.second.get_endsNode().first-1, item.second.get_endsNode().second-1);
    }
    fprintf(xml, "  </Residue>\n");
    fprintf(xml, " </Residues>\n\n");



//
// HarmonicBondForce
    fprintf(xml, " <HarmonicBondForce>\n");
    double length = 0, k_bond;
    Edge edge;
    for (const auto & item :_graph.get_edges().member()) {
        edge = item.second;
        if (edge.length())
            length = RISE_PER_BP*(edge.length()+1);
        else {
            if (edge.is_crossover()) length = CROSSOVER_DIS;
            else length = RISE_PER_BP;
        }
        k_bond = edge.is_ds() ? STRETCH_DS : STRETCH_SS;
        fprintf(xml, "  <Bond type1=\"%d\" type2=\"%d\" length=\"%lf\" k=\"%lf\"/>\n",
                edge.get_endsNode().first, edge.get_endsNode().second, length, k_bond);

    }
    fprintf(xml, " </HarmonicBondForce>\n\n");


// HarmonicAngleForce
    fprintf(xml, " <HarmonicAngleForce>\n");
    int a, b, a1, a2, b1, b2;
    for (const auto item1: _graph.get_HJ()) {
        a = item1.get_endsNode().first;
        b = item1.get_endsNode().second;
        a1 = _graph.connectFrom(a)[0];
        a2 = _graph.connectFrom(a)[1];
        if (a1 == b) a1 = _graph.connectFrom(a)[2];
        else if (a2 == b) a2 = _graph.connectFrom(a)[2];
        b1 = _graph.connectFrom(b)[0];
        b2 = _graph.connectFrom(b)[1];
        if (b1 == a) b1 = _graph.connectFrom(b)[2];
        else if (b2 == a) b2 = _graph.connectFrom(b)[2];
        fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
                a1, a, a2, 3.14, ANGLE_S);
        fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
                b1, b, b2, 3.14, ANGLE_S);
        fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
                b1, b, a, 1.57, ANGLE_S);
        fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
                b2, b, a, 1.57, ANGLE_S);
        fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
                a1, a, b, 1.57, ANGLE_S);
        fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
                a2, a, b, 1.57, ANGLE_S);

    }

// Dihedral Force
//    fprintf(xml, " <PeriodicTorsionForce>\n");
//    for (const auto item1: _graph.get_HJ()) {
//        a = item1.get_endsNode().first;
//        b = item1.get_endsNode().second;
//        a1 = _graph.connectFrom(a)[0];
//        a2 = _graph.connectFrom(a)[1];
//        if (a1 == b) a1 = _graph.connectFrom(a)[2];
//        else if (a2 == b) a2 = _graph.connectFrom(a)[2];
//        b1 = _graph.connectFrom(b)[0];
//        b2 = _graph.connectFrom(b)[1];
//        if (b1 == a) b1 = _graph.connectFrom(b)[2];
//        else if (b2 == a) b2 = _graph.connectFrom(b)[2];
//        fprintf(xml, "  <Proper type1=\"%d\" type2=\"%d\" type3=\"%d\" type4=\"%d\" periodicity1=\"1\" phase1=\"0.0\""
//                        " k=\"%lf\"/>\n", a1, a, b, b1, 3.14);
//
//        fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
//                b1, b, b2, 3.14, ANGLE_S);
//        fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
//                b1, b, a, 1.57, ANGLE_S);
//        fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
//                b1, b, a, 1.57, ANGLE_S);
//
//    }
//    fprintf(xml, " </PeriodicTorsionForce>\n\n");

//    for (const auto & item2 : _graph.get_nodes().member()) {
//        if (item2.second.)
//    }



    int i;
    Edge edge1, edge2;
    int edge1index, edge2index;
    for (const auto &item2 : _graph.get_nodes().member()) {
        if (item2.second.is_inHJ()) continue;
        if (item2.second.get_type()==3) continue;
        i = item2.second.get_num();
        if (_graph.showGraph()[i].size()==1) continue;


        for (int a = 0; a < _graph.showGraph()[i].size(); ++a)
            for (int b = a; b < _graph.showGraph()[i].size(); ++b) {
                if (a==b) continue;
                edge1index = _graph.showGraph()[i][a];
                edge2index = _graph.showGraph()[i][b];
                edge1 = _graph.findEdgeFromEnds(i, edge1index);
                edge2 = _graph.findEdgeFromEnds(i, edge2index);
                if (edge1.is_ds()&&edge2.is_ds())
                    fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
                            edge1index, i, edge2index, 3.14, ANGLE_W);
                else if ((edge1.is_crossover()&&edge2.is_ds())||(edge2.is_crossover()&&edge1.is_ds()))
                    fprintf(xml, "  <Angle type1=\"%d\" type2=\"%d\" type3=\"%d\" angle=\"%lf\" k=\"%lf\"/>\n",
                            edge1index, i, edge2index, 1.57, ANGLE_W);
            }
    }

    fprintf(xml, " </HarmonicAngleForce>\n\n");



    fprintf(xml, "</ForceField>\n\n");



    fclose(xml);

}