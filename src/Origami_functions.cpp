//
// Created by Mingqiu Wang on 1/23/16.
//

#include "../include/Origami.h"


using namespace std;

bool Origami::isDiscontinuity(int strand, int base, int pairStrand, int pairBase,
                              int strandOld, int baseOld, int pairStrandOld, int pairBaseOld) {

    if (base == 2 || base == 1) return true; // starting and ending residues are breaks

    // with the judge above, below baseOld is at least 2.
    ID oldoldCom = _nucleotide[ID{strandOld, baseOld - 1}].pairID();
    if (pairStrand != pairStrandOld) return true; // changed strand
    if (oldoldCom.strandID() != pairStrandOld) return true; // changed strand

    // if double strand
    if (pairBaseOld != -1) {
        if (!((pairBase - pairBaseOld) == 1 || (pairBase - pairBaseOld) == -1)) return true; // changed base
        if (!((oldoldCom.baseID() - pairBaseOld) == 1 || (oldoldCom.baseID() - pairBaseOld) == -1)) return true; // changed base

    }
        // if single strand
    else {
        if (_nucleotide[ID{strandOld, baseOld - 1}].withPair()) return true;
        if (pairBase != -1 ) return true;

    }
    return false;
}


void Origami::makeHB(int strand, int base) {
    if (_breaksOnEachStand.find(strand) == _breaksOnEachStand.end()) {
        vector<int> thisStrand;
        thisStrand.push_back(base);
        _breaksOnEachStand[strand] = thisStrand;
    }
    else
        _breaksOnEachStand[strand].push_back(base);
}



Vector3Dd Origami::helicalCenter(ID id1) {
    ID id2 = _nucleotide[id1].pairID();

    if (id1.baseID()!=-1&&id2.baseID()!=-1)
        return (_nucleotide[id1].coordinate()+_nucleotide[id2].coordinate())*0.5;
    else if (id1.baseID()==-1)
        return (_nucleotide[id2].coordinate());
    else if (id2.baseID()==-1)
        return (_nucleotide[id1].coordinate());
}




Edge Origami::makeEdgeCrossover(ID id1, ID id2) {
    int a = _graph.findNodeNum(id1);
    int b = _graph.findNodeNum(id2);
    int c = _graph.findNodeType(id1);
    int d = _graph.findNodeType(id2);
    return Edge({a, b}, {c, d}, {id1, id2}, true, {});
}

Node Origami::makeNode(AxialDiscontinuity pair1, AxialDiscontinuity pair2) {
    int type = 1;
    std::vector<std::pair<ID,ID>> ids;
    ids.push_back(pair1.get_ids());
    ids.push_back(pair2.get_ids());
    double mass = 2 * MASS;
    Vector3Dd position = (helicalCenter(pair1.get_ids().first)
                          + helicalCenter(pair2.get_ids().first)) * 0.5;
    double vdWradii =VDWRADII_1;

    return Node{type, ids, mass, position, vdWradii};

}

Node Origami::makeNode(AxialDiscontinuity pair) {
    int type = pair.is_single() ? 3 : 2;
    std::vector<std::pair<ID,ID>> ids;

    ids.push_back(pair.get_ids());
//    if (ids.at(0).first.baseID()==9467)
//        cout << endl;
    double mass = 1.0/(type - 1) * MASS;
    Vector3Dd position = helicalCenter(pair.get_ids().first);
    double vdWradii = pair.is_single() ? VDWRADII_3 : VDWRADII_2;

    return Node{type, ids, mass, position, vdWradii};
}


Edge Origami::makeEdge(ID id1, ID id2) {


    int a = _graph.findNodeNum(id1);
    int b = _graph.findNodeNum(id2);
    int c = _graph.findNodeType(id1);
    int d = _graph.findNodeType(id2);

    vector<pair<ID, ID>> ids;

    int strand = id1.strandID();
    if (id1.baseID() < id2.baseID())
        for (int item = id1.baseID()+1; item < id2.baseID(); ++item)
            ids.push_back({{strand, item}, _nucleotide[{strand, item}].pairID()});
    else
        for (int item = id2.baseID()+1; item < id1.baseID(); ++item)
            ids.push_back({{strand, item}, _nucleotide[{strand, item}].pairID()});

    return Edge({a, b}, {c, d}, {id1, id2}, false, ids);
}

