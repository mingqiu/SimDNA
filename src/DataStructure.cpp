//
// Created by Mingqiu Wang on 1/22/16.
//

#include "../include/DataStructure.h"
using namespace std;

std::ostream& operator<<(std::ostream& os, const ID & idnum) {
    os << idnum.strandID() << "\t" << idnum.baseID();
    return os;
}

std::ostream &operator<<(std::ostream &os, const AxialDiscontinuity &hbp) {
    os << hbp.get_ids().first << "\t" << hbp.get_ids().second << "\t"
    << (hbp.is_typeA()&&hbp.is_typeB()) << "\t" << hbp.is_single() << endl;
    return os;

}


void OrigamiGraph::insertEdge(Edge eg, int c1, int c2) {
    if (!findEdge(c1, c2)) {
        _edges.insert(eg);
        _origamiGraph[c1].push_back(c2);
        _origamiGraph[c2].push_back(c1);
    }
}

bool OrigamiGraph::findEdge(int id1, int id2) const {
    for (const auto &it : _origamiGraph[id1]) {
        if (it==id2) { return true;}
    }
    return false;
}

int OrigamiGraph::findNodeType(ID id) {
    return _nodes.findTypeFromID(id).get_type();
}

const vector<pair<ID, ID>> &OrigamiGraph::findIDsfromIndex(int a) {
    return _nodes.findTypeFromIndex(a).get_ids();
}

int OrigamiGraph::findNodeNum(ID id) {
    return _nodes.findIndexFromID(id);
}

int OrigamiGraph::findNodeType(int a) {
    return _nodes.findTypeFromIndex(a).get_type();
}

void OrigamiGraph::findHollidayJ() {

    Edge edge;
    int sum = 0;
    int a, b;
    for (auto &&item : _edges.updateMember()) {
        edge = item.second;
        if (edge.is_crossover()) {
            _crossovers.push_back(item.second);
            if (edge.get_types().first == 1 && edge.get_types().second == 1) {
                ++sum;
                item.second.set_isHJ(true);
                a = edge.get_endsNode().first;
                b = edge.get_endsNode().second;
                _nodes.findTypeFromIndex(a).set_inHJ(true);
                _nodes.findTypeFromIndex(b).set_inHJ(true);
                _HJ.push_back(item.second);
            }
        }
    }
    _numHJ = sum;

}