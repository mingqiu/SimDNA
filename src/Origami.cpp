//
// Created by Mingqiu Wang on 1/22/16.
//

#include "../include/Origami.h"

using namespace std;

Origami::Origami(std::string s) : _fileName(s) {
        std::cout << "Step 1: Input file reading \n";
        if (!input()) {
            std::cout << "Error occured in input file\n" << s <<
            "\nNo file found" << std::endl;
            exit(1);
        };
        std::cout << "................ DONE\n";
}

bool Origami::input() {
    ifstream input{_fileName};
    if (!input.is_open()) {
        return false;
    }

    int strand, base; double x, y, z; int pairStrand, pairBase;
    int strandOld, baseOld; double xOld, yOld, zOld; int pairStrandOld, pairBaseOld;
    bool isABreak;

    input >> _totalResidueNum >> strandOld
        >> baseOld >> xOld >> yOld >> zOld >> pairStrandOld >> pairBaseOld;
//    if (baseOld!=-1)strandOld = 0; if (pairBaseOld!=-1)pairStrandOld = 0;
//    if (baseOld!=-1)++baseOld; if (pairBaseOld!=-1)++pairBaseOld;
    // read two lines(residues), and judge the first line (old residue)
    for (int i = 1; i < _totalResidueNum; ++i) {
        input >> strand >> base >> x >> y >> z >> pairStrand >> pairBase;

//        if (base!=-1)strand = 0; if (pairBase!=-1)pairStrand = 0;

//        if (base!=-1)++base; if (pairBase!=-1)++pairBase;
        isABreak = isDiscontinuity(strand, base, pairStrand, pairBase,
                                   strandOld, baseOld, pairStrandOld, pairBaseOld);

        _nucleotide[ ID{strandOld, baseOld} ] =
                Nucleotide{ ID{strandOld, baseOld}, Vector3Dd{xOld, yOld, zOld},
                            ID{pairStrandOld, pairBaseOld}, isABreak};
        if (isABreak) makeHB(strandOld, baseOld);

        if (_resNumInEachStrand.find(strandOld) == _resNumInEachStrand.end())
            _resNumInEachStrand[strandOld] = 1;
        else
            ++_resNumInEachStrand[strandOld];

        strandOld = strand;
        baseOld = base;
        xOld = x;
        yOld = y;
        zOld = z;
        pairStrandOld = pairStrand;
        pairBaseOld = pairBase;
    }

    // assign last residue
    _nucleotide[ ID{strand, base} ] =
            Nucleotide{ ID{strand, base}, Vector3Dd{x, y, z}, ID{pairStrand, pairBase}, true};
    makeHB(strand, base);
    ++_resNumInEachStrand[strand];

    _strandNum = _resNumInEachStrand.size();

    input.close();
//    testInput();
    return true;
}

void Origami::identifyDiscontinuity() {
    std::cout << "Step 2: Axial discontinuities generating \n";

    ID id, idleft, idright;

    for (int i = 0; i < _strandNum; i++) {
        vector<int> strand = _breaksOnEachStand[i];
        for (auto item = 1; item < strand.size()-1; ++item) {
            id = {i, strand[item]};


//            if (i==0&&strand[item]==11949)
//                cout << endl;

// if belongs to single strand
            if (!_nucleotide[id].withPair()) {
                _axialDiscons.insert({id, {-1, -1}}, false, false, true);
                continue;
            }

            idleft = {i, strand[item-1]};
            idright = {i, strand[item+1]};

// if belongs to double strand
            // if it has a left neighbor discontinuity
            if ((id.baseID() - idleft.baseID()) == 1) {
                double i = dist(helicalCenter(idright), helicalCenter(id));
                if (dist(helicalCenter(idleft), helicalCenter(id)) > THRES_CROSSOVER_DIS) {
                    _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                         false, true, false);
                    continue;
                }
                if ((idright.baseID() - id.baseID()) == 1) {
                    if (dist(helicalCenter(idright), helicalCenter(id)) > THRES_CROSSOVER_DIS) {
                        _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                             false, true, false);
                        _crossovers.push_back({id, idright});
//                        cout << "test" << endl;
                    }
                    else
                        _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                             _nucleotide[idright].withPair(), false, false);
                    continue;
                }
                else _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                          true, false, false);
                continue;
            }

            // if it has a right neighbor discontinuity
            if ((idright.baseID() - id.baseID()) == 1) {
                if (dist(helicalCenter(idright), helicalCenter(id)) > THRES_CROSSOVER_DIS) {
                    _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                         false, true, false);
                    _crossovers.push_back({id, idright});
//                    cout << "test" << endl;
                }
                else
                    _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                         _nucleotide[idright].withPair(), false, false);
                continue;
            }
            _axialDiscons.insert({id, _nucleotide[id].pairID()}, false, false, false);
        }

// start residue in a strand
        id = {i, strand[0]};
        idright = {i, strand[1]};

//        if (i == 42 && strand[0] == 1)
//            cout << endl;
        if (!_nucleotide[id].withPair())
            _axialDiscons.insert({id, {-1, -1}}, false, false, true);

        else if ((idright.baseID() - id.baseID()) == 1) {

            if (dist(helicalCenter(idright), helicalCenter(id)) > THRES_CROSSOVER_DIS) {
                _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                     false, true, false);
                _crossovers.push_back({id, idright});
                if (TEST) cout << "crossover at start" << endl;
            }
            else _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                      _nucleotide[idright].withPair(), false, false);
        }
        else _axialDiscons.insert({id, _nucleotide[id].pairID()}, false, false, false);

// end residue in a strand
        id = {i, strand[strand.size()-1]};
        idleft = {i, strand[strand.size()-2]};

        if (!_nucleotide[id].withPair()) _axialDiscons.insert({id, {-1, -1}}, false, false, true);
        else if ((id.baseID() - idleft.baseID()) == 1) {
            if (dist(helicalCenter(idleft), helicalCenter(id)) > THRES_CROSSOVER_DIS) {
                _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                     false, true, false);
                _crossovers.push_back({idleft, id});
                if (TEST) cout << "crossover at end" << endl;

            }
            else _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                      true, false, false);
        }
        else
            _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                 false, false, false);
    }

//    for (const auto & item : _crossovers)
//        cout << item.first << "\t" << item.second << endl;
    std::cout << "................ DONE\n";


}

void Origami::processNodes() {

    std::cout << "Step 3: Nodes generating \n";

    int hb, hbold;
    AxialDiscontinuity hbp, hbpold;
    for (int i = 0; i < _strandNum; i++) {
        vector<int> strand = _breaksOnEachStand[i];
        hb = strand[0];
        hbp = _axialDiscons.findTypeFromID({i, hb});
        for (int item = 1; item < strand.size(); ++item) {

            hbold = hb;
            hb = strand[item];
            hbpold = hbp;
            hbp = _axialDiscons.findTypeFromID({i, hb});
            if ((hb - hbold)==1) {
                if (AxialDiscontinuity::similar(hbpold, hbp) && hbp.is_typeA() ) {
                    if (dist(helicalCenter({i, hb}), helicalCenter({i, hbold})) > THRES_CROSSOVER_DIS && hbp.is_typeB()) {
                        if (item == strand.size() - 1) break;
                        hb = strand[item + 1];
                        hbp = _axialDiscons.findTypeFromID({i, hb});
                        ++item;
                        continue;
                    }

                    _graph.insertNode(makeNode(hbpold, hbp));
                    if (item == strand.size() - 1) break;
                    hb = strand[item + 1];
                    hbp = _axialDiscons.findTypeFromID({i, hb});
                    ++item;
                    continue;
                }
            }
            _graph.insertNode(makeNode(hbpold));
        }
        _graph.insertNode(makeNode(hbp));
    }
    _graph.resizeGraph();
    std::cout << "................ DONE\n";

}


void Origami::processCrossovers() {
    std::cout << "Step 4: Crossovers generating \n";

    ID id1, id2;
    Edge edge;
    for (const auto &cross : _crossovers) {
        id1 = cross.first;
        id2 = cross.second;
        edge = makeEdgeCrossover(id1, id2);
        _graph.insertEdge(edge, edge.get_endsNode().first, edge.get_endsNode().second);
    }
//    _graph.howManyFourWays();
    std::cout << "................ DONE\n";

}


void Origami::connecting() {
    std::cout << "Step 5: Graph connecting \n";

    int node1, node2;
    ID id1, id2;
    for (int i = 0; i < _strandNum; i++) {
        vector<int> strand = _breaksOnEachStand[i];
        id2 = {i, strand[0]};
        node2 = _graph.findNodeNum(id2);
        for (int item = 1; item < strand.size(); ++item) {
            id1 = id2;
            node1 = node2;
            id2 = {i, strand[item]};
            node2 = _graph.findNodeNum(id2);
//            if (id2.baseID()==1546)
//                cout << endl;
            if (node1 == node2) continue;
//            if (node1==135||node2==135)
//                cout << endl;
            _graph.insertEdge(makeEdge(id1, id2), node1, node2);
        }
    }

    std::cout << "................ DONE\n";


    std::cout << "Step 6: Holliday junctions identifying \n";
    _graph.findHollidayJ();

    std::cout << "................ DONE\n";

}




void Origami::processStackedJuncs() {


    queue<int> tips;

    for (int i = 0; i < _graph.showGraph().size(); ++i)
        if (_graph.showGraph().at(i).size() == 2 && _graph.findNodeType(i) == 1)
            tips.push(i);

    vector<int> visited;
    int a, b;

    while (!tips.empty()) {
        a = tips.front();
        tips.pop();
        if (find(visited, a)) continue;
        queue<int> toTrace;
        vector<int> st;
        toTrace.push(a);

        for (int i = 0; !toTrace.empty() && i < 100; ++i) { // don't search a stack with more than 100 nodes
            b = toTrace.front();
            toTrace.pop();
            st.push_back(b);
            visited.push_back(b);
            for (auto item : _graph.connectFrom(b))
                if (!find(visited, item)) {
                    toTrace.push(item);
                }
        }

        if (st.at(st.size() - 1) == st.at(st.size() - 2)) {
            st.pop_back();
            _stacks.push_back(st);
        }
    }

    pair<int , int> endsNode;
    double dis;
    for (auto item3: _stacks)
        for (int i = 1; i < item3.size() - 1; i += 2) {
            endsNode = {item3.at(i), item3.at(i + 1)};
            dis = dist(_graph.findNodeFromNum(item3.at(i)).get_position(),
                              _graph.findNodeFromNum(item3.at(i+1)).get_position());
            _graph.insertEdge(Edge{endsNode, dis}, item3.at(i), item3.at(i + 1));
//            cout << dist (_graph.findNodeFromNum(item3.at(i)).get_position(),
//                  _graph.findNodeFromNum(item3.at(i+1)).get_position()) << endl;
        }
}
