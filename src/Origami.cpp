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

    // read two lines(residues), and judge the first line (old residue)
    for (int i = 1; i < _totalResidueNum; ++i) {
        input >> strand >> base >> x >> y >> z >> pairStrand >> pairBase;

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

vector<pair<ID, ID>> Origami::makeDiscontinuity() {
    std::cout << "Step 2: Axial discountinuities generating \n";

    ID id, idleft, idright;
    vector<pair<ID, ID>> crossovers;

    for (int i = 0; i < _strandNum; i++) {
        vector<int> strand = _breaksOnEachStand[i];
        for (auto item = 1; item < strand.size()-1; ++item) {
            id = {i, strand[item]};

            // if belongs to single strand
            if (!_nucleotide[id].withPair()) {
                _axialDiscons.insert({id, {-1, -1}}, false, false, true);
                continue;
            }

            idleft = {i, strand[item-1]};
            idright = {i, strand[item+1]};

            if ((id.baseID() - idleft.baseID()) == 1) {
                if (dist(helicalCenter(idleft), helicalCenter(id))>CROSSOVER_DIS) {
                    _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                         false, true, false);
                    crossovers.push_back({idleft, id});
                    continue;
                }

                if ((idright.baseID() - id.baseID()) == 1) {
                    if (dist(helicalCenter(idright), helicalCenter(id))>CROSSOVER_DIS) {
                        _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                             false, true, false);
                        crossovers.push_back({idright, id});

                    }
                    else {
                        _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                             _nucleotide[idright].withPair(), false, false);

                    }
                    continue;
                }


                else _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                          true, false, false);
                continue;
            }
            if ((idright.baseID() - id.baseID()) == 1) {
                if (dist(helicalCenter(idright), helicalCenter(id))>10) {
                    _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                         false, true, false);
                    crossovers.push_back({idright, id});

                }
                else {
                    _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                         _nucleotide[idright].withPair(), false, false);

                }
                continue;
            }
            _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                 false, false, false);

        }

        id = {i, strand[0]};
        idright = {i, strand[1]};

        if (!_nucleotide[id].withPair()) {
            _axialDiscons.insert({id, {-1, -1}}, false, false, true);
        }

        else if ((idright.baseID() - id.baseID()) == 1) {
            if (dist(helicalCenter(idright), helicalCenter(id))>10) {
                _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                     false, true, false);
                crossovers.push_back({idright, id});

            }

            else _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                      _nucleotide[idright].withPair(), false, false);
        }
        else
        _axialDiscons.insert({id, _nucleotide[id].pairID()},
                             false, false, false);


        id = {i, strand[strand.size()-1]};
        idleft = {i, strand[strand.size()-2]};

        if (!_nucleotide[id].withPair()) {
            _axialDiscons.insert({id, {-1, -1}}, false, false, true);
        }

        else if ((id.baseID() - idleft.baseID()) == 1) {
            if (dist(helicalCenter(idleft), helicalCenter(id))>10) {
                _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                     false, true, false);
                crossovers.push_back({idleft, id});

            }
            else _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                      true, false, false);
        }
        else
            _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                 false, false, false);
    }

    std::cout << "................ DONE\n";

    return crossovers;

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
                    if (dist(helicalCenter({i, hb}), helicalCenter({i, hbold}))>10&& hbp.is_typeB()) {
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


void Origami::testhbpAssign() {
    int hb, hbold, hbp;
    ID id, idleft, idright;

    std::ofstream myfile;
    myfile.open ("test2.txt");
    for (int i = 0; i < _strandNum; i++) {
        vector<int> strand = _breaksOnEachStand[i];

        for (int item = strand.size()-1; item >=0; --item) {
            id = {i,strand[item]};
            hbp = _axialDiscons.findIndexFromID(id);
            myfile << id << "\t" << _nucleotide[id].pairID() << "\t" <<
            (_axialDiscons.findTypeFromIndex(hbp).is_typeA() &&
             _axialDiscons.findTypeFromIndex(hbp).is_typeB()) << endl;
        }

    }
    myfile.close();

}

Node Origami::makeNode(AxialDiscontinuity pair1, AxialDiscontinuity pair2) {
    int type = 1;
    std::vector<std::pair<ID,ID>> ids;
    ids.push_back(pair1.get_ids());
    ids.push_back(pair2.get_ids());
    double mass = 2 * MASS;
    Vector3Dd position = (helicalCenter(pair1.get_ids().first)
                         + helicalCenter(pair2.get_ids().first)) * 0.5;
    double vdWradii = 1;

    return Node{type, ids, mass, position, vdWradii};

}

Node Origami::makeNode(AxialDiscontinuity pair) {
    int type = pair.is_single() ? 3 : 2;
    std::vector<std::pair<ID,ID>> ids;
    ids.push_back(pair.get_ids());

    double mass = 1.0/(type - 1) * MASS;
    Vector3Dd position = helicalCenter(pair.get_ids().first);
    double vdWradii = 0.2;

    return Node{type, ids, mass, position, vdWradii};
}




Edge Origami::makeEdgeCrossover(ID id1, ID id2) {
    int a = _graph.findNodeNum(id1);
    int b = _graph.findNodeNum(id2);
    int c = _graph.findNodeType(id1);
    int d = _graph.findNodeType(id2);
    return Edge({a, b}, {c, d}, {id1, id2}, true, {});
}

void Origami::processCrossovers(std::vector<std::pair<ID, ID>> crossovers) { // insert edges into _graph
    std::cout << "Step 4: Crossovers generating \n";

    ID id1, id2;
    Edge edge;
    for (const auto & cross : crossovers) {
        id1 = cross.first;
        id2 = cross.second;
        edge = makeEdgeCrossover(id1, id2);
        _graph.insertEdge(edge, edge.get_endsNode().first, edge.get_endsNode().second);
    }
    _graph.howManyFourWays();
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
            if (node1 == node2) continue;
            _graph.insertEdge(makeEdge(id1, id2), node1, node2);
        }
    }

    std::cout << "................ DONE\n";

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

