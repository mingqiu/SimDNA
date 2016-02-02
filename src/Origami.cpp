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

vector<pair<ID, ID>> Origami::identifyDiscontinuity() {
    std::cout << "Step 2: Axial discontinuities generating \n";

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

// if belongs to double strand
            // if it has a left neighbor discontinuity
            if ((id.baseID() - idleft.baseID()) == 1) {
                if (dist(helicalCenter(idleft), helicalCenter(id)) > THRES_CROSSOVER_DIS) {
                    _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                         false, true, false);
                    crossovers.push_back({idleft, id});
                    continue;
                }
                if ((idright.baseID() - id.baseID()) == 1) {
                    if (dist(helicalCenter(idright), helicalCenter(id)) > THRES_CROSSOVER_DIS) {
                        _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                             false, true, false);
                        crossovers.push_back({idright, id});
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
                    crossovers.push_back({idright, id});
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

        if (!_nucleotide[id].withPair())
            _axialDiscons.insert({id, {-1, -1}}, false, false, true);

        else if ((idright.baseID() - id.baseID()) == 1) {
            if (dist(helicalCenter(idright), helicalCenter(id)) > THRES_CROSSOVER_DIS) {
                _axialDiscons.insert({id, _nucleotide[id].pairID()},
                                     false, true, false);
                crossovers.push_back({idright, id});
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

//    for (const auto &item1: crossovers)
//        cout << item1.first << "\t" << item1.second << endl;
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


void Origami::processCrossovers(std::vector<std::pair<ID, ID>> crossovers) {
    std::cout << "Step 4: Crossovers generating \n";

    ID id1, id2;
    Edge edge;
    for (const auto &cross : crossovers) {
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
            if (node1 == node2) continue;
            _graph.insertEdge(makeEdge(id1, id2), node1, node2);
        }
    }

    std::cout << "................ DONE\n";


    std::cout << "Step 6: Holliday junctions identifying \n";
    _graph.findHollidayJ();

    std::cout << "................ DONE\n";

}




