//
// Created by Mingqiu Wang on 1/22/16.
//

#include "../include/Origami.h"

using namespace std;

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

        isABreak = isBreak(strand, base, pairStrand, pairBase,
                           strandOld, baseOld, pairStrandOld, pairBaseOld);
        _nucleotide[ ID{strandOld, baseOld} ] =
                Nucleotide{ ID{strandOld, baseOld}, Vector3Dd{xOld, yOld, zOld},
                            ID{pairStrandOld, pairBaseOld}, isABreak};
        if (isABreak) makeHP(strandOld, baseOld);

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
    makeHP(strand, base);
    ++_resNumInEachStrand[strand];

    _strandNum = _resNumInEachStrand.size();

    input.close();
    testInput();
    return true;
}

/*
 * Tell if ID{strandOld, baseOld} belongs to helical break
 *
 * the connectivity is
 * 5' --> {oldold}    --      {strandOld, baseOld}    --     {strand, base}     --> 3'
 * 3' <-- {oldoldCom} -- {pairStrandOld, pairPaseOld} -- {pairStrand, pairBase} <-- 5'
 *
 * The detecting process recognizes an ID to belong to helical break if it doesn't belong to
 *   in ds
 *      1) it has complementary base,
 *      & 2) it has both 3' and 5' neighbors,
 *      & 3) both its 3' and 5' neighbors have complementary, and
 *      & 4) all three complementarity are on the same strand
 *   in ss
 *      1) it has both 3' and 5' neighbors
 *      & 2) both of its 3' and 5' neighbors are single nucleotides
 *
 * */

bool Origami::isBreak(int strand, int base, int pairStrand, int pairBase,
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


void Origami::makeHP(int strand, int base) {
    if (_breaksOnEachStand.find(strand) == _breaksOnEachStand.end()) {
        vector<int> thisStrand;
        thisStrand.push_back(base);
        _breaksOnEachStand[strand] = thisStrand;
    }
    else
        _breaksOnEachStand[strand].push_back(base);
}


vector<pair<ID, ID>> Origami::makeHBPs() {

    ID id, idleft, idright;

    vector<pair<ID, ID>> crossovers;

    for (int i = 0; i < _strandNum; i++) {


        vector<int> strand = _breaksOnEachStand[i];

        for (auto item = 1; item < strand.size()-1; ++item) {
            if (strand[item]==46&&i==148)
                cout << endl;
            id = {i, strand[item]};

            if (!_nucleotide[id].withPair()) {
                _helicalBreakPairs.insert({id, {-1, -1}}, false, false, true);
                continue;
            }

            idleft = {i, strand[item-1]};
            idright = {i, strand[item+1]};

            if ((id.baseID() - idleft.baseID()) == 1) {
                if (dist(helicalCenter(idleft), helicalCenter(id))>18) {
                    _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                              false, true, false);
                    crossovers.push_back({idleft, id});
                }

                if ((idright.baseID() - id.baseID()) == 1) {
                    if (dist(helicalCenter(idright), helicalCenter(id))>10) {
                        _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                                  false, true, false);
                        crossovers.push_back({idright, id});

                    }
                    else {
                        _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                                  _nucleotide[idright].withPair(), false, false);

                    }
                    continue;
                }


                else _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                               true, false, false);
                continue;
            }
            if ((idright.baseID() - id.baseID()) == 1) {
                if (dist(helicalCenter(idright), helicalCenter(id))>10) {
                    _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                              false, true, false);
                    crossovers.push_back({idright, id});

                }
                else {
                    _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                              _nucleotide[idright].withPair(), false, false);

                }
                continue;
            }
            _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                      false, false, false);

        }

        id = {i, strand[0]};
        idright = {i, strand[1]};

        if (!_nucleotide[id].withPair()) {
            _helicalBreakPairs.insert({id, {-1, -1}}, false, false, true);
        }

        else if ((idright.baseID() - id.baseID()) == 1) {
            if (dist(helicalCenter(idright), helicalCenter(id))>10) {
                _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                          false, true, false);
                crossovers.push_back({idright, id});

            }

            else _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                           _nucleotide[idright].withPair(), false, false);
        }
        else
        _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                  false, false, false);


        id = {i, strand[strand.size()-1]};
        idleft = {i, strand[strand.size()-2]};

        if (!_nucleotide[id].withPair()) {
            _helicalBreakPairs.insert({id, {-1, -1}}, false, false, true);
        }

        else if ((id.baseID() - idleft.baseID()) == 1) {
            if (dist(helicalCenter(idleft), helicalCenter(id))>10) {
                _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                          false, true, false);
                crossovers.push_back({idleft, id});

            }
            else _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                           true, false, false);
        }
        else
            _helicalBreakPairs.insert({id, _nucleotide[id].pairID()},
                                      false, false, false);
    }


    return crossovers;
//    _helicalBreakPairs.print();

}

void Origami::processNodes() {
    int hb, hbold;
    HelicalBreakPair hbp, hbpold;
    for (int i = 0; i < _strandNum; i++) {
        vector<int> strand = _breaksOnEachStand[i];
        hb = strand[0];
        if (i==148)
            cout << endl;

        for (int item = 1; item < strand.size(); ++item) {

            hbold = hb;

            hb = strand[item];
            hbpold = hbp;
            hbp = _helicalBreakPairs.findTypeFromID({i,hb});
            if (hb==46&&i==148)
                cout << endl;
            if ((hb - hbold)==1) {
                if (HelicalBreakPair::similar(hbpold, hbp) && hbp.is_typeA() ) {
                    if (dist(helicalCenter({i, hb}), helicalCenter({i, hbold}))>10&& hbp.is_typeB()) {
                        if (item == strand.size() - 1) break;
                        hb = strand[item + 1];
                        hbp = _helicalBreakPairs.findTypeFromID({i, hb});
                        ++item;
                        continue;
                    }
                    _graph.insertNode(makeNode(hbpold, hbp));
                    if (item == strand.size() - 1) break;
                    hb = strand[item + 1];
                    hbp = _helicalBreakPairs.findTypeFromID({i, hb});
                    ++item;
                    continue;
                }
            }
            _graph.insertNode(makeNode(hbpold));
        }
        _graph.insertNode(makeNode(hbpold));
    }
    _graph.resizeGraph();
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
            hbp = _helicalBreakPairs.findIndexFromID(id);
            myfile << id << "\t" << _nucleotide[id].pairID() << "\t" <<
            (_helicalBreakPairs.findTypeFromIndex(hbp).is_typeA()&&
                    _helicalBreakPairs.findTypeFromIndex(hbp).is_typeB()) << endl;
        }

    }
    myfile.close();

}

Node Origami::makeNode(HelicalBreakPair pair1, HelicalBreakPair pair2) {
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

Node Origami::makeNode(HelicalBreakPair pair) {
    int type = pair.is_single() ? 3 : 2;
    std::vector<std::pair<ID,ID>> ids;
    ids.push_back(pair.get_ids());

    double mass = 1.0/(type - 1) * MASS;
    Vector3Dd position = helicalCenter(pair.get_ids().first);
    double vdWradii = 0.2;

    return Node{type, ids, mass, position, vdWradii};
}

Edge Origami::makeEdge(Node node) {

    node.get_num();

//    Edge(const std::pair<int, int> &_endsNode, const std::pair<int, int> &_types,
//    bool _crossover, const std::vector<std::pair<ID, ID>> &_ids)
//    : _endsNode(_endsNode),_types(_types), _crossover(_crossover), _ids(_ids)

    return Edge();
}


Edge Origami::makeEdgeCrossover(ID id1, ID id2) {
    int a = _graph.findNodeNum(id1);
    int b = _graph.findNodeNum(id2);
    int c = _graph.findNodeType(id1);
    int d = _graph.findNodeType(id2);
    if (id1.baseID()==2477||id2.baseID()==2477)
        cout << a << "\t"<< b << endl;
    return Edge({a, b}, {c, d}, {id1, id2}, true, {});
}




