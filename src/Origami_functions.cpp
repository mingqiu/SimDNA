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


// return coordinate of helical center given a helical break pair
// if only one resiude, return coor. of this residue
// otherwise return the average of two residues
// input:
// output: coordinate
Vector3Dd Origami::helicalCenter(ID id1) {
    ID id2 = _nucleotide[id1].pairID();

    if (id1.baseID()!=-1&&id2.baseID()!=-1)
        return (_nucleotide[id1].coordinate()+_nucleotide[id2].coordinate())*0.5;
    else if (id1.baseID()==-1)
        return (_nucleotide[id2].coordinate());
    else if (id2.baseID()==-1)
        return (_nucleotide[id1].coordinate());
}
