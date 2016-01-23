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
    ++_resNumInEachStrand[strand];

    _strandNum = _resNumInEachStrand.size();

    input.close();
//    testInput();
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

