//
// Created by Mingqiu Wang on 1/23/16.
//

#include "../include/Origami.h"


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
