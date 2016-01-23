//
// Created by Mingqiu Wang on 1/22/16.
//

#ifndef SIMDNA_DATASTRUCTURE_H
#define SIMDNA_DATASTRUCTURE_H

#include <stdio.h>
#include <iostream>

#include "../lib/Vector3D.h"


/*
 * ID of each nucleotide residue {strand No., base No.}
 *
 *     default constructor assigns both strand and base number as -1
 *     therefore if there exists a base number as -1, meaning there is not such base or it hasn't defined
 * */
class ID {
    int _strandID;
    int _baseID;
public:

    ID() : _strandID{-1}, _baseID{-1} {}
    ID(int s, int b) : _strandID{s}, _baseID{b} {}

    int strandID() const { return _strandID; }
    int baseID() const { return _baseID; }

    bool operator<(const ID& rhs) const {
        return _strandID < rhs._strandID || (rhs._strandID >= _strandID && _baseID < rhs._baseID);
    }
    bool operator==(const ID& other) const {
        return (_strandID == other._strandID && _baseID == other._baseID);
    }

};

std::ostream& operator<<(std::ostream& os, const ID &idnum);

struct IDHasher {
    int operator()(const ID& k) const {
        return ((std::hash<int>()(k.strandID())
                 ^ (std::hash<int>()(k.baseID()) << 1)^(std::hash<int>()(k.baseID()) << 1)) >> 1);
    }
};

/*
 * Information (pair, coordinate, etc..) of a nucleotide
 *
 * */
class Nucleotide {
    ID _id;
    Vector3Dd _coordinate;
    ID _pairID;
    bool _isSS; // true if pair residue doesn't exist, _pairID == {-1, -1}
    bool _isBreak; // whether this residue belongs to helical break

public:
    Nucleotide() {}
    Nucleotide(ID IDnum, Vector3Dd coor, ID comID, bool isBreak) :
            _id{IDnum}, _pairID{comID}, _coordinate{coor},
            _isSS{ !((comID.strandID() == -1) && (comID.baseID() == -1))},
            _isBreak{ isBreak } {}


    ID id() const {return _id;}
    ID pairID() const {return _pairID;}
    Vector3Dd coordinate() const {return _coordinate;}
    bool withPair() const {return _isSS;}
    bool isBreak() const { return _isBreak;}

};



#endif //SIMDNA_DATASTRUCTURE_H
