//
// Created by Mingqiu Wang on 1/22/16.
//

#ifndef SIMDNA_ORIGAMI_H
#define SIMDNA_ORIGAMI_H

#include <string>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <set>
#include <queue>

#include "../include/DataStructure.h"
#include "../include/Configure.h"

class Origami {

    Origami() = delete;


// Initialized as reading input, when function input() is invoked
    std::string _fileName;
    // name of input *.pairs file

    int _totalResidueNum;
    // total numbers of nucleotide residues in this origami

    int _strandNum;
    // how many strands

    std::unordered_map<int, int> _resNumInEachStrand;
    // store number of nucleotide residue in each strand, key = strand ID, value = number of residues in this strand

    std::unordered_map<ID, Nucleotide, IDHasher> _nucleotide;
    // given on ID, find all its information

    std::unordered_map<int, std::vector<int>> _breaksOnEachStand;
    // helical breaks on each strand



// Initialized when function makeHBPs() is invoked
    HelicalBreakPairs _helicalBreakPairs;


    OrigamiGraph _graph;



public:
    Origami(std::string s);
    std::vector<std::pair<ID, ID>> makeHBPs();
    void processNodes(); // insert nodes into _graph
    void processCrossovers(std::vector<std::pair<ID, ID>> crossovers);
    void connecting();



private:

/****** Involved in class initialization  *****/

/*
 * Read input *.pairs file to initialize an Origami class
 */
    bool input();
/*
 * Examine whether a residue is at helical break
 */
    bool isBreak(int, int, int, int, int, int, int, int);
/*
 * Collect all ID belonging to helical break into _breaksOnEachStand
 */
    void makeHB(int strand, int base);


/****** Useful functions  *****/
    Vector3Dd helicalCenter(ID id1);
    Node makeNode(HelicalBreakPair pair1, HelicalBreakPair pair2);
    Node makeNode(HelicalBreakPair pair1);
    Edge makeEdgeCrossover(ID id1, ID id2);
    Edge makeEdge(int node1, int node2);
    Edge makeEdge(ID id1, ID id2);


/***********  Test functions below  *******************/

    void testInput();
    void testhbpAssign();


public:

    void toPDB(std::string str);

};




#endif //SIMDNA_ORIGAMI_H
