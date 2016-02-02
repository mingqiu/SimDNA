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



// Initialized when function identifyDiscontinuity() is invoked
    AxialDiscons _axialDiscons;

// Initialized when function processNodes, processCrossovers, connecting are invoked
    OrigamiGraph _graph;




public:
    Origami(std::string s);
/*
 * Identify the type of all axial discontinuities and initialize _axialDiscons
 */
    std::vector<std::pair<ID, ID>> identifyDiscontinuity();
/*
 * Insert nodes into _graph
 */
    void processNodes();
/*
 * Insert edges(crossovers) into _graph
 */
    void processCrossovers(std::vector<std::pair<ID, ID>> crossovers);
/*
 * Insert edges(other than crossovers) into _graph
 */
    void connecting();
/*
 * Output PDB and XML file
 */
    void toPDB(std::string str);
    void toXML(std::string str);






private:

/****** Involved in class initialization  *****/

/*
 * Read input *.pairs file to initialize an Origami class
 */
    bool input();
/*
 * Examine whether a residue is at helical break
 *
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
 *
 */
    bool isDiscontinuity(int, int, int, int, int, int, int, int);
/*
 * Collect all ID belonging to helical break into _breaksOnEachStand
 */
    void makeHB(int strand, int base);



    /****** Useful functions  *****/


/* return coordinate of helical center given a helical break pair
 * if only one resiude, return coor. of this residue
 * otherwise return the average of two residues
 * input:
 * output: coordinate
 */
    Vector3Dd helicalCenter(ID id1);

    Node makeNode(AxialDiscontinuity pair1, AxialDiscontinuity pair2);

    Node makeNode(AxialDiscontinuity pair1);

    Edge makeEdgeCrossover(ID id1, ID id2);

    Edge makeEdge(ID id1, ID id2);



/***********  Test functions below  *******************/

    void testInput();
    void testhbpAssign();

public:
    void toIDs(int i) {
        for (const auto & item : _graph.findIDsfromIndex(i))
            std::cout << item.first << "\t" << item.second << std::endl;
    }

};




#endif //SIMDNA_ORIGAMI_H
