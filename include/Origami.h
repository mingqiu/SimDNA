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



// Initialized as reading input, when function input() is invoked

    HelicalBreakPairs _helicalBreakPairs;


    OrigamiGraph _graph;



public:
    Origami(std::string s) : _fileName(s) {
        std::cout << "Step 1: Input file reading ... ";

        if (!input()) {
            std::cout << "Error occured in input file\n" << s <<
            "\nNo file found" << std::endl;
            exit(1);
        };
        std::cout << "DONE\n";

    }

    std::vector<std::pair<ID, ID>> makeHBPs();
    Node makeNode(HelicalBreakPair pair1, HelicalBreakPair pair2);
    Node makeNode(HelicalBreakPair pair1);
    Edge makeEdge(Node node);
    Edge makeEdgeCrossover(ID id1, ID id2);

    void processNodes(); // insert nodes into _graph
    void processCrossovers(std::vector<std::pair<ID, ID>> crossovers) { // insert edges into _graph
        ID id1, id2;
        Edge edge;
        for (const auto & cross : crossovers) {
            id1 = cross.first;
            id2 = cross.second;
            edge = makeEdgeCrossover(id1, id2);
//            if (id1.baseID()==2501||id2.baseID()==2501)
//                std::cout << //id1 << "\t" << id2 <<
//                        edge.get_endsHB().first << "\t" <<
//                        edge.get_endsHB().second << "\t" <<
//                        edge.is_ds() <<
//                        std::endl;
            _graph.insertEdge(edge, edge.get_endsNode().first, edge.get_endsNode().second);
        }
        _graph.howManyFourWays();
    }



private:

/*
 * Read input *.pairs file to initialize an Origami class
 * Involved in class initialization
 */
    bool input();
/*
 * Examine whether a residue is at helical break
 * Involved in class initialization
 */
    bool isBreak(int, int, int, int, int, int, int, int);
/*
 * Collect all ID belonging to helical break into _breaksOnEachStand
 * */
    void makeHP(int strand, int base);
    Vector3Dd helicalCenter(ID id1);



/***********  Test functions below  *******************/

    void testInput();
    void testhbpAssign();

};




#endif //SIMDNA_ORIGAMI_H
