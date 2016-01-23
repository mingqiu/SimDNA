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
    size_t _strandNum;
    // how many strands
    std::unordered_map<int, int> _resNumInEachStrand;
    // store number of nucleotide residue in each strand, key = strand ID, value = number of residues in this strand
    std::unordered_map<ID, Nucleotide, IDHasher> _nucleotide;
    // given on ID, find all its information


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




/***********  Test functions below  *******************/

    void testInput();

};




#endif //SIMDNA_ORIGAMI_H
