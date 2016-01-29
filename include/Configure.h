//
// Created by Mingqiu Wang on 1/22/16.
//

#ifndef SIMDNA_CONFIGURE_H
#define SIMDNA_CONFIGURE_H


#define TEST 1 // if 1: turn on test; 0: off

//const std::string inputJsonFileName = "fourhelix";
const std::string inputJsonFileName = "A_2D";
//const std::string inputJsonFileName = "pointer_v1_12_no_deletion";
//const std::string inputJsonFileName = "aNANO_3D_7_14";

#define MASS 649
#define STRETCH_DS 1
#define STRETCH_SS 0.1
#define CROSSOVER_DIS 18 // if a neighbor bp is seperated by this distance, it's a crossover
#endif //SIMDNA_CONFIGURE_H
