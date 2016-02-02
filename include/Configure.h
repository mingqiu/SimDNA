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

#define MASS 649 // (g/mol)
#define STRETCH_DS 1 // (pN/nm)
#define STRETCH_SS 0.1 // (pN/nm)
#define ANGLE_W 0.1
#define ANGLE_S 1.0
#define THRES_CROSSOVER_DIS 18 // (A) if two neighbor bps are seperated by this distance, it's a crossover
#define CROSSOVER_DIS 2.2 // (nm) equilibrium length of a crossover
#define RISE_PER_BP 0.332 // (nm) helical rise per bp
#define VDWRADII_1 1 // (nm) vdw radii of type 1# node
#define VDWRADII_2 0.2 // (nm) vdw radii of type 2# node
#define VDWRADII_3 1.2 // (nm) vdw radii of type 3# node

#endif //SIMDNA_CONFIGURE_H
