#ifndef MESH_GEN
#define MESH_GEn

#include <petscksp.h>
#include <iostream>
#include <limits>
// generate equidistant grid
/*--------------------------------------------------------------------------*/
extern Vec mesh_generation_1D_VX(const double &xmin, const double &xmax, const unsigned int &Number_Of_Elements);
/*--------------------------------------------------------------------------*/
extern Mat mesh_generation_1D_EtoV(const double &xmin, const double &xmax, const unsigned int &Number_Of_Elements);
/*--------------------------------------------------------------------------*/
Mat FaceToFace_1D(const unsigned int &Number_Of_Elements, const Mat &EtoV);
/*--------------------------------------------------------------------------*/
#endif
