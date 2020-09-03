#ifndef MESH_GEN
#define MESH_GEn

#include <petscksp.h>
#include <iostream>
#include <limits>

#include "Elements.hpp"
/*--------------------------------------------------------------------------*/
void Compute_Vertex_Coordinates_Uniform_Rectangle_2D(const double &xmin, const double &xmax, const double &ymin, const double &ymax, const unsigned int &Number_Of_Elements_X, const unsigned int &Number_Of_Elements_Y, std::vector<VertexCoordinates2D> &List_Of_Vertices, std::vector<Boundaries2D> &List_Of_Boundaries, std::vector<Elements2D> &List_Of_Elements);
/*--------------------------------------------------------------------------*/
// generate equidistant grid
/*--------------------------------------------------------------------------*/
extern Vec mesh_generation_1D_VX(const double &xmin, const double &xmax, const unsigned int &Number_Of_Elements);
/*--------------------------------------------------------------------------*/
extern Mat mesh_generation_1D_EtoV(const double &xmin, const double &xmax, const unsigned int &Number_Of_Elements);
/*--------------------------------------------------------------------------*/
Mat FaceToFace_1D(const unsigned int &Number_Of_Elements, const Mat &EtoV);
/*--------------------------------------------------------------------------*/
#endif
