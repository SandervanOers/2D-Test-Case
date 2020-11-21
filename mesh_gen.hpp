#ifndef MESH_GEN
#define MESH_GEn

#include <petscksp.h>
#include <iostream>
#include <limits>
#include <algorithm>
#include <set>

#include "Elements.hpp"
#include "gmsh_io.hpp"
/*--------------------------------------------------------------------------*/
#define NODETOL 1e-12
/*--------------------------------------------------------------------------*/
typedef std::pair<int, int> pairs;
inline auto compare = [](pairs lhs, pairs rhs) //custom compare lambda function
   {
      if(lhs.first > lhs.second ) lhs = pairs{lhs.second, lhs.first };
      if(rhs.first > rhs.second ) rhs = pairs{rhs.second, rhs.first };
      return lhs< rhs;
   };
/*--------------------------------------------------------------------------*/
struct BoundaryInfo {
    PetscInt left;
    PetscInt right;
    PetscInt faceleft;
    PetscInt faceright;
    BoundaryInfo(PetscInt left_v, PetscInt right_v, PetscInt faceleft_v, PetscInt faceright_v) : left(left_v), right(right_v), faceleft(faceleft_v), faceright(faceright_v) {}
};
inline auto compare_BoundaryInfo = [](BoundaryInfo lhs, BoundaryInfo rhs)
{
      if(lhs.left > lhs.right ) lhs = BoundaryInfo(lhs.right, lhs.left, lhs.faceright, lhs.faceleft );
      if(rhs.left > rhs.right ) rhs = BoundaryInfo(rhs.right, rhs.left, rhs.faceright, rhs.faceleft );
      return lhs.left<rhs.left || (!(rhs.left<lhs.left) && lhs.right<rhs.right);
};
/*--------------------------------------------------------------------------*/
void Compute_Vertex_Coordinates_Uniform_Rectangle_2D(const double &xmin, const double &xmax, const double &ymin, const double &ymax, const unsigned int &Number_Of_Elements_X, const unsigned int &Number_Of_Elements_Y, std::vector<VertexCoordinates2D> &List_Of_Vertices, std::vector<Boundaries2D> &List_Of_Boundaries, std::vector<Elements2D> &List_Of_Elements);
/*--------------------------------------------------------------------------*/
// generate equidistant grid
/*--------------------------------------------------------------------------*/
extern Vec mesh_generation_1D_VX(const double &xmin, const double &xmax, const unsigned int &Number_Of_Elements);
/*--------------------------------------------------------------------------*/
extern Mat mesh_generation_1D_EtoV(const double &xmin, const double &xmax, const unsigned int &Number_Of_Elements);
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian(std::vector<Elements2D> &List_Of_Elements2D, const std::vector<VertexCoordinates2D> &List_Of_Vertices);
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_boundaries(std::vector<Boundaries2D> &List_Of_Boundaries2D, const std::vector<VertexCoordinates2D> &List_Of_Vertices);
/*--------------------------------------------------------------------------*/
void set_Order_Polynomials_Uniform(std::vector<Elements2D> &List_Of_Elements2D, const unsigned int &N);
/*--------------------------------------------------------------------------*/
void set_Order_Polynomials_Uniform(std::vector<Squares2D> &List_Of_Elements2D, const unsigned int &N);
/*--------------------------------------------------------------------------*/
extern unsigned int get_Number_Of_Nodes(std::vector<Elements2D> &List_Of_Elements2D);
/*--------------------------------------------------------------------------*/
extern unsigned int get_Number_Of_Nodes(std::vector<Squares2D> &List_Of_Elements2D);
/*--------------------------------------------------------------------------*/
void set_theta_Uniform(std::vector<Boundaries2D> &List_Of_Boundaries2D, const double &theta);
/*--------------------------------------------------------------------------*/
void set_theta_Uniform(std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries2D, const double &theta);
/*--------------------------------------------------------------------------*/
void load_msh_mesh2D(const std::string &mesh_name, Vec &VX, Vec &VY, Mat &EToV, std::vector<VertexCoordinates2D> &List_Of_Vertices, std::vector<Squares2D> &List_Of_Elements, int &element_num, int &node_num);
/*--------------------------------------------------------------------------*/
void Connect2D(const Mat &EToV, const unsigned int &Number_Of_Elements, const unsigned int &Number_Of_Vertices, Mat &EToE, Mat &EToF, std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries);
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Square(std::vector<Squares2D> &List_Of_Elements, const std::vector<VertexCoordinates2D> &List_Of_Vertices);
/*--------------------------------------------------------------------------*/
//void Calculate_Jacobian_Quadrilateral(const Squares2D &Quad, const std::vector<VertexCoordinates2D> &List_Of_Vertices, const double &r_p, const double &s_p, double &Jacobian, double &drdx, double &drdy, double &dsdx, double &dsdy);
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Boundaries_Square(const std::vector<Squares2D> &List_Of_Elements, std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<VertexCoordinates2D> &List_Of_Vertices);
/*--------------------------------------------------------------------------*/
std::string mesh_name_trapezoid(const unsigned int n);
/*--------------------------------------------------------------------------*/
#endif
