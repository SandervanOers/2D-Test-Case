#ifndef MESH_GEN
#define MESH_GEn

#include <petscksp.h>
#include <iostream>
#include <limits>
#include <algorithm>
#include <set>

#include <memory>

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
inline void sortrows(std::vector<std::vector<unsigned int>>& matrix, int col) { // can not sort by last column //stable_sort
    std::sort(matrix.begin(),
              matrix.end(),
              [col](const std::vector<unsigned int>& lhs, const std::vector<unsigned int>& rhs) {
                  return lhs[col] < rhs[col] || ((lhs[col]==rhs[col]) && lhs[col+1] < rhs[col+1]);
              });
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
//void load_msh_mesh(const std::string &mesh_name, Vec &VX, Vec &VY, Vec &VZ, Mat &EToV, std::vector<VertexCoordinates3D> &List_Of_Vertices, std::vector<Cuboid> &List_Of_Elements, int &element_num,  int &node_num);
/*--------------------------------------------------------------------------*/
//void Compute_Vertex_Coordinates_Uniform_Rectangle_2D(const double &xmin, const double &xmax, const double &ymin, const double &ymax, const unsigned int &Number_Of_Elements_X, const unsigned int &Number_Of_Elements_Y, std::vector<VertexCoordinates2D> &List_Of_Vertices, std::vector<Boundaries2D> &List_Of_Boundaries, std::vector<Elements2D> &List_Of_Elements);
/*--------------------------------------------------------------------------*/
void set_Order_Polynomials_Uniform(std::vector<Squares2D> &List_Of_Elements2D, const unsigned int &N);
/*--------------------------------------------------------------------------*/
void set_Order_Polynomials_Uniform(std::vector<Cuboid> &List_Of_Elements, const unsigned int &Nx, const unsigned int &Ny, const unsigned int &Nz);
/*--------------------------------------------------------------------------*/
extern unsigned int get_Number_Of_Nodes(std::vector<Cuboid> &List_Of_Elements);
/*--------------------------------------------------------------------------*/
extern unsigned int get_Number_Of_Nodes(std::vector<Squares2D> &List_Of_Elements2D);
/*--------------------------------------------------------------------------*/
void set_theta_Uniform(std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries2D, const double &theta);
/*--------------------------------------------------------------------------*/
void load_msh_mesh2D(const std::string &mesh_name, Vec &VX, Vec &VY, Mat &EToV, std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, std::vector<Squares2D> &List_Of_Elements, int &element_num, int &node_num);
/*--------------------------------------------------------------------------*/
//void load_msh_mesh3D(const std::string &mesh_name, Vec &VX, Vec &VY, Vec &VZS, Mat &EToV, std::vector<VertexCoordinates3D> &List_Of_Vertices, std::vector<Cuboid> &List_Of_Elements, int &element_num,  int &node_num);
/*--------------------------------------------------------------------------*/
void Connect2D(const Mat &EToV, const unsigned int &Number_Of_Elements, const unsigned int &Number_Of_Vertices, Mat &EToE, Mat &EToF, std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries);
/*--------------------------------------------------------------------------*/
//void Connect3D(const Mat &EToV, const unsigned int &Number_Of_Elements, const unsigned int &Number_Of_Vertices, Mat &EToE, Mat &EToF, std::vector<InternalBoundariesCuboid> &List_Of_Boundaries);
/*--------------------------------------------------------------------------*/
//void Calculate_CuboidFaceNormals(const std::vector<Cuboid> &List_Of_Elements, std::vector<InternalBoundariesCuboid> &List_Of_Boundaries, const std::vector<VertexCoordinates3D> &List_Of_Vertices);
/*--------------------------------------------------------------------------*/
//void create_ListInternalBoundaries(const unsigned int &Number_Of_Elements, Mat &EToE, Mat &EToF, std::vector<InternalBoundariesCuboid> &List_Of_Boundaries);
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Square(std::vector<Squares2D> &List_Of_Elements, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices);
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Boundaries_Square(const std::vector<Squares2D> &List_Of_Elements, std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices);
/*--------------------------------------------------------------------------*/
std::string mesh_name_trapezoid(const unsigned int n);
/*--------------------------------------------------------------------------*/
void Calculate_Area_Square(std::vector<Squares2D> &List_Of_Elements, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices);
/*--------------------------------------------------------------------------*/
#endif
