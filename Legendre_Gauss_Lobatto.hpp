#ifndef LGL_GRID
#define LGL_GRID

#include <petscksp.h>
#include <slepceps.h>
#include <slepcsys.h>
#include <iostream>
#include <limits>
#include <array>
#include <vector>
#include <chrono>
#include <memory>
#include "Elements.hpp"
#include "Cubature2D.hpp"
/*--------------------------------------------------------------------------*/
#define NODETOL 1e-12
/*--------------------------------------------------------------------------*/
extern Vec JacobiGL(const double &alpha, const double &beta, const unsigned int &N);
/*--------------------------------------------------------------------------*/
extern Vec JacobiGL_withWeights(const double &alpha, const double &beta, const unsigned int &N, Vec &Weights);
/*--------------------------------------------------------------------------*/
Vec JacobiGQ(const double &alpha, const double &beta, const unsigned int &N);
/*--------------------------------------------------------------------------*/
extern Vec JacobiGQ_withWeights(const double &alpha, const double &beta, const unsigned int &N, Vec &Weights);
/*--------------------------------------------------------------------------*/
extern Vec JacobiP(const Vec &x, const double &alpha, const double &beta, const unsigned int &N);
// Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
//          (alpha+beta <> -1) at points x for order N and returns P[1:length(xp))]
// Note   : They are normalized to be orthonormal.
/*--------------------------------------------------------------------------*/
Vec GradJacobiP(const Vec &r, const double &alpha, const double &beta, const unsigned int &N);
// Purpose: Evaluate the derivative of the Jacobi polynomial of type (alpha,beta)>-1,
//	       at points r for order N and returns dP[1:length(r))]
/*--------------------------------------------------------------------------*/
extern Mat Vandermonde1D(const Vec &r, const unsigned int &N);
// Initialize the 1D Vandermonde matrix, V_{ij} = phi_j(r_i);
/*--------------------------------------------------------------------------*/
Mat GradVandermonde1D(const Vec &r, const unsigned int &N);
// Initialize  the gradient of modal basis (i) at (r) at order N
/*--------------------------------------------------------------------------*/
extern Mat DMatrix1D(const Vec &r, const unsigned int &N, const Mat &V);
// Initialize the (r) differentiation matrices on the interval, evaluated at (r) at order N
/*--------------------------------------------------------------------------*/
extern Mat Lift1D(const unsigned int &N, const Mat &V);
/*--------------------------------------------------------------------------*/
extern Mat GeometricFactors1D(const Mat &x, const Mat &Dr);
/*--------------------------------------------------------------------------*/
extern Mat normals1D(const unsigned int &N, const unsigned int &Number_Of_Elements);
/*--------------------------------------------------------------------------*/
extern double LagrangePolynomial(const Vec &r, const double &x, const unsigned int &i);
/*--------------------------------------------------------------------------*/
extern double LagrangePolynomialDeriv(const Vec &r, const double &x, const unsigned int &i);
/*--------------------------------------------------------------------------*/
// 2D //
/*--------------------------------------------------------------------------*/
extern void Nodes2D(unsigned int N, Vec &X, Vec &Y);
/*--------------------------------------------------------------------------*/
extern void XYtoRS(const Vec &X, const Vec &Y, Vec &R, Vec &S);
/*--------------------------------------------------------------------------*/
Vec Warpfactor(const unsigned int &N, const Vec &rout);
/*--------------------------------------------------------------------------*/
extern Mat Vandermonde2D(const unsigned int &N, const Vec &r, const Vec &s);
/*--------------------------------------------------------------------------*/
void RStoAB(const Vec &R, const Vec &S, Vec &A, Vec &B);
/*--------------------------------------------------------------------------*/
Vec Simplex2DP(const Vec &A, const Vec &B, const unsigned int &i, const unsigned int &j);
/*--------------------------------------------------------------------------*/
void store_Nodes_Reference_Triangle();
/*--------------------------------------------------------------------------*/
extern void Read_RS_Coordinates_Reference_Triangle(const unsigned int &N, Vec &R, Vec &S);
/*--------------------------------------------------------------------------*/
void GradSimplex2D(const Vec &a, const Vec &b, const unsigned int &id, const unsigned int &jd, Vec &dmodedr, Vec &dmodeds);
/*--------------------------------------------------------------------------*/
void GradVandermonde2D(const unsigned int &N, const Vec &R, const Vec &S, Mat &V2Dr, Mat &V2Ds);
/*--------------------------------------------------------------------------*/
extern void DMatrices2D(const unsigned int &N, const Vec &R, const Vec &S, const Mat &V, Mat &Dr, Mat &Ds);
/*--------------------------------------------------------------------------*/
extern Mat InterpMatrix2D(const unsigned int &N, const Vec &R, const Vec &S);
/*--------------------------------------------------------------------------*/
extern Mat Inverse_Matrix(const Mat &V);
/*--------------------------------------------------------------------------*/
extern Mat load_VandermondeMatrix(const unsigned int N);
/*--------------------------------------------------------------------------*/
extern Mat load_InverseVandermondeMatrix(const unsigned int N);
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix2D(const unsigned int &N);
/*--------------------------------------------------------------------------*/
extern Mat InverseMassMatrix2D(const unsigned int &N);
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix2D_Cubature(const unsigned int &N, const Mat &cubV, const Vec &cubW, const unsigned int &Ncub);
/*--------------------------------------------------------------------------*/
//extern Vec LagrangePolynomial2D(const Mat &V, const Vec &P, const unsigned int &Np);
/*--------------------------------------------------------------------------*/
extern void store_LagrangePolynomial_Cubature();
/*--------------------------------------------------------------------------*/
extern void store_DerivativeLagrangePolynomial_Cubature();
/*--------------------------------------------------------------------------*/
extern Mat load_LagrangePolynomial_Cubature(const unsigned int &Order_Polynomials, const unsigned int &Cubature_Order);
/*--------------------------------------------------------------------------*/
extern Mat load_DerivativeLagrangePolynomial_Cubature(const unsigned int &Order_Polynomials, const unsigned int &Cubature_Order, const bool &deriv);
/*--------------------------------------------------------------------------*/
extern Mat VandermondeMultiply(const Mat &V2DInv, const Mat &Matrix);
/*--------------------------------------------------------------------------*/
void GeometricFactors2D(const Vec &x, const Vec &y, const Mat &Dr, const Mat &Ds, const unsigned int &Np, Vec &rx, Vec &ry, Vec &sx, Vec &sy, double &J);
/*--------------------------------------------------------------------------*/
extern void stdVectorToPetscVec(const std::vector<double> VecIn, Vec &PetscVec);
/*--------------------------------------------------------------------------*/
extern void NodesSquares2D(unsigned int N, Vec &XX, Vec &YY);
/*--------------------------------------------------------------------------*/
void set_Node_Coordinates_Uniform_Square2D(std::vector<Squares2D> &List_Of_Elements, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, const unsigned int &N);
/*--------------------------------------------------------------------------*/
extern double LagrangePolynomial_Test(const Vec &r, const double &x, const unsigned int &i);
/*--------------------------------------------------------------------------*/
extern double LagrangePolynomialDeriv_Test(const Vec &r, const double &x, const unsigned int &i);
/*--------------------------------------------------------------------------*/
Vec SecondDerivJacobiP(const Vec &r, const double &alpha, const double &beta, const unsigned int &N);
/*--------------------------------------------------------------------------*/
//void set_Node_Coordinates_Cuboid(std::vector<Cuboid> &List_Of_Elements, const std::vector<VertexCoordinates3D> &List_Of_Vertices, const unsigned int &Nx, const unsigned int &Ny, const unsigned int &Nz);
/*--------------------------------------------------------------------------*/
// 3D //
/*--------------------------------------------------------------------------*/
void NodesCuboid(unsigned int Nx, unsigned int Ny, unsigned int Nz, Vec &XX, Vec &YY, Vec &ZZ);
/*--------------------------------------------------------------------------*/
void set_Node_Coordinates_Uniform(std::vector<std::unique_ptr<Element>> &List_Of_Elements, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, const unsigned int &Nx, const unsigned int &Ny, const unsigned int &Nz);
/*--------------------------------------------------------------------------*/
#endif
