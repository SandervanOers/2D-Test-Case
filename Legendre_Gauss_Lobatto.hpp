#ifndef LGL_GRID
#define LGL_GRID

#include <petscksp.h>
#include <slepceps.h>
#include <slepcsys.h>
#include <iostream>
#include <limits>
#include <array>
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
#endif
