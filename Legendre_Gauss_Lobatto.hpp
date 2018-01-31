#ifndef LGL_GRID
#define LGL_GRID

#include <petscksp.h>
#include <slepceps.h>
#include <slepcsys.h>
#include <iostream>
#include <limits>
/*--------------------------------------------------------------------------*/
extern Vec JacobiGL(const double &alpha, const double &beta, const unsigned int &N);
/*--------------------------------------------------------------------------*/
Vec JacobiGQ(const double &alpha, const double &beta, const unsigned int &N);
/*--------------------------------------------------------------------------*/
Vec JacobiP(const Vec &x, const double &alpha, const double &beta, const unsigned int &N);
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
#endif
