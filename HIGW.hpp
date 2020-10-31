#ifndef HIGW
#define HIGW

#include <petscksp.h>
#include <iostream>
#include <iomanip>

#include "Elements.hpp"
#include "Legendre_Gauss_Lobatto.hpp"
#include "initial_cond.hpp"
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix_local(const Mat &V);
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix_inverse_local(const Mat &V);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian(const Mat &M1, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian_pres(const Mat &M1, const Mat &M2, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian_comp(const Mat &M1, const Mat &M2, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np);
/*--------------------------------------------------------------------------*/
extern double calculate_Error(const Vec &Exact, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np, const double &DeltaX);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian2D(const Mat &M1, const Vec &Solution, const std::vector<Elements2D> &List_Of_Elements2D, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian2D(const Mat &M1, const Vec &Solution, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const double &DeltaX, const double &DeltaY, const unsigned int &Np);
/*--------------------------------------------------------------------------*/
extern void create_Matrices(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &NMat, Mat &NDerivMat);
/*--------------------------------------------------------------------------*/
extern void create_Compressible_System_MidPoint(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, Mat &A, Mat &B);
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecU, Vec &VecW, Vec &VecR, Vec &VecP);
/*--------------------------------------------------------------------------*/
extern void Simulate(const Mat &A, const Mat &B, const Mat &M1_small, const Vec &Initial_Condition, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol);
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
#endif
