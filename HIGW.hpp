#ifndef HIGW
#define HIGW

#include <petscksp.h>
#include <iostream>
#include <iomanip>

#include "Elements.hpp"
#include "Legendre_Gauss_Lobatto.hpp"
#include "initial_cond.hpp"
//#include "mesh_gen.hpp"
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
extern double calculate_Hamiltonian2D_full(const Mat &M1, const Mat &M2, const Vec &Solution, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian2D_IC(const Mat &M1, const Mat &M2, const Vec &Solution, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const double &DeltaX, const double &DeltaY, const unsigned int &Np);
/*--------------------------------------------------------------------------*/
extern void create_Matrices(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &NMat, Mat &NDerivMat);
/*--------------------------------------------------------------------------*/
void create_Matrices_Quads(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &NMat, Mat &NDerivMat);
/*--------------------------------------------------------------------------*/
void create_Matrices_Quads_Full(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, const double &rho_0_Deriv, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &M2_small, Mat &NMat, Mat &NDerivMat);
/*--------------------------------------------------------------------------*/
extern void create_Compressible_System_MidPoint(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, Mat &A, Mat &B);
/*--------------------------------------------------------------------------*/
extern void create_Compressible_System_MidPoint_Full(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B);
/*--------------------------------------------------------------------------*/
extern void create_Incompressible_System_MidPoint_Full(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B, Mat &DIV);
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecU, Vec &VecW, Vec &VecR, Vec &VecP);
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition_system2(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecU, Vec &VecW, Vec &VecR, Vec &VecP, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc);
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition_Incompressible(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &kxmode, const double &kzmode, Vec &Initial_Condition, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc, const Mat &DIV);
/*--------------------------------------------------------------------------*/
extern void Simulate(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Vec &Initial_Condition, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables);
/*--------------------------------------------------------------------------*/
extern void Simulate_IC(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables);
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D_Quad(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const std::vector<Squares2D> &List_Of_Elements, const std::vector<VertexCoordinates2D> &List_Of_Vertices, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Quadrilateral(const Squares2D &Quad, const std::vector<VertexCoordinates2D> &List_Of_Vertices, const double &r_p, const double &s_p, double &Jacobian, double &drdx, double &drdy, double &dsdx, double &dsdy, double &x, double &y);
/*--------------------------------------------------------------------------*/
void printFullMatrixInfo(Mat& matrix, const std::string& name);
/*--------------------------------------------------------------------------*/
extern void correctInitialProjectionOfVelocity(const unsigned int &N_Nodes, Vec &UInit, const Mat &DIV);
/*--------------------------------------------------------------------------*/
void compute_Divergence_Velocity(const Vec &Initial_Condition, const double &N_Nodes, const Mat &DIV);
/*--------------------------------------------------------------------------*/
void create_Matrices_Quads_Full_IC(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, const double &rho_0_Deriv, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &M2_small, Mat &NMat, Mat &NDerivMat);
/*--------------------------------------------------------------------------*/
#endif
