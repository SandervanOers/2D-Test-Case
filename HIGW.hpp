#ifndef HIGW
#define HIGW

#include <petscksp.h>
#include <iostream>
#include <iomanip>
#include <memory>

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
extern double calculate_Hamiltonian2D(const Mat &M1, const Vec &Solution, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian2D_full(const Mat &M1, const Mat &M2, const Vec &Solution, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian2D_IC(const Mat &M1, const Mat &M2, const Vec &Solution, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, double &M);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian3D_IC(const Mat &M1, const Mat &M2, const Vec &Solution, const std::vector<std::unique_ptr<Element>> &List_Of_Elements, const unsigned int &N_Nodes, double &M);
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const double &DeltaX, const double &DeltaY, const unsigned int &Np);
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
extern void compute_InitialCondition_Incompressible(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecNodes, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc, const Mat &DIV);
/*--------------------------------------------------------------------------*/
extern void Simulate(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Vec &Initial_Condition, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables);
/*--------------------------------------------------------------------------*/
extern void Simulate_IC(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const Vec &VecNodes, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega);
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D_Quad(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const std::vector<Squares2D> &List_Of_Elements, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
extern double calculate_Error3D_FV(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const double &DeltaX, const double &DeltaY, const double &DeltaZ, const unsigned int &Np);
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Quadrilateral(const Squares2D &Quad, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, const double &r_p, const double &s_p, double &Jacobian, double &drdx, double &drdy, double &dsdx, double &dsdy, double &x, double &y);
/*--------------------------------------------------------------------------*/
//void Calculate_Jacobian_Cuboid(const Cuboid &Element, const std::vector<VertexCoordinates3D> &List_Of_Vertices, const double &r_p, const double &s_p, const double &t_p, double &det_J, double &drdx, double &drdy, double &drdz, double &dsdx, double &dsdy, double &dsdz, double &dtdx, double &dtdy, double &dtdz, double &x, double &y, double &z);
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Cuboid(const std::unique_ptr<Element> &Element, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, const double &r_p, const double &s_p, const double &t_p, double &det_J, double &drdx, double &drdy, double &drdz, double &dsdx, double &dsdy, double &dsdz, double &dtdx, double &dtdy, double &dtdz, double &x, double &y, double &z);
/*--------------------------------------------------------------------------*/
void printFullMatrixInfo(Mat& matrix, const std::string& name);
/*--------------------------------------------------------------------------*/
extern void correctInitialProjectionOfVelocity(const unsigned int &N_Nodes, Vec &UInit, const Mat &DIV);
/*--------------------------------------------------------------------------*/
void compute_Divergence_Velocity(const Vec &Initial_Condition, const double &N_Nodes, const Mat &DIV);
/*--------------------------------------------------------------------------*/
extern void create_WA_System_MidPoint(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B, Mat &DIV, const double &Re, const double &Fr);
/*--------------------------------------------------------------------------*/
extern void Simulate_WA(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const Vec &VecNodes, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega);
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition_WA(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &Fr, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecNodes, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc, const Mat &DIV);
/*--------------------------------------------------------------------------*/
// EB //
/*--------------------------------------------------------------------------*/
void create_Matrices_Quads_EB(const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, const std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, const double &rho_0_Deriv, const double &Fr, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &M2_small, Mat &NMat, Mat &NDerivMat);
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition_EB(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &Fr, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecNodes, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc, const Mat &DIV);
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition_3DEB(const std::vector<std::unique_ptr<Element>> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &Fr, const double &kxmode, const double &kymode, const double &kzmode, Vec &Initial_Condition, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc, const Mat &DIV);
/*--------------------------------------------------------------------------*/
extern void create_EB_System_MidPoint(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B, Mat &DIV, const double &Re, const double &Fr);
/*--------------------------------------------------------------------------*/
extern void Simulate_EB(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const Vec &VecNodes, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega);
/*--------------------------------------------------------------------------*/
extern void ComputeNodesForcing(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, Vec &VecNodes);
/*--------------------------------------------------------------------------*/
extern void create_WA_System_Forced_MidPoint(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const Vec &Forcing_a, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B, Mat &DIV, const double &Re, const double &Fr, const double &gamma);
/*--------------------------------------------------------------------------*/
extern void create_WA_System_Forced_MidPoint3D(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const Vec &Forcing_a, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B, Mat &DIV, const double &Re, const double &Fr, const double &gamma);
/*--------------------------------------------------------------------------*/
extern void ComputeForcing(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, Vec &Forcing_a);
/*--------------------------------------------------------------------------*/
extern void ComputeForcing(const std::vector<std::unique_ptr<Element>> &List_Of_Elements, const unsigned int &N_Nodes, Vec &Forcing_a);
/*--------------------------------------------------------------------------*/
extern void Simulate_WA_Forced(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega);
/*--------------------------------------------------------------------------*/
extern void Simulate_WA_Forced3D(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const std::vector<std::unique_ptr<Element>> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega);
/*--------------------------------------------------------------------------*/
extern void Simulate_WA_Forced_Continuous(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega);
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition_EB_Bucket(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &Fr, const double &kxmode, const double &kzmode, Vec &Initial_Condition, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc, const Mat &DIV);
/*--------------------------------------------------------------------------*/
extern void create_rot_WA_System_Forced_MidPoint(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const Vec &Forcing_a, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B, Mat &DIV, const double &Re, const double &Fr, const double &gamma);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian2D_kin(const Mat &M1, const Vec &Solution, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes);
/*--------------------------------------------------------------------------*/
extern void Simulate_rot(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const Vec &VecNodes, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega);
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition_rot_WA(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &Fr, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecNodes, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc, const Mat &DIV);
/*--------------------------------------------------------------------------*/
extern void create_Matrices_Cuboids(const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, const std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries, const std::vector<std::unique_ptr<Element>> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Nx_e, const unsigned int &Ny_e, const unsigned int &Nz_e, const unsigned int &N_Q, const double &Fr, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &M2_small, Mat &NMat, Mat &NDerivMat);
/*--------------------------------------------------------------------------*/
#endif
