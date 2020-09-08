#ifndef INITIAL_CONDITION
#define INITIAL_CONDITION

#include <petscksp.h>
/*--------------------------------------------------------------------------*/
// Stratified Test Problem
/*--------------------------------------------------------------------------*/
extern double calculate_sigma(const double &N2, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double rho_0(const double &x, const double &N2);
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv(const double &x, const double &N2);
/*--------------------------------------------------------------------------*/
extern double N_2(const double &x, const double &N2);
/*--------------------------------------------------------------------------*/
// Viscous Homogeneous Test Problem
/*--------------------------------------------------------------------------*/
extern double viscous_calculate_sigma(const double &nu, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double viscous_Exact_Solution_m(const double &x, const double &t, const double &nu, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double viscous_Exact_Solution_p(const double &x, const double &t, const double &nu, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double viscous_rho_0(const double &x, const double &N2);
/*--------------------------------------------------------------------------*/
extern double viscous_rho_0_deriv(const double &x, const double &N2);
/*--------------------------------------------------------------------------*/
// System 1 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_system1(const double &N2, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m_system1(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_system1(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
// System 3 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_system3(const double &beta, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m_system3(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_system3(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double rho_0_compressible(const double &x, const double &beta);
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv_compressible(const double &x, const double &beta);
/*--------------------------------------------------------------------------*/
extern double N_2_compressible(const double &x, const double &beta);
/*--------------------------------------------------------------------------*/
// System 4 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_system4(const double &beta, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m_system4(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_system4(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_system4(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
// System 5 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_system5(const double &beta, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m_system5(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_system5(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_system5(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
// System 6: Verification: 1D Compressible Nonviscous Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_system1dcom(const double &beta, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m_system1dcom(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_system1dcom(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_system1dcom(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
// 2D
/*--------------------------------------------------------------------------*/
// System 1 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_2D_system1(const double &beta, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mx_2D_system1(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mz_2D_system1(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_2D_system1(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_2D_system1(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double rho_0_2D_system1(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv_2D_system1(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern double N_2_2D_system1(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/

#endif
