#ifndef INITIAL_CONDITION
#define INITIAL_CONDITION

#include <iostream>
#include <petscksp.h>
#include <vector>
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
extern std::vector<double> rho_0_2D_system1(const std::vector<double> &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv_2D_system1(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern std::vector<double> rho_0_deriv_2D_system1(const std::vector<double> &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern double N_2_2D_system1(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern std::vector<double> N_2_2D_system1(const std::vector<double> &z, const double &beta);
/*--------------------------------------------------------------------------*/
// System 2 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_2D_system2(const double &beta, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mx_2D_system2(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mz_2D_system2(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_2D_system2(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_2D_system2(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double rho_0_2D_system2(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern std::vector<double> rho_0_2D_system2(const std::vector<double> &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv_2D_system2(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern std::vector<double> rho_0_deriv_2D_system2(const std::vector<double> &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern double N_2_2D_system2(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern std::vector<double> N_2_2D_system2(const std::vector<double> &z, const double &beta);
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
//  Incompresisble Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_2DIC(const double &beta, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mx_2DIC(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mz_2DIC(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_2DIC(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_2DIC(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz);
/*--------------------------------------------------------------------------*/
extern double rho_0_2DIC(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/
//extern std::vector<double> rho_0_2DIC(const std::vector<double> &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv_2DIC(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/
//extern std::vector<double> rho_0_deriv_2DIC(const std::vector<double> &z, const double &beta);
/*--------------------------------------------------------------------------*/
extern double N_2_2DIC(const double &z, const double &beta);
/*--------------------------------------------------------------------------*/
//extern std::vector<double> N_2_2DIC(const std::vector<double> &z, const double &beta);
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
//  WA Case - Experiments
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_2DWA(const double &beta, const unsigned int &kx, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_mx_2DWA(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_mz_2DWA(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_r_2DWA(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_p_2DWA(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr);
extern double rho_0_2DWA(const double &z, const double &beta, const double &Fr);
extern double rho_0_deriv_2DWA(const double &z, const double &beta, const double &Fr);
extern double N_2_2DWA(const double &z, const double &beta, const double &Fr);
/*--------------------------------------------------------------------------*/
//  EB Case - Experiments
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_2DEB(const double &beta, const unsigned int &kx, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_mx_2DEB(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_mz_2DEB(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_r_2DEB(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_p_2DEB(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr);
extern double rho_0_2DEB(const double &z, const double &beta, const double &Fr);
extern double rho_0_deriv_2DEB(const double &z, const double &beta, const double &Fr);
extern double N_2_2DEB(const double &z, const double &beta, const double &Fr);
/*--------------------------------------------------------------------------*/
// EB - Bucket Test Case
extern double calculate_sigma_2DEB_Bucket();
extern double period_2DEB_Bucket();
extern double Exact_Solution_mx_2DEB_Bucket(const double &x, const double &z, const double &t);
extern double Exact_Solution_mz_2DEB_Bucket(const double &x, const double &z, const double &t);
extern double Exact_Solution_r_2DEB_Bucket(const double &x, const double &z, const double &t);
extern double Exact_Solution_p_2DEB_Bucket(const double &x, const double &z, const double &t);
extern double rho_0_2DEB_Bucket(const double &z, const double &beta, const double &Fr);
extern double rho_0_deriv_2DEB_Bucket(const double &z, const double &beta, const double &Fr);
extern double N_2_2DEB_Bucket(const double &z, const double &beta, const double &Fr);
/*--------------------------------------------------------------------------*/
// EB - Exact Solution 3D
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_3DEB(const double &beta, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_mx_3DEB(const double &x, const double &y, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_my_3DEB(const double &x, const double &y, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_mz_3DEB(const double &x, const double &y, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_r_3DEB(const double &x, const double &y, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr);
extern double Exact_Solution_p_3DEB(const double &x, const double &y, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr);
extern double rho_0_3DEB(const double &z, const double &beta, const double &Fr);
extern double rho_0_deriv_3DEB(const double &z, const double &beta, const double &Fr);
extern double N_2_3DEB(const double &z, const double &beta, const double &Fr);
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif
