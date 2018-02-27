# include "initial_cond.hpp"
/*--------------------------------------------------------------------------*/
extern double calculate_sigma(const double &N2, const unsigned int &k)
{
    return sqrt(1/4.0*N2*N2+2.0*PETSC_PI*2.0*PETSC_PI*k*k)/2.0/PETSC_PI;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k)
{
    return exp(-0.5*N2*x)*sin(2.0*PETSC_PI*k*x)*sin(2.0*PETSC_PI*sigma*(t+0.0));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k)
{
    return exp(-0.5*N2*x)*(N2/4.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.0));
}
/*--------------------------------------------------------------------------*/
extern double rho_0(const double &x, const double &N2)
{
    return exp(-N2*x);
}
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv(const double &x, const double &N2)
{
    return -N2*exp(-N2*x);
}
/*--------------------------------------------------------------------------*/
