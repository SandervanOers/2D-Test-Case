# include "initial_cond.hpp"
/*--------------------------------------------------------------------------*/
// Stratified Test Problem
/*--------------------------------------------------------------------------*/
extern double calculate_sigma(const double &N2, const unsigned int &k)
{
    return sqrt(1/4.0*N2*N2+2.0*PETSC_PI*2.0*PETSC_PI*k*k)/2.0/PETSC_PI;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k)
{
    return exp(-0.5*N2*x)*sin(2.0*PETSC_PI*k*x)*sin(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k)
{
    return exp(-0.5*N2*x)*(N2/4.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.1));
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
extern double N_2(const double &x, const double &N2)
{
    // Can generalize to nonconst N2;
    return N2;//
}
/*--------------------------------------------------------------------------*/
// Viscous Homogeneous Test Problem
/*--------------------------------------------------------------------------*/
extern double viscous_calculate_sigma(const double &nu, const unsigned int &k)
{
    return k*sqrt(1.0-PETSC_PI*PETSC_PI*k*k*nu*nu);
}
/*--------------------------------------------------------------------------*/
extern double viscous_Exact_Solution_m(const double &x, const double &t, const double &nu, const double &sigma, const unsigned int &k)
{
    return sin(2.0*PETSC_PI*k*x)*exp(-2.0*k*k*PETSC_PI*PETSC_PI*nu*t)*sin(2*PETSC_PI*sigma*t);
}
/*--------------------------------------------------------------------------*/
extern double viscous_Exact_Solution_p(const double &x, const double &t, const double &nu, const double &sigma, const unsigned int &k)
{
    return cos(2.0*PETSC_PI*k*x)*exp(-2.0*k*k*PETSC_PI*PETSC_PI*nu*t)*(k*nu*PETSC_PI*sin(2*PETSC_PI*sigma*t)+sigma/((double)k)*cos(2.0*PETSC_PI*sigma*t));
}
/*--------------------------------------------------------------------------*/
extern double viscous_rho_0(const double &x, const double &N2)
{
    return 1.0;
}
/*--------------------------------------------------------------------------*/
extern double viscous_rho_0_deriv(const double &x, const double &N2)
{
    return 0.0;
}
/*--------------------------------------------------------------------------*/
// System 1 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_system1(const double &N2, const unsigned int &k)
{
    return sqrt((N2+1.0)/N2*(0.25*N2*N2-N2+1.0+4.0*PETSC_PI*PETSC_PI*k*k))/2.0/PETSC_PI;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m_system1(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k)
{
    return exp(-0.5*N2*x)*sin(2.0*PETSC_PI*k*x)*sin(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_system1(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k)
{
    return exp(-0.5*N2*x)*((N2/2.0-1.0)/2.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
// System 3 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_system3(const double &beta, const unsigned int &k)
{
    double N2 = beta-1.0;
    return sqrt((N2+1.0)/N2*(0.25*N2*N2-0.5*N2+0.25+4.0*PETSC_PI*PETSC_PI*k*k))/2.0/PETSC_PI;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m_system3(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*sin(2.0*PETSC_PI*k*x)*sin(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_system3(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*((N2-1.0)/4.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
extern double rho_0_compressible(const double &x, const double &beta)
{
    return exp(-beta*x);
}
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv_compressible(const double &x, const double &beta)
{
    return -beta*exp(-beta*x);
}
/*--------------------------------------------------------------------------*/
extern double N_2_compressible(const double &x, const double &beta)
{
    return beta-1.0;
}
/*--------------------------------------------------------------------------*/
// System 4 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_system4(const double &beta, const unsigned int &k)
{
    double N2 = beta-1.0;
    return sqrt(4.0*PETSC_PI*PETSC_PI*k*k*(1.0+2.0/N2)+0.25*N2*N2+0.25*(1.0+2.0/N2))/2.0/PETSC_PI;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m_system4(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*sin(2.0*PETSC_PI*k*x)*sin(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_system4(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*((N2+1.0)/4.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_system4(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*((N2-1.0)/4.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
// System 5 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_system5(const double &beta, const unsigned int &k)
{
    double N2 = beta-1.0;
    return sqrt(4.0*PETSC_PI*PETSC_PI*k*k*(1.0+2.0/N2)+0.25*N2*N2+0.25*(1.0+2.0/N2))/2.0/PETSC_PI;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m_system5(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*sin(2.0*PETSC_PI*k*x)*sin(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_system5(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*(-(N2+1.0)/4.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_system5(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*((N2-1.0)/4.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
// System 6: Verification: 1D Compressible Nonviscous Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_system1dcom(const double &beta, const unsigned int &k)
{
    double N2 = beta-1.0;
    return sqrt(4.0*PETSC_PI*PETSC_PI*k*k+0.25*N2*N2+0.5*N2+0.25)/2.0/PETSC_PI;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m_system1dcom(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*sin(2.0*PETSC_PI*k*x)*sin(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_system1dcom(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*(-(N2+1.0)/4.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_system1dcom(const double &x, const double &t, const double &beta, const double &sigma, const unsigned int &k)
{
    double N2 = beta-1.0;
    return exp(-0.5*beta*x)*((N2-1.0)/4.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.1));
}
/*--------------------------------------------------------------------------*/
// 2D
/*--------------------------------------------------------------------------*/
// System 1 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_2D_system1(const double &beta, const unsigned int &kx, const unsigned int &kz)
{
    return sqrt(beta*beta/16.0/PETSC_PI/PETSC_PI+kx*kx+kz*kz);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mx_2D_system1(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    //return exp(-0.5*beta*z)*kx/(pow(sigma,2.0)-pow(kx,2.0))*(beta/4.0/PETSC_PI*sin(2.0*PETSC_PI*kz*z)+kz*cos(2.0*PETSC_PI*kz*z))*sin(2.0*PETSC_PI*kx*x)*sin(2*PETSC_PI*sigma*t);
    return 2.0*PETSC_PI*kx*cos(2.0*PETSC_PI*kz*z)*sin(2.0*PETSC_PI*kx*x)*sin(2*PETSC_PI*sigma*t);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mz_2D_system1(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    //return sin(2.0*PETSC_PI*kz*z)*cos(2.0*PETSC_PI*kx*x)*sin(2*PETSC_PI*sigma*t);
    return 2.0*PETSC_PI*kz*sin(2.0*PETSC_PI*kz*z)*cos(2.0*PETSC_PI*kx*x)*sin(2*PETSC_PI*sigma*t);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_2D_system1(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    return 0;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_2D_system1(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    //return exp(-0.5*beta*z)*sigma/(pow(sigma,2.0)-pow(kx,2.0))*(beta/4.0/PETSC_PI*sin(2.0*PETSC_PI*kz*z)+kz*cos(2.0*PETSC_PI*kz*z))*cos(2.0*PETSC_PI*kx*x)*cos(2*PETSC_PI*sigma*t);
    return cos(2.0*PETSC_PI*kz*z)*cos(2.0*PETSC_PI*kx*x)*cos(2*PETSC_PI*sigma*t);
}
/*--------------------------------------------------------------------------*/
extern double rho_0_2D_system1(const double &z, const double &beta)
{
    return exp(-beta*z);
}
/*--------------------------------------------------------------------------*/
extern std::vector<double> rho_0_2D_system1(const std::vector<double> &z, const double &beta)
{
    std::vector<double> returnvector;
    for (const double& i : z)
    {
        returnvector.push_back(exp(-beta*(i)));
        //returnvector.push_back(1.0+i+i*i);
    }
    return returnvector;
}
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv_2D_system1(const double &z, const double &beta)
{
    return -beta*exp(-beta*z);
}
/*--------------------------------------------------------------------------*/
extern std::vector<double> rho_0_deriv_2D_system1(const std::vector<double> &z, const double &beta)
{
    std::vector<double> returnvector;
    for (const double& i : z)
    {
        returnvector.push_back(-beta*exp(-beta*(i)));
        //returnvector.push_back(1.0+2.0*i);
    }
    return returnvector;
}
/*--------------------------------------------------------------------------*/
extern double N_2_2D_system1(const double &z, const double &beta)
{
    return beta-1.0;
}
/*--------------------------------------------------------------------------*/
extern std::vector<double> N_2_2D_system1(const std::vector<double> &z, const double &beta)
{
    std::vector<double> returnvector;
    for (const double& i : z)
    {
        returnvector.push_back(beta-1.0);
    }
    return returnvector;
}
/*--------------------------------------------------------------------------*/
// System 2 Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_2D_system2(const double &beta, const unsigned int &kx, const unsigned int &kz)
{
    //return sqrt(beta*beta/16.0/PETSC_PI/PETSC_PI+kx*kx+kz*kz);
    double N2 = beta - 1;
    double kxd = (double)kx;
    double kzd = (double)kz;
    //return sqrt(beta*beta*beta/N2/16.0/PETSC_PI/PETSC_PI+beta/N2*(kx*kx+kz*kz));
    //return sqrt((beta-2.0)*(beta-2.0)/16.0/PETSC_PI/PETSC_PI + kx*kx + kz * kz);
    //return sqrt((N2+1.0)/N2*(kx*kx+kz*kz+(N2-1.0)*(N2-1.0)/16.0/PETSC_PI/PETSC_PI));
    double b = -kxd*kxd-kzd*kzd-(N2+1.0)*(N2+1.0)/16.0/PETSC_PI/PETSC_PI;
    double c = kxd*kxd*N2/4.0/PETSC_PI/PETSC_PI;
    return sqrt(-b/2.0+0.5*sqrt(b*b-4.0*c));
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mx_2D_system2(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    //return exp(-0.5*beta*z)*kx/(pow(sigma,2.0)-pow(kx,2.0))*(beta/4.0/PETSC_PI*sin(2.0*PETSC_PI*kz*z)+kz*cos(2.0*PETSC_PI*kz*z))*sin(2.0*PETSC_PI*kx*x)*sin(2*PETSC_PI*sigma*t);
    //return exp(-0.5*beta*z)*(4.0*PETSC_PI*kx)/(16.0*PETSC_PI*PETSC_PI*kz*kz+beta*beta)*(beta*sin(2.0*PETSC_PI*kz*z)+4.0*PETSC_PI*kz*cos(2.0*PETSC_PI*kz*z))*sin(2.0*PETSC_PI*kx*x)*sin(2*PETSC_PI*sigma*t);
    //return exp(-0.5*beta*z)*(4.0*PETSC_PI*kx)/(16.0*PETSC_PI*PETSC_PI*kz*kz+(beta-2.0)*(beta-2.0))*((beta-2.0)*sin(2.0*PETSC_PI*kz*z)+4.0*PETSC_PI*kz*cos(2.0*PETSC_PI*kz*z))*sin(2.0*PETSC_PI*kx*x)*sin(2*PETSC_PI*sigma*t);
    return 0.0;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mz_2D_system2(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    return exp(-0.5*beta*z)*sin(2.0*PETSC_PI*kz*z)*cos(2.0*PETSC_PI*kx*x)*sin(2*PETSC_PI*sigma*t);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_2D_system2(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    double N2 = beta - 1;
    //return 0;
    return exp(-0.5*beta*z)*sigma/(kx*kx-sigma*sigma)*(((N2+1.0)/4.0/PETSC_PI-kx*kx*N2/2.0/PETSC_PI/sigma/sigma)*sin(2.0*PETSC_PI*kz*z)-kz*cos(2.0*PETSC_PI*kz*z))*cos(2.0*PETSC_PI*kx*x)*cos(2.0*PETSC_PI*sigma*t);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_2D_system2(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    double N2 = beta - 1;
    //return exp(-0.5*beta*z)*sigma/(pow(sigma,2.0)-pow(kx,2.0))*(beta/4.0/PETSC_PI*sin(2.0*PETSC_PI*kz*z)+kz*cos(2.0*PETSC_PI*kz*z))*cos(2.0*PETSC_PI*kx*x)*cos(2*PETSC_PI*sigma*t);
    //return exp(-0.5*beta*z)*(4.0*PETSC_PI*sigma*N2)/(beta*(16.0*PETSC_PI*PETSC_PI*kz*kz+beta*beta))*(beta*sin(2.0*PETSC_PI*kz*z)+4.0*PETSC_PI*kz*cos(2.0*PETSC_PI*kz*z))*cos(2.0*PETSC_PI*kx*x)*cos(2*PETSC_PI*sigma*t);
    //return exp(-0.5*beta*z)*(1.0)/(beta*(16.0*PETSC_PI*PETSC_PI*kz*kz+(beta-2.0)*(beta-2.0)))*(4.0*PETSC_PI*sigma*(beta-2.0)*sin(2.0*PETSC_PI*kz*z)+16.0*PETSC_PI*PETSC_PI*sigma*kz*cos(2.0*PETSC_PI*kz*z))*cos(2.0*PETSC_PI*kx*x)*cos(2*PETSC_PI*sigma*t);

    ///return exp(-0.5*beta*z)*(4.0*PETSC_PI*sigma*sqrt(beta))/(16.0*PETSC_PI*PETSC_PI*kz*kz+(beta-2.0)*(beta-2.0))*N2/sqrt(N2+1.0)*((beta-2.0)*sin(2.0*PETSC_PI*kz*z)+4.0*PETSC_PI*kz*cos(2.0*PETSC_PI*kz*z))*cos(2.0*PETSC_PI*kx*x)*cos(2*PETSC_PI*sigma*t);
    //return exp(-0.5*beta*z)*1.0/(4.0*sigma*PETSC_PI*(16.0*PETSC_PI*PETSC_PI*kz*kz+(beta-2.0)*(beta-2.0))) *((16.0*PETSC_PI*PETSC_PI*(kx*kx+kz*kz)+(beta-2.0)*(beta-2.0))*(beta-2.0)*sin(2.0*PETSC_PI*kz*z)+(16.0*PETSC_PI*PETSC_PI*(kx*kx+kz*kz)+(beta-2.0)*(beta-2.0))* 4.0*PETSC_PI*kz*cos(2.0*PETSC_PI*kz*z))*cos(2.0*PETSC_PI*kx*x)*cos(2*PETSC_PI*sigma*t);
    return exp(-0.5*beta*z)*sigma/(kx*kx-sigma*sigma)*(-(N2-1.0)/4.0/PETSC_PI*sin(2.0*PETSC_PI*kz*z)-kz*cos(2.0*PETSC_PI*kz*z))*cos(2.0*PETSC_PI*kx*x)*cos(2.0*PETSC_PI*sigma*t);

}
/*--------------------------------------------------------------------------*/
extern double rho_0_2D_system2(const double &z, const double &beta)
{
    return exp(-beta*z);
}
/*--------------------------------------------------------------------------*/
extern std::vector<double> rho_0_2D_system2(const std::vector<double> &z, const double &beta)
{
    std::vector<double> returnvector;
    for (const double& i : z)
    {
        returnvector.push_back(exp(-beta*(i)));
    }
    return returnvector;
}
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv_2D_system2(const double &z, const double &beta)
{
    return -beta*exp(-beta*z);
}
/*--------------------------------------------------------------------------*/
extern std::vector<double> rho_0_deriv_2D_system2(const std::vector<double> &z, const double &beta)
{
    std::vector<double> returnvector;
    for (const double& i : z)
    {
        returnvector.push_back(-beta*exp(-beta*(i)));
    }
    return returnvector;
}
/*--------------------------------------------------------------------------*/
extern double N_2_2D_system2(const double &z, const double &beta)
{
    return beta-1.0;
}
/*--------------------------------------------------------------------------*/
extern std::vector<double> N_2_2D_system2(const std::vector<double> &z, const double &beta)
{
    std::vector<double> returnvector;
    for (const double& i : z)
    {
        returnvector.push_back(beta-1.0);
    }
    return returnvector;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
//  Incompresisble Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_2DIC(const double &beta, const unsigned int &kx, const unsigned int &kz)
{
    double N2 = beta;
    double kxd = (double)kx;
    double kzd = (double)kz;
    double b = 16.0*PETSC_PI*PETSC_PI*(kxd*kxd+kzd*kzd)+N2*N2;
    double c = 4.0*N2*kxd*kxd;
    return sqrt(c/b);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mx_2DIC(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    double N2 = beta;
    double kxd = (double)kx;
    double kzd = (double)kz;
    return exp(-0.5*beta*z)*(-N2/(2.0*2.0*PETSC_PI*kxd)*sin(2.0*PETSC_PI*kzd*z)-kzd/kxd*cos(2.0*PETSC_PI*kzd*z))*sin(2.0*PETSC_PI*kxd*x)*sin(sigma*t+0.1);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mz_2DIC(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    double kxd = (double)kx;
    double kzd = (double)kz;
    return exp(-0.5*beta*z)*sin(2.0*PETSC_PI*kzd*z)*cos(2.0*PETSC_PI*kxd*x)*sin(2*PETSC_PI*sigma*t+0.1);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_2DIC(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    double N2 = beta;
    double kxd = (double)kx;
    double kzd = (double)kz;
    return -N2/(2.0*PETSC_PI*sigma)*exp(-0.5*beta*z)*sin(2.0*PETSC_PI*kzd*z)*cos(2.0*PETSC_PI*kxd*x)*cos(2*PETSC_PI*sigma*t+0.1);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_2DIC(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz)
{
    double N2 = beta;
    double kxd = (double)kx;
    double kzd = (double)kz;
    return exp(-0.5*beta*z)*(-N2*sigma/(2.0*2.0*PETSC_PI*kxd*kxd)*sin(2.0*PETSC_PI*kzd*z)-kz*sigma/kxd/kxd*cos(2.0*PETSC_PI*kzd*z))*cos(2.0*PETSC_PI*kxd*x)*cos(2*PETSC_PI*sigma*t+0.1);
}
/*--------------------------------------------------------------------------*/
extern double rho_0_2DIC(const double &z, const double &beta)
{
    return exp(-beta*z);
}
/*--------------------------------------------------------------------------*/
/*extern std::vector<double> rho_0_2DIC(const std::vector<double> &z, const double &beta)
{
    std::vector<double> returnvector;
    for (const double& i : z)
    {
        returnvector.push_back(exp(-beta*(i)));
    }
    return returnvector;
}*/
/*--------------------------------------------------------------------------*/
extern double rho_0_deriv_2DIC(const double &z, const double &beta)
{
    return -beta*exp(-beta*z);
}
/*--------------------------------------------------------------------------*/
/*extern std::vector<double> rho_0_deriv_2DIC(const std::vector<double> &z, const double &beta)
{
    std::vector<double> returnvector;
    for (const double& i : z)
    {
        returnvector.push_back(-beta*exp(-beta*(i)));
    }
    return returnvector;
}*/
/*--------------------------------------------------------------------------*/
extern double N_2_2DIC(const double &z, const double &beta)
{
    return beta;
}
/*--------------------------------------------------------------------------*/
/*extern std::vector<double> N_2_2DIC(const std::vector<double> &z, const double &beta)
{
    std::vector<double> returnvector;
    for (const double& i : z)
    {
        returnvector.push_back(beta);
    }
    return returnvector;
}*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
//  WA Case - Experiments
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_2DWA(const double &beta, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    return 0.0;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mx_2DWA(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    return -PETSC_PI*sin(PETSC_PI*x)*cos(PETSC_PI*z);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_mz_2DWA(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    return PETSC_PI*cos(PETSC_PI*x)*sin(PETSC_PI*z);
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_r_2DWA(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    return 0.0;
}
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p_2DWA(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    return 0.0;
}
/*--------------------------------------------------------------------------*/
extern double rho_0_2DWA(const double &z, const double &beta, const double &Fr)
{
    return 1.0;
}
extern double rho_0_deriv_2DWA(const double &z, const double &beta, const double &Fr)
{
    return Fr*Fr*(-beta);
}
extern double N_2_2DWA(const double &z, const double &beta, const double &Fr)
{
    return beta;
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
//  EB Case - Test Case
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_2DEB(const double &beta, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kzd = (double)kz;
    double N2 = N_2_2DEB(0.0, beta, Fr);
    return sqrt(Fr*Fr*N2*kxd/(4.0*PETSC_PI*PETSC_PI*(kxd*kxd+kzd*kzd)));
}
extern double period_2DEB(const double &beta, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    double sigma = calculate_sigma_2DEB(beta, kx, kz, Fr);
    return 1.0/sigma;
}
extern double Exact_Solution_mx_2DEB(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kzd = (double)kz;
    return -kzd/kxd*cos(2.0*PETSC_PI*kzd*z)*sin(2.0*PETSC_PI*kxd*x)*sin(2.0*PETSC_PI*sigma*t);
}
extern double Exact_Solution_mz_2DEB(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kzd = (double)kz;
    return sin(2.0*PETSC_PI*kzd*z)*cos(2.0*PETSC_PI*kxd*x)*sin(2.0*PETSC_PI*sigma*t);
}
extern double Exact_Solution_r_2DEB(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kzd = (double)kz;
    double N2 = N_2_2DEB(z, beta, Fr);
    return -Fr*Fr*N2/sigma/2.0/PETSC_PI*sin(2.0*PETSC_PI*kzd*z)*cos(2.0*PETSC_PI*kxd*x)*cos(2.0*PETSC_PI*sigma*t);
}
extern double Exact_Solution_p_2DEB(const double &x, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kzd = (double)kz;
    return -kzd*sigma/kxd/kxd*cos(2.0*PETSC_PI*kzd*z)*cos(2.0*PETSC_PI*kxd*x)*cos(2.0*PETSC_PI*sigma*t);
}
extern double rho_0_2DEB(const double &z, const double &beta, const double &Fr)
{
    return 1.0;
}
extern double rho_0_deriv_2DEB(const double &z, const double &beta, const double &Fr)
{
    return Fr*Fr*(-beta);
}
extern double N_2_2DEB(const double &z, const double &beta, const double &Fr)
{
    return beta;
}
/*--------------------------------------------------------------------------*/
// EB - Bucket Test Case
extern double calculate_sigma_2DEB_Bucket()
{
    double N2 = 1.0;
    return 0.5*sqrt(2.0);
}
extern double period_2DEB_Bucket()
{
    double sigma = calculate_sigma_2DEB_Bucket();
    return 2.0*PETSC_PI/sigma;
}
extern double Exact_Solution_mx_2DEB_Bucket(const double &x, const double &z, const double &t)
{
    double u;
    if (x-0.5 <= z && -x+1.5 <= z) // region I
    {
        u = 0.2e1 * ((double) (int) pow((double) (x - 1), (double) 2) + pow(z - 0.1e1 / 0.2e1, 0.2e1) - 0.2e1 * z + 0.1e1) * sin(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 > z && -x+1.5 < z) // region II
    {
        u = (-(double) (7 * (int) pow((double) (x - 1), (double) 2)) + 0.18e2 * (double) (x - 1) * (z - 0.1e1 / 0.2e1) - 0.7e1 * pow(z - 0.1e1 / 0.2e1, 0.2e1) + (double) (6 * x) - 0.1e1 - 0.10e2 * z) * sin(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 < z && -x+1.5 > z) // region III
    {
        u = (-(double) (7 * (int) pow((double) (x - 1), (double) 2)) - 0.18e2 * (double) (x - 1) * (z - 0.1e1 / 0.2e1) - 0.7e1 * pow(z - 0.1e1 / 0.2e1, 0.2e1) - (double) (6 * x) + 0.11e2 - 0.10e2 * z) * sin(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 >= z && -x+1.5 >= z) // region IV
    {
        u = -0.16e2 * ((double) (int) pow((double) (x - 1), (double) 2) + pow(z - 0.1e1 / 0.2e1, 0.2e1) + z - 0.1e1 / 0.2e1) * sin(sqrt(0.2e1) * t / 0.2e1);
    }
    else
    {
        std::cout << "Error in Initial Condition u" << std::endl;
        std::cout << "x = " << x << ", z = " << z << std::endl;
        u = 0.0;
    }
    return u;
}
extern double Exact_Solution_mz_2DEB_Bucket(const double &x, const double &z, const double &t)
{
    double w;
    if (x-0.5 <= z && -x+1.5 <= z) // region I
    {
        w = -0.4e1 * (double) (x - 1) * (z - 0.3e1 / 0.2e1) * sin(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 > z && -x+1.5 < z) // region II
    {
        w = (-(double) (9 * (int) pow((double) (x - 1), (double) 2)) + 0.14e2 * (double) (x - 1) * (z - 0.1e1 / 0.2e1) - 0.9e1 * pow(z - 0.1e1 / 0.2e1, 0.2e1) + (double) (10 * x) - 0.7e1 - 0.6e1 * z) * sin(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 < z && -x+1.5 > z) // region III
    {
        w = ((double) (9 * (int) pow((double) (x - 1), (double) 2)) + 0.14e2 * (double) (x - 1) * (z - 0.1e1 / 0.2e1) + 0.9e1 * pow(z - 0.1e1 / 0.2e1, 0.2e1) + (double) (10 * x) - 0.13e2 + 0.6e1 * z) * sin(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 >= z && -x+1.5 >= z) // region IV
    {
        w = 0.32e2 * z * (double) (x - 1) * sin(sqrt(0.2e1) * t / 0.2e1);
    }
    else
    {
        std::cout << "Error in Initial Condition w" << std::endl;
        w = 0.0;
    }
    return w;
}
extern double Exact_Solution_r_2DEB_Bucket(const double &x, const double &z, const double &t)
{
    double r;
    if (x-0.5 <= z && -x+1.5 <= z) // region I
    {
        r = 0.4e1 * sqrt(0.2e1) * (double) (x - 1) * (z - 0.3e1 / 0.2e1) * cos(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 > z && -x+1.5 < z) // region II
    {
        //r = -sqrt(0.2e1) * (-(double) (9 * (int) pow((double) (x - 1), (double) 2)) + 0.14e2 * (double) (x - 1) * (z - 0.1e1 / 0.2e1) - 0.9e1 * pow(z - 0.1e1 / 0.2e1, 0.2e1) + (double) (10 * x) - 0.7e1 - 0.6e1 * z) * cos(sqrt(0.2e1) * t / 0.2e1);
        r = -sqrt(2.0) * (-(double) (9.0 *  pow((double) (x - 1.0), (double) 2.0)) + 0.14e2 * (double) (x - 1.0) * (z - 0.5) - 0.9e1 * pow(z - 0.5, 2.0) + (double) (10.0 * x) - 0.7e1 - 0.6e1 * z) * cos(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 < z && -x+1.5 > z) // region III
    {
        //r = -sqrt(0.2e1) * ((double) (9 * (int) pow((double) (x - 1), (double) 2)) + 0.14e2 * (double) (x - 1) * (z - 0.1e1 / 0.2e1) + 0.9e1 * pow(z - 0.1e1 / 0.2e1, 0.2e1) + (double) (10 * x) - 0.13e2 + 0.6e1 * z) * cos(sqrt(0.2e1) * t / 0.2e1);
        r = -sqrt(0.2e1) * ((double) (9.0 * pow((double) (x - 1), (double) 2)) + 0.14e2 * (double) (x - 1) * (z - 0.1e1 / 0.2e1) + 0.9e1 * pow(z - 0.1e1 / 0.2e1, 0.2e1) + (double) (10 * x) - 0.13e2 + 0.6e1 * z) * cos(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 >= z && -x+1.5 >= z) // region IV
    {
        r = -0.32e2 * sqrt(0.2e1) * z * (double) (x - 1) * cos(sqrt(0.2e1) * t / 0.2e1);
    }
    else
    {
        std::cout << "Error in Initial Condition r" << std::endl;
        r = 0.0;
    }
    return r;
}
extern double Exact_Solution_p_2DEB_Bucket(const double &x, const double &z, const double &t)
{
    double p;
    if (x-0.5 <= z && -x+1.5 <= z) // region I
    {
        p = -sqrt(0.2e1) * cos(sqrt(0.2e1) * t / 0.2e1) * (double) (x - 1.0) * ((double)  pow((double) (x - 1.0), (double) 2.0) + 0.3e1 * pow(z - 0.1e1 / 0.2e1, 0.2e1) - 0.6e1 * z + 0.3e1) / 0.3e1;
    }
    else if (x-0.5 > z && -x+1.5 < z) // region II
    {
        p = sqrt(0.2e1) * cos(sqrt(0.2e1) * t / 0.2e1) * (double) (x - 1.0) * ((double) (7.0 * pow((double) (x - 1.0), (double) 2.0)) - 0.27e2 * (double) (x - 1) * (z - 0.1e1 / 0.2e1) + 0.21e2 * pow(z - 0.1e1 / 0.2e1, 0.2e1) - (double) (9 * x) - 0.6e1 + 0.30e2 * z) / 0.6e1 - 0.3e1 / 0.2e1 * sqrt(0.2e1) * (pow(z - 0.1e1 / 0.2e1, 0.3e1) + pow(z - 0.1e1 / 0.2e1, 0.2e1)) * cos(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 < z && -x+1.5 > z) // region III
    {
        p = sqrt(0.2e1) * cos(sqrt(0.2e1) * t / 0.2e1) * (double) (x - 1.0) * ((double) (7.0 * pow((double) (x - 1.0), (double) 2.0)) + 0.27e2 * (double) (x - 1) * (z - 0.1e1 / 0.2e1) + 0.21e2 * pow(z - 0.1e1 / 0.2e1, 0.2e1) + (double) (9 * x) - 0.24e2 + 0.30e2 * z) / 0.6e1 + 0.3e1 / 0.2e1 * sqrt(0.2e1) * (pow(z - 0.1e1 / 0.2e1, 0.3e1) + pow(z - 0.1e1 / 0.2e1, 0.2e1)) * cos(sqrt(0.2e1) * t / 0.2e1);
    }
    else if (x-0.5 >= z && -x+1.5 >= z) // region IV
    {
        p = 0.8e1 / 0.3e1 * sqrt(0.2e1) * cos(sqrt(0.2e1) * t / 0.2e1) * (double) (x - 1.0) * ((double) pow((double) (x - 1.0), (double) 2.0) + 0.3e1 * pow(z - 0.1e1 / 0.2e1, 0.2e1) + 0.3e1 * z - 0.3e1 / 0.2e1);
    }
    else
    {
        std::cout << "Error in Initial Condition p" << std::endl;
        p = 0.0;
    }
    return p;
}
extern double rho_0_2DEB_Bucket(const double &z, const double &beta, const double &Fr)
{
    return 1.0;
}
extern double rho_0_deriv_2DEB_Bucket(const double &z, const double &beta, const double &Fr)
{
    return -1.0;
}
extern double N_2_2DEB_Bucket(const double &z, const double &beta, const double &Fr)
{
    return 1.0;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// EB - Exact Solution 3D
/*--------------------------------------------------------------------------*/
extern double calculate_sigma_3DEB(const double &beta, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kyd = (double)ky;
    double kzd = (double)kz;
    double N2 = N_2_3DEB(0.0, beta, Fr);
    return sqrt(Fr*Fr*N2*(kxd*kxd+kyd*kyd)/(4.0*PETSC_PI*PETSC_PI*(kxd*kxd+kyd*kyd+kzd*kzd)));
}
extern double Exact_Solution_mx_3DEB(const double &x, const double &y, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kyd = (double)ky;
    double kzd = (double)kz;
    double N2 = N_2_3DEB(0.0, beta, Fr);
    return (sigma*sigma-N2/(4.0*PETSC_PI*PETSC_PI))*kxd/kzd*sin(2.0*PETSC_PI*kxd*x)*cos(2.0*PETSC_PI*kyd*y)*cos(2.0*PETSC_PI*kzd*z)*sin(2.0*PETSC_PI*sigma*t);
}
extern double Exact_Solution_my_3DEB(const double &x, const double &y, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kyd = (double)ky;
    double kzd = (double)kz;
    double N2 = N_2_3DEB(0.0, beta, Fr);
    return (sigma*sigma-N2/(4.0*PETSC_PI*PETSC_PI))*kyd/kzd*cos(2.0*PETSC_PI*kxd*x)*sin(2.0*PETSC_PI*kyd*y)*cos(2.0*PETSC_PI*kzd*z)*sin(2.0*PETSC_PI*sigma*t);
}
extern double Exact_Solution_mz_3DEB(const double &x, const double &y, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kyd = (double)ky;
    double kzd = (double)kz;
    double N2 = N_2_3DEB(0.0, beta, Fr);
    return cos(2.0*PETSC_PI*kxd*x)*cos(2.0*PETSC_PI*kyd*y)*sin(2.0*PETSC_PI*kzd*z)*sin(2.0*PETSC_PI*sigma*t);
}
extern double Exact_Solution_r_3DEB(const double &x, const double &y, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kyd = (double)ky;
    double kzd = (double)kz;
    double N2 = N_2_3DEB(0.0, beta, Fr);
    return -N2/(2.0*PETSC_PI*sigma)*cos(2.0*PETSC_PI*kxd*x)*cos(2.0*PETSC_PI*kyd*y)*sin(2.0*PETSC_PI*kzd*z)*cos(2.0*PETSC_PI*sigma*t);
}
extern double Exact_Solution_p_3DEB(const double &x, const double &y, const double &z, const double &t, const double &beta, const double &sigma, const unsigned int &kx, const unsigned int &ky, const unsigned int &kz, const double &Fr)
{
    double kxd = (double)kx;
    double kyd = (double)ky;
    double kzd = (double)kz;
    double N2 = N_2_3DEB(0.0, beta, Fr);
    return (sigma*sigma-N2/(4.0*PETSC_PI*PETSC_PI))*1.0/(kzd*sigma)*cos(2.0*PETSC_PI*kxd*x)*cos(2.0*PETSC_PI*kyd*y)*cos(2.0*PETSC_PI*kzd*z)*cos(2.0*PETSC_PI*sigma*t);
}
extern double rho_0_3DEB(const double &z, const double &beta, const double &Fr)
{
    return 1.0;
}
extern double rho_0_deriv_3DEB(const double &z, const double &beta, const double &Fr)
{
    return Fr*Fr*(-beta);
}
extern double N_2_3DEB(const double &z, const double &beta, const double &Fr)
{
    return beta;
}
/*--------------------------------------------------------------------------*/
