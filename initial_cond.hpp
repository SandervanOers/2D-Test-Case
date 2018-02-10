#ifndef INITIAL_CONDITION
#define INITIAL_CONDITION

#include <petscksp.h>
/*--------------------------------------------------------------------------*/
extern double calculate_sigma(const double &N2, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_m(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
extern double Exact_Solution_p(const double &x, const double &t, const double &N2, const double &sigma, const unsigned int &k);
/*--------------------------------------------------------------------------*/
#endif
