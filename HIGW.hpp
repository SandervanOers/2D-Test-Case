#ifndef HIGW
#define HIGW

#include <petscksp.h>
#include <iostream>
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix_local(const Mat &V, const double &deltaX, const unsigned int &N);
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix_inverse_local(const Mat &V, const double &deltaX);
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian(const Mat &M1, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np);
/*--------------------------------------------------------------------------*/
#endif
