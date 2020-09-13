#ifndef CubatureData2D_H__INCLUDED
#define CubatureData2D_H__INCLUDED

#include "CubatureData2D.hpp"


#include <petscksp.h>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
/*--------------------------------------------------------------------------*/
extern void Cubature2D(const unsigned int &Corder, Vec &R, Vec &S, Vec &W, unsigned int &Ncub);
/*--------------------------------------------------------------------------*/










#endif  // CubatureData2D_H__INCLUDED
