// GNU General Public License Agreement
// Copyright (C) 2004-2010 CodeCogs, Zyba Ltd, Broadwood, Holford, TA5 1DU, England.
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by CodeCogs. 
// You must retain a copy of this licence in all copies. 
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
// ---------------------------------------------------------------------------------
//! Calculates the Airy functions Ai and Bi and their derivatives.

#ifndef MATHS_SPECIAL_AIRY_AIRY_H
#define MATHS_SPECIAL_AIRY_AIRY_H

//#include <math.h>
#include <cmath>
#include "poly_eval.hpp"

namespace Maths
{

namespace Special
{

namespace Airy
{

//! Calculates the Airy functions Ai and Bi and their derivatives.

int airy(double x, double &Ai, double &Aip, double &Bi, double &Bip );

}

}

}


#endif

