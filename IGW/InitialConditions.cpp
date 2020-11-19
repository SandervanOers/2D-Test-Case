/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include"InitialConditions.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/FaceCacheData.hpp"
#include "LinearAlgebra/Matrix.hpp"



EBSingleMode::EBSingleMode():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 2.0*Pi_;
	l_ = 2.0*Pi_;
	m_ = 2.0*Pi_;
	N2_ = 2.0;
	sigma_ = sqrt(((k_*k_+l_*l_)*N2_)/(k_*k_+l_*l_+m_*m_));
}

double
EBSingleMode::getU(const PointPhysicalT& pPhys, double t)const
{
	return -m_*k_/(k_*k_+l_*l_)*cos(m_*pPhys[2])*sin(k_*pPhys[0])*cos(l_*pPhys[1])*sin(sigma_*t+0.1);
	//cos(k_*pPhys[0]+l_*pPhys[1]-sigma_*t);
}

double
EBSingleMode::getV(const PointPhysicalT& pPhys, double t)const
{
	return -m_*l_/(k_*k_+l_*l_)*cos(m_*pPhys[2])*cos(k_*pPhys[0])*sin(l_*pPhys[1])*sin(sigma_*t+0.1);
	//cos(k_*pPhys[0]+l_*pPhys[1]-sigma_*t);
}

double
EBSingleMode::getW(const PointPhysicalT& pPhys, double t)const
{
	return sin(m_*pPhys[2])*cos(k_*pPhys[0])*cos(l_*pPhys[1])*sin(sigma_*t+0.1);
	//sin(k_*pPhys[0]+l_*pPhys[1]-sigma_*t);
}
		
double
EBSingleMode::getR(const PointPhysicalT& pPhys, double t)const
{
	return -N2_/sigma_*sin(m_*pPhys[2])*cos(k_*pPhys[0])*cos(l_*pPhys[1])*cos(sigma_*t+0.1);
	//cos(k_*pPhys[0]+l_*pPhys[1]-sigma_*t);
}
		
double
EBSingleMode::getP(const PointPhysicalT& pPhys, double t)const
{
	return	-m_*sigma_/(k_*k_+l_*l_)*cos(m_*pPhys[2])*cos(k_*pPhys[0])*cos(l_*pPhys[1])*cos(sigma_*t+0.1);
	//cos(k_*pPhys[0]+l_*pPhys[1]-sigma_*t);
}

double
EBSingleMode::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-N2_*pPhys[1]);
}

double
EBSingleMode::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-N2_*exp(-N2_*pPhys[1]);
}

double
EBSingleMode::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	N2_;
}

WA2D::WA2D():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	N2_ = 1.0;
	Fn_ = 1.0;
	Fm_ = 1.0;
}

double
WA2D::getU(const PointPhysicalT& pPhys, double t)const
{
	return -Fm_*Pi_*sin(Fn_*Pi_*pPhys[0])*cos(Fm_*Pi_*pPhys[1]);
}

double
WA2D::getV(const PointPhysicalT& pPhys, double t)const
{
	return Fn_*Pi_*cos(Fn_*Pi_*pPhys[0])*sin(Fm_*Pi_*pPhys[1]);
}

double
WA2D::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
WA2D::getR(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
WA2D::getP(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;
}

double
WA2D::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	1.0;
}

double
WA2D::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;
}

double
WA2D::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;
}

EBsemiellipse11mode::EBsemiellipse11mode():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 0.0;
	N2_ = 2.0;
	sigma_ = sqrt(0.5*N2_);
}

double
EBsemiellipse11mode::getU(const PointPhysicalT& pPhys, double t)const
{
	return (-3.0*(pPhys[0]*pPhys[0]+pPhys[1]*pPhys[1])+3.0)*sin(sigma_*t+0.1);
}

double
EBsemiellipse11mode::getV(const PointPhysicalT& pPhys, double t)const
{
	return 6.0*pPhys[0]*pPhys[1]*sin(sigma_*t+0.1);
}

double
EBsemiellipse11mode::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
EBsemiellipse11mode::getR(const PointPhysicalT& pPhys, double t)const
{
	return -N2_/sigma_*6.0*pPhys[0]*pPhys[1]*cos(sigma_*t+0.1);
}
		
double
EBsemiellipse11mode::getP(const PointPhysicalT& pPhys, double t)const
{
	return	sigma_*(3.0*pPhys[0]*pPhys[1]*pPhys[1]+pPhys[0]*pPhys[0]*pPhys[0]-3.0*pPhys[0])*cos(sigma_*t+0.1);
}

double
EBsemiellipse11mode::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-N2_*pPhys[1]);
}

double
EBsemiellipse11mode::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-N2_*exp(-N2_*pPhys[1]);
}

double
EBsemiellipse11mode::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	N2_;
}

EBbucketb::EBbucketb():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 0.0;
	N2_ = 2.0;
	sigma_ = sqrt(0.5*N2_);
}

double
EBbucketb::getU(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0)
			return ((pPhys[0]-pPhys[1]+1.0)*(pPhys[0]-pPhys[1]+1.0)-2.0+(pPhys[0]+pPhys[1]-1.0)*(pPhys[0]+pPhys[1]-1.0))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0)
			return (-2.0*(1.0-2.0*pPhys[0]+2.0*pPhys[1])*(1.0-2.0*pPhys[0]+2.0*pPhys[1])+4.0-2.0*(-1.0-2.0*pPhys[0]-2.0*pPhys[1])*(-1.0-2.0*pPhys[0]-2.0*pPhys[1]))*sin(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
		{
		if (pPhys[0] >= 0.0)
			return (-2.0*(1.0-2.0*pPhys[0]+2.0*pPhys[1])*(1.0-2.0*pPhys[0]+2.0*pPhys[1])+1.0+(pPhys[0]+pPhys[1]-1.0)*(pPhys[0]+pPhys[1]-1.0))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0)
			return ((pPhys[0]-pPhys[1]+1.0)*(pPhys[0]-pPhys[1]+1.0)+1.0-2.0*(-1.0-2.0*pPhys[0]-2.0*pPhys[1])*(-1.0-2.0*pPhys[0]-2.0*pPhys[1]))*sin(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for U" << endl;
}

double
EBbucketb::getV(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0)
			return ((pPhys[0]-pPhys[1]+1.0)*(pPhys[0]-pPhys[1]+1.0)-(pPhys[0]+pPhys[1]-1.0)*(pPhys[0]+pPhys[1]-1.0))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0)
			return (-2.0*(1.0-2.0*pPhys[0]+2.0*pPhys[1])*(1.0-2.0*pPhys[0]+2.0*pPhys[1])+2.0*(-1.0-2.0*pPhys[0]-2.0*pPhys[1])*(-1.0-2.0*pPhys[0]-2.0*pPhys[1]))*sin(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
		{
		if (pPhys[0] >= 0.0)
			return (-2.0*(1.0-2.0*pPhys[0]+2.0*pPhys[1])*(1.0-2.0*pPhys[0]+2.0*pPhys[1])+3.0-(pPhys[0]+pPhys[1]-1.0)*(pPhys[0]+pPhys[1]-1.0))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0)
			return ((pPhys[0]-pPhys[1]+1.0)*(pPhys[0]-pPhys[1]+1.0)-3.0+2.0*(-1.0-2.0*pPhys[0]-2.0*pPhys[1])*(-1.0-2.0*pPhys[0]-2.0*pPhys[1]))*sin(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for V (W)" << endl;
}

double
EBbucketb::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
EBbucketb::getR(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0)
			return -N2_/sigma_*((pPhys[0]-pPhys[1]+1.0)*(pPhys[0]-pPhys[1]+1.0)-(pPhys[0]+pPhys[1]-1.0)*(pPhys[0]+pPhys[1]-1.0))*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0)
			return -N2_/sigma_*(-2.0*(1.0-2.0*pPhys[0]+2.0*pPhys[1])*(1.0-2.0*pPhys[0]+2.0*pPhys[1])+2.0*(-1.0-2.0*pPhys[0]-2.0*pPhys[1])*(-1.0-2.0*pPhys[0]-2.0*pPhys[1]))*cos(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
		{
		if (pPhys[0] >= 0.0)
			return -N2_/sigma_*(-2.0*(1.0-2.0*pPhys[0]+2.0*pPhys[1])*(1.0-2.0*pPhys[0]+2.0*pPhys[1])+3.0-(pPhys[0]+pPhys[1]-1.0)*(pPhys[0]+pPhys[1]-1.0))*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0)
			return -N2_/sigma_*((pPhys[0]-pPhys[1]+1.0)*(pPhys[0]-pPhys[1]+1.0)-3.0+2.0*(-1.0-2.0*pPhys[0]-2.0*pPhys[1])*(-1.0-2.0*pPhys[0]-2.0*pPhys[1]))*cos(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for R" << endl;
}

double
EBbucketb::getP(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0)
			return -sigma_*((pPhys[0]-pPhys[1]+1.0)*(pPhys[0]-pPhys[1]+1.0)*(pPhys[0]-pPhys[1]+1.0)/3.0+(pPhys[0]+pPhys[1]-1.0)*(pPhys[0]+pPhys[1]-1.0)*(pPhys[0]+pPhys[1]-1.0)/3.0-2.0*pPhys[0])*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0)
			return -sigma_*((1.0-2.0*pPhys[0]+2.0*pPhys[1])*(1.0-2.0*pPhys[0]+2.0*pPhys[1])*(1.0-2.0*pPhys[0]+2.0*pPhys[1])/3.0+(-1.0-2.0*pPhys[0]-2.0*pPhys[1])*(-1.0-2.0*pPhys[0]-2.0*pPhys[1])*(-1.0-2.0*pPhys[0]-2.0*pPhys[1])/3.0 +4.0*pPhys[0])*cos(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
		{
		if (pPhys[0] >= 0.0)
			return -sigma_*((1.0-2.0*pPhys[0]+2.0*pPhys[1])*(1.0-2.0*pPhys[0]+2.0*pPhys[1])*(1.0-2.0*pPhys[0]+2.0*pPhys[1])/3.0+(pPhys[0]+pPhys[1]-1.0)*(pPhys[0]+pPhys[1]-1.0)*(pPhys[0]+pPhys[1]-1.0)/3.0+pPhys[0]-3.0*pPhys[1])*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0)
			return -sigma_*((pPhys[0]-pPhys[1]+1.0)*(pPhys[0]-pPhys[1]+1.0)*(pPhys[0]-pPhys[1]+1.0)/3.0+(-1.0-2.0*pPhys[0]-2.0*pPhys[1])*(-1.0-2.0*pPhys[0]-2.0*pPhys[1])*(-1.0-2.0*pPhys[0]-2.0*pPhys[1])/3.0+pPhys[0]+3.0*pPhys[1])*cos(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for P" << endl;
}

double
EBbucketb::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-N2_*pPhys[1]);
}

double
EBbucketb::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-N2_*exp(-N2_*pPhys[1]);
}

double
EBbucketb::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	N2_;
}

EBbucketd::EBbucketd(): // EBbucketd not confirmed yet
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 0.0;
	N2_ = 2.0;
	sigma_ = sqrt(0.5*N2_);
}

double
EBbucketd::getU(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0) // region I
			return ((-12.5*pPhys[0]+12.5*pPhys[1]-8.75)*exp(-(2.5*pPhys[0]-2.5*pPhys[1]+1.75)*(2.5*pPhys[0]-2.5*pPhys[1]+1.75))+(-12.5*pPhys[0]-12.5*pPhys[1]+16.25)*exp(-(2.5*pPhys[0]+2.5*pPhys[1]-3.25)*(2.5*pPhys[0]+2.5*pPhys[1]-3.25)))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0) // region IV
			return ((17.5-50.0*pPhys[0]+50.0*pPhys[1])*exp(-(-5.0*pPhys[0]+5.0*pPhys[1]+1.75)*(-5.0*pPhys[0]+5.0*pPhys[1]+1.75))+(-32.5-50.0*pPhys[0]-50.0*pPhys[1])*exp(-(5.0*pPhys[0]+5.0*pPhys[1]+3.25)*(5.0*pPhys[0]+5.0*pPhys[1]+3.25)))*sin(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
		{
		if (pPhys[0] >= 0.0) // region II
			return ((17.5-50.0*pPhys[0]+50.0*pPhys[1])*exp(-(-5.0*pPhys[0]+5.0*pPhys[1]+1.75)*(-5.0*pPhys[0]+5.0*pPhys[1]+1.75))+(16.25-12.5*pPhys[0]-12.5*pPhys[1])*exp(-(2.5*pPhys[0]+2.5*pPhys[1]-3.25)*(2.5*pPhys[0]+2.5*pPhys[1]-3.25)))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0) // region III
			return ((-8.75-12.5*pPhys[0]+12.5*pPhys[1])*exp(-(2.5*pPhys[0]-2.5*pPhys[1]+1.75)*(2.5*pPhys[0]-2.5*pPhys[1]+1.75))+(-32.5-50.0*pPhys[0]-50.0*pPhys[1])*exp(-(5.0*pPhys[0]+5.0*pPhys[1]+3.25)*(5.0*pPhys[0]+5.0*pPhys[1]+3.25)))*sin(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for U" << endl;
}

double
EBbucketd::getV(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0)
			return ((-12.5*pPhys[0]+12.5*pPhys[1]-8.75)*exp(-(2.5*pPhys[0]-2.5*pPhys[1]+1.75)*(2.5*pPhys[0]-2.5*pPhys[1]+1.75))+(12.5*pPhys[0]+12.5*pPhys[1]-16.25)*exp(-(2.5*pPhys[0]+2.5*pPhys[1]-3.25)*(2.5*pPhys[0]+2.5*pPhys[1]-3.25)))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0)
			return ((17.5-50.0*pPhys[0]+50.0*pPhys[1])*exp(-(-5.0*pPhys[0]+5.0*pPhys[1]+1.75)*(-5.0*pPhys[0]+5.0*pPhys[1]+1.75))+(32.5+50.0*pPhys[0]+50.0*pPhys[1])*exp(-(5.0*pPhys[0]+5.0*pPhys[1]+3.25)*(5.0*pPhys[0]+5.0*pPhys[1]+3.25)))*sin(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
		{
		if (pPhys[0] >= 0.0)
			return ((17.5-50.0*pPhys[0]+50.0*pPhys[1])*exp(-(-5.0*pPhys[0]+5.0*pPhys[1]+1.75)*(-5.0*pPhys[0]+5.0*pPhys[1]+1.75))+(-16.25+12.5*pPhys[0]+12.5*pPhys[1])*exp(-(2.5*pPhys[0]+2.5*pPhys[1]-3.25)*(2.5*pPhys[0]+2.5*pPhys[1]-3.25)))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0)
			return ((-8.75-12.5*pPhys[0]+12.5*pPhys[1])*exp(-(2.5*pPhys[0]-2.5*pPhys[1]+1.75)*(2.5*pPhys[0]-2.5*pPhys[1]+1.75))+(32.5+50.0*pPhys[0]+50.0*pPhys[1])*exp(-(5.0*pPhys[0]+5.0*pPhys[1]+3.25)*(5.0*pPhys[0]+5.0*pPhys[1]+3.25)))*sin(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for V (W)" << endl;
}

double
EBbucketd::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
EBbucketd::getR(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0)
			return -N2_/sigma_*((-12.5*pPhys[0]+12.5*pPhys[1]-8.75)*exp(-(2.5*pPhys[0]-2.5*pPhys[1]+1.75)*(2.5*pPhys[0]-2.5*pPhys[1]+1.75))+(12.5*pPhys[0]+12.5*pPhys[1]-16.25)*exp(-(2.5*pPhys[0]+2.5*pPhys[1]-3.25)*(2.5*pPhys[0]+2.5*pPhys[1]-3.25)))*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0)
			return -N2_/sigma_*((17.5-50.0*pPhys[0]+50.0*pPhys[1])*exp(-(-5.0*pPhys[0]+5.0*pPhys[1]+1.75)*(-5.0*pPhys[0]+5.0*pPhys[1]+1.75))+(32.5+50.0*pPhys[0]+50.0*pPhys[1])*exp(-(5.0*pPhys[0]+5.0*pPhys[1]+3.25)*(5.0*pPhys[0]+5.0*pPhys[1]+3.25)))*cos(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
	{
		if (pPhys[0] >= 0.0)
			return -N2_/sigma_*((17.5-50.0*pPhys[0]+50.0*pPhys[1])*exp(-(-5.0*pPhys[0]+5.0*pPhys[1]+1.75)*(-5.0*pPhys[0]+5.0*pPhys[1]+1.75))+(-16.25+12.5*pPhys[0]+12.5*pPhys[1])*exp(-(2.5*pPhys[0]+2.5*pPhys[1]-3.25)*(2.5*pPhys[0]+2.5*pPhys[1]-3.25)))*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0)
			return -N2_/sigma_*((-8.75-12.5*pPhys[0]+12.5*pPhys[1])*exp(-(2.5*pPhys[0]-2.5*pPhys[1]+1.75)*(2.5*pPhys[0]-2.5*pPhys[1]+1.75))+(32.5+50.0*pPhys[0]+50.0*pPhys[1])*exp(-(5.0*pPhys[0]+5.0*pPhys[1]+3.25)*(5.0*pPhys[0]+5.0*pPhys[1]+3.25)))*cos(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for R" << endl;
}

double
EBbucketd::getP(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0)
			return -sigma_*(exp(-(2.5*pPhys[0]-2.5*pPhys[1]+1.75)*(2.5*pPhys[0]-2.5*pPhys[1]+1.75))+exp(-(2.5*pPhys[0]+2.5*pPhys[1]-3.25)*(2.5*pPhys[0]+2.5*pPhys[1]-3.25)))*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0)
			return -sigma_*(exp(-(-5.0*pPhys[0]+5.0*pPhys[1]+1.75)*(-5.0*pPhys[0]+5.0*pPhys[1]+1.75))+exp(-(5.0*pPhys[0]+5.0*pPhys[1]+3.25)*(5.0*pPhys[0]+5.0*pPhys[1]+3.25)))*cos(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
		{
		if (pPhys[0] >= 0.0)
			return -sigma_*(exp(-(-5.0*pPhys[0]+5.0*pPhys[1]+1.75)*(-5.0*pPhys[0]+5.0*pPhys[1]+1.75))+exp(-(2.5*pPhys[0]+2.5*pPhys[1]-3.25)*(2.5*pPhys[0]+2.5*pPhys[1]-3.25)))*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0)
			return -sigma_*(exp(-(2.5*pPhys[0]-2.5*pPhys[1]+1.75)*(2.5*pPhys[0]-2.5*pPhys[1]+1.75))+exp(-(5.0*pPhys[0]+5.0*pPhys[1]+3.25)*(5.0*pPhys[0]+5.0*pPhys[1]+3.25)))*cos(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for P" << endl;
}

double
EBbucketd::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-N2_*pPhys[1]);
}

double
EBbucketd::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-N2_*exp(-N2_*pPhys[1]);
}

double
EBbucketd::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	N2_;
}

EBbuckete::EBbuckete():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 0.0;
	N2_ = 2.0;
	sigma_ = sqrt(0.5*N2_);
}

double
EBbuckete::getU(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0) // region I
			return (cosh(pPhys[0]-pPhys[1]+1.0)-2.0*cosh(1.0)+cosh(pPhys[0]+pPhys[1]-1.0))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0) // region IV
			return (-2.0*cosh(2.0*pPhys[0]-2.0*pPhys[1]-1.0)+4.0*cosh(1.0)-2.0*cosh(2.0*pPhys[0]+2.0*pPhys[1]+1.0))*sin(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
		{
		if (pPhys[0] >= 0.0) // region II
			return (-2.0*cosh(-1.0+2.0*pPhys[0]-2.0*pPhys[1])+cosh(1.0)+cosh(pPhys[0]+pPhys[1]-1.0))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0) // region III
			return (cosh(pPhys[0]-pPhys[1]+1.0)+cosh(1.0)-2.0*cosh(2.0*pPhys[0]+2.0*pPhys[1]+1.0))*sin(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for U" << endl;
}

double
EBbuckete::getV(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0)
			return (cosh(pPhys[0]-pPhys[1]+1.0)-cosh(pPhys[0]+pPhys[1]-1.0))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0)
			return (-2.0*cosh(2.0*pPhys[0]-2.0*pPhys[1]-1.0)+2.0*cosh(2.0*pPhys[0]+2.0*pPhys[1]+1.0))*sin(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
		{
		if (pPhys[0] >= 0.0)
			return (-2.0*cosh(-1.0+2.0*pPhys[0]-2.0*pPhys[1])+3.0*cosh(1.0)-cosh(pPhys[0]+pPhys[1]-1.0))*sin(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0)
			return (cosh(pPhys[0]-pPhys[1]+1.0)-3.0*cosh(1.0)+2.0*cosh(2.0*pPhys[0]+2.0*pPhys[1]+1.0))*sin(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for V (W)" << endl;
}

double
EBbuckete::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
EBbuckete::getR(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0)
			return -N2_/sigma_*(cosh(pPhys[0]-pPhys[1]+1.0)-cosh(pPhys[0]+pPhys[1]-1.0))*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0)
			return -N2_/sigma_*(-2.0*cosh(2.0*pPhys[0]-2.0*pPhys[1]-1.0)+2.0*cosh(2.0*pPhys[0]+2.0*pPhys[1]+1.0))*cos(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
	{
		if (pPhys[0] >= 0.0)
			return -N2_/sigma_*(-2.0*cosh(-1.0+2.0*pPhys[0]-2.0*pPhys[1])+3.0*cosh(1.0)-cosh(pPhys[0]+pPhys[1]-1.0))*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0)
			return -N2_/sigma_*(cosh(pPhys[0]-pPhys[1]+1.0)-3.0*cosh(1.0)+2.0*cosh(2.0*pPhys[0]+2.0*pPhys[1]+1.0))*cos(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for R" << endl;
}

double
EBbuckete::getP(const PointPhysicalT& pPhys, double t)const
{
	if (abs(pPhys[1]) >= abs(pPhys[0]))
	{
		if (pPhys[1] >= 0.0)
			return -sigma_*(sinh(pPhys[0]-pPhys[1]+1.0)-cosh(1.0)*(pPhys[0]-pPhys[1]+1.0)+sinh(pPhys[0]+pPhys[1]-1.0)-cosh(1.0)*(pPhys[0]+pPhys[1]-1.0))*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[1] < 0.0)
			return -sigma_*(sinh(1.0-2.0*pPhys[0]+2.0*pPhys[1])-cosh(1.0)*(1.0-2.0*pPhys[0]+2.0*pPhys[1])+sinh(-2.0*pPhys[0]-2.0*pPhys[1]-1.0)-cosh(1.0)*(-2.0*pPhys[0]-2.0*pPhys[1]-1.0))*cos(sigma_*t+Pi_/4.0);
	}
	else if (abs(pPhys[1]) < abs(pPhys[0]))
		{
		if (pPhys[0] >= 0.0)
			return -sigma_*(sinh(-2.0*pPhys[0]+2.0*pPhys[1]+1.0)-cosh(1.0)*(-2.0*pPhys[0]+2.0*pPhys[1]+1.0)+sinh(pPhys[0]+pPhys[1]-1.0)-cosh(1.0)*(pPhys[0]+pPhys[1]-1.0))*cos(sigma_*t+Pi_/4.0);
		else if (pPhys[0] < 0.0)
			return -sigma_*(sinh(pPhys[0]-pPhys[1]+1.0)-cosh(1.0)*(pPhys[0]-pPhys[1]+1.0)+sinh(-2.0*pPhys[0]-2.0*pPhys[1]-1.0)-cosh(1.0)*(-2.0*pPhys[0]-2.0*pPhys[1]-1.0))*cos(sigma_*t+Pi_/4.0);
	}
	else
		cout << "Error in Initial Condition for P" << endl;
}

double
EBbuckete::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-N2_*pPhys[1]);
}

double
EBbuckete::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-N2_*exp(-N2_*pPhys[1]);
}

double
EBbuckete::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	N2_;
}

EBTrapezoidWAtau3Over2::EBTrapezoidWAtau3Over2():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 0.0;
	N2_ = 2.0;
	sigma_ = sqrt(0.5*N2_);
	mmax_ = 1500;
	nmax_ = 10;
}

double
EBTrapezoidWAtau3Over2::getU(const PointPhysicalT& pPhys, double t)const
{
	double sumn;
	double bn;
	double a2;
	double ret;
	ret = 0.0;
	for (double m = 1;m<mmax_+1;++m) 
	{
		sumn = 0.0;
		for (double n = 1;n<nmax_+1;++n) 
		{
			bn = 1.5*pow(5.0,n);
			sumn = sumn + sin(m*Pi_/(2.0*bn))*(1.0/(bn*bn-m*m)-1.0/(25.0*bn*bn-m*m));
		}
		a2 = 2.0*m*pow(-1.0,m+1.0)/Pi_*sumn;
		ret = ret + a2*4.0/3.0*m*Pi_*sin(4.0/3.0*m*Pi_*(pPhys[0]+1.0))*cos(4.0/3.0*m*Pi_*pPhys[1]);
	}
	return -ret*sin(sigma_*t+Pi_/4.0);
}
double
EBTrapezoidWAtau3Over2::getV(const PointPhysicalT& pPhys, double t)const
{
	double sumn;
	double bn;
	double a2;
	double ret;
	ret = 0.0;
	for (double m = 1;m<mmax_+1;++m) 
	{
		sumn = 0.0;
		for (double n = 1;n<nmax_+1;++n) 
		{
			bn = 1.5*pow(5.0,n);
			sumn = sumn + sin(m*Pi_/(2.0*bn))*(1.0/(bn*bn-m*m)-1.0/(25.0*bn*bn-m*m));
		}
		a2 = 2.0*m*pow(-1.0,m+1.0)/Pi_*sumn;
		ret = ret + a2*4.0/3.0*m*Pi_*cos(4.0/3.0*m*Pi_*(pPhys[0]+1.0))*sin(4.0/3.0*m*Pi_*pPhys[1]);
	}
	return ret*sin(sigma_*t+Pi_/4.0);
}
double
EBTrapezoidWAtau3Over2::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}		
double
EBTrapezoidWAtau3Over2::getR(const PointPhysicalT& pPhys, double t)const
{
	double sumn;
	double bn;
	double a2;
	double ret;
	ret = 0.0;
	for (double m = 1;m<mmax_+1;++m) 
	{
		sumn = 0.0;
		for (double n = 1;n<nmax_+1;++n) 
		{
			bn = 1.5*pow(5.0,n);
			sumn = sumn + sin(m*Pi_/(2.0*bn))*(1.0/(bn*bn-m*m)-1.0/(25.0*bn*bn-m*m));
		}
		a2 = 2.0*m*pow(-1.0,m+1.0)/Pi_*sumn;
		ret = ret + a2*4.0/3.0*m*Pi_*cos(4.0/3.0*m*Pi_*(pPhys[0]+1.0))*sin(4.0/3.0*m*Pi_*pPhys[1]);
	}
	return -N2_/sigma_*ret*cos(sigma_*t+Pi_/4.0);
}
double
EBTrapezoidWAtau3Over2::getP(const PointPhysicalT& pPhys, double t)const
{
	double sumn;
	double bn;
	double a2;
	double ret;
	ret = 0.0;
	for (double m = 1;m<mmax_+1;++m) 
	{
		sumn = 0.0;
		for (double n = 1;n<nmax_+1;++n) 
		{
			bn = 1.5*pow(5.0,n);
			sumn = sumn + sin(m*Pi_/(2.0*bn))*(1.0/(bn*bn-m*m)-1.0/(25.0*bn*bn-m*m));
		}
		a2 = 2.0*m*pow(-1.0,m+1.0)/Pi_*sumn;
		ret = ret + a2*cos(4.0/3.0*m*Pi_*(pPhys[0]+1.0))*cos(4.0/3.0*m*Pi_*pPhys[1]);
	}
	return -sigma_*ret*cos(sigma_*t+Pi_/4.0);
}
double
EBTrapezoidWAtau3Over2::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-N2_*pPhys[1]);
}
double
EBTrapezoidWAtau3Over2::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-N2_*exp(-N2_*pPhys[1]);
}
double
EBTrapezoidWAtau3Over2::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	N2_;
}

EB2DSingleMode::EB2DSingleMode():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 2.0*Pi_;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = 2.0;
	sigma_ = sqrt(((k_*k_+l_*l_)*N2_)/(k_*k_+l_*l_+m_*m_));
}

double
EB2DSingleMode::getU(const PointPhysicalT& pPhys, double t)const
{
	return -m_*k_/(k_*k_+l_*l_)*cos(m_*pPhys[1])*sin(k_*pPhys[0])*sin(sigma_*t+0.1);
}

double
EB2DSingleMode::getV(const PointPhysicalT& pPhys, double t)const
{
	return sin(m_*pPhys[1])*cos(k_*pPhys[0])*sin(sigma_*t+0.1);
}

double
EB2DSingleMode::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
EB2DSingleMode::getR(const PointPhysicalT& pPhys, double t)const
{
	return -N2_/sigma_*sin(m_*pPhys[1])*cos(k_*pPhys[0])*cos(sigma_*t+0.1);
}
		
double
EB2DSingleMode::getP(const PointPhysicalT& pPhys, double t)const
{
	return	-m_*sigma_/(k_*k_+l_*l_)*cos(m_*pPhys[1])*cos(k_*pPhys[0])*cos(sigma_*t+0.1);
}

double
EB2DSingleMode::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-N2_*pPhys[1]);
}

double
EB2DSingleMode::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-N2_*exp(-N2_*pPhys[1]);
}

double
EB2DSingleMode::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	N2_;
}

RE0::RE0():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = -1.0;
	beta_	= N2_+1.0;
	sigma_ = sqrt(m_*m_);
	//sqrt(beta_*beta_+m_*m_);
}

double
RE0::getU(const PointPhysicalT& pPhys, double t)const
{
	return sin(m_*pPhys[0])*sin(sigma_*t+0.1);
	//exp(-0.5*beta_*pPhys[0])*sin(m_*pPhys[0])*sin(sigma_*t+0.1);
}

double
RE0::getV(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}

double
RE0::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
RE0::getR(const PointPhysicalT& pPhys, double t)const
{
	//cout << "m_ = " << m_ << ", sigma_ = " << sigma_ << endl;
	return m_/sigma_*cos(m_*pPhys[0])*cos(sigma_*t+0.1);
	
	//exp(-0.5*beta_*pPhys[0])*(0.5*beta_/sigma_*sin(m_*pPhys[0])+m_/sigma_*cos(m_*pPhys[0]))*cos(sigma_*t+0.1);
}
		
double
RE0::getP(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;
}

double
RE0::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	1.0;//exp(-beta_*pPhys[0]);
}

double
RE0::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

double
RE0::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

RES0::RES0():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = 2.0;
	beta_	= N2_+1.0;
	sigma_ = sqrt(0.25*beta_*beta_+m_*m_);
}

double
RES0::getU(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*sin(m_*pPhys[0])*sin(sigma_*t+0.1);
}

double
RES0::getV(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}

double
RES0::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
RES0::getR(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*(0.5*beta_/sigma_*sin(m_*pPhys[0])+m_/sigma_*cos(m_*pPhys[0]))*cos(sigma_*t+0.1);
}
		
double
RES0::getP(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;
}

double
RES0::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-beta_*pPhys[0]);
}

double
RES0::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-beta_*exp(-beta_*pPhys[0]);
}

double
RES0::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

RES1::RES1():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = 2.0;
	beta_	= N2_+1.0;
	sigma_ = sqrt(beta_/N2_*(0.25*beta_*beta_+m_*m_));
}
double
RES1::getU(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*sin(m_*pPhys[0])*sin(sigma_*t+0.1);
}
double
RES1::getV(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
double
RES1::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
double
RES1::getR(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*(0.5*beta_/sigma_*sin(m_*pPhys[0])+m_/sigma_*cos(m_*pPhys[0]))*cos(sigma_*t+0.1);
}
double
RES1::getP(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;
}
double
RES1::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-beta_*pPhys[0]);
}
double
RES1::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-beta_*exp(-beta_*pPhys[0]);
}

double
RES1::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

RES2::RES2():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = 2.0;
	beta_	= N2_+1.0;
	sigma_ = sqrt(beta_/N2_*(0.25*(N2_-1.0)*(N2_-1.0)+m_*m_));
}
double
RES2::getU(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*sin(m_*pPhys[0])*sin(sigma_*t+0.1);
}
double
RES2::getV(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
double
RES2::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
double
RES2::getR(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*(0.5*(N2_-1.0)/sigma_*sin(m_*pPhys[0])+m_/sigma_*cos(m_*pPhys[0]))*cos(sigma_*t+0.1);
}
double
RES2::getP(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;
}
double
RES2::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-beta_*pPhys[0]);
}
double
RES2::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-beta_*exp(-beta_*pPhys[0]);
}
double
RES2::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

RESDO0::RESDO0():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = 1.0;
	beta_	= N2_+1.0;
	sigma_ = sqrt((beta_*beta_+4.0*m_*m_)/(4.0*N2_));
}

double
RESDO0::getU(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*sin(m_*pPhys[0])*sin(sigma_*t+0.1);
}

double
RESDO0::getV(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}

double
RESDO0::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
RESDO0::getR(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*(-0.5*beta_/sigma_*sin(m_*pPhys[0])+m_/sigma_*cos(m_*pPhys[0]))*cos(sigma_*t+0.1);
}
		
double
RESDO0::getP(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;
}

double
RESDO0::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-beta_*pPhys[0]);
}

double
RESDO0::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-beta_*exp(-beta_*pPhys[0]);
} 
double
RESDO0::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

REMTC0::REMTC0():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = 13.645697703860776;//-1.0;//
	beta_	= N2_+1.0;
	sigma_ = sqrt(((N2_+1.0)*(N2_-1.0)-4.0*m_*m_)/(2.0*N2_));
}

double
REMTC0::getU(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*sin(m_*pPhys[0])*sin(sigma_*t+0.1);
}

double
REMTC0::getV(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}

double
REMTC0::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
REMTC0::getR(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*(-0.5*beta_/sigma_*sin(m_*pPhys[0])+m_/sigma_*cos(m_*pPhys[0]))*cos(sigma_*t+0.1);
}
		
double
REMTC0::getP(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-0.5*beta_*pPhys[0])*(0.5*(N2_-1.0)/sigma_*sin(m_*pPhys[0])+m_/sigma_*cos(m_*pPhys[0]))*cos(sigma_*t+0.1);
}

double
REMTC0::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-beta_*pPhys[0]);
}

double
REMTC0::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-beta_*exp(-beta_*pPhys[0]);
}
double
REMTC0::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
} 

CS1D::CS1D():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 0.0;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = 2.0;//-1.0;//
	beta_	= N2_+1.0;
	sigma_ = 0.5*sqrt(beta_*beta_+4.0*m_*m_);
}

double
CS1D::getU(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*sin(m_*pPhys[0])*sin(sigma_*t+0.1);
}

double
CS1D::getV(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}

double
CS1D::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
CS1D::getR(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[0])*(-0.5*beta_/sigma_*sin(m_*pPhys[0])+m_/sigma_*cos(m_*pPhys[0]))*cos(sigma_*t+0.1);
}
		
double
CS1D::getP(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-0.5*beta_*pPhys[0])*(0.5*(N2_-1.0)/sigma_*sin(m_*pPhys[0])+m_/sigma_*cos(m_*pPhys[0]))*cos(sigma_*t+0.1);
}

double
CS1D::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-beta_*pPhys[0]);
}

double
CS1D::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-beta_*exp(-beta_*pPhys[0]);
}
double
CS1D::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

RES02D::RES02D():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 2.0*Pi_;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = 2.0;
	beta_	= N2_+1.0;
	sigma_ = sqrt(0.25*beta_*beta_+m_*m_+k_*k_);
}

double
RES02D::getU(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[1])*k_/(sigma_*sigma_-k_*k_)*(0.5*beta_*sin(m_*pPhys[1])+m_*cos(m_*pPhys[1]))*sin(k_*pPhys[0])*sin(sigma_*t+0.1);
}

double
RES02D::getV(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[1])*sin(m_*pPhys[1])*cos(k_*pPhys[0])*sin(sigma_*t+0.1);
}

double
RES02D::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
RES02D::getR(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
RES02D::getP(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-0.5*beta_*pPhys[1])*sigma_/(sigma_*sigma_-k_*k_)*(0.5*beta_*sin(m_*pPhys[1])+m_*cos(m_*pPhys[1]))*cos(k_*pPhys[0])*cos(sigma_*t+0.1);
}

double
RES02D::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-beta_*pPhys[1]);
}

double
RES02D::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-beta_*exp(-beta_*pPhys[1]);
}
double
RES02D::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

CS2D::CS2D():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 2.0*Pi_;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = 2.0;
	beta_	= N2_+1.0;
	sigma_ = sqrt(0.5*(k_*k_+0.25*beta_*beta_+m_*m_+sqrt((k_*k_+0.25*beta_*beta_+m_*m_)*(k_*k_+0.25*beta_*beta_+m_*m_)-4.0*k_*k_*N2_)));
}

double
CS2D::getU(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[1])*k_/(sigma_*sigma_-k_*k_)*(0.5*(N2_-1.0)*sin(m_*pPhys[1])+m_*cos(m_*pPhys[1]))*sin(k_*pPhys[0])*sin(sigma_*t+0.1);
}

double
CS2D::getV(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[1])*sin(m_*pPhys[1])*cos(k_*pPhys[0])*sin(sigma_*t+0.1);
}

double
CS2D::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
CS2D::getR(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[1])*sigma_/(sigma_*sigma_-k_*k_)*((k_*k_*N2_/(sigma_*sigma_)-0.5*beta_)*sin(m_*pPhys[1])+m_*cos(m_*pPhys[1]))*cos(k_*pPhys[0])*cos(sigma_*t+0.1);
}
		
double
CS2D::getP(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-0.5*beta_*pPhys[1])*sigma_/(sigma_*sigma_-k_*k_)*(0.5*(N2_-1.0)*sin(m_*pPhys[1])+m_*cos(m_*pPhys[1]))*cos(k_*pPhys[0])*cos(sigma_*t+0.1);
}

double
CS2D::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-beta_*pPhys[1]);
}

double
CS2D::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-beta_*exp(-beta_*pPhys[1]);
}
double
CS2D::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

Lambd3DMixed::Lambd3DMixed():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 2.0*Pi_;
	l_ = 2.0*Pi_;
	m_ = 0.0;
	N2_ = 2.0;
	beta_	= 3.0;
	sigma_ = sqrt(k_*k_+l_*l_);
}

double
Lambd3DMixed::getU(const PointPhysicalT& pPhys, double t)const
{
	return -k_/sigma_*exp(-pPhys[2])*cos(k_*pPhys[0]+l_*pPhys[1]+sigma_*t);
}	 

double
Lambd3DMixed::getV(const PointPhysicalT& pPhys, double t)const
{
	return -l_/sigma_*exp(-pPhys[2])*cos(k_*pPhys[0]+l_*pPhys[1]+sigma_*t);
}	

double
Lambd3DMixed::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}	
		
double
Lambd3DMixed::getR(const PointPhysicalT& pPhys, double t)const
{
	return exp(-pPhys[2])*cos(k_*pPhys[0]+l_*pPhys[1]+sigma_*t);
}	
		
double
Lambd3DMixed::getP(const PointPhysicalT& pPhys, double t)const
{
	return exp(-pPhys[2])*cos(k_*pPhys[0]+l_*pPhys[1]+sigma_*t);
}
	
double
Lambd3DMixed::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-beta_*pPhys[2]);;//1.0;
}

double
Lambd3DMixed::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-beta_*exp(-beta_*pPhys[2]);//0.0;
}
double
Lambd3DMixed::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

CS3D::CS3D():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 2.0*Pi_;
	l_ = 2.0*Pi_;
	m_ = 2.0*Pi_;
	N2_ = 2.0;
	beta_	= N2_+1.0;
	sigma_ = sqrt(0.5*(k_*k_+l_*l_+0.25*beta_*beta_+m_*m_+sqrt((k_*k_+l_*l_+0.25*beta_*beta_+m_*m_)*(k_*k_+l_*l_+0.25*beta_*beta_+m_*m_)-4.0*(k_*k_+l_*l_)*N2_)));
}

double
CS3D::getU(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[2])*k_/(sigma_*sigma_-k_*k_-l_*l_)*(0.5*(N2_-1.0)*sin(m_*pPhys[2])+m_*cos(m_*pPhys[2]))*sin(k_*pPhys[0])*cos(l_*pPhys[1])*sin(sigma_*t+0.1);
}	 

double
CS3D::getV(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[2])*l_/(sigma_*sigma_-k_*k_-l_*l_)*(0.5*(N2_-1.0)*sin(m_*pPhys[2])+m_*cos(m_*pPhys[2]))*cos(k_*pPhys[0])*sin(l_*pPhys[1])*sin(sigma_*t+0.1);
}	

double
CS3D::getW(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[2])*sin(m_*pPhys[2])*cos(k_*pPhys[0])*cos(l_*pPhys[1])*sin(sigma_*t+0.1);
}	
		
double
CS3D::getR(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[2])*sigma_/(sigma_*sigma_-k_*k_-l_*l_)*(((k_*k_+l_*l_)*N2_/(sigma_*sigma_)-0.5*beta_)*sin(m_*pPhys[2])+m_*cos(m_*pPhys[2]))*cos(k_*pPhys[0])*cos(l_*pPhys[1])*cos(sigma_*t+0.1);
}	
		
double
CS3D::getP(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*beta_*pPhys[2])*sigma_/(sigma_*sigma_-k_*k_-l_*l_)*(0.5*(N2_-1.0)*sin(m_*pPhys[2])+m_*cos(m_*pPhys[2]))*cos(k_*pPhys[0])*cos(l_*pPhys[1])*cos(sigma_*t+0.1);
}
	
double
CS3D::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-beta_*pPhys[2]);
}

double
CS3D::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-beta_*exp(-beta_*pPhys[2]);
}
double
CS3D::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//-beta_*exp(-beta_*pPhys[0]);
}

ICS2D::ICS2D():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	k_ = 2.0*Pi_;
	l_ = 0.0;
	m_ = 2.0*Pi_;
	N2_ = 2.0;
	sigma_ = sqrt((4.0*k_*k_*N2_)/(N2_*N2_+4.0*k_*k_+4.0*m_*m_));
}

double
ICS2D::getU(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*N2_*pPhys[1])*(-0.5*N2_/k_*sin(m_*pPhys[1])-m_/k_*cos(m_*pPhys[1]))*sin(k_*pPhys[0])*sin(sigma_*t+0.1);
}

double
ICS2D::getV(const PointPhysicalT& pPhys, double t)const
{
	return exp(-0.5*N2_*pPhys[1])*sin(m_*pPhys[1])*cos(k_*pPhys[0])*sin(sigma_*t+0.1);
}

double
ICS2D::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
ICS2D::getR(const PointPhysicalT& pPhys, double t)const
{
	return -exp(-0.5*N2_*pPhys[1])*N2_/sigma_*sin(m_*pPhys[1])*cos(k_*pPhys[0])*cos(sigma_*t+0.1);
}
		
double
ICS2D::getP(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-0.5*N2_*pPhys[1])*(-0.5*N2_*sigma_/k_/k_*sin(m_*pPhys[1])-m_*sigma_/k_/k_*cos(m_*pPhys[1]))*cos(k_*pPhys[0])*cos(sigma_*t+0.1);
}

double
ICS2D::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-N2_*pPhys[1]);
}

double
ICS2D::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-N2_*exp(-N2_*pPhys[1]);
}
double
ICS2D::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	2.0;
}

IGW2D::IGW2D():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	N2_ = 2.0;
	sigma_ = sqrt(N2_/2.0);
	nmax_ = 10;
	
}

double
IGW2D::getU(const PointPhysicalT& pPhys, double t)const
{
	double ret;
	ret = 0.0;
	for (double n = 1;n<nmax_+1;++n)
	{
		ret = ret + cos(n*Pi_*pPhys[1])*cos(n*Pi_*pPhys[0]-sigma_*t);
	}
	return ret;
}

double
IGW2D::getV(const PointPhysicalT& pPhys, double t)const
{
	double ret;
	ret = 0.0;
	for (double n = 1;n<nmax_+1;++n)
	{
		ret = ret + sin(n*Pi_*pPhys[1])*sin(n*Pi_*pPhys[0]-sigma_*t);
	}
	return ret;
}

double
IGW2D::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
IGW2D::getR(const PointPhysicalT& pPhys, double t)const
{
	double ret;
	ret = 0.0;
	for (double n = 1;n<nmax_+1;++n)
	{
		ret = ret + N2_/sigma_*sin(n*Pi_*pPhys[1])*cos(n*Pi_*pPhys[0]-sigma_*t);
	}
	return ret;
}
		
double
IGW2D::getP(const PointPhysicalT& pPhys, double t)const
{
	double ret;
	ret = 0.0;
	for (double n = 1;n<nmax_+1;++n)
	{
		ret = ret + sigma_/n/Pi_*cos(n*Pi_*pPhys[1])*cos(n*Pi_*pPhys[0]-sigma_*t);
	}
	return ret;
}

double
IGW2D::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	exp(-N2_*pPhys[1]);
}

double
IGW2D::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	-N2_*exp(-N2_*pPhys[1]);
}
double
IGW2D::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	N2_;
}

IGWN2LinearAiry::IGWN2LinearAiry():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	N_0_ = 1.0;
	lambda_ = 0.5;
	sigma_ = sqrt(2.0/3.0);
	k1_ = 7.82220337369810;
	
	k1hat_ = sqrt(k1_*k1_*lambda_/sigma_/sigma_);
	A0_	= (N_0_*N_0_-sigma_*sigma_)/lambda_;
	
}

double
IGWN2LinearAiry::getU(const PointPhysicalT& pPhys, double t)const
{
	double Ai, Bi, Aip, Bip;
	Maths::Special::Airy::airy(-pow(k1hat_,2.0/3.0)*(pPhys[1]-1.0+A0_), Ai, Aip, Bi, Bip);
	double Aic, Bic, Aipc, Bipc;
	Maths::Special::Airy::airy(-pow(k1hat_,2.0/3.0)*(A0_), Aic, Aipc, Bic, Bipc);	
	return pow(k1hat_,2.0/3.0)/k1_*(-Aip+Aic/Bic*Bip)*cos(k1_*pPhys[0]-sigma_*t);
}

double
IGWN2LinearAiry::getV(const PointPhysicalT& pPhys, double t)const
{
	double Ai, Bi, Aip, Bip;
	Maths::Special::Airy::airy(-pow(k1hat_,2.0/3.0)*(pPhys[1]-1.0+A0_), Ai, Aip, Bi, Bip);
	double Aic, Bic, Aipc, Bipc;
	Maths::Special::Airy::airy(-pow(k1hat_,2.0/3.0)*(A0_), Aic, Aipc, Bic, Bipc);
	return (Ai-Aic/Bic*Bi)*sin(k1_*pPhys[0]-sigma_*t);
}

double
IGWN2LinearAiry::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
IGWN2LinearAiry::getR(const PointPhysicalT& pPhys, double t)const
{
	double Ai, Bi, Aip, Bip;
	Maths::Special::Airy::airy(-pow(k1hat_,2.0/3.0)*(pPhys[1]-1.0+A0_), Ai, Aip, Bi, Bip);
	double Aic, Bic, Aipc, Bipc;
	Maths::Special::Airy::airy(-pow(k1hat_,2.0/3.0)*(A0_), Aic, Aipc, Bic, Bipc);
	return getN2(pPhys, t)/sigma_*(Ai-Aic/Bic*Bi)*cos(k1_*pPhys[0]-sigma_*t);
	
	//-getN2(pPhys, t)/sigma_*(Ai-Aic/Bic*Bi)*cos(k1_*pPhys[0]-sigma_*t);
	
	//Ai-Aic/Bic*Bi*cos(k1_*pPhys[0]-sigma_*t);
}
		
double
IGWN2LinearAiry::getP(const PointPhysicalT& pPhys, double t)const
{
	double Ai, Bi, Aip, Bip;
	Maths::Special::Airy::airy(-pow(k1hat_,2.0/3.0)*(pPhys[1]-1.0+A0_), Ai, Aip, Bi, Bip);
	double Aic, Bic, Aipc, Bipc;
	Maths::Special::Airy::airy(-pow(k1hat_,2.0/3.0)*(A0_), Aic, Aipc, Bic, Bipc);	
	return pow(k1hat_,2.0/3.0)/k1_*sigma_/k1_*(-Aip+Aic/Bic*Bip)*cos(k1_*pPhys[0]-sigma_*t);
}

double
IGWN2LinearAiry::getrho0(const PointPhysicalT& pPhys, double t)const
{
	//cout << exp(-(N_0_*N_0_-lambda_)*pPhys[1]-0.5*lambda_*pPhys[1]*pPhys[1]) << endl;
	return	1.0; //exp(-(N_0_*N_0_-lambda_)*pPhys[1]-0.5*lambda_*pPhys[1]*pPhys[1]);
}

double
IGWN2LinearAiry::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;//(-(N_0_*N_0_-lambda_)-lambda_*pPhys[1])*exp(-(N_0_*N_0_-lambda_)*pPhys[1]-0.5*lambda_*pPhys[1]*pPhys[1]);
}

double
IGWN2LinearAiry::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	N_0_*N_0_+lambda_*(pPhys[1]-1.0);
}

IGW2Dprolonged::IGW2Dprolonged():
ExactSolutionBase()
{
	Pi_ = 3.141592653589793;
	N2_ = 2.0;
	sigma_ = sqrt(N2_/2.0);
	nmax_ = 10;
	
}

double
IGW2Dprolonged::getU(const PointPhysicalT& pPhys, double t)const
{
	double ret;
	ret = 0.0;
	//if (pPhys[0] >= 2.0 && pPhys[0] <= 4.0) 
	{
	for (double n = 1;n<nmax_+1;++n)
	{
		ret = ret + cos(n*Pi_*pPhys[1])*cos(n*Pi_*pPhys[0]-sigma_*t);
	}
	}
	return ret;
}

double
IGW2Dprolonged::getV(const PointPhysicalT& pPhys, double t)const
{
	double ret;
	ret = 0.0;
	//if (pPhys[0] >= 2.0 && pPhys[0] <= 4.0)
	{
	for (double n = 1;n<nmax_+1;++n)
	{
		ret = ret + sin(n*Pi_*pPhys[1])*sin(n*Pi_*pPhys[0]-sigma_*t);
	}
	}
	return ret;
}

double
IGW2Dprolonged::getW(const PointPhysicalT& pPhys, double t)const
{
	return 0.0;
}
		
double
IGW2Dprolonged::getR(const PointPhysicalT& pPhys, double t)const
{
	double ret;
	ret = 0.0;
	//if (pPhys[0] >= 2.0 && pPhys[0] <= 4.0) 
	{
	for (double n = 1;n<nmax_+1;++n)
	{
		ret = ret + N2_/sigma_*sin(n*Pi_*pPhys[1])*cos(n*Pi_*pPhys[0]-sigma_*t);
	}
	}
	return ret;
}
		
double
IGW2Dprolonged::getP(const PointPhysicalT& pPhys, double t)const
{
	double ret;
	ret = 0.0;
	//if (pPhys[0] >= 2.0 && pPhys[0] <= 4.0) 
	{
	for (double n = 1;n<nmax_+1;++n)
	{
		ret = ret + sigma_/n/Pi_*cos(n*Pi_*pPhys[1])*cos(n*Pi_*pPhys[0]-sigma_*t);
	}
	}
	return ret;
}

double
IGW2Dprolonged::getrho0(const PointPhysicalT& pPhys, double t)const
{
	return	1.0;
}

double
IGW2Dprolonged::getrho0Deriv(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;
}

double
IGW2Dprolonged::getN2(const PointPhysicalT& pPhys, double t)const
{
	return	0.0;
}
/*
Compressible3DOneThirdPeriodic::Compressible3DOneThirdPeriodic(double lx, double ly, double lz):
    ExactSolutionBase(),
    iComplex_(0,-1),
    lx_(lx),
    ly_(ly),
    lz_(lz)
{
	k_ 			= 1*2*Pi/lx_;
	l_ 			= 1*2*Pi/ly_;
	sigma_ 		= 1/(sqrt(k_*k_+l_*l_+1));

	cout<< "k="<<k_<<", l="<<l_<<", sigma="<<sigma_<<endl;
}

double
Compressible3DOneThirdPeriodic::getU(const PointPhysicalT& pPhys, double t)const
{
// 	ComplexNumber u = 0;
// 	double xx = pPhys[0];
// 	double yy = pPhys[1];
// 	double zz = pPhys[2];
// 	
// 	u = iComplex_*k_*(sigma_*sigma_*sigma_)/(1-sigma_*sigma_)*(1+(l_/(sigma_*k_))*(l_/(sigma_*k_)))*sin(k_*xx);
// 	u = u*(cos(l_*yy)-iComplex_*sin(l_*yy))*cos(Pi*zz/lz_)*(-cos(sigma_*t)+iComplex_*sin(sigma_*t));
// 	return real(u);

	double u = 0;
	double xx = pPhys[0];
	double yy = pPhys[1];
	double zz = pPhys[2];
	
	u = k_*(sigma_*sigma_*sigma_)/(1-sigma_*sigma_)*(1+(l_*l_)/(sigma_*k_*sigma_*k_))*sin(k_*xx)*(-sin(l_*yy-sigma_*t))*cos(2*Pi*zz/lz_);
	
	return u;
}

double
Compressible3DOneThirdPeriodic::getV(const PointPhysicalT& pPhys, double t)const
{
// 	ComplexNumber v = 0;
// 	double xx = pPhys[0];
// 	double yy = pPhys[1];
// 	double zz = pPhys[2];
// 	
// 	v =(-l_*sigma_*cos(k_*xx)+1/k_*sin(k_*xx))*(cos(l_*yy)-iComplex_*sin(l_*yy))*cos(Pi*zz/lz_)*(-cos(sigma_*t)+iComplex_*sin(sigma_*t));
// 	return real(v);


	double v = 0;
	double xx = pPhys[0];
	double yy = pPhys[1];
	double zz = pPhys[2];
	
	v =(-l_*sigma_*cos(k_*xx)+1/k_*sin(k_*xx))*(cos(l_*yy-sigma_*t))*cos(2*Pi*zz/lz_);
	return v;
}
		
double
Compressible3DOneThirdPeriodic::getW(const PointPhysicalT& pPhys, double t)const
{
// 	ComplexNumber w = 0;
// 	double xx = pPhys[0];
// 	double yy = pPhys[1];
// 	double zz = pPhys[2];
// 	
// 	w = (-iComplex_*sigma_)*(cos(k_*xx)+l_/(k_*sigma_)*sin(k_*xx))*(cos(l_*yy)-iComplex_*sin(l_*yy))*sin(Pi*zz/lz_)*(-cos(sigma_*t)+iComplex_*sin(sigma_*t));
// 	return real(w);

	double w = 0;
	double xx = pPhys[0];
	double yy = pPhys[1];
	double zz = pPhys[2];
	
	w = (sigma_)*(cos(k_*xx)+l_/(k_*sigma_)*sin(k_*xx))*(sin(l_*yy-sigma_*t))*sin(2*Pi*zz/lz_);
	return w;
}
		
double
Compressible3DOneThirdPeriodic::getP(const PointPhysicalT& pPhys, double t)const
{

// 	ComplexNumber p = 0;
// 	double xx = pPhys[0];
// 	double yy = pPhys[1];
// 	double zz = pPhys[2];
// 	
// 	p = (-iComplex_*sigma_)*(-iComplex_*sigma_)*(cos(k_*xx)+l_/(k_*sigma_)*sin(k_*xx))*(cos(l_*yy)-iComplex_*sin(l_*yy))*cos(Pi*zz/lz_)*(-cos(sigma_*t)+iComplex_*sin(sigma_*t));
// 	return real(p);

	double p = 0;
	double xx = pPhys[0];
	double yy = pPhys[1];
	double zz = pPhys[2];
	
	p = -(sigma_)*(sigma_)*(cos(k_*xx)+l_/(k_*sigma_)*sin(k_*xx))*(cos(l_*yy-sigma_*t))*cos(2*Pi*zz/lz_);
	return p;
}	

Compressible3DPeriodic::Compressible3DPeriodic():
ExactSolutionBase()
{
}

double
Compressible3DPeriodic::getU(const PointPhysicalT& pPhys, double t)const
{
//	return 1;
    return std::cos(2*Pi*(t-pPhys[0])); 
}

double
Compressible3DPeriodic::getV(const PointPhysicalT& pPhys, double t)const
{
//    return 2;
	return std::cos(2*Pi*(t-pPhys[1]));
}
		
double
Compressible3DPeriodic::getW(const PointPhysicalT& pPhys, double t)const
{
//    return 3;
	return std::cos(2*Pi*(t-pPhys[2]));
}
		
double
Compressible3DPeriodic::getP(const PointPhysicalT& pPhys, double t)const
{
//    return 4;
      return std::cos(2*Pi*(t-pPhys[0]))+ std::cos(2*Pi*(t-pPhys[1]))+std::cos(2*Pi*(t-pPhys[2]));
}

Compressible3DWalls::Compressible3DWalls():
ExactSolutionBase()
{
}

double
Compressible3DWalls::getU(const PointPhysicalT& pPhys, double t)const
{
	return -cos(2*Pi*(t))*std::sin(2*Pi*(pPhys[0]));
}

double
Compressible3DWalls::getV(const PointPhysicalT& pPhys, double t)const
{
	return -cos(2*Pi*(t))*std::sin(2*Pi*(pPhys[1]));
}
		
double
Compressible3DWalls::getW(const PointPhysicalT& pPhys, double t)const
{
	return -cos(2*Pi*(t))*std::sin(2*Pi*(pPhys[2]));
}
		
double
Compressible3DWalls::getP(const PointPhysicalT& pPhys, double t)const
{
	return	cos(2*Pi*(pPhys[0]))*sin(2*Pi*(t)) + cos(2*Pi*(pPhys[1]))*sin(2*Pi*(t))+cos(2*Pi*(pPhys[2]))*sin(2*Pi*(t));
}




Incompressible3DPeriodic::Incompressible3DPeriodic():
	ExactSolutionBase()
{
}

double
Incompressible3DPeriodic::getU(const PointPhysicalT& pPhys, double t)const
{
	return (sqrt(3)*cos(2*Pi*(pPhys[0] + pPhys[1] + pPhys[2]) + sqrt(3)*t/3)+3*sin(2*Pi*(pPhys[0] + pPhys[1] + pPhys[2])  + sqrt(3)*t/3))/Pi/2;
}

double
Incompressible3DPeriodic::getV(const PointPhysicalT& pPhys, double t)const
{
	return (sqrt(3)*cos(2*Pi*(pPhys[0] + pPhys[1] + pPhys[2]) + sqrt(3)*t/3)-3*sin(2*Pi*(pPhys[0] + pPhys[1] + pPhys[2])  + sqrt(3)*t/3))/Pi/2;
}
		
double
Incompressible3DPeriodic::getW(const PointPhysicalT& pPhys, double t)const
{
	return -(2*sqrt(3)*cos(2*Pi*(pPhys[0] + pPhys[1] + pPhys[2]) + sqrt(3)*t/3))/Pi/2;
}
		
double
Incompressible3DPeriodic::getP(const PointPhysicalT& pPhys, double t)const
{
	return 2*cos(2*Pi*(pPhys[0] + pPhys[1] + pPhys[2]) + sqrt(3)*t/3)/(4*Pi*Pi);
}

InitialVelocityConstructorTaylor::InitialVelocityConstructorTaylor():
    ExactSolutionBase()
{
    
}

InitialVelocityConstructorTaylor::InitialVelocityConstructorTaylor(double lx, double ly, double lz):
	ExactSolutionBase(),
	lx_(lx),
	ly_(ly),
	lz_(lz),
	sigma_(0.6570),//sigma eigenfrequency calculated in Matlab code
	a_(lx_/(Pi*sigma_)),
	size_(50),//N=10 in Matlab code
	muj_(size_+1),
	vCoeff_(size_+1),
	iComplex_(0,-1)
{
		try
		{
			
			
			cout <<"Lx= "<<lx_<<" Ly= "<<ly_<<" Lz= "<<lz_<<endl;
			ifstream bCoeff("bCoeff.txt");// do some checking with existence!!!
			
			std::copy(IstreamIterator(bCoeff), IstreamIterator(), std::back_inserter(b_));
			
//************creation of Symmetric modes*******************
			ComplexNumber a;
			for (int j=1;j<=size_;++j)
			{
				//a = Pi/lx_*sqrt((lx_/Pi)*(lx_/Pi) + j*j - a_*a_);
				muj_[j] = Pi/lx_*std::sqrt(static_cast<ComplexNumber>((lx_/Pi)*(lx_/Pi) + j*j - a_*a_));
				
				std::cerr<<endl;
				std::cerr<<"muj[j]"<<muj_[j]<<endl;
				std::cerr<<"b[j]"<<b_[j-1]<<endl;
				
				if (j%2==0)
				{
					vCoeff_[j] = -iComplex_*sinh(ly_/2)/sinh(muj_[j]*(ly_/2))*a_*Pi/lx_*cos(a_*Pi/2)*b_[j-1];
				}
				else
				{
					vCoeff_[j] = -(cosh(ly_/2)/cosh(muj_[j]*(ly_/2)))*(a_*Pi)/lx_*sin(a_*Pi/2)*b_[j-1];
				}
				std::cerr<<"vCoeff_[j]"<<vCoeff_[j]<<endl;
			}
		}
		catch(const char* c)
		{cout<<c;}
}

double 
InitialVelocityConstructorTaylor::getU(const PointPhysicalT& pPhys, double t)const
{
	ComplexNumber u = 0;
	double xx = pPhys[0];
	double yy = pPhys[1];
	double zz = pPhys[2];

	for(int j=1;j<=size_;++j)
	{
		if (j%2==0)
		{
			u = u + (iComplex_*vCoeff_[j]*(a_/j - j/a_)*pow(-1, (double)j/2)*sin(j*(Pi/lx_*xx-Pi/2))*sinh(muj_[j]*(yy-ly_/2)));

		}
		else
		{
			u = u + (iComplex_*vCoeff_[j]*(a_/j - j/a_)*pow(-1,(double)(j - 1)/2)*cos(j*(Pi/lx_*xx-Pi/2))*cosh(muj_[j]*(yy-ly_/2)));
		}
	}
	u = u*cos(Pi*zz/lz_)*(cos(sigma_*t)-iComplex_*sin(sigma_*t));
	//cout <<"Constructed u" <<real(u)<<endl;
	return real(u);
}

double 
InitialVelocityConstructorTaylor::getV(const PointPhysicalT& pPhys, double t)const
{
	ComplexNumber v = 0;
	ComplexNumber bla=0;
	double xx = pPhys[0];
	double yy = pPhys[1];
	double zz = pPhys[2];
	v         = cosh(yy-ly_/2)*cos(a_*(Pi/lx_*xx-Pi/2)) - iComplex_*sinh(yy-ly_/2)*sin(a_*(Pi/lx_*xx-Pi/2));
	
	for(int j=1;j<=size_;++j)
	{
		if (j%2==0)
		{
			v = v + (vCoeff_[j]*lx_/Pi*pow(-1, (double)j/2)*(iComplex_*muj_[j]/a_*cos(j*(Pi/lx_*xx-Pi/2))*cosh(muj_[j]*(yy-ly_/2)) + sin(j*(Pi/lx_*xx-Pi/2))*sinh(muj_[j]*(yy-ly_/2))/(ComplexNumber)j));

		}
		else
		{  /*if (j==3)
			{
				v = v+ (vCoeff_[j]*lx_/Pi*pow(-1,(j - 1)/2)*(-iComplex_*muj_[j]/a_*sin(j*(Pi/lx_*xx-Pi/2))*sinh(muj_[j]*(yy-ly_/2)) + (cos(j*(Pi/lx_*xx-Pi/2))*cosh(muj_[j]*(yy-ly_/2)))/(ComplexNumber)j));
			}* /
 			v = v + (vCoeff_[j]*lx_/Pi*pow(-1,(double)(j - 1)/2)*(-iComplex_*muj_[j]/a_*sin(j*(Pi/lx_*xx-Pi/2))*sinh(muj_[j]*(yy-ly_/2)) + (cos(j*(Pi/lx_*xx-Pi/2))*cosh(muj_[j]*(yy-ly_/2)))/(ComplexNumber)j));
		}
	}
	v = real(v)*cos(Pi*zz/lz_)*(cos(sigma_*t)-iComplex_*sin(sigma_*t));
	return real(v);
}

double 
InitialVelocityConstructorTaylor::getW(const PointPhysicalT& pPhys, double t)const
{
	ComplexNumber w = 0;
	double xx = pPhys[0];
	double yy = pPhys[1];
	double zz = pPhys[2];
	
	w= -iComplex_*sigma_*(-iComplex_*(1/sigma_)*(sinh(yy-ly_/2)*cos(a_*(Pi/lx_*xx-Pi/2))-iComplex_*cosh(yy-ly_/2)*sin(a_*(Pi/lx_*xx-Pi/2))));
	
	for(int j=1;j<=size_;++j)
	{
		if (j%2==0)
		{
			w = w - iComplex_*sigma_*(-iComplex_)*(1/sigma_)*vCoeff_[j]*lx_/Pi*pow(-1, (double)j/2)*(iComplex_*(muj_[j]/a_)*cos(j*(Pi/lx_*xx-Pi/2))*muj_[j]*sinh(muj_[j]*(yy-ly_/2)) + (sin(j*(Pi/lx_*xx-Pi/2))*muj_[j]*cosh(muj_[j]*(yy-ly_/2)))/(ComplexNumber)j);
			w = w - iComplex_*sigma_*(-iComplex_)*(1/sigma_)*iComplex_*vCoeff_[j]*(a_/j - j/a_)*pow(-1, (double)j/2)*(double)j*(Pi/lx_)*cos(j*(Pi/lx_*xx-Pi/2))*sinh(muj_[j]*(yy-ly_/2));

		}
		else
		{
			w = w - iComplex_*sigma_* (-iComplex_*(1/sigma_)*vCoeff_[j]*lx_/Pi*pow(-1,(double)(j - 1)/2)*(-iComplex_*muj_[j]/a_*sin(j*(Pi/lx_*xx-Pi/2))*muj_[j]*cosh(muj_[j]*(yy-ly_/2)) + (cos(j*(Pi/lx_*xx-Pi/2))*muj_[j]*sinh(muj_[j]*(yy-ly_/2)))/(ComplexNumber)j));
			w = w - iComplex_*sigma_*(-iComplex_*(1/sigma_)*iComplex_*vCoeff_[j]*(a_/j - j/a_)*pow(-1,(double)(j - 1)/2)*(-j*Pi/lx_)*sin(j*(Pi/lx_*xx-Pi/2))*cosh(muj_[j]*(yy-ly_/2)));
		}
	}
	w = w*sin(Pi*zz/lz_)*(cos(sigma_*t)-iComplex_*sin(sigma_*t));
	return real(w);
}

double 
InitialVelocityConstructorTaylor::getP(const PointPhysicalT& pPhys, double t)const
{
	ComplexNumber p = 0;
	double xx = pPhys[0];
	double yy = pPhys[1];
	double zz = pPhys[2];
	
	p= (-iComplex_*sigma_)*(-iComplex_*sigma_)*(-iComplex_*(1/sigma_)*(sinh(yy-ly_/2)*cos(a_*(Pi/lx_*xx-Pi/2))-iComplex_*cosh(yy-ly_/2)*sin(a_*(Pi/lx_*xx-Pi/2))));
	
	for(int j=1;j<=size_;++j)
	{
		if (j%2==0)
		{
			p = p - (-iComplex_*sigma_)*iComplex_*sigma_*(-iComplex_)*(1/sigma_)*vCoeff_[j]*lx_/Pi*pow(-1, (double)j/2)*(iComplex_*(muj_[j]/a_)*cos(j*(Pi/lx_*xx-Pi/2))*muj_[j]*sinh(muj_[j]*(yy-ly_/2)) + (sin(j*(Pi/lx_*xx-Pi/2))*muj_[j]*cosh(muj_[j]*(yy-ly_/2)))/(ComplexNumber)(ComplexNumber)j);
			p = p - (-iComplex_*sigma_)*iComplex_*sigma_*(-iComplex_)*(1/sigma_)*iComplex_*vCoeff_[j]*(a_/j - j/a_)*pow(-1, (double)j/2)*(double)j*(Pi/lx_)*cos(j*(Pi/lx_*xx-Pi/2))*sinh(muj_[j]*(yy-ly_/2));

		}
		else
		{
			p = p - (-iComplex_*sigma_)*iComplex_*sigma_* (-iComplex_*(1/sigma_)*vCoeff_[j]*lx_/Pi*pow(-1,(double)(j - 1)/2)*(-iComplex_*muj_[j]/a_*sin(j*(Pi/lx_*xx-Pi/2))*muj_[j]*cosh(muj_[j]*(yy-ly_/2)) + (cos(j*(Pi/lx_*xx-Pi/2))*muj_[j]*sinh(muj_[j]*(yy-ly_/2)))/(ComplexNumber)j));
			p = p - (-iComplex_*sigma_)*iComplex_*sigma_*(-iComplex_*(1/sigma_)*iComplex_*vCoeff_[j]*(a_/j - j/a_)*pow(-1,(double)(j - 1)/2)*(-j*Pi/lx_)*sin(j*(Pi/lx_*xx-Pi/2))*cosh(muj_[j]*(yy-ly_/2)));
		}
	}
	p = p*cos(Pi*zz/lz_)*(cos(sigma_*t)-iComplex_*sin(sigma_*t));
	return real(p);
}
*/
InitCondU::InitCondU(const ExactSolutionBase* init) :
    InitCond(init)
{
}

void 
InitCondU::operator()(const ElementT* element, const PointReferenceT& pRef, LinearAlgebra::NumericalVector& r) const
{

    PointPhysicalT                    pPhys(3);  // Declare and...
    
    element->referenceToPhysical(pRef, pPhys); // ...transform the point.
//
    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();
//
    double uVal = InitCond::velocity_->getU(pPhys, 0);
//
    
        //cout << "**************************"<<numberOfDegreesOfFreedom<<endl;
    for (unsigned int j=0; j < numberOfDegreesOfFreedom; ++j)
    {
            //
        r[j] = element->basisFunction(j,pRef) * uVal;
    }

}

// initial condition for v
InitCondV::InitCondV(const ExactSolutionBase* init) :
    InitCond(init)
{
}

void 
InitCondV::operator()(const ElementT* element, const PointReferenceT& pRef, LinearAlgebra::NumericalVector& r) const
{
    PointPhysicalT                    pPhys(3);  // Declare and...
    
    element->referenceToPhysical(pRef, pPhys); // ...transform the point.
    
    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();
    
    double vVal = InitCond::velocity_->getV(pPhys, 0);
     
    for (unsigned int j=0; j < numberOfDegreesOfFreedom; ++j)
    {
        r(j) = element->basisFunction(j,pRef) * vVal;
    }
}

// initial condition for w
InitCondW::InitCondW(const ExactSolutionBase* init) :
    InitCond(init)
{
}
 
void 
InitCondW::operator()(const ElementT* element, const PointReferenceT& pRef, LinearAlgebra::NumericalVector& r) const
{
    
    PointPhysicalT                    pPhys(3);  // Declare and...
    
    element->referenceToPhysical(pRef, pPhys); // ...transform the point.
    
    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();
    
    double wVal = InitCond::velocity_->getW(pPhys, 0);
    
    for (unsigned int j=0; j < numberOfDegreesOfFreedom; ++j)
    {
        r(j) = element->basisFunction(j,pRef) * wVal;
    }
}

InitCondR::InitCondR(const ExactSolutionBase* init) :
    InitCond(init)
{
}

void
InitCondR::operator()(const ElementT* element, const PointReferenceT& pRef, LinearAlgebra::NumericalVector& r) const
{
    PointPhysicalT                    pPhys(3);  // Declare and...
    
    element->referenceToPhysical(pRef, pPhys); // ...transform the point.
    
    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();
    
    double pVal = InitCond::velocity_->getR(pPhys, 0);
    
    for (unsigned int j=0; j < numberOfDegreesOfFreedom; ++j)
    {
        r(j) = element->basisFunction(j,pRef) * pVal;
    }
}


InitCondP::InitCondP(const ExactSolutionBase* init) :
    InitCond(init)
{
}

void
InitCondP::operator()(const ElementT* element, const PointReferenceT& pRef, LinearAlgebra::NumericalVector& r) const
{
    PointPhysicalT                    pPhys(3);  // Declare and...
    
    element->referenceToPhysical(pRef, pPhys); // ...transform the point.
    
    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();
    
    double pVal = InitCond::velocity_->getP(pPhys, 0);
    
    for (unsigned int j=0; j < numberOfDegreesOfFreedom; ++j)
    {
        r(j) = element->basisFunction(j,pRef) * pVal;
    }
}

