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

#include <iostream>
#include <string>
using namespace std;

#include <petscmat.h>
#include <petscksp.h>


#include "HEuler.hpp"
    // ============================================================
static char help[] = "Petsc Help";

int main(int argc, char **argv)
{
	try {
		int n, p;
		double P;
		string ProblemName;
		if (argc > 4) {
			ProblemName = string(argv[1]);
			n = std::atoi(argv[2]);
			p = std::atoi(argv[3]);
			P = std::atof(argv[4]); 
			if (n < 0) throw "Number of elements must be positive \n";
			if (p < 0) throw "Polynomial Order must be positive \n";
			if (P < 0) throw "Number of periods must be positive \n";
			argv[4] = argv[0];
			argc -= 4;
			argv += 4;
		} 
		else {
			throw "usage: IGW.out ProblemName n p P [petsc-args] \n";
		}    
    
	PetscInitialize(&argc,&argv,PETSC_NULL,help); 	
	
	PetscMPIInt    rank,size;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank );
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Hello World from %d\n",rank);
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);	

	//MPI_Barrier(PETSC_COMM_WORLD);
	//if (rank == 0)
	//{
    clock_t start, finish;
    	
	unsigned int DIM;
	
	cout << "Number of Periods = " << P << endl;
	//PR(P);

	DIM = 2; 
	HEulerConfigurationData config(DIM, 7, p, 1, HEulerConfigurationData::EBbucketb);//EBTrapezoidWAtau3Over2);

    HEulerGlobalVariables globalVar;
	
	//if (config.solutionType_ == HEulerConfigurationData::EBSingleMode || config.solutionType_ == HEulerConfigurationData::EB2DSingleMode || config.solutionType_ == HEulerConfigurationData::RE0 || config.solutionType_ == HEulerConfigurationData::RES0)
	//{
        config.lx_ = 1;
        config.ly_ = 1;
        config.lz_ = 1;
		
		config.nx_ = n;
		config.ny_ = n;
		config.nz_ = n;	
		
		config.WA_ = 0;
		config.gamma_ = 0.0;
	//}
	
	config.EB_ = 0;
	// add exceptions (problems here)
	
	/* No Longer works:
	 * 
	if (config.solutionType_ == HEulerConfigurationData::RE0)
	{
		config.onePeriod_ = 1.0;//1.804827278170969;
        config.stratified_ = 1;
		config.N2_ = -1.0;//2.0;
		config.incompressible_ = 0;		
	}
	
	if (config.solutionType_ == HEulerConfigurationData::RES0)
	{
		config.onePeriod_ = 0.972666392338948;
        config.stratified_ = 1;
		config.N2_ = 2.0;
		config.incompressible_ = 0;		
	}	
	if (config.solutionType_ == HEulerConfigurationData::RES1)
	{
		config.onePeriod_ = 0.794178783728057;
        config.stratified_ = 1;
		config.N2_ = 2.0;
		config.incompressible_ = 0;		
	}
	if (config.solutionType_ == HEulerConfigurationData::RES2)
	{
		config.onePeriod_ = 0.813923532513087;
        config.stratified_ = 1;
		config.N2_ = 2.0;
		config.incompressible_ = 0;		
	}			
	if (config.solutionType_ == HEulerConfigurationData::RESDO0)
	{
		config.onePeriod_ = 1.375558003710250;
        config.stratified_ = 1;
		config.N2_ = 1.0;
		config.incompressible_ = 0;		
	}	
	if (config.solutionType_ == HEulerConfigurationData::REMTC0)
	{
		config.onePeriod_ = 2.0*3.141592653589793;//1.414213562373095;//
        config.stratified_ = 1;
		config.N2_ = 13.645697703860776;//-1.0;//
		config.incompressible_ = 0;		
	}		
	*/
	if (config.solutionType_ == HEulerConfigurationData::CS1D)
	{
		config.onePeriod_ = 0.972666392338948;//1.0;
        config.stratified_ = 1;
		config.N2_ = 2.0;//-1.0;
		config.incompressible_ = 0;		
	}

	if (config.solutionType_ == HEulerConfigurationData::RES02D)
	{
		config.onePeriod_ = 0.697242057215755;
        config.stratified_ = 1;
		config.N2_ = 2.0;
		config.incompressible_ = 0;		
	}
	
	if (config.solutionType_ == HEulerConfigurationData::CS2D)
	{
		config.onePeriod_ = 0.701506121432065;
        config.stratified_ = 1;
		config.N2_ = 2.0;
		config.incompressible_ = 0;		
	}	
	if (config.solutionType_ == HEulerConfigurationData::CS3D)
	{
		config.onePeriod_ = 0.575103916194890;
        config.stratified_ = 1;
		config.N2_ = 2.0;
		config.incompressible_ = 0;		
	}	
	if (config.solutionType_ == HEulerConfigurationData::ICS2D)
	{
		config.onePeriod_ = 6.322848851930389;
        config.stratified_ = 1; 
		config.N2_ = 2.0;
		config.incompressible_ = 1;		
	}		
	if (config.solutionType_ == HEulerConfigurationData::EB2DSingleMode)
	{
		config.onePeriod_ = 6.283185307179586;
        config.stratified_ = 0;
		config.N2_ = 2.0;//1.0;//
		config.incompressible_ = 1;	
		config.EB_ = 1;
	}
	if (config.solutionType_ == HEulerConfigurationData::EBsemiellipse11mode)
	{
		config.onePeriod_ = 6.283185307179586;
        config.stratified_ = 0;
		config.N2_ = 2.0;
		config.incompressible_ = 1;	
		config.EB_ = 1;
	}	
	if (config.solutionType_ == HEulerConfigurationData::EBbucketb)
	{
		config.onePeriod_ = 6.283185307179586;
        config.stratified_ = 0;
		config.N2_ = 2.0;
		config.incompressible_ = 1;	
		config.EB_ = 1;
	}	
	if (config.solutionType_ == HEulerConfigurationData::EBbucketd)
	{
		config.onePeriod_ = 6.283185307179586;
        config.stratified_ = 0;
		config.N2_ = 2.0;
		config.incompressible_ = 1;	
		config.EB_ = 1;
	}
	if (config.solutionType_ == HEulerConfigurationData::EBbuckete)
	{
		config.onePeriod_ = 6.283185307179586;
        config.stratified_ = 0;
		config.N2_ = 2.0;
		config.incompressible_ = 1;	
		config.EB_ = 1;
	}	
	if (config.solutionType_ == HEulerConfigurationData::EBTrapezoidWAtau3Over2)
	{
		config.onePeriod_ = 6.283185307179586;
        config.stratified_ = 0;
		config.N2_ = 2.0;
		config.incompressible_ = 1;	
		config.EB_ = 1;
	}	
	if (config.solutionType_ == HEulerConfigurationData::EBSingleMode)
	{
		config.onePeriod_ = 5.441398092702654;  
        config.stratified_ = 0;
		config.N2_ = 2.0;
		config.incompressible_ = 1;	
		config.EB_ = 1;
	}	
	if (config.solutionType_ == HEulerConfigurationData::IGW2D)
	{
		config.onePeriod_ = 6.283185307179586; //8.885765876316732;//
        config.stratified_ = 0;
		config.N2_ = 2.0;
		config.incompressible_ = 1;	
		config.EB_ = 1;
		
        config.lx_ = 2;
        config.ly_ = 1;
        config.lz_ = 1;
		
		config.nx_ = 2*n;
		config.ny_ = n;
		config.nz_ = n;				
	}
	
	if (config.solutionType_ == HEulerConfigurationData::IGWN2LinearAiry)
	{
		config.onePeriod_ = 7.695298980971185;
        config.stratified_ = 1;
		config.N2_ = 1.0;
		config.incompressible_ = 1;	
		config.EB_ = 1;
		
        config.lx_ = 0.803250057177826;
        config.ly_ = 1;
        config.lz_ = 1;
		
		config.nx_ = ceil(0.803250057177826*n);
		cout << config.nx_  << endl;
		config.ny_ = n;
		config.nz_ = n;				
	}	
	if (config.solutionType_ == HEulerConfigurationData::IGW2Dprolonged)
	{
		config.onePeriod_ = 6.283185307179586; //8.885765876316732;//
        config.stratified_ = 0;
		config.N2_ = 2.0;
		config.incompressible_ = 1;	
		
        config.lx_ = 6;
        config.ly_ = 1;
        config.lz_ = 1;
		
		config.nx_ = 6*n;
		config.ny_ = n;
		config.nz_ = n;				
	}			
	
	if (config.solutionType_ == HEulerConfigurationData::WA2D)
	{
		config.onePeriod_ = 1; 
        config.stratified_ = 0;
		config.N2_ = 1.0;
		config.incompressible_ = 1;	
		config.WA_ = 1;
		config.gamma_ = 3.141592653589793/20.0;//3.141592653589793/4.0;//0.0;//
	}	
	
	if (config.solutionType_ == HEulerConfigurationData::Lambd3DMixed)
	{
		config.onePeriod_ = 0.707106781186548; 
        config.stratified_ = 1.0;//0;
		config.N2_ = 2.0;//-1.0;
		config.incompressible_ = 0;	
	}
	
	if (config.solutionType_ == HEulerConfigurationData::CS3D)
	{
		config.onePeriod_ = 0.575103916194890;
        config.stratified_ = 1;
		config.N2_ = 2.0;
		config.incompressible_ = 0;		
	}		

	config.numOfPeriods_= static_cast<double>(P)-1.0+0.9999999999;
	config.numOfPeriodsInOnePlotStep_ = static_cast<double>(10*n*n);	//0000////(10*n*n);////config.lx_/config.onePeriod_ //0.5*10*n*n
	config.numOfTimeStepInOnePeriod_ = static_cast<double>(10*n*n); //(10*n*n);//// //1024
	config.numOfHamInOnePeriod_ = static_cast<double>(1);

	int msec;
	
	HEuler eulerProblem(&globalVar, &config);
	
	start = clock();
    eulerProblem.initialiseMesh();
	finish = clock();
	msec = (finish-start) * 1000 / CLOCKS_PER_SEC;
	printf("Time taken for Mesh Initialization: %d ticks  %d milliticks \n", msec/1000, msec%1000);	
	
	start = clock();
	eulerProblem.initialConditions();
	finish = clock();
	msec = (finish-start) * 1000 / CLOCKS_PER_SEC;
	printf("Time taken for Mesh IC: %d ticks  %d milliticks \n", msec/1000, msec%1000);	
	
	eulerProblem.output(0);
   
	if (config.incompressible_ == 0)
	{
		eulerProblem.createCompressibleSystem(); //new
	}
	else if (config.incompressible_ == 1)
	{
		eulerProblem.createIncompressibleSystem();
	}
	
    start = clock();
	//eulerProblem.solve();    	//new
	finish = clock();
	
	eulerProblem.output(config.onePeriod_*config.numOfPeriods_);
	
    cout << "Time = " << finish-start << endl;
    msec = (finish-start) * 1000 / CLOCKS_PER_SEC;
    printf("Time taken for solve: %d ticks  %d milliticks \n", msec/1000, msec%1000);	
	printf(" ");
   // }
	
	//if (rank != 0)
	//{
	//	MPI_Finalize();
	//}
		return 0;
		
	} catch (const char* e) {
		std::cout << e;
	}
}