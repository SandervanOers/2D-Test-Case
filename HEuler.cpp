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

#include "HEuler.hpp"
#include "Base/RectangularMeshDescriptor.hpp"

HEuler::HEuler(HEulerGlobalVariables* global, const HEulerConfigurationData* config):
    Base::HpgemUI(global, config),
    P_(),
    Q_(),
	B_(),
	C_(),
	R_()
{
    switch(config->solutionType_)
	{
		case HEulerConfigurationData::WA2D:
			exactSolution_ = new WA2D();
            break;			
		case HEulerConfigurationData::IGW2D:
			exactSolution_ = new IGW2D();
            break;			
		case HEulerConfigurationData::ICS2D:
			exactSolution_ = new ICS2D();
            break;			
		case HEulerConfigurationData::CS3D:
			exactSolution_ = new CS3D();
            break;			
		case HEulerConfigurationData::CS2D:
			exactSolution_ = new CS2D();
            break;			
		case HEulerConfigurationData::RES02D:
			exactSolution_ = new RES02D();
            break;			
		case HEulerConfigurationData::CS1D:
			exactSolution_ = new CS1D();
            break;			
		case HEulerConfigurationData::REMTC0:
			exactSolution_ = new REMTC0();
            break;				
		case HEulerConfigurationData::RES2:
			exactSolution_ = new RES2();
            break;			
		case HEulerConfigurationData::RES1:
			exactSolution_ = new RES1();
            break;			
		case HEulerConfigurationData::RES0:
			exactSolution_ = new RES0();
            break;			
		case HEulerConfigurationData::RE0:
			exactSolution_ = new RE0();
            break;	
		case HEulerConfigurationData::RESDO0:
			exactSolution_ = new RESDO0();
            break;						
		case HEulerConfigurationData::EB2DSingleMode:
			exactSolution_ = new EB2DSingleMode();
            break;			
		case HEulerConfigurationData::EBSingleMode:
			exactSolution_ = new EBSingleMode();
            break;		
		case HEulerConfigurationData::Lambd3DMixed:
			exactSolution_ = new Lambd3DMixed();
			break;
		case HEulerConfigurationData::IGW2Dprolonged:
			exactSolution_ = new IGW2Dprolonged();
			break;	
		case HEulerConfigurationData::EBsemiellipse11mode:
			exactSolution_ = new EBsemiellipse11mode();
			break;
		case HEulerConfigurationData::IGWN2LinearAiry:
			exactSolution_ = new IGWN2LinearAiry();
			break;
		case HEulerConfigurationData::EBbucketb:
			exactSolution_ = new EBbucketb();
			break;
		case HEulerConfigurationData::EBbucketd:
			exactSolution_ = new EBbucketd();
			break;
		case HEulerConfigurationData::EBbuckete:
			exactSolution_ = new EBbuckete();
			break;
		case HEulerConfigurationData::EBTrapezoidWAtau3Over2:
			exactSolution_ = new EBTrapezoidWAtau3Over2();
			break;			
		default:
			cout<<"Can't happen. Major error!"<<endl;
            break; 
	}


}

HEuler::~HEuler()
{
    
        MatDestroy(&P_);
        MatDestroy(&Q_);
		MatDestroy(&R_);
		MatDestroy(&B_);
		MatDestroy(&C_);
        PetscFinalize();
}

void
HEuler::outputMatrix(Mat& matrix, const string& name)
{
    cout << "Mat Assembly for "<< name<<endl;
    MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY);
    
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_SELF,name.c_str(),&viewer);
    MatView(matrix,viewer);
}
void
HEuler::outputMatrix(const Mat& matrix, const string& name)const 
{
    cout << "Mat Assembly for "<< name<<endl;
    MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY);
    
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_SELF,name.c_str(),&viewer);
    MatView(matrix,viewer);
}
void
HEuler::outputVectorMatlab(const Vec& vec, const string& name)const
{
    cout << "Vec Assembly for "<< name<<endl;
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
    
    PetscViewer viewer;
        //PetscViewerBinaryOpen(PETSC_COMM_WORLD, name.c_str(), FILE_MODE_WRITE, &viewer);
     PetscViewerASCIIOpen(PETSC_COMM_SELF,name.c_str(),&viewer);
    VecView(vec, viewer);
}

void
HEuler::outputVectorMatlab(Vec& vec, const string& name)
{
    cout << "Vec Assembly for "<< name<<endl;
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
    
    PetscViewer viewer;
        //    PetscViewerBinaryOpen(PETSC_COMM_WORLD, name.c_str(), FILE_MODE_WRITE, &viewer);
    PetscViewerASCIIOpen(PETSC_COMM_SELF,name.c_str(),&viewer);

    VecView(vec, viewer);
}

void
HEuler::printFullMatrixInfo(Mat& matrix, const string& name)
{
    PetscInt m=0,n=0;
    MatGetSize(matrix,&m,&n);
    
    MatInfo info;
    MatGetInfo(matrix,MAT_LOCAL, &info);
    cout<<name<<endl;
    //printf("N = %d, N = %d, block_size = %d, memory = %d, assemblies = %d, mallocs = %d, matrix nonzeros (SeqBAIJ format) = %d, allocated nonzeros= %d\n", m, n, (int)info.block_size, (int)info.memory, (int)info.assemblies, (int)info.mallocs,(int)info.nz_used,(int)info.nz_allocated);
	printf("N = %d, N = %d, block_size = %d, memory = %u, assemblies = %d, mallocs = %d, matrix nonzeros (SeqBAIJ format) = %d, allocated nonzeros= %d\n", m, n, (int)info.block_size, (unsigned int)info.memory, (int)info.assemblies, (int)info.mallocs,(int)info.nz_used,(int)info.nz_allocated);
}

bool
HEuler::initialiseMesh()
{
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
	MeshId id;
	if (config->solutionType_== HEulerConfigurationData::EBsemiellipse11mode)//HEulerConfigurationData::EB2DSingleMode)
	{
		std::stringstream filename;
		filename << "../../../../mesh/semiellipse2n"<<config->nx_<<".hyb"; 
		std::string myString = filename.str();
		const char* myPointer = myString.c_str();
		MeshId id = addMesh(myPointer);	
	}
	else if (config->solutionType_ == HEulerConfigurationData::EBbucketb ||config->solutionType_ == HEulerConfigurationData::EBbucketd || config->solutionType_ == HEulerConfigurationData::EBbuckete)
	{
		std::stringstream filename;
		filename << "../../../../mesh/bucket2n"<<config->nx_<<".hyb"; 
		std::string myString = filename.str();
		const char* myPointer = myString.c_str();
		MeshId id = addMesh(myPointer);	
	}
	else if (config->solutionType_ == HEulerConfigurationData::EBTrapezoidWAtau3Over2)
	{
		std::stringstream filename;
		filename << "../../../../mesh/trapezoidn"<<config->nx_<<".hyb"; 
		std::string myString = filename.str();
		const char* myPointer = myString.c_str();
		MeshId id = addMesh(myPointer);	
	}
	else
	{
	RectangularMeshDescriptor rectangularMesh(config->DIM_);
	
    if (config->solutionType_==HEulerConfigurationData::EBSingleMode)
    {
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor::SOLID_WALL;//PERIODIC;
        rectangularMesh.boundaryConditions_[1] = RectangularMeshDescriptor::SOLID_WALL;//PERIODIC;
		rectangularMesh.boundaryConditions_[2] = RectangularMeshDescriptor::SOLID_WALL;
		rectangularMesh.bottomLeft_[0]       = 0;
		rectangularMesh.bottomLeft_[1]       = 0;
		rectangularMesh.bottomLeft_[2]       = 0;
		rectangularMesh.topRight_[0]          = config->lx_;
		rectangularMesh.topRight_[1]          = config->ly_;
		rectangularMesh.topRight_[2]          = config->lz_;
		rectangularMesh.numElementsInDIM_[0] = config->nx_;
		rectangularMesh.numElementsInDIM_[1] = config->ny_;
		rectangularMesh.numElementsInDIM_[2] = config->nz_;
    }
	
    if (config->solutionType_==HEulerConfigurationData::EB2DSingleMode || config->solutionType_==HEulerConfigurationData::WA2D)
    {
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor::SOLID_WALL;//PERIODIC;
        rectangularMesh.boundaryConditions_[1] = RectangularMeshDescriptor::SOLID_WALL;
		rectangularMesh.bottomLeft_[0]       = 0;
		rectangularMesh.bottomLeft_[1]       = 0;
		rectangularMesh.topRight_[0]          = config->lx_;
		rectangularMesh.topRight_[1]          = config->ly_;
		rectangularMesh.numElementsInDIM_[0] = config->nx_;
		rectangularMesh.numElementsInDIM_[1] = config->ny_;
    }	
	
    if (config->solutionType_==HEulerConfigurationData::IGW2D || config->solutionType_==HEulerConfigurationData::IGW2Dprolonged || config->solutionType_==HEulerConfigurationData::IGWN2LinearAiry)
    {
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor::PERIODIC;
        rectangularMesh.boundaryConditions_[1] = RectangularMeshDescriptor::SOLID_WALL;
		rectangularMesh.bottomLeft_[0]       = 0;
		rectangularMesh.bottomLeft_[1]       = 0;
		rectangularMesh.topRight_[0]          = config->lx_;
		rectangularMesh.topRight_[1]          = config->ly_;
		rectangularMesh.numElementsInDIM_[0] = config->nx_;
		rectangularMesh.numElementsInDIM_[1] = config->ny_;
    }		
	
    if (config->solutionType_==HEulerConfigurationData::RES02D || config->solutionType_==HEulerConfigurationData::CS2D || config->solutionType_==HEulerConfigurationData::ICS2D)
    {
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor::SOLID_WALL;
        rectangularMesh.boundaryConditions_[1] = RectangularMeshDescriptor::SOLID_WALL;
		rectangularMesh.bottomLeft_[0]       = 0;
		rectangularMesh.bottomLeft_[1]       = 0;
		rectangularMesh.topRight_[0]          = config->lx_;
		rectangularMesh.topRight_[1]          = config->ly_;
		rectangularMesh.numElementsInDIM_[0] = config->nx_;
		rectangularMesh.numElementsInDIM_[1] = config->ny_;
    }	

    if (config->solutionType_==HEulerConfigurationData::CS3D)
    {
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor::SOLID_WALL;
        rectangularMesh.boundaryConditions_[1] = RectangularMeshDescriptor::SOLID_WALL;
		rectangularMesh.boundaryConditions_[2] = RectangularMeshDescriptor::SOLID_WALL;
		rectangularMesh.bottomLeft_[0]       = 0;
		rectangularMesh.bottomLeft_[1]       = 0;
		rectangularMesh.bottomLeft_[2]       = 0;
		rectangularMesh.topRight_[0]          = config->lx_;
		rectangularMesh.topRight_[1]          = config->ly_;
		rectangularMesh.topRight_[2]          = config->lz_;
		rectangularMesh.numElementsInDIM_[0] = config->nx_;
		rectangularMesh.numElementsInDIM_[1] = config->ny_;
		rectangularMesh.numElementsInDIM_[2] = config->nz_;
    }	
	
	if (config->solutionType_==HEulerConfigurationData::Lambd3DMixed)
	{
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor::PERIODIC;
        rectangularMesh.boundaryConditions_[1] = RectangularMeshDescriptor::PERIODIC;
		rectangularMesh.boundaryConditions_[2] = RectangularMeshDescriptor::SOLID_WALL;
		rectangularMesh.bottomLeft_[0]       = 0;
		rectangularMesh.bottomLeft_[1]       = 0;
		rectangularMesh.bottomLeft_[2]       = 0;
		rectangularMesh.topRight_[0]          = config->lx_;
		rectangularMesh.topRight_[1]          = config->ly_;
		rectangularMesh.topRight_[2]          = config->lz_;
		rectangularMesh.numElementsInDIM_[0] = config->nx_;
		rectangularMesh.numElementsInDIM_[1] = config->ny_;
		rectangularMesh.numElementsInDIM_[2] = config->nz_;		
	}

    if (config->solutionType_==HEulerConfigurationData::RE0 || config->solutionType_==HEulerConfigurationData::RES0 || config->solutionType_==HEulerConfigurationData::RES1 || config->solutionType_==HEulerConfigurationData::RES2 || config->solutionType_==HEulerConfigurationData::RESDO0 || config->solutionType_==HEulerConfigurationData::REMTC0 || config->solutionType_==HEulerConfigurationData::CS1D) 
    {
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor::SOLID_WALL; //PERIODIC;//
		rectangularMesh.bottomLeft_[0]       = 0;
		rectangularMesh.topRight_[0]          = config->lx_;
		rectangularMesh.numElementsInDIM_[0] = config->nx_;
    }		
	/*
    rectangularMesh.bottomLeft_[0]       = 0;
    rectangularMesh.bottomLeft_[1]       = 0;
    rectangularMesh.bottomLeft_[2]       = 0;
    rectangularMesh.topRight_[0]          = config->lx_;
    rectangularMesh.topRight_[1]          = config->ly_;
    rectangularMesh.topRight_[2]          = config->lz_;
    rectangularMesh.numElementsInDIM_[0] = config->nx_;
    rectangularMesh.numElementsInDIM_[1] = config->ny_;
    rectangularMesh.numElementsInDIM_[2] = config->nz_;
    */
    MeshId id = addMesh(rectangularMesh, RECTANGULAR);//TRIANGULAR);//
    }
	
    HEulerGlobalVariables* globalData = static_cast<HEulerGlobalVariables*>(globalData_);
    globalData->nElements_ = getNumberOfElements(id);
    globalData->dt_         = config->onePeriod_/config->numOfTimeStepInOnePeriod_;	
	
    return true;
}

    /// create Mass Matrices, store them as a User Defined Element Data
    /// calculate projection of the every unknown on the FEM spaces.
void
HEuler::initialConditions()
{
    unsigned int ldof = static_cast<const HEulerConfigurationData*>(configData_)->numberOfBasisFunctions_;
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);

    HEulerElementData* elemData;
    ElementIntegralBackgroundDensity gradMass;
    bool useCache = false;
    
    ElementIntegralT   elIntegral(useCache);
    
    typedef void  (HEuler::*Integrand)(const ElementT* , const PointReferenceT&, ElementIntegralBackgroundDensity&);
    //Integrand massMatrixIntegrand = &HEuler::calculateMassMatrix;
    
    LinearAlgebra::NumericalVector rightHand(ldof);
    
    
    ElementT::SolutionVector       numerical(ldof);
    
    LinearAlgebra::Matrix          invMassM(ldof, ldof);
    
    cout << "ldof="<<ldof<<endl;
    
    InitCondU       uEx(exactSolution_);
    InitCondV       vEx(exactSolution_);
    InitCondW       wEx(exactSolution_);
	InitCondR       rhoEx(exactSolution_);
	InitCondP       pEx(exactSolution_);
   
    cout <<"start calculations of inital condition!"<<endl;
    unsigned int count =0;
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {
        ElementT* elem =(*cit);
        
        elemData = new HEulerElementData(ldof);
        
        LinearAlgebra::Matrix& massMatrix = elemData->massMatrix_;
        LinearAlgebra::Matrix& invMassM   = elemData->invMassMatrix_;
		LinearAlgebra::Matrix& massMatrixOverBackgroundDensity = elemData->massMatrixOverBackgroundDensity_;
		LinearAlgebra::Matrix& massMatrixTimesBackgroundDensity = elemData->massMatrixTimesBackgroundDensity_;
		LinearAlgebra::Matrix& massMatrixTimesBackgroundDensityDeriv = elemData->massMatrixTimesBackgroundDensityDeriv_;
		LinearAlgebra::Matrix& massMatrixOverBackgroundDensityDeriv = elemData->massMatrixOverBackgroundDensityDeriv_;
		LinearAlgebra::Matrix& massMatrixTimesN2 = elemData->massMatrixTimesN2_;
		LinearAlgebra::Matrix& massMatrixOverN2 = elemData->massMatrixOverN2_;
        
        elIntegral.integrate<ElementIntegralBackgroundDensity>(elem, this, gradMass);
		massMatrix								= gradMass.massMatrixR_;
        massMatrixOverBackgroundDensity 		= gradMass.massMatrixOverBackgroundDensityR_;		
		massMatrixTimesBackgroundDensity	 	= gradMass.massMatrixTimesBackgroundDensityR_;		
		massMatrixTimesBackgroundDensityDeriv 	= gradMass.massMatrixTimesBackgroundDensityDerivR_;
		massMatrixOverBackgroundDensityDeriv 	= gradMass.massMatrixOverBackgroundDensityDerivR_;
		massMatrixTimesN2					 	= gradMass.massMatrixTimesN2R_;
		massMatrixOverN2					 	= gradMass.massMatrixOverN2R_;		
		
        massMatrix.inverse(invMassM);
        
        elem->setUserData(elemData);
        
		if (config->DIM_ == 3) 
		{
			elIntegral.integrate(elem, &wEx, rightHand);
			numerical = invMassM*rightHand;// projection of W
			elem->setTimeLevelData(0,2,numerical);			

		}
		if (config->DIM_ >= 2)
		{
			elIntegral.integrate(elem, &vEx, rightHand);
			numerical = invMassM*rightHand;// projection of V
			elem->setTimeLevelData(0,1,numerical);			
		}
		
			elIntegral.integrate(elem, &uEx, rightHand);
			numerical = invMassM*rightHand;// projection of U
			elem->setTimeLevelData(0,0,numerical);		
					
			elIntegral.integrate(elem, &rhoEx, rightHand);
			numerical = invMassM*rightHand;	// projection of rho
			elem->setTimeLevelData(0,3,numerical);
			
			elIntegral.integrate(elem, &pEx, rightHand);
			numerical = invMassM*rightHand;// projection of P
			elem->setTimeLevelData(0,4,numerical);	
    }
    cout <<"finish calculations of inital condition!"<<endl;
}

void //computes the mass matrix
HEuler::elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& massMatrix)
{
}
void //computes the mass matrix and 1Mij
HEuler::elementIntegrand(const ElementT* element, const PointReferenceT& p, ElementIntegralBackgroundDensity& returnObject)
{
    unsigned int numOfDegreesOfFreedom = element->getNrOfBasisFunctions();
    LinearAlgebra::Matrix& massMatrix = returnObject.massMatrixR_;
    LinearAlgebra::Matrix& massMatrixOverBackgroundDensity = returnObject.massMatrixOverBackgroundDensityR_;
	LinearAlgebra::Matrix& massMatrixTimesBackgroundDensity = returnObject.massMatrixTimesBackgroundDensityR_;
	LinearAlgebra::Matrix& massMatrixTimesBackgroundDensityDeriv = returnObject.massMatrixTimesBackgroundDensityDerivR_;
	LinearAlgebra::Matrix& massMatrixOverBackgroundDensityDeriv = returnObject.massMatrixOverBackgroundDensityDerivR_;
	LinearAlgebra::Matrix& massMatrixTimesN2					 = returnObject.massMatrixTimesN2R_;
	LinearAlgebra::Matrix& massMatrixOverN2					 = returnObject.massMatrixOverN2R_;
	
    massMatrix.resize(numOfDegreesOfFreedom, numOfDegreesOfFreedom);
    massMatrixOverBackgroundDensity.resize(numOfDegreesOfFreedom, numOfDegreesOfFreedom);
	massMatrixTimesBackgroundDensity.resize(numOfDegreesOfFreedom, numOfDegreesOfFreedom);
	massMatrixTimesBackgroundDensityDeriv.resize(numOfDegreesOfFreedom, numOfDegreesOfFreedom);
	massMatrixOverBackgroundDensityDeriv.resize(numOfDegreesOfFreedom, numOfDegreesOfFreedom);
	massMatrixTimesN2.resize(numOfDegreesOfFreedom, numOfDegreesOfFreedom);
	massMatrixOverN2.resize(numOfDegreesOfFreedom, numOfDegreesOfFreedom);

	double rho0 = 1.0;
	double rho0Deriv = 0.0;
	double N2 = 1.0;
	const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
	if (config->stratified_==1 && config->EB_ == 0)
	{
		PointPhysical pPhys(config->DIM_);
		element->referenceToPhysical(p, pPhys);
		rho0		= exactSolution_->getrho0(pPhys, 0.0);
		rho0Deriv	= exactSolution_->getrho0Deriv(pPhys, 0.0);
	}	
	if (config->stratified_ == 1 || config->EB_ == 1) //EB no longer counts as stratified, while we need it (N2 not always equal to one)
	{
		PointPhysical pPhys(config->DIM_);
		element->referenceToPhysical(p, pPhys);
		N2 	= exactSolution_->getN2(pPhys,0.0);
	}
    for (unsigned int i=0; i < numOfDegreesOfFreedom; ++i)
    {
        for (unsigned int j=0; j < numOfDegreesOfFreedom; ++j)
        {
            massMatrix(i,j) = element->basisFunction(j,p) * element->basisFunction(i,p);
			massMatrixOverBackgroundDensity(i,j)  = element->basisFunction(j,p) * element->basisFunction(i,p)/rho0;
			massMatrixTimesBackgroundDensity(i,j)  = element->basisFunction(j,p) * element->basisFunction(i,p)*rho0;
			massMatrixTimesBackgroundDensityDeriv(i,j)  = element->basisFunction(j,p) * element->basisFunction(i,p)*rho0Deriv;
			massMatrixOverBackgroundDensityDeriv(i,j)  = element->basisFunction(j,p) * element->basisFunction(i,p)/rho0Deriv;
			massMatrixTimesN2(i,j)  = element->basisFunction(j,p) * element->basisFunction(i,p)*N2*rho0;
			massMatrixOverN2(i,j)  = element->basisFunction(j,p) * element->basisFunction(i,p)/N2/rho0;			
        }
		
        /*for (unsigned int j=0; j <= i; ++j)
        {
            massMatrix(i,j) = massMatrix(j,i) = element->basisFunction(j,p) * element->basisFunction(i,p);
			massMatrixOverBackgroundDensity(i,j) = massMatrixOverBackgroundDensity(j,i) = element->basisFunction(j,p) * element->basisFunction(i,p)/rho0;
			massMatrixTimesBackgroundDensity(i,j) = massMatrixTimesBackgroundDensity(j,i) = element->basisFunction(j,p) * element->basisFunction(i,p)*rho0;
			massMatrixTimesBackgroundDensityDeriv(i,j) = massMatrixTimesBackgroundDensityDeriv(j,i) = element->basisFunction(j,p) * element->basisFunction(i,p)*rho0Deriv;
        }*/		
    }

}

void
HEuler::elementIntegrand(const ElementT* element, const PointReferenceT& p, ElementIntegralData& returnObject)
{
    
    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();
    
    LinearAlgebra::Matrix& xDerReturnData = returnObject.xGrad_;
    LinearAlgebra::Matrix& yDerReturnData = returnObject.yGrad_;
    LinearAlgebra::Matrix& zDerReturnData = returnObject.zGrad_;
    
    xDerReturnData.resize(numberOfDegreesOfFreedom, numberOfDegreesOfFreedom);
    yDerReturnData.resize(numberOfDegreesOfFreedom, numberOfDegreesOfFreedom);
    zDerReturnData.resize(numberOfDegreesOfFreedom, numberOfDegreesOfFreedom);
    
	double rho0 = 1.0;
	double rho0Derivx = 0.0;
	double rho0Derivy = 0.0;
	double rho0Derivz = 0.0;
	const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
	if (config->stratified_==1 && config->EB_ == 0)
	{
		PointPhysical pPhys(config->DIM_);
		element->referenceToPhysical(p, pPhys);
		rho0		= exactSolution_->getrho0(pPhys, 0.0);	
		if (config->DIM_ == 1)
		{
			rho0Derivx	= exactSolution_->getrho0Deriv(pPhys, 0.0);	
		}	
		if (config->DIM_ == 2)
		{
			rho0Derivy	= exactSolution_->getrho0Deriv(pPhys, 0.0);	
		}	
		if (config->DIM_ == 3)
		{
			rho0Derivz	= exactSolution_->getrho0Deriv(pPhys, 0.0);	
		}	
	}	
    
    LinearAlgebra::NumericalVector grads(config->DIM_);
    
    for (unsigned int i=0; i < numberOfDegreesOfFreedom; ++i)
    {
    	element->basisFunctionDeriv(i,p,grads);
        //Utilities::PhysGradientOfBasisFunction obj(element, i);
        //obj(p, grads);
        
        for (unsigned int j=0; j < numberOfDegreesOfFreedom; ++j)
        {
            xDerReturnData(i,j) = element->basisFunction(j,p) * grads[0]*rho0 + element->basisFunction(j,p) * element->basisFunction(i,p)*rho0Derivx;
            yDerReturnData(i,j) = element->basisFunction(j,p) * grads[1]*rho0 + element->basisFunction(j,p) * element->basisFunction(i,p)*rho0Derivy;
            zDerReturnData(i,j) = element->basisFunction(j,p) * grads[2]*rho0 + element->basisFunction(j,p) * element->basisFunction(i,p)*rho0Derivz;
        }
    }
}

void
HEuler::faceIntegrand(const FaceT* face,          const LinearAlgebra::NumericalVector& normal,
                   const PointReferenceOnTheFaceT& p,  FluxData& ret)
{
    if (face->isInternal())
    {
		//cout << "Internal" << endl;
		const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
		//const ExactSolutionBase*  		     velocity_;
		
        const double magn                     = Utilities::norm2(normal);
        unsigned int numberOfDegreesOfFreedom = face->getPtrElementLeft()->getNrOfBasisFunctions();
        
        PointReferenceT 	pL(3), pR(3);
        double              bFL, bFR, BFevalL, BFevalR;
        
        double theta = 0.5;// = 0.357;//((double) rand() / (RAND_MAX));//static_cast<const HEulerConfigurationData*>(configData_)->theta_; //
        //cout << theta << endl; //config->theta_;//
        double nx    = normal[0]/magn;
        double ny    = normal[1]/magn;
        double nz    = normal[2]/magn;
        
        const ElementT* const left   = face->getPtrElementLeft();
        const ElementT* const right  = face->getPtrElementRight();
        
        face->mapRefFaceToRefElemL(p, pL);
        face->mapRefFaceToRefElemR(p, pR);
		
		double rho0L, rho0R;
		rho0L = 1.0;
		rho0R = 1.0;
		if (config->stratified_==1 && config->EB_ == 0)
		{
			PointPhysical pPhysL(config->DIM_);
			PointPhysical pPhysR(config->DIM_);
			left->referenceToPhysical(pL, pPhysL);
			right->referenceToPhysical(pR, pPhysR);
			rho0L		= exactSolution_->getrho0(pPhysL, 0.0);		
			rho0R		= exactSolution_->getrho0(pPhysR, 0.0);	
		}
		
		//if (rho0L > 0.223130160148430)
		//{
			//theta = 0.75;//25;
        //}
		//else 
		//{
		//	theta = 0.5;
		//}
        for (int j = 0; j< numberOfDegreesOfFreedom; ++j)
        {
            bFL 	= 	left->basisFunction(j,pL);
            bFR 	=   right->basisFunction(j, pR);
            
            LinearAlgebra::Matrix& leftReturnData   = ret.left_[j];
            
            LinearAlgebra::Matrix& rightReturnData  = ret.right_[j];
            

            for (unsigned int i=0; i < numberOfDegreesOfFreedom; ++i)
            {
                
                BFevalL = left->basisFunction(i,pL);
                
                BFevalR = right->basisFunction(i,pR);
                
                leftReturnData(0,i) = 	BFevalL*nx*(1-theta)*bFL*rho0L;
                leftReturnData(1,i) = 	BFevalL*ny*(1-theta)*bFL*rho0L;
                leftReturnData(2,i) = 	BFevalL*nz*(1-theta)*bFL*rho0L;
                leftReturnData(3,i) = 	BFevalR*nx*(1-theta)*bFL*rho0L;
                leftReturnData(4,i) = 	BFevalR*ny*(1-theta)*bFL*rho0L;
                leftReturnData(5,i) = 	BFevalR*nz*(1-theta)*bFL*rho0L;
				
                rightReturnData(0,i) = 	BFevalL*nx*(theta)*bFR*rho0R;
                rightReturnData(1,i) = 	BFevalL*ny*(theta)*bFR*rho0R;
                rightReturnData(2,i) = 	BFevalL*nz*(theta)*bFR*rho0R;
                rightReturnData(3,i) = 	BFevalR*nx*(theta)*bFR*rho0R;
                rightReturnData(4,i) = 	BFevalR*ny*(theta)*bFR*rho0R;
                rightReturnData(5,i) = 	BFevalR*nz*(theta)*bFR*rho0R;
//456                
                leftReturnData(6,i) = 	bFL*nx*(1-theta)*BFevalL*rho0L;
                leftReturnData(7,i) = 	bFL*ny*(1-theta)*BFevalL*rho0L;
                leftReturnData(8,i) = 	bFL*nz*(1-theta)*BFevalL*rho0L;
                leftReturnData(9,i) = 	bFL*nx*(theta)*BFevalR*rho0L;
                leftReturnData(10,i) = 	bFL*ny*(theta)*BFevalR*rho0L;
                leftReturnData(11,i) = 	bFL*nz*(theta)*BFevalR*rho0L;
                   
                rightReturnData(6,i) = 	bFR*nx*(1-theta)*BFevalL*rho0R;
                rightReturnData(7,i) = 	bFR*ny*(1-theta)*BFevalL*rho0R;
                rightReturnData(8,i) = 	bFR*nz*(1-theta)*BFevalL*rho0R;
                rightReturnData(9,i) = 	bFR*nx*(theta)*BFevalR*rho0R;
                rightReturnData(10,i) = bFR*ny*(theta)*BFevalR*rho0R;
                rightReturnData(11,i) = bFR*nz*(theta)*BFevalR*rho0R;       
            }
        }
    }
	else {
		//cout << "External" << endl;
	}
    
}

void
HEuler::correctInitialProjectionOfVelocity(const Vec& UInit, Vec& UCorrected)const
{
	const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
	unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nu     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
    unsigned int Nv     = 0;
	unsigned int Nw     = 0;
	if (config->DIM_ == 3)
	{
		Nv = Nu;
	}
	if (config->DIM_ >= 2)
	{
		Nw = Nu;
	}	
    
    unsigned int Nl     = Nu;
    unsigned int N      = Nu+Nv+Nw;

	// globalDIV is DIV in createIncompressibleSystem
	Mat& globalDIV      = static_cast<HEulerGlobalVariables*>(globalData_)->DivergenceFreeMatrix_;
	
	PetscInt a(10);
	PetscReal max;
    
	double reltol = 1.e-16;
	double abstol = 1.e-16;
	
    Vec RHS, MU;
	
	VecCreateSeq(PETSC_COMM_SELF, Nl, &RHS);
	VecCreateSeq(PETSC_COMM_SELF, N, &UCorrected);
	
	MatMult(globalDIV, UInit, RHS);
        // 	outputVectorMatlab(RHS,"rhs.dat");
	
    
	VecScale(RHS,-1);
	VecAssemblyBegin(RHS);
	VecAssemblyEnd(RHS);
	VecMax(RHS, &a, &max);
	cout << "DIV U Error ="<<max<<endl;
	
	if (max > abstol)
	{
			//[A=DIV*DIV'/6;]
            //    [mu=gmres(A,b,10, 1.e-15, 1.e+8);]
		KSP ksp;
		
            // Preconditioner
		PC pc;
		KSPCreate(PETSC_COMM_SELF, &ksp);
		//KSPSetOperators(ksp, globalDIV, globalDIV, SAME_NONZERO_PATTERN);
		KSPSetOperators(ksp, globalDIV, globalDIV);
		cout << "2"<<endl;
        
            //	        KSPSetType(ksp, KSPPREONLY);
		KSPSetType(ksp, KSPLSQR);
		cout << "3"<<endl;
		KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
		cout << "4"<<endl;
		
		
		KSPGetPC(ksp,&pc);
		cout << "5"<<endl;
        
		
		cout << "6"<<endl;
		
		KSPSetFromOptions(ksp);
            //PCFactorSetLevels(pc, 1);
            //PCFactorSetFill(pc, 2.945);
		
		cout << "7"<<endl;
		
		PCSetType(pc, PCNONE);
		
        
		KSPSetTolerances(ksp, reltol, abstol, PETSC_DEFAULT, PETSC_DEFAULT);
		KSPSetUp(ksp);
		
		/*! Solving linear system.*/
		cout << "8"<<endl;
		KSPSolve(ksp, RHS, UCorrected);
            //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
		
		KSPConvergedReason reason;
		KSPGetConvergedReason(ksp, &reason);
			//CHKERRQ(ierr);
		if(reason < 0)
		{
			PetscInt its;
			KSPGetIterationNumber(ksp, &its);
			PetscReal rnorm;
			KSPGetResidualNorm(ksp, &rnorm);
			cout << "\t\tPetSc: Solving the linear system has failed. Reason code: "
            << reason << endl << "Check KSPConvergedReason for the reason" << endl
            << "\t\tPetSc: Residual after " << int(its) << " iterations : ";
			cout.setf(ios::scientific, ios::floatfield);
			cout.precision(4);
			cout << rnorm << endl;
			cout.setf(ios::fixed, ios::floatfield);
			cout.precision(5);
		}
        
        else
		{
			if(reason > 0)
			{
				PetscInt its;
				KSPGetIterationNumber(ksp, &its);
				PetscReal rnorm;
				KSPGetResidualNorm(ksp, &rnorm);
				cout << "\t\tPetSc: Solving the linear system has succeeded. Reason code: "
                << reason << endl << "Check KSPConvergedReason for the reason" << endl
                << "\t\tPetsc: Convergence in " << int(its) << " iterations : ";
				cout.setf(ios::scientific, ios::floatfield);
				cout.precision(4);
				cout << rnorm << endl;
				cout.setf(ios::fixed, ios::floatfield);
				cout.precision(5);
			}
			else
			{
				cout << "\t\tPetSc: Solving the linear system is still under way" << endl;
			}
		}
			//[u=velocity+DIV'*mu;]
            //outputVectorMatlab(UCorrected,"mu.dat");
        
        cout << "UC is ready!"<<endl;
        
        a	= 10;
        max	= 0;
        VecMax(UCorrected, &a, &max);
        cout << "Difference between UInit and UCorrected "<<max<<endl;
			//[maxE=max(u+velocity);]
        VecAXPY(UCorrected, 1, UInit);
            //   [error=DIV*u;]
        KSPDestroy(&ksp);
	}
	else
	{
		VecCopy(UInit,UCorrected); 
	}
	
    cout << "*************************************************************************"<<endl;
	
    VecDestroy(&RHS);

	PetscScalar* XTEMP = new PetscScalar [3*Nl];
	VecGetArray(UCorrected, &XTEMP);

    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {
        Base::Element* element= *cit;
        int k = element->getID();

        unsigned int pos = k*nb;

        for (int i =0; i<nb;++i)
        {
            element->setData(0,0,i, XTEMP[pos+i]);   //set U
			if (config->DIM_ == 3)
			{
				element->setData(0,2,i, XTEMP[Nu+Nv+pos+i]);//set W
			}
            
			element->setData(0,1,i, XTEMP[Nu+pos+i]);//set V or W 
        }
    }
    VecRestoreArray(UCorrected, &XTEMP);

}

void
HEuler::createCompressibleSystem()
{
 	PetscInfoAllow(PETSC_TRUE, "history.txt");
	
	//Mat C;
	Mat BF;
    
	//Mat C1;
	Mat BF1;
	Mat MInv;
	Mat sMInv;
	//Mat A;
	//Mat Ah;
    Mat globalDIV;
	Mat DIV;
	
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);	
	
    unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nw     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
	unsigned int Nv		= 0;
	unsigned int Nu     = 0;
	if (config->DIM_ == 3)
	{
         Nv= Nw;
	}
	if (config->DIM_ >= 2)
	{
         Nu = Nw;
	}	
    
    unsigned int Nl     = Nw;
    unsigned int N      = Nu+Nv+Nw;	
	
	double N2			= config->N2_;
	
    double dt           = static_cast<HEulerGlobalVariables*>(globalData_)->dt_;
    cout << "dt=" <<dt<<endl;	
	
	cout << "**************Starting create Matrices*************"<<endl;
	Mat M1;
	Mat sM1;
	Mat VecNMat;
	Mat NMat;
 	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    1*nb, 				PETSC_NULL, &M1);
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &sM1);
	
 	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    1*nb, 				PETSC_NULL, &MInv);
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &sMInv);	
	
	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nw, 	Nl,          1*nb, 				PETSC_NULL, &VecNMat);
	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nw, 	Nl,          1*nb, 				PETSC_NULL, &NMat);		
	
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nl,         (7)*nb, 			PETSC_NULL, &BF);//number of possible nonzero blocks are 7: element and his 6 neighbours
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 	Nu+Nv+Nw,    (7*3)*nb, 				PETSC_NULL, &globalDIV);//number of possible nonzero blocks are 7: element and his 6 neighbours	
	
	cout << "**************Mass Matrix Calculated*************"<<endl;

//	timer t0;
//	std::cout << t0.elapsed() << " - : Element Integration started\n";
    ElementIntegralData gradMass;
    bool useCache = false;
    ElementIntegralT   elIntegral(useCache);
    typedef void  (HEuler::*Integrand)(const ElementT* , const PointReferenceT&, ElementIntegralData&);
    Integrand gradMassInteg = &HEuler::elementIntegrand;
    
    std::cout << " - : Element Integration started\n";
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {

        int m  = (*cit)->getID();
        unsigned int pos = (*cit)->getID()*nb;
        
        HEulerElementData* elementData = static_cast<HEulerElementData*>((*cit)->getUserData());
        
        ElementT* element =(*cit);
        
        elIntegral.integrate<ElementIntegralData>(element, this, gradMass);
		
		for (unsigned int j = 0; j<nb;++j)
		{
			for (unsigned int i=0; i<nb;++i)//can be optimized later!!!
 			{
                //MatSetValue(NMat, 	Nu+Nv+pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensity_(j,i), 		ADD_VALUES);
				//MatSetValue(VecNMat, 	Nu+Nv+pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensityDeriv_(j,i), 	ADD_VALUES);
                
				//MatSetValue(NMat, 		pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensity_(j,i), 		ADD_VALUES);
				//MatSetValue(VecNMat, 	pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensityDeriv_(j,i), 	ADD_VALUES);	

                //MatSetValue(NMatT, 		pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensity_(i,j), 		ADD_VALUES);
				//MatSetValue(VecNMatT, 	pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensityDeriv_(i,j), 	ADD_VALUES);	
                
                MatSetValue(M1,   		pos+j, 			pos+i,     	elementData->massMatrixOverBackgroundDensity_(j,i), 		ADD_VALUES);
				MatSetValue(sM1,   		pos+j, 			pos+i,      elementData->massMatrixOverBackgroundDensity_(j,i), 		ADD_VALUES);
            
                MatSetValue(MInv, 		pos+j,       	pos+i,    	elementData->invMassMatrix_(j,i) , 							ADD_VALUES);
                MatSetValue(sMInv, 		pos+j, 			pos+i,      elementData->invMassMatrix_(j,i) , 							ADD_VALUES);
                
                MatSetValue(BF, 		pos+j, 		 	pos+i,   	-gradMass.xGrad_(i,j),  									ADD_VALUES);
				MatSetValue(globalDIV, 	pos+j, 			pos+i,   	gradMass.xGrad_(j,i),  										ADD_VALUES);

				if (config->DIM_ >= 2)
				{
					MatSetValue(MInv, 		Nu+pos+j, 		Nu+pos+i, 	elementData->invMassMatrix_(j,i) , 							ADD_VALUES);
					MatSetValue(M1,   		Nu+pos+j, 		Nu+pos+i,   elementData->massMatrixOverBackgroundDensity_(j,i), 		ADD_VALUES);
					MatSetValue(BF, 		Nu+pos+j, 	 	pos+i,   	-gradMass.yGrad_(i,j),  									ADD_VALUES);
					MatSetValue(globalDIV, 	pos+j, 			Nu+pos+i,   gradMass.yGrad_(j,i),  										ADD_VALUES);					
				}				
				if (config->DIM_ == 3)
				{
					MatSetValue(MInv, 		Nu+Nv+pos+j,    Nu+Nv+pos+i,    elementData->invMassMatrix_(j,i), 					ADD_VALUES);	
					MatSetValue(M1,   		Nu+Nv+pos+j, 	Nu+Nv+pos+i,    elementData->massMatrixOverBackgroundDensity_(j,i), ADD_VALUES);					
					MatSetValue(BF, 		Nu+Nv+pos+j, 	pos+i,   		-gradMass.zGrad_(i,j),  							ADD_VALUES);
					MatSetValue(globalDIV, 	pos+j, 			Nu+Nv+pos+i,	gradMass.zGrad_(j,i),  								ADD_VALUES);
				}			
                                                                                          
 			}
		}
    }
	cout << "Mat Assembly for "<< "MINV"<<endl;
    MatAssemblyBegin(MInv,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MInv,MAT_FINAL_ASSEMBLY);
	printFullMatrixInfo(MInv, "MInv");
    cout << "Mat Assembly for "<< "MINVSm"<<endl;
    MatAssemblyBegin(sMInv,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sMInv,MAT_FINAL_ASSEMBLY);
	printFullMatrixInfo(sMInv, "sMInv");
    cout << "Mat Assembly for "<< "M1"<<endl;
    MatAssemblyBegin(M1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for "<< "sM1"<<endl;
    MatAssemblyBegin(sM1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sM1,MAT_FINAL_ASSEMBLY);
    //cout << "Mat Assembly for "<< "NMat"<<endl;
    //MatAssemblyBegin(NMat, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(NMat,MAT_FINAL_ASSEMBLY);
    //cout << "Mat Assembly for "<< "VecNMat"<<endl;
    //MatAssemblyBegin(VecNMat,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(VecNMat,MAT_FINAL_ASSEMBLY);
	
    //cout << "Mat Assembly for "<< "NMatT"<<endl;
    //MatAssemblyBegin(NMatT, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(NMatT,MAT_FINAL_ASSEMBLY);
    //cout << "Mat Assembly for "<< "VecNMatT"<<endl;
    //MatAssemblyBegin(VecNMatT,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(VecNMatT,MAT_FINAL_ASSEMBLY);	
		
    unsigned int posR=0;
    unsigned int posL=0;
    unsigned int eR=0;
    unsigned int eL=0;

    std::cout << " - : Face Integration started\n";
        //
        //
    FluxData fData(nb);
    //typedef void  (HEuler::*FaceIntegrand)(const FaceT*, const PointPhysicalT& normal , const PointReferenceOnTheFaceT&, FluxData&);
    //FaceIntegrand faceInteg = &HEuler::faceIntegrand;
    FaceIntegralT   faceIntegral(useCache);
    
    for (ConstFaceIterator citFe = faceColBegin(); citFe != faceColEnd(); ++citFe)
    {
        
        if ((*citFe)->getPtrElementRight()== NULL) // boundary face
        {
        }
        else
        {
            eR = (*citFe)->getPtrElementRight()->getID();
            eL = (*citFe)->getPtrElementLeft()->getID();
            
            posR  = eR*nb;
            
            posL  = eL*nb;
            
            faceIntegral.integrate((*citFe), this, fData);

            for (unsigned int j=0;j<nb;++j)
            {
                for (unsigned int i=0;i<nb;++i)
                {
                    MatSetValue(BF, posR+j, 		posL+i,     fData.right_[j](0,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posR+j, 		posR+i,    -fData.right_[j](3,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posL+j, 		posL+i,   	fData.left_[j](0,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posL+j, 		posR+i,    -fData.left_[j](3,i),  	ADD_VALUES);//U coefficient	

                    MatSetValue(globalDIV, posR+j, posL+i,    	  fData.right_[j](6,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posL+j, posL+i,   	 -fData.left_[j](6,i),  	ADD_VALUES);//U coefficient		
                    MatSetValue(globalDIV, posL+j, posR+i,   	 -fData.left_[j](9,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posR+j, posR+i,    	  fData.right_[j](9,i),  	ADD_VALUES);//U coefficient		
						
					if (config->DIM_ >= 2)
					{					
                    MatSetValue(BF, Nu+posR+j,      posL+i,     fData.right_[j](1,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+posR+j,      posR+i,    -fData.right_[j](4,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+posL+j,      posL+i,   	fData.left_[j](1,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+posL+j,      posR+i,    -fData.left_[j](4,i),  	ADD_VALUES);//V coefficient
				
                    MatSetValue(globalDIV, posL+j, Nu+posL+i,    -fData.left_[j](7,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posR+j, Nu+posL+i,     fData.right_[j](7,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posR+j, Nu+posR+i,     fData.right_[j](10,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posL+j, Nu+posR+i,    -fData.left_[j](10,i),  	ADD_VALUES);//V coefficient
					}
					
					if (config->DIM_ >= 3) // changed 2 to 3
					{						
                    MatSetValue(BF, Nu+Nv+posL+j, 	posR+i,    -fData.left_[j](5,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posL+j, 	posL+i,   	fData.left_[j](2,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posR+j, 	posR+i,    -fData.right_[j](5,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posR+j, 	posL+i,     fData.right_[j](2,i),  	ADD_VALUES);//W coefficient
					
                    MatSetValue(globalDIV, posL+j, Nu+Nv+posL+i, -fData.left_[j](8,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posR+j, Nu+Nv+posL+i,  fData.right_[j](8,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posL+j, Nu+Nv+posR+i, -fData.left_[j](11,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posR+j, Nu+Nv+posR+i,  fData.right_[j](11,i),  	ADD_VALUES);//W coefficient
					}
                }
            }
                //                cout <<"***********************************"<<endl;
        }
    }
	cout << "Mat Assembly for " << "globalDIV"<<endl;
    MatAssemblyBegin(globalDIV,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalDIV,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for " << "BF"<<endl;
    MatAssemblyBegin(BF,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(BF,MAT_FINAL_ASSEMBLY);
	
    printFullMatrixInfo(globalDIV, "globalDIV");
    printFullMatrixInfo(BF, "BF");	
	
    double fillBF=1;//(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
    //MatMatMult(MInv, BF, MAT_INITIAL_MATRIX, fillBF, &BF1);
	//MatMatMult(M1, BF1, MAT_INITIAL_MATRIX, fillBF, &BF1);
	if (config->stratified_ == 1)
	{	
    //MatMatMult(MInv, BF, MAT_INITIAL_MATRIX, fillBF, &BF1);
	//MatMatMult(M1, BF1, MAT_INITIAL_MATRIX, fillBF, &BF1);		
	
		MatMatMult(BF, sMInv, MAT_INITIAL_MATRIX, fillBF, &BF1);
		MatMatMult(BF1, sM1, MAT_INITIAL_MATRIX, fillBF, &BF1);
		MatMatMult(MInv, BF1, MAT_INITIAL_MATRIX, fillBF, &BF1);
	}
	else
	{
		MatMatMult(MInv, BF, MAT_INITIAL_MATRIX, fillBF, &BF1);	
	}

	double fillDIV=1;//(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
    //MatMatMult(sMInv, globalDIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
	//MatMatMult(sM1, DIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
	
	if (config->stratified_ == 1)
	{
    //MatMatMult(sMInv, globalDIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
	//MatMatMult(sM1, DIV, MAT_INITIAL_MATRIX, fillBF, &DIV);		
	
		MatMatMult(globalDIV, MInv, MAT_INITIAL_MATRIX, fillBF, &DIV);
		MatMatMult(DIV, M1, MAT_INITIAL_MATRIX, fillBF, &DIV);	
		MatMatMult(sMInv, DIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
	}
	else
	{		
		MatMatMult(sMInv, globalDIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
	}
	/*
	if (config->stratified_ == 1)
	{	// 1D Case:
		MatMatMult(sMInv, NMat, MAT_INITIAL_MATRIX, 1., &C);
		MatMatMult(sM1, C, MAT_INITIAL_MATRIX, 1., &C);		
		MatMatMult(sMInv, C, MAT_INITIAL_MATRIX, 1., &C);

		MatMatMult(sMInv, VecNMat, MAT_INITIAL_MATRIX, 1., &C1);
		MatMatMult(sM1, C1, MAT_INITIAL_MATRIX, 1., &C1);		
		MatMatMult(sMInv, C1, MAT_INITIAL_MATRIX, 1., &C1);	
	}
	else	
	{
		MatMatMult(sMInv, NMat, MAT_INITIAL_MATRIX, 1., &C);
		MatMatMult(sMInv, VecNMat, MAT_INITIAL_MATRIX, 1., &C1);
	}
	*/
	
		//MatMatMult(sMInv, NMatT, MAT_INITIAL_MATRIX, 1., &CT);
		//MatMatMult(sM1, CT, MAT_INITIAL_MATRIX, 1., &CT);		
		//MatMatMult(sMInv, CT, MAT_INITIAL_MATRIX, 1., &CT);

		//MatMatMult(sMInv, VecNMatT, MAT_INITIAL_MATRIX, 1., &C1T);
		//MatMatMult(sM1, C1T, MAT_INITIAL_MATRIX, 1., &C1T);		
		//MatMatMult(sMInv, C1T, MAT_INITIAL_MATRIX, 1., &C1T);			
		
		// 2D Case:
		//MatMatMult(NMat, sMInv, MAT_INITIAL_MATRIX, fillBF, &C);
		//MatMatMult(C, sM1, MAT_INITIAL_MATRIX, fillBF, &C);
		//MatMatMult(MInv, C, MAT_INITIAL_MATRIX, fillBF, &C);

		//MatMatMult(VecNMat, sMInv, MAT_INITIAL_MATRIX, fillBF, &C1);
		//MatMatMult(C1, sM1, MAT_INITIAL_MATRIX, fillBF, &C1);
		//MatMatMult(MInv, C1, MAT_INITIAL_MATRIX, fillBF, &C1);		
	
	//outputMatrix(sMInv, "sMInv.dat");
	//outputMatrix(sM1, "sM1.dat");
	
	//outputMatrix(NMat, "NMat.dat");
	//outputMatrix(NMatT, "NMatT.dat");	
	
	//outputMatrix(VecNMat, "VecNMat.dat");
	//outputMatrix(VecNMatT, "VecNMatT.dat");
	
	MatDestroy(&BF);
	//MatDestroy(&NMat);
	//MatDestroy(&VecNMat);
	MatDestroy(&globalDIV);
	//MatDestroy(&C1);
	//MatDestroy(&C);
	MatDestroy(&sMInv);
	MatDestroy(&MInv);
	MatDestroy(&sM1);
	MatDestroy(&M1);
	
    cout << "Mat Assembly for "<< "BF1"<<endl;
    MatAssemblyBegin(BF1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(BF1,MAT_FINAL_ASSEMBLY);
    
    cout << "Mat Assembly for "<< "DIV"<<endl;
    MatAssemblyBegin(DIV,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(DIV,MAT_FINAL_ASSEMBLY);
	
    //cout << "Mat Assembly for "<< "C"<<endl;
    //MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);

	//cout << "Mat Assembly for "<< "C1"<<endl;
    //MatAssemblyBegin(C1,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(C1,MAT_FINAL_ASSEMBLY);	
    
    printFullMatrixInfo(DIV, "DIV");
    printFullMatrixInfo(BF1, "BF1");	
	//printFullMatrixInfo(C, "C");
	//printFullMatrixInfo(C1, "C1");
	
    //cout << "Mat Assembly for "<< "CT"<<endl;
    //MatAssemblyBegin(CT,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(CT,MAT_FINAL_ASSEMBLY);

	//cout << "Mat Assembly for "<< "C1T"<<endl;
    //MatAssemblyBegin(C1T,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(C1T,MAT_FINAL_ASSEMBLY);	
	
	//outputMatrix(C, "C.dat");
	//outputMatrix(C1, "C1.dat");	
	
	//outputMatrix(C, "CT.dat");
	//outputMatrix(C1, "C1T.dat");

			
	
    double dummy=0;
    const PetscInt* cols;
    const PetscScalar* values;
    PetscInt 	numberOfNonZeros;
    std::cout << " - : Started Creating P matrix!\n";
    PetscInt* nnzP= new PetscInt[Nu+Nv+Nw+Nl+Nl];
    PetscInt* nnzQ= new PetscInt[Nu+Nv+Nw+Nl+Nl];
	cout << Nu << " " << Nv << " " << Nw << " " << Nl << " " << N << " " << N+Nl << " " << nb << endl;
    for (int i =0;i<Nl;++i)
    {
        nnzP[i]=(7*nb+1);
        nnzQ[i]=(7*nb+1);
		//if (Nu > 0)
		//{
			nnzP[Nu+i]=(7*nb+1);
			nnzQ[Nu+i]=(7*nb+1);
		//}
		//if (Nv > 0 )
		//{
			nnzP[Nu+Nv+i]=(7*nb+1+1-1);
			nnzQ[Nu+Nv+i]=(7*nb+1+1-1);
		//}
        nnzP[Nu+Nv+Nw+i]=(3*7*nb+1-1); 
        nnzQ[Nu+Nv+Nw+i]=(3*7*nb+1-1); 
        nnzP[Nu+Nv+Nw+Nl+i]=(3*7*nb+1-1); 
        nnzQ[Nu+Nv+Nw+Nl+i]=(3*7*nb+1-1); 
    }

    std::cout << " stage1"<<endl;
    
    MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl+Nl, Nu+Nv+Nw+Nl+Nl, PETSC_NULL, nnzP,&P_);
    std::cout << " stage2"<<endl;

    
    MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl+Nl, Nu+Nv+Nw+Nl+Nl, PETSC_NULL, nnzQ,&Q_);
    delete [] nnzP;
    delete [] nnzQ;
	
	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	N+Nl+Nl, N+Nl+Nl,       (15)*nb+1, 			PETSC_NULL, &P_);
    //MatCreateSeqAIJ(PETSC_COMM_SELF, 	N+Nl+Nl, N+Nl+Nl,       (15)*nb+1, 			PETSC_NULL, &Q_);
    
    std::cout << " stage3"<<endl;
    cout << endl << "N = " << Nu+Nv+Nw << endl;
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nl+Nl,         (7)*nb, 			PETSC_NULL, &B_);//number of possible nonzero blocks are 7: element and his 6 neighbours
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl+Nl, 	Nu+Nv+Nw,    (7*3)*nb, 				PETSC_NULL, &C_);
	//Mat& B = static_cast<HEulerGlobalVariables*>(globalData_)->B_;
	//Mat& C = static_cast<HEulerGlobalVariables*>(globalData_)->C_;	
	
    for (unsigned int i =0;i<N+Nl+Nl;++i)
    {
        MatSetValue(P_, i,     i,     1.0, ADD_VALUES); //creating 2*I matrix in the upper part
        
        MatSetValue(Q_, i,     i,     1.0, ADD_VALUES); //creating 2*I matrix in the upper part

		// replacing C and C1:
		if (i<Nl)	
		{
			MatSetValue(P_, 	Nu+Nv+i,     N+i,     -0.5*dt*1.0/N2, 	ADD_VALUES);
			MatSetValue(Q_, 	Nu+Nv+i,     N+i,     0.5*dt*1.0/N2, 	ADD_VALUES);
			
			MatSetValue(B_, 	Nu+Nv+i, 		i, -0.5*dt*1.0/N2, 	ADD_VALUES);
			
			MatSetValue(P_, 	Nu+Nv+i,     N+Nl+i,     0.5*dt*(1.0/N2+1.0), 	ADD_VALUES);
			MatSetValue(Q_, 	Nu+Nv+i,     N+Nl+i,     -0.5*dt*(1.0/N2+1.0), 	ADD_VALUES);
			
			MatSetValue(B_, 	Nu+Nv+i, 		Nl+i, 0.5*dt*(1.0/N2+1.0), 	ADD_VALUES);
			
			MatSetValue(P_, 	N+Nl+i,    Nu+Nv+i,     -0.5*dt, 	ADD_VALUES);
			MatSetValue(Q_, 	N+Nl+i,    Nu+Nv+i,     0.5*dt, 	ADD_VALUES);	
			
			MatSetValue(C_, 	Nl+i, 		Nu+Nv+i, -0.5*dt*1.0/N2, 	ADD_VALUES);
			
			MatSetValue(P_, 	Nu+Nv+i,     N+i,     -0.5*dt*1.0/N2*-3.0, 	ADD_VALUES);
			MatSetValue(Q_, 	Nu+Nv+i,     N+i,     0.5*dt*1.0/N2*-3.0, 	ADD_VALUES);
			
			MatSetValue(B_, 	Nu+Nv+i, 		i, -0.5*dt*1.0/N2*-3.0, 	ADD_VALUES);
			
			MatSetValue(P_, 	Nu+Nv+i,     N+Nl+i,     0.5*dt*(1.0/N2)*-3.0, 	ADD_VALUES);
			MatSetValue(Q_, 	Nu+Nv+i,     N+Nl+i,     -0.5*dt*(1.0/N2)*-3.0, 	ADD_VALUES);
			
			MatSetValue(B_, 	Nu+Nv+i, 		Nl+i, 0.5*dt*(1.0/N2)*-3.0, 	ADD_VALUES);
			
			MatSetValue(P_, 	N+i,    Nu+Nv+i,     0.5*dt*-3.0, 	ADD_VALUES);
			MatSetValue(Q_, 	N+i,    Nu+Nv+i,     -0.5*dt*-3.0, 	ADD_VALUES);	

			MatSetValue(C_, 	i, 		Nu+Nv+i, -0.5*dt*1.0/N2, 	ADD_VALUES);		
						
		}	

    }
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    
    cout <<endl<< "N = " << N<<endl;
	
    for (unsigned int i=0;i<N;++i)
    {
        MatGetRow(BF1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {

				//MatSetValue(P_, 	i,     N+cols[j],     -0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
				//MatSetValue(Q_, 	i,     N+cols[j],     0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
				
				//MatSetValue(P_, 	i,     N+cols[j],     0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
				//MatSetValue(Q_, 	i,     N+cols[j],     -0.5*dt*1.0/N2*dummy, 	ADD_VALUES);

				//MatSetValue(P_, 	i,     N+Nl+cols[j],     0.5*dt*(1.0/N2)*dummy, 	ADD_VALUES);
                //MatSetValue(Q_, 	i,     N+Nl+cols[j],     -0.5*dt*(1.0/N2)*dummy, 		ADD_VALUES);	
				
				//MatSetValue(P_, 	i,     N+Nl+cols[j],     -0.5*dt*(1.0/N2+1.0)*dummy, 	ADD_VALUES);
                //MatSetValue(Q_, 	i,     N+Nl+cols[j],     0.5*dt*(1.0/N2+1.0)*dummy, 		ADD_VALUES);	
				
				MatSetValue(P_, 	i,     N+Nl+cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
                MatSetValue(Q_, 	i,     N+Nl+cols[j],     0.5*dt*dummy, 		ADD_VALUES);	

				MatSetValue(B_, 	i, 		Nl+cols[j], -0.5*dt*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);

        if (i<Nl)
        {
			// C and C1 now added directly: NMat and VecNMat no longer needed
			/*
			MatGetRow(C, i, &numberOfNonZeros, &cols, &values);
			for (int j=0;j<numberOfNonZeros;++j)
			{
				dummy = (values[j]);
				
				if (dummy!=0.0)
				{
					MatSetValue(P_, 	Nu+Nv+i,     N+cols[j],     -0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	Nu+Nv+i,     N+cols[j],     0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
					
					MatSetValue(P_, 	Nu+Nv+i,     N+Nl+cols[j],     0.5*dt*(1.0/N2+1.0)*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	Nu+Nv+i,     N+Nl+cols[j],     -0.5*dt*(1.0/N2+1.0)*dummy, 	ADD_VALUES);
					
					MatSetValue(P_, 	N+Nl+i,    Nu+Nv+cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	N+Nl+i,    Nu+Nv+cols[j],     0.5*dt*dummy, 	ADD_VALUES);	
				}
			}
			MatRestoreRow(C, i, &numberOfNonZeros, &cols, &values);		

			MatGetRow(C1, i, &numberOfNonZeros, &cols, &values);
			for (int j=0;j<numberOfNonZeros;++j)
			{
				cout << dummy << endl;
				dummy = (values[j]);
				if (dummy!=0.0)
				{
					
					MatSetValue(P_, 	Nu+Nv+i,     N+cols[j],     -0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	Nu+Nv+i,     N+cols[j],     0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
					
					MatSetValue(P_, 	Nu+Nv+i,     N+Nl+cols[j],     0.5*dt*(1.0/N2)*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	Nu+Nv+i,     N+Nl+cols[j],     -0.5*dt*(1.0/N2)*dummy, 	ADD_VALUES);
					
					MatSetValue(P_, 	N+i,    Nu+Nv+cols[j],     0.5*dt*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	N+i,    Nu+Nv+cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
				}
			}
			MatRestoreRow(C1, i, &numberOfNonZeros, &cols, &values);	
			*/ 
            MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
                dummy = (values[j]);
                if (dummy!=0.0)
                {
                    MatSetValue(P_, 	N+i,     cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
                    MatSetValue(Q_, 	N+i,     cols[j],     0.5*dt*dummy, 	ADD_VALUES);
                    MatSetValue(P_, 	N+Nl+i,     cols[j],     -0.5*dt*dummy, ADD_VALUES);
                    MatSetValue(Q_, 	N+Nl+i,     cols[j],     0.5*dt*dummy, 	ADD_VALUES);					
					
					MatSetValue(C_, 	i, 		cols[j], -0.5*dt*dummy, 	ADD_VALUES);
					MatSetValue(C_, 	Nl+i, 		cols[j], -0.5*dt*dummy, 	ADD_VALUES);
                }
            }
            MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);
			
        }
    }
    cout << "Mat Assembly for "<< "P_"<<endl;
    MatAssemblyBegin(P_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(P_,MAT_FINAL_ASSEMBLY);
    
    cout << "Mat Assembly for "<< "Q_"<<endl;
    MatAssemblyBegin(Q_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q_,MAT_FINAL_ASSEMBLY);
    
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
	
    cout << "Mat Assembly for "<< "B_"<<endl;
    MatAssemblyBegin(B_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B_,MAT_FINAL_ASSEMBLY);	
	printFullMatrixInfo(B_, "B_");
    cout << "Mat Assembly for "<< "C_"<<endl;
    MatAssemblyBegin(C_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(C_,MAT_FINAL_ASSEMBLY);	
	printFullMatrixInfo(C_, "C_");	
    
    
	outputMatrix(P_, "P.dat");
	outputMatrix(Q_, "Q.dat");
	
	MatDestroy(&BF1);
	MatDestroy(&DIV);
	
}

void
HEuler::createIncompressibleSystem()
{
 	PetscInfoAllow(PETSC_TRUE, "history.txt");
	
	Mat C;
	Mat BF;
    
	Mat C1;
	Mat BF1;
	Mat MInv;
	Mat sMInv;
	//Mat globalDIV;
  	//Mat Ml;
	Mat A;
	Mat Ah;	
	
    Mat& globalDIV = static_cast<HEulerGlobalVariables*>(globalData_)->DivergenceFreeMatrix_;
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
	
    unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nu     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
	unsigned int Nv		= 0;
	if (config->DIM_ == 3)
	{
         Nv= Nu;
	}	
    unsigned int Nw     = Nu;
    unsigned int Nl     = Nu;
    unsigned int N      = Nu+Nv+Nw;
    
	double N2			= config->N2_;
	double gamma		= config->gamma_;
    double dt           = static_cast<HEulerGlobalVariables*>(globalData_)->dt_;
    cout << "dt=" <<dt<<endl;
    
	Vec 												UInit, RHS, UCorrected;
	PetscViewer 										viewer;
	
	MatInfo info;

	// to pass initial condition to correctInitialProjectionOfVelocity
	double u, v, w, rho; 	

	cout << "**************Starting create Matrices*************"<<endl;
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    1*nb, 				PETSC_NULL, &MInv);
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &sMInv);
	Mat M1;
	Mat sM1;
	Mat VecNMat;
	Mat NMat;
	Mat Identity3;
 	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    1*nb, 				PETSC_NULL, &M1);
	
	Mat M2;
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    1*nb, 				PETSC_NULL, &M2);
	Mat sM1GradRho0;
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu, 	Nu,    1*nb, 				PETSC_NULL, &sM1GradRho0);
	Mat M1GradRho0;
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    1*nb, 				PETSC_NULL, &M1GradRho0);	
	
	Mat sM2;
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &sM2);	
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &sM1);	
	Mat sEBM1GradRho0;
	Mat sEBM2;
	if (config->EB_==1){
		MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &sEBM2);
		MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu, 	Nu,    1*nb, 				PETSC_NULL, &sEBM1GradRho0);	
	}
	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nw, 	Nl,          1*nb, 				PETSC_NULL, &VecNMat);
	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nw, 	Nl,          1*nb, 				PETSC_NULL, &NMat);		
 	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &Ml);
    
	    //MatCreateSeqAIJ(PETSC_COMM_SELF, 	N, 		N,    	1*nb, 			PETSC_NULL, &M);
        //MatCreateSeqAIJ(PETSC_COMM_SELF, 	N1, 	N1,    	(3*7)*nb, 			PETSC_NULL, &A);
        //MatCreateSeqAIJ(PETSC_COMM_SELF, 	N1, 	3*N,    (3*7+7)*nb, 			PETSC_NULL, &Ah);
	
	
	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    12*nb, 			PETSC_NULL, &C);
	
    MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nl,         nb, 			PETSC_NULL, &Identity3);
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nl,         (7)*nb, 			PETSC_NULL, &BF);//number of possible nonzero blocks are 7: element and his 6 neighbours
	
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 	Nu+Nv+Nw,    (3*7)*nb, 			PETSC_NULL, &globalDIV);//number of possible nonzero blocks are 7: element and his 6 neighbours
    
	VecCreateSeq(PETSC_COMM_SELF, Nu+Nv+Nw, &UInit);
	VecCreateSeq(PETSC_COMM_SELF, Nu+Nv+Nw, &UCorrected);
	VecCreateSeq(PETSC_COMM_SELF, Nl, &RHS);
//	VecCreateSeq(PETSC_COMM_SELF, Nl, &gv.LambdaConstraint_);
	
	cout << "**************Mass Matrix Calculated*************"<<endl;

//	timer t0;
//	std::cout << t0.elapsed() << " - : Element Integration started\n";
    ElementIntegralData gradMass;
    bool useCache = false;
    ElementIntegralT   elIntegral(useCache);
    typedef void  (HEuler::*Integrand)(const ElementT* , const PointReferenceT&, ElementIntegralData&);
    Integrand gradMassInteg = &HEuler::elementIntegrand;
    
    std::cout << " - : Element Integration started\n";
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {

        int m  = (*cit)->getID();
        unsigned int pos = (*cit)->getID()*nb;
        
        HEulerElementData* elementData = static_cast<HEulerElementData*>((*cit)->getUserData());
        
        ElementT* element =(*cit);
        
        elIntegral.integrate<ElementIntegralData>(element, this, gradMass);
		
		for (unsigned int i = 0; i<nb;++i)
 		{
			u = element->getData(0,0,i);
			VecSetValue(UInit, pos+i,     u, INSERT_VALUES);
			if (config->DIM_ == 3)
			{
				v = element->getData(0,1,i);
				VecSetValue(UInit, Nu+pos+i,     v, INSERT_VALUES);
			}			
            w = element->getData(0,config->DIM_-1,i);
			//cout << u << " " << config->DIM_-1 << " " <<w << endl;
			VecSetValue(UInit, Nu+Nv+pos+i, w, INSERT_VALUES);
			
			if (config->WA_ == 0)
			{
				//MatSetValue(Identity3, Nu+Nv+pos+i, pos+i, 1.0, ADD_VALUES);
			}
			else
			{
				MatSetValue(Identity3, pos+i, pos+i, sin(gamma), ADD_VALUES);
				MatSetValue(Identity3, Nu+Nv+pos+i, pos+i, cos(gamma), ADD_VALUES);
			}
				
		}
        
		for (unsigned int j = 0; j<nb;++j)
		{
			for (unsigned int i=0; i<nb;++i)//can be optimized later!!!
 			{
                    //integrateOverElement(*itM, scalarProduct<dim>(ExpansionServer::basisFunction(i),ExpansionServer::basisFunction(j)), 					massVij=calculatePhiIDotPhiJForUWithV(i,j, elementM, itM);
                ///massij=elementData->massMatrix_(j,i);
                //MatSetValue(C, pos+j,      Nu+pos+i,            massij*(omega3),  ADD_VALUES);//V coefficient
				///MatSetValue(C, Nu+Nv+pos+j,      Nu+Nv+pos+i,            massij*0.5*dt,  ADD_VALUES);
//                MatSetValue(C, pos+j,       Nu+Nv+pos+i,   		massij*(-omega2), ADD_VALUES);   //W coefficient  omega2=0 done for memory shortage
                
                //MatSetValue(C, Nu+pos+j,   pos+i,               massij*(-omega3), ADD_VALUES);//U coefficient
//                MatSetValue(C, Nu+pos+j,    Nu+Nv+pos+i,        massij*(omega1),  ADD_VALUES);//W coefficient  omega1=0 done for memory shortage
//                MatSetValue(C, Nu+Nv+pos+j, pos+i,              massij*(omega2),  ADD_VALUES);//U coefficient    omega2=0 done for memory shortage
//                MatSetValue(C, Nu+Nv+pos+j, Nu+pos+i,           massij*(-omega1), ADD_VALUES);   //V coefficient omega1=0 done for memory shortage
                
                ///MatSetValue(Ml,    pos+j, 		pos+i,     massij , ADD_VALUES);
				MatSetValue(M1,   		pos+j, 			pos+i,     	elementData->massMatrixOverBackgroundDensity_(j,i), 		ADD_VALUES);
				
				MatSetValue(sM1,   		pos+j, 			pos+i,      elementData->massMatrixOverBackgroundDensity_(j,i), 		ADD_VALUES);
				MatSetValue(M1,   		Nu+pos+j, 		Nu+pos+i,   elementData->massMatrixOverBackgroundDensity_(j,i), 		ADD_VALUES);

                MatSetValue(MInv, pos+j,       pos+i,    	elementData->invMassMatrix_(j,i) , ADD_VALUES);
                MatSetValue(MInv, Nu+Nv+pos+j, Nu+Nv+pos+i, elementData->invMassMatrix_(j,i) , ADD_VALUES);
                
                MatSetValue(sMInv, pos+j, pos+i,           elementData->invMassMatrix_(j,i) , ADD_VALUES);

                MatSetValue(BF, pos+j, 		 pos+i,   	-gradMass.xGrad_(i,j),  ADD_VALUES);//U coefficient
				MatSetValue(globalDIV, pos+j, pos+i,   	    gradMass.xGrad_(j,i),  ADD_VALUES);//U coefficient

				MatSetValue(BF, Nu+pos+j, 	 	pos+i,   	-gradMass.yGrad_(i,j),  ADD_VALUES);
				MatSetValue(globalDIV, pos+j, Nu+pos+i,   	gradMass.yGrad_(j,i),  ADD_VALUES);

				if (config->EB_==0)
				{
					MatSetValue(sM2,   		pos+j, 			pos+i,      elementData->massMatrixOverBackgroundDensityDeriv_(j,i), 		ADD_VALUES);
					MatSetValue(M2,   		pos+j, 			pos+i,      elementData->massMatrixOverBackgroundDensityDeriv_(j,i), 		ADD_VALUES);
					MatSetValue(M2,   		Nu+pos+j, 		Nu+pos+i,      elementData->massMatrixOverBackgroundDensityDeriv_(j,i), 	ADD_VALUES);
				}				
				MatSetValue(sM1GradRho0,   		pos+j, 			pos+i,     	elementData->massMatrixTimesBackgroundDensityDeriv_(j,i), 		ADD_VALUES);
				
				if (config->EB_==1){
					MatSetValue(sEBM1GradRho0,   		pos+j, 			pos+i,     	elementData->massMatrixTimesN2_(j,i), 		ADD_VALUES);
					MatSetValue(sEBM2,   		pos+j, 			pos+i,     	elementData->massMatrixOverN2_(j,i), 		ADD_VALUES);
				}		
				

				MatSetValue(M1GradRho0,   		Nu+pos+j, 			Nu+pos+i,     	elementData->massMatrixTimesBackgroundDensityDeriv_(j,i), 		ADD_VALUES);				
				//MatSetValue(M1GradRho0,   		pos+j, 			pos+i,     	elementData->massMatrixTimesBackgroundDensityDeriv_(j,i), 		ADD_VALUES);
				//MatSetValue(M1GradRho0,   		pos+j, 			pos+i,     	elementData->massMatrixTimesBackgroundDensityDeriv_(j,i), 		ADD_VALUES);
				
				if (config->DIM_ == 3)
				{
					MatSetValue(MInv, Nu+pos+j,    Nu+pos+i,    elementData->invMassMatrix_(j,i) , ADD_VALUES);	

					MatSetValue(M1, Nu+Nv+pos+j, 	Nu+Nv+pos+i,    elementData->massMatrixOverBackgroundDensity_(j,i), ADD_VALUES);
					MatSetValue(BF, Nu+Nv+pos+j, pos+i,   	-gradMass.zGrad_(i,j),  ADD_VALUES);//W coefficient
					MatSetValue(globalDIV, pos+j, Nu+Nv+pos+i,	gradMass.zGrad_(j,i),  ADD_VALUES);//W coefficient
				}			
                                                                                          
 			}
		}
    }
	

        //        outputMatrix(globalDIV, "globalDIV.dat");
        //        outputMatrix(BF, "BF.dat");
    
    printFullMatrixInfo(globalDIV, "globalDIV");
    printFullMatrixInfo(BF, "BF");
    
	cout << "Mat Assembly for "<< "MINV"<<endl;
    MatAssemblyBegin(MInv,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MInv,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for "<< "MINVSm"<<endl;
    MatAssemblyBegin(sMInv,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sMInv,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for "<< "M1"<<endl;
    MatAssemblyBegin(M1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for "<< "sM1"<<endl;
    MatAssemblyBegin(sM1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sM1,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for "<< "sM2"<<endl;
    MatAssemblyBegin(sM2,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sM2,MAT_FINAL_ASSEMBLY);	
    cout << "Mat Assembly for "<< "M2"<<endl;
    MatAssemblyBegin(M2,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M2,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for "<< "M1GradRho0"<<endl;
    MatAssemblyBegin(M1GradRho0,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1GradRho0,MAT_FINAL_ASSEMBLY);	
    //cout << "Mat Assembly for "<< "Identity3"<<endl;
    //MatAssemblyBegin(Identity3,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(Identity3,MAT_FINAL_ASSEMBLY);	
    
	cout << "Vec Assembly for "<< "UInit"<<endl;
	VecAssemblyBegin(UInit);
    VecAssemblyEnd(UInit);	
	
    cout << "Mat Assembly for "<< "sM1GradRho0"<<endl;
    MatAssemblyBegin(sM1GradRho0,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sM1GradRho0,MAT_FINAL_ASSEMBLY);	
	Mat EBCWinR;
	if (config->EB_==1)
	{
    cout << "Mat Assembly for "<< "sEBM1GradRho0"<<endl;
    MatAssemblyBegin(sEBM1GradRho0,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sEBM1GradRho0,MAT_FINAL_ASSEMBLY);
	
	MatMatMult(sMInv,sEBM1GradRho0,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&EBCWinR);
    cout << "Mat Assembly for "<< "EBCWinR"<<endl;
    MatAssemblyBegin(EBCWinR,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EBCWinR,MAT_FINAL_ASSEMBLY);
	printFullMatrixInfo(EBCWinR, "EBCWinR");
	outputMatrix(EBCWinR, "EBCWinR.dat");	
	}	
	
	
    //cout << "Mat Assembly for "<< "M1GradRho0"<<endl;
    //MatAssemblyBegin(M1GradRho0,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(M1GradRho0,MAT_FINAL_ASSEMBLY);		
	
	Mat CWinR; // Matrix multiplying the W coefficients in the equations for R ( = approx N^2)
	MatMatMult(sM1GradRho0,sMInv,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CWinR);
	if (config->stratified_==1 && config->EB_ == 0)
	{
	MatMatMult(CWinR,sM1,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CWinR);
	MatMatMult(sMInv,CWinR,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CWinR);
	}
    cout << "Mat Assembly for "<< "CWinR"<<endl;
    MatAssemblyBegin(CWinR,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(CWinR,MAT_FINAL_ASSEMBLY);
	printFullMatrixInfo(CWinR, "CWinR");
	outputMatrix(CWinR, "CWinR.dat");
	
	Mat sCRinW; // Matrix multiplying the R coefficients in the equations for W ( = approx 1)
	if (config->EB_ == 0)
	{
	MatMatMult(sM1GradRho0,sMInv,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&sCRinW);	
	MatMatMult(sCRinW,sM2,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&sCRinW);
	MatMatMult(sMInv,sCRinW,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&sCRinW);
	}
	else 
	{
    cout << "Mat Assembly for "<< "sEBM2"<<endl;
    MatAssemblyBegin(sEBM2,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sEBM2,MAT_FINAL_ASSEMBLY);			
	MatMatMult(sEBM1GradRho0,sMInv,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&sCRinW);	
	MatMatMult(sCRinW,sEBM2,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&sCRinW);
	MatMatMult(sMInv,sCRinW,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&sCRinW);
	}

    cout << "Mat Assembly for "<< "sCRinW"<<endl;
    MatAssemblyBegin(sCRinW,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sCRinW,MAT_FINAL_ASSEMBLY);
	printFullMatrixInfo(sCRinW, "sCRinW");
	outputMatrix(sCRinW, "sCRinW.dat");	
	
	
	
	/*
	Mat CRinW; // Matrix multiplying the R coefficients in the equations for W ( = approx 1) (in the incompressibility condition)
	MatMatMult(M1GradRho0,MInv,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CRinW);	
	MatMatMult(CRinW,M2,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CRinW);
	MatMatMult(MInv,CRinW,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CRinW);
    cout << "Mat Assembly for "<< "CRinW"<<endl;
    MatAssemblyBegin(CRinW,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(CRinW,MAT_FINAL_ASSEMBLY);
	printFullMatrixInfo(CRinW, "CRinW");
	outputMatrix(CRinW, "CRinW.dat");	
	 * */

 	double dummy=0;
	const PetscInt* cols;
	const PetscScalar* values;
	PetscInt 	numberOfNonZeros;
	for (unsigned int i = 0; i<Nu;++i)
	{
		//if (config->EB_ == 1)
		//{
			//MatSetValue(Identity3, 	Nu+i,     i,     1.0, 	ADD_VALUES);		
		//}
		//else
		{		
        MatGetRow(sCRinW, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
				MatSetValue(Identity3, 	Nu+Nv+i,     cols[j],     dummy, 	ADD_VALUES); // No 3D?
            }
        }
		}
        MatRestoreRow(sCRinW, i, &numberOfNonZeros, &cols, &values);	
	}
    cout << "Mat Assembly for "<< "Identity3"<<endl;
    MatAssemblyBegin(Identity3,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Identity3,MAT_FINAL_ASSEMBLY);		
	/*
	Mat CCC;
	MatMatMult(M1GradRho0,MInv,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CCC);	
	MatMatMult(CCC,M1,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CCC);
	MatMatMult(CCC,MInv,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CCC);
	

    cout << "Mat Assembly for "<< "CCC"<<endl;
    MatAssemblyBegin(CCC,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(CCC,MAT_FINAL_ASSEMBLY);
	printFullMatrixInfo(CCC, "CCC");
	outputMatrix(CCC, "CCC.dat");
	*/
	
	
    unsigned int posR=0;
    unsigned int posL=0;
    unsigned int eR=0;
    unsigned int eL=0;
    
    std::cout << " - : Face Integration started\n";
        //
        //
    FluxData fData(nb);
    //typedef void  (HEuler::*FaceIntegrand)(const FaceT*, const PointPhysicalT& normal , const PointReferenceOnTheFaceT&, FluxData&);
    //FaceIntegrand faceInteg = &HEuler::faceIntegrand;
    FaceIntegralT   faceIntegral(useCache);
    

	for (ConstFaceIterator citFe = faceColBegin(); citFe != faceColEnd(); ++citFe)
    {
        
        if ((*citFe)->getPtrElementRight()== NULL) // boundary face
        {
        }
        else
        {
            eR = (*citFe)->getPtrElementRight()->getID();
            eL = (*citFe)->getPtrElementLeft()->getID();
            
            posR  = eR*nb;
            
            posL  = eL*nb;
            
            faceIntegral.integrate((*citFe), this, fData);
            
            for (unsigned int j=0;j<nb;++j)
            {
                for (unsigned int i=0;i<nb;++i)
                {
                    MatSetValue(BF, posR+j, 		posL+i,     fData.right_[j](0,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posR+j, 		posR+i,    -fData.right_[j](3,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posL+j, 		posL+i,   	fData.left_[j](0,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posL+j, 		posR+i,    -fData.left_[j](3,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posR+j, posL+i,    	  fData.right_[j](6,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posL+j, posL+i,   	 -fData.left_[j](6,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posL+j, posR+i,   	 -fData.left_[j](9,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posR+j, posR+i,    	  fData.right_[j](9,i),  	ADD_VALUES);//U coefficient
					
					MatSetValue(BF, Nu+posR+j,      posL+i,     fData.right_[j](1,i),  	ADD_VALUES);//V coefficient
					MatSetValue(BF, Nu+posR+j,      posR+i,    -fData.right_[j](4,i),  	ADD_VALUES);//V coefficient	
					MatSetValue(BF, Nu+posL+j,      posL+i,   	fData.left_[j](1,i),  	ADD_VALUES);//V coefficient
					MatSetValue(BF, Nu+posL+j,      posR+i,    -fData.left_[j](4,i),  	ADD_VALUES);//V coefficient	
					MatSetValue(globalDIV, posL+j, Nu+posL+i,    -fData.left_[j](7,i),  	ADD_VALUES);//V coefficient
					MatSetValue(globalDIV, posR+j, Nu+posL+i,     fData.right_[j](7,i),  	ADD_VALUES);//V coefficient
					MatSetValue(globalDIV, posR+j, Nu+posR+i,     fData.right_[j](10,i),  	ADD_VALUES);//V coefficient
					MatSetValue(globalDIV, posL+j, Nu+posR+i,    -fData.left_[j](10,i),  	ADD_VALUES);//V coefficient	
						
					if (config->DIM_ == 3)
					{
						MatSetValue(BF, Nu+Nv+posL+j, 	posR+i,    -fData.left_[j](5,i),  	ADD_VALUES);//W coefficient
						MatSetValue(BF, Nu+Nv+posL+j, 	posL+i,   	fData.left_[j](2,i),  	ADD_VALUES);//W coefficient
						MatSetValue(BF, Nu+Nv+posR+j, 	posR+i,    -fData.right_[j](5,i),  	ADD_VALUES);//W coefficient
						MatSetValue(BF, Nu+Nv+posR+j, 	posL+i,     fData.right_[j](2,i),  	ADD_VALUES);//W coefficient
						MatSetValue(globalDIV, posL+j, Nu+Nv+posL+i, -fData.left_[j](8,i),  	ADD_VALUES);//W coefficient
						MatSetValue(globalDIV, posR+j, Nu+Nv+posL+i,  fData.right_[j](8,i),  	ADD_VALUES);//W coefficient
						MatSetValue(globalDIV, posL+j, Nu+Nv+posR+i, -fData.left_[j](11,i),  	ADD_VALUES);//W coefficient
						MatSetValue(globalDIV, posR+j, Nu+Nv+posR+i,  fData.right_[j](11,i),  	ADD_VALUES);//W coefficient		
					}					
				}
                //                cout <<"***********************************"<<endl;
			}	
		}
	}
    
	cout << "Mat Assembly for " << "globalDIV"<<endl;
    MatAssemblyBegin(globalDIV,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalDIV,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for " << "BF"<<endl;
    MatAssemblyBegin(BF,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(BF,MAT_FINAL_ASSEMBLY);
	
    printFullMatrixInfo(globalDIV, "globalDIV");
    printFullMatrixInfo(BF, "BF");	

	
	double fillBF=1; //(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
	//if (config->stratified_ == 1)
	if (config->stratified_==1 && config->EB_ == 0)
	{	
		MatMatMult(BF, sMInv, MAT_INITIAL_MATRIX, fillBF, &BF1);
		MatMatMult(BF1, sM1, MAT_INITIAL_MATRIX, fillBF, &BF1);
		MatMatMult(MInv, BF1, MAT_INITIAL_MATRIX, fillBF, &BF1);
	}
	else
	{
		MatMatMult(MInv, BF, MAT_INITIAL_MATRIX, fillBF, &BF1);	
	}

	double fillDIV=1; //(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
	//if (config->stratified_ == 1)
	if (config->stratified_==1 && config->EB_ == 0)
	{
		MatMatMult(globalDIV, MInv, MAT_INITIAL_MATRIX, fillDIV, &globalDIV);
		MatMatMult(globalDIV, M1, MAT_INITIAL_MATRIX, fillDIV, &globalDIV);	
		// Not needed since we scale both sides of the equation (= constructed incompressibility)
		//MatMatMult(sMInv, DIV, MAT_INITIAL_MATRIX, fillDIV, &DIV);
		//MatMatMult(sMInv, globalDIV, MAT_INITIAL_MATRIX, fillDIV, &globalDIV);
	}
	else
	{		
		// Not needed since we scale both sides of the equation (= constructed incompressibility)
		//MatMatMult(sMInv, globalDIV, MAT_INITIAL_MATRIX, fillDIV, &globalDIV);
	}

	MatDestroy(&BF);
	MatDestroy(&NMat);
	MatDestroy(&VecNMat);
	//MatDestroy(&globalDIV);
	
    cout << "Mat Assembly for "<< "BF1"<<endl;
    MatAssemblyBegin(BF1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(BF1,MAT_FINAL_ASSEMBLY);
    
    cout << "Mat Assembly for "<< "DIV"<<endl;
    MatAssemblyBegin(globalDIV,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalDIV,MAT_FINAL_ASSEMBLY);

    printFullMatrixInfo(globalDIV, "globalDIV");
    printFullMatrixInfo(BF1, "BF1");	
	
	//MatMatMult(globalDIV, BF1, MAT_INITIAL_MATRIX, 1, &A);
	MatMatMult(globalDIV, BF1, MAT_INITIAL_MATRIX, 1, &A);
	std::cout << " - : Create rotational matrix3!\n";
	//MatMatMult(globalDIV, C1, MAT_INITIAL_MATRIX, 1, &Ah);
	//MatMatMult(globalDIV, Identity3, MAT_INITIAL_MATRIX, 1, &Ah);
	MatMatMult(globalDIV, Identity3, MAT_INITIAL_MATRIX, 1, &Ah);
	
	//MatMatMult(globalDIV, CRinW, MAT_INITIAL_MATRIX, 1, &Ah);
	std::cout << " - : Finished rotational matrix4!\n";
	
	MatDestroy(&Identity3);
	printFullMatrixInfo(Ah, "Ah");
	printFullMatrixInfo(A, "A");
	


	cout << "Mat Assembly for "<< "A"<<endl;
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	cout << "Mat Assembly for "<< "Ah"<<endl;
	MatAssemblyBegin(Ah,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Ah,MAT_FINAL_ASSEMBLY);
	
	printFullMatrixInfo(Ah, "Ah");
	printFullMatrixInfo(A, "A");

	outputMatrix(Ah, "Ah.dat");
	
    //correctInitialProjectionOfVelocity(UInit, UCorrected);
        //UCorrected=UInit;
	cout << "Correction is done!"<<endl;
	
	cout << "Correction is done!"<<endl;
	MatMult(globalDIV, UInit, RHS);

    PetscInt a(10);
    PetscReal max;
    VecMax(RHS, &a, &max);
    cout << "Discrete divergence OLD = " << std::setprecision (20) <<max<<endl;
    MatMult(globalDIV, UCorrected, RHS);
    VecMax(RHS, &a, &max);
	cout << "Discrete divergence NEW = " << std::setprecision (20) <<max<<endl;


	std::cout<< " - : Started Creating P matrix!\n";
	/*
	PetscInt* nnzP= new PetscInt[Nu+Nv+Nw+Nl];
	PetscInt* nnzQ= new PetscInt[Nu+Nv+Nw+Nl];
	PetscInt* nnzR= new PetscInt[Nu+Nv+Nw+Nl];
   
    for (int i =0;i<Nu;++i)
    {
		// Needs to be optimized for 2D/3D 
		// polynomial orders = 0, 1, 2, 3
		// 3 gives allocation errors, 4 also
		// nnzQ does not allocate enough
		// use Matlab script
        nnzP[i]=(3*nb+1+10);
        nnzQ[i]=(3*nb+1);
        nnzP[Nu+i]=(3*nb+1+10);
        nnzQ[Nu+i]=(3*nb+1);
        nnzP[Nu+Nv+i]=(3*nb+1+10);
        nnzQ[Nu+Nv+i]=(3*nb+1);
        nnzP[Nu+Nv+Nw+i]=(21*nb+1);
        nnzQ[Nu+Nv+Nw+i]=(21*nb+1);
        nnzR[i]=(3*nb+1);
        nnzR[Nu+i]=(3*nb+1);
        nnzR[Nu+Nv+i]=(3*nb+1);
        nnzR[Nu+Nv+Nw+i]=(21*nb+1);		
    }
	 * */
//        // 	nnzP[Nu+Nv+Nw+Nl]=Nl;
//        // 	nnzQ[Nu+Nv+Nw+Nl]=Nl;
//	
	/*if (config->solutionType_ == HEulerConfigurationData::IGW2D) // 
	{
		cout << "yeah" << endl;
		PetscInt* nnzP= new PetscInt[Nu+Nv+Nw+Nl+Nl];
		PetscInt* nnzQ= new PetscInt[Nu+Nv+Nw+Nl+Nl];
		cout << Nu << " " << Nv << " " << Nw << " " << Nl << " " << N << " " << N+Nl << " " << nb << endl;
		for (int i =0;i<Nl;++i)
		{
			nnzP[i]=(3); //9*nb+1
			nnzQ[i]=(1);
			if (Nu > 0)
			{
				nnzP[Nu+i]=(6);
				nnzQ[Nu+i]=(2);
			}
			//if (Nv > 0 )
			//{
			//	nnzP[Nu+Nv+i]=(7*nb+1);
			//	nnzQ[Nu+Nv+i]=(7*nb+1);
			//}
			nnzP[Nu+Nv+Nw+i]=(2);
			nnzQ[Nu+Nv+Nw+i]=(2);
			nnzP[Nu+Nv+Nw+Nl+i]=(15);
			nnzQ[Nu+Nv+Nw+Nl+i]=(4*nb);				
		}
		
		std::cout<< " stage1"<<endl;
		
		//MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl+Nl, Nu+Nv+Nw+Nl+Nl, PETSC_NULL, nnzP,&P_);
		
		std::cout<< " stage2"<<endl;
		cout <<endl<< "nb = " << nb <<endl;
		MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl+Nl, Nu+Nv+Nw+Nl+Nl, PETSC_NULL, nnzQ,&Q_);
		MatCreateSeqAIJ(PETSC_COMM_SELF, 	N+Nl+Nl, N+Nl+Nl,       (3*7)*nb+1, 			PETSC_NULL, &P_);
		delete [] nnzP;
		delete [] nnzQ;
	}
	else
	{
		
		
	//MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl, Nu+Nv+Nw+Nl, PETSC_NULL, nnzR,&R_);		
	

	//delete [] nnzR;*/
	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	N+Nl+Nl, N+Nl+Nl,       (3*7)*nb+1, 			PETSC_NULL, &P_);
    //MatCreateSeqAIJ(PETSC_COMM_SELF, 	N+Nl+Nl, N+Nl+Nl,       (3*7)*nb+1, 			PETSC_NULL, &Q_);
	
	MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, N+Nl+Nl, N+Nl+Nl, (3*7)*nb+1, PETSC_NULL, (3*7)*nb+1, PETSC_NULL,&P_);
	MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, N+Nl+Nl, N+Nl+Nl, (3*7)*nb+1, PETSC_NULL, (3*7)*nb+1, PETSC_NULL,&Q_);
	
	PetscInt m, Istart, Iend;
	
	MatGetOwnershipRange(P_,&Istart,&Iend);
	m    = Iend-Istart;
	
	PR(m);
	
	
	//}

	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	N+Nl+Nl, N+Nl+Nl,       (3*7)*nb+1, 			PETSC_NULL, &P_);
    //cout <<endl<<Nu+Nv+Nw<<endl;
    
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");

						

    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    
    cout <<endl<< "N = " << N<<endl;
    
    for (unsigned int i =0;i<N+Nl;++i)
    {
		// P(ressure) variable does not evolve with modified midpoint rule. So P does not get a 1 here
        MatSetValue(P_, i,     i,     1.0, ADD_VALUES);//creating 2*I matrix in the upper part
        MatSetValue(Q_, i,     i,     1.0, ADD_VALUES);//creating 2*I matrix in the upper part

		if (i<Nl)	
		{
			// from C (from NMat)
			//MatSetValue(P_, 	Nu+Nv+i,     N+i,     -0.5*dt*1.0/N2, 	ADD_VALUES);
			//MatSetValue(Q_, 	Nu+Nv+i,     N+i,     0.5*dt*1.0/N2, 	ADD_VALUES);
			
			//MatSetValue(P_, 	Nu+Nv+i,     N+Nl+i,     0.5*dt*(1.0/N2+1.0), 	ADD_VALUES);
			//MatSetValue(Q_, 	Nu+Nv+i,     N+Nl+i,     -0.5*dt*(1.0/N2+1.0), 	ADD_VALUES);
			
			//MatSetValue(P_, 	N+Nl+i,    Nu+Nv+i,     -0.5*dt, 	ADD_VALUES);
			//MatSetValue(Q_, 	N+Nl+i,    Nu+Nv+i,     0.5*dt, 	ADD_VALUES);	
			
			if (config->WA_ == 0)
			{
			// 3 comes from beta_ (=N2+1)
			// from C1 (from VecNMat)
			//MatSetValue(P_, 	Nu+Nv+i,     N+i,     -0.5*dt*-1.0, 	ADD_VALUES);
			//MatSetValue(Q_, 	Nu+Nv+i,     N+i,     0.5*dt*-1.0, 	ADD_VALUES);
			
			//MatSetValue(P_, 	Nu+Nv+i,     N+Nl+i,     0.5*dt*(1.0/N2)*-3.0, 	ADD_VALUES);
			//MatSetValue(Q_, 	Nu+Nv+i,     N+Nl+i,     -0.5*dt*(1.0/N2)*-3.0, 	ADD_VALUES);
			
			//MatSetValue(P_, 	N+i,    Nu+Nv+i,     0.5*dt*-N2, 	ADD_VALUES); //-2.0
			//MatSetValue(Q_, 	N+i,    Nu+Nv+i,     -0.5*dt*-N2, 	ADD_VALUES);		
			}
			else
			{
			MatSetValue(P_, 	Nu+Nv+i,     N+i,     -0.5*dt*-1.0*cos(gamma), 	ADD_VALUES);
			MatSetValue(Q_, 	Nu+Nv+i,     N+i,     0.5*dt*-1.0*cos(gamma), 	ADD_VALUES);
			MatSetValue(P_, 	i,     N+i,     -0.5*dt*-1.0*sin(gamma), 	ADD_VALUES);
			MatSetValue(Q_, 	i,     N+i,     0.5*dt*-1.0*sin(gamma), 	ADD_VALUES);			

			
			MatSetValue(P_, 	N+i,    Nu+Nv+i,     0.5*dt*-N2*cos(gamma), 	ADD_VALUES); //-2.0
			MatSetValue(Q_, 	N+i,    Nu+Nv+i,     -0.5*dt*-N2*cos(gamma), 	ADD_VALUES);		
			MatSetValue(P_, 	N+i,   i,     0.5*dt*-N2*sin(gamma), 	ADD_VALUES); //-2.0
			MatSetValue(Q_, 	N+i,    i,     -0.5*dt*-N2*sin(gamma), 	ADD_VALUES);	
			}
		}		
    }
    
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    
    cout <<endl<<N<<endl;

    numberOfNonZeros=0;
	for (unsigned int i=0;i<N;++i)
	{	
		/*
        MatGetRow(C1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
				dummy = -dummy;
                MatSetValue(P_, 	i,     cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
                MatSetValue(Q_, 	i,     cols[j],     0.5*dt*dummy, 	ADD_VALUES);
				MatSetValue(R_, 	i,     cols[j],     dummy*2.0, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C1, i, &numberOfNonZeros, &cols, &values);
        */		
        MatGetRow(BF1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {	
				MatSetValue(P_, 	i,     N+Nl+cols[j],     -dt*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);
		
        if (i < Nl)
        {
			if (config->EB_==1)
			{
			MatGetRow(EBCWinR, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
				
                dummy = (values[j]);
                if (dummy!=0)
                {
					dummy = -dummy; 
					//cout << i << " " << cols[j] << " " <<  dummy << endl;
                    MatSetValue(P_, 	N+i,     Nu+Nv+cols[j],     -0.5*dt*dummy*-1.0, 	ADD_VALUES); // minus for definition N2=-drho/dz
                    MatSetValue(Q_, 	N+i,     Nu+Nv+cols[j],     0.5*dt*dummy*-1.0, 	ADD_VALUES);
                }
            }
            MatRestoreRow(EBCWinR, i, &numberOfNonZeros, &cols, &values);
			////MatSetValue(P_, 	N+i,     Nu+Nv+i,     -0.5*dt*N2, 	ADD_VALUES);
			////MatSetValue(Q_, 	N+i,     Nu+Nv+i,     0.5*dt*N2, 	ADD_VALUES);	
			
			//MatSetValue(P_, 	Nu+Nv+i,     N+i,     -0.5*dt*-1.0, 	ADD_VALUES);
			//MatSetValue(Q_, 	Nu+Nv+i,     N+i,     0.5*dt*-1.0, 	ADD_VALUES);			
			}
			else 
			{
            MatGetRow(CWinR, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
                dummy = (values[j]);
                if (dummy!=0)
                {
					dummy = -dummy; 
					//cout << i << " " << cols[j] << " " <<  dummy << endl;
                    MatSetValue(P_, 	N+i,     Nu+Nv+cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
                    MatSetValue(Q_, 	N+i,     Nu+Nv+cols[j],     0.5*dt*dummy, 	ADD_VALUES);
                }
            }
			MatRestoreRow(CWinR, i, &numberOfNonZeros, &cols, &values);
			}
            
            MatGetRow(sCRinW, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
                dummy = (values[j]);
                if (dummy!=0)
                {
					dummy = -dummy; 
                    MatSetValue(P_, 	Nu+Nv+i,     N+cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
                    MatSetValue(Q_, 	Nu+Nv+i,     N+cols[j],     0.5*dt*dummy, 	ADD_VALUES);
                }
            }
            MatRestoreRow(sCRinW, i, &numberOfNonZeros, &cols, &values);			
			
		
			
            MatGetRow(Ah, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
                dummy = (values[j]);
                if (dummy!=0)
                {
					//dummy = -dummy;
                    MatSetValue(P_, 	N+Nl+i,     N+cols[j],     -0.5*dummy, 	ADD_VALUES);
                    MatSetValue(Q_, 	N+Nl+i,     N+cols[j],     0.5*dummy, 	ADD_VALUES);
					//MatSetValue(R_, 	N+Nl+i,     Nu+cols[j],  -2.0/dt*dummy, 	ADD_VALUES);
                }
            }
            MatRestoreRow(Ah, i, &numberOfNonZeros, &cols, &values);
            
            MatGetRow(A, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
                dummy = (values[j]);

                if (dummy!=0)
                {
                    MatSetValue(P_, 	N+Nl+i,     	N+Nl+cols[j],     dummy, 		ADD_VALUES);
                }
            }
            MatRestoreRow(A, i, &numberOfNonZeros, &cols, &values);
        }
        
    }
					
    cout << "Mat Assembly for "<< "P"<<endl;
    MatAssemblyBegin(P_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(P_,MAT_FINAL_ASSEMBLY);
    
    cout << "Mat Assembly for "<< "Q"<<endl;
    MatAssemblyBegin(Q_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q_,MAT_FINAL_ASSEMBLY);
	
    //cout << "Mat Assembly for "<< "R"<<endl;
    //MatAssemblyBegin(R_,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(R_,MAT_FINAL_ASSEMBLY);	
    
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
	//printFullMatrixInfo(R_, "R_");
    
    
	outputMatrix(P_, "P.dat");
	//outputMatrix(Q_, "Q.dat");
	//outputMatrix(R_, "R.dat");
		
		
    std::cout<< " - : Finished Creating matrix P !\n";

	/*if (config->DIM_==3)
	{
		Mat& DX      = static_cast<HEulerGlobalVariables*>(globalData_)->DX_;
		Mat& DY      = static_cast<HEulerGlobalVariables*>(globalData_)->DY_;
		Mat& DZ      = static_cast<HEulerGlobalVariables*>(globalData_)->DZ_;
		
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu, 	Nl,         7*3*nb, 			PETSC_NULL, &DX);
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu, 	Nl,         7*3*nb, 			PETSC_NULL, &DY);
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu, 	Nl,         7*3*nb, 			PETSC_NULL, &DZ);

 	double dummy=0;
	const PetscInt* cols;
	const PetscScalar* values;
	PetscInt 	numberOfNonZeros;
	for (unsigned int i = 0; i<Nu;++i)
	{
		{		
        MatGetRow(BF1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
				MatSetValue(DX, 	i,     cols[j],     dummy, 	ADD_VALUES); 
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);	
        MatGetRow(BF1, Nu+i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
				MatSetValue(DY, 	i,     cols[j],     dummy, 	ADD_VALUES); 
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);	
		MatGetRow(BF1, Nu+Nv+i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
				MatSetValue(DZ, 	i,     cols[j],     dummy, 	ADD_VALUES); 
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);	
		}
		
	}
    cout << "Mat Assembly for "<< "DX"<<endl;
    MatAssemblyBegin(DX,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(DX,MAT_FINAL_ASSEMBLY);	
	cout << "Mat Assembly for "<< "DY"<<endl;
    MatAssemblyBegin(DY,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(DY,MAT_FINAL_ASSEMBLY);	
	cout << "Mat Assembly for "<< "DZ"<<endl;
    MatAssemblyBegin(DZ,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(DZ,MAT_FINAL_ASSEMBLY);		
	}
	*/
	MatDestroy(&BF1);
	MatDestroy(&A);
	MatDestroy(&Ah);
	//MatDestroy(&C1);
	
	MatDestroy(&CWinR);
	MatDestroy(&sCRinW);
	
	VecDestroy(&UInit);
	VecDestroy(&UCorrected);
	VecDestroy(&RHS);
    cout << "**************END of contructing matrices*************"<<endl;
}

void
HEuler::solve()
{
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
    HEulerGlobalVariables* globalData = static_cast<HEulerGlobalVariables*>(globalData_);
    
	double dt           = static_cast<HEulerGlobalVariables*>(globalData_)->dt_;
    
    unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nw     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
	unsigned int Nv		= 0;
	unsigned int Nu		= 0;
	if (config->DIM_ == 3)
	{
         Nv = Nw;
	}	
	if (config->DIM_ >= 2)
	{
         Nu = Nw;
	}	
	
    unsigned int Nl     = Nw;
    unsigned int N      = Nu+Nv+Nw;
    
    double endTime      = config->numOfPeriods_*config->onePeriod_ ;
    int nplotPerPeriod  = config->numOfPeriodsInOnePlotStep_;
    
	double Ham, TotalMass, xMom, yMom, zMom;
	double Kin, Pot, Int;
    Mat& globalDIV      = static_cast<HEulerGlobalVariables*>(globalData_)->DivergenceFreeMatrix_;
	//Mat& DX      = static_cast<HEulerGlobalVariables*>(globalData_)->DX_;
	//Mat& DY      = static_cast<HEulerGlobalVariables*>(globalData_)->DY_;
	//Mat& DZ      = static_cast<HEulerGlobalVariables*>(globalData_)->DZ_;
	
	/*
	Vec UVelU, UVelV, UVelW, URho;
	 VecCreateSeq(PETSC_COMM_SELF, Nu, &UVelU);
	 VecCreateSeq(PETSC_COMM_SELF, Nu, &UVelV);
	 VecCreateSeq(PETSC_COMM_SELF, Nu, &UVelW);
	 VecCreateSeq(PETSC_COMM_SELF, Nu, &URho);
	*/
    cout << "endTime="<<endTime<<endl;
    Vec UExact, Lambda, RHS, RH, UVel; //RHS1,RHO,
    
	//VecCreateSeq(PETSC_COMM_SELF, N+Nl+Nl, &RHO);
    VecCreateSeq(PETSC_COMM_SELF, N+Nl+Nl, &Lambda);
    VecCreateSeq(PETSC_COMM_SELF, N+Nl+Nl, &RHS);
	//VecCreateSeq(PETSC_COMM_SELF, N+Nl+Nl, &RHS1);
    
	 VecCreateSeq(PETSC_COMM_SELF, N+Nl+Nl, &UExact);
	 VecCreateSeq(PETSC_COMM_SELF, N, &UVel);
	 
	 VecCreateSeq(PETSC_COMM_SELF, Nl, &RH);
	 
    int Num =  config->numOfTimeStepInOnePeriod_*config->numOfPeriods_+2+1;
	
    Vec Hamiltonian, TotalMassVector, Divergence, xMomentum, yMomentum, zMomentum;
	
	Vec KineticEnergy, PotentialEnergy, InternalEnergy;

	
	VecCreateSeq(PETSC_COMM_SELF, Num, &KineticEnergy);
	VecCreateSeq(PETSC_COMM_SELF, Num, &PotentialEnergy);
	VecCreateSeq(PETSC_COMM_SELF, Num, &InternalEnergy);
	
    VecCreateSeq(PETSC_COMM_SELF, Num, &Hamiltonian);
    VecCreateSeq(PETSC_COMM_SELF, Num, &TotalMassVector);
    VecCreateSeq(PETSC_COMM_SELF, Num, &Divergence);
	VecCreateSeq(PETSC_COMM_SELF, Num, &xMomentum);
	VecCreateSeq(PETSC_COMM_SELF, Num, &yMomentum);
	VecCreateSeq(PETSC_COMM_SELF, Num, &zMomentum);

    double u, v, w, rho, p;
    
    //cout<<"Initializing the solver with velocity..."<<endl;
    double sum;
	double index;
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {
        Base::Element* element= *cit;
        int m = element->getID();
        
        unsigned int pos = m*nb;
        
        
        for (int i =0; i<nb;++i)
        {
            u = element->getData(0,0,i);
            VecSetValue(Lambda, pos+i,     u, INSERT_VALUES);
            VecSetValue(UExact, pos+i,     u, INSERT_VALUES);
			VecSetValue(UVel, pos+i,     u, INSERT_VALUES);
            //VecSetValue(UVelU, pos+i,     u, INSERT_VALUES);
            
			if (config->DIM_ >= 2)
			{
				v = element->getData(0,1,i);
				VecSetValue(Lambda, Nu+pos+i,   v, INSERT_VALUES);
				VecSetValue(UExact, Nu+pos+i,   v, INSERT_VALUES);
				VecSetValue(UVel, Nu+pos+i,   v, INSERT_VALUES);
				//VecSetValue(UVelV, pos+i,   v, INSERT_VALUES);
			}
            if (config->DIM_ == 3)
			{
				w = element->getData(0,2,i);
				VecSetValue(UExact, Nu+Nv+pos+i, w, INSERT_VALUES);
				VecSetValue(Lambda, Nu+Nv+pos+i, w, INSERT_VALUES);
				VecSetValue(UVel, Nu+Nv+pos+i, w, INSERT_VALUES);
				//VecSetValue(UVelW, pos+i, w, INSERT_VALUES);
			}
            
            rho = element->getData(0,3,i);
            VecSetValue(UExact, N+pos+i, rho, INSERT_VALUES);
            VecSetValue(Lambda, N+pos+i, rho, INSERT_VALUES);	
			//VecSetValue(URho, pos+i, rho, INSERT_VALUES);	
			
			if (config->incompressible_==0)
			{
            p = element->getData(0,4,i);
            VecSetValue(UExact, N+Nl+pos+i, p, INSERT_VALUES);
            VecSetValue(Lambda, N+Nl+pos+i, p, INSERT_VALUES);
			// incompressible case: old value of p is never used
			}

		
        }
    }
    VecAssemblyBegin(UExact);
    VecAssemblyEnd(UExact);
    //VecAssemblyBegin(RHO);
    //VecAssemblyEnd(RHO);	
    VecAssemblyBegin(Lambda);
    VecAssemblyEnd(Lambda);		
    VecAssemblyBegin(UVel);
    VecAssemblyEnd(UVel);	
	/*
	VecAssemblyBegin(UVelU);
    VecAssemblyEnd(UVelU);
    VecAssemblyBegin(UVelV);
    VecAssemblyEnd(UVelV);
	VecAssemblyBegin(UVelW);
    VecAssemblyEnd(UVelW);
    VecAssemblyBegin(URho);
    VecAssemblyEnd(URho);	
	*/
    PetscScalar* XTEMP = new PetscScalar [Nu+Nv+Nw+Nl+Nl];
	//PetscScalar* XTEMP1 = new PetscScalar [Nu+Nv+Nw+Nl+Nl];
	//PetscScalar* XTEMPRHO = new PetscScalar [Nu+Nv+Nw+Nl+Nl];
	HEuler::calculateEnergy(Lambda,Ham, Kin, Pot, Int, xMom, yMom, zMom, TotalMass);
	
	PRp(Ham);	
		
    VecSetValue(TotalMassVector, 0, 	TotalMass, INSERT_VALUES);
    VecSetValue(KineticEnergy,   0, 	Kin, INSERT_VALUES);
    VecSetValue(PotentialEnergy, 0, 	Pot, INSERT_VALUES);	
	VecSetValue(InternalEnergy,  0, 	Int, INSERT_VALUES);	
    VecSetValue(Hamiltonian,     0, 	Ham, INSERT_VALUES);
	VecSetValue(xMomentum,     	 0, 	xMom, INSERT_VALUES);
	VecSetValue(yMomentum,     	 0, 	yMom, INSERT_VALUES);
	VecSetValue(zMomentum,     	 0, 	zMom, INSERT_VALUES);

	if (config->incompressible_ == 1)
	{
		MatMult(globalDIV, UVel, RH);
		PetscInt a(10);
		PetscReal maxDiv;
		VecMax(RH, &a, &maxDiv);
		PRp(maxDiv);
		//cout << "Divergence = " << std::setprecision (20) << maxDiv<<endl;	
		VecSetValue(Divergence, 0, 	maxDiv, INSERT_VALUES);
		/*
		Vec RH1, RH2, RH3;
		VecCreateSeq(PETSC_COMM_SELF, Nl, &RH1);
		VecCreateSeq(PETSC_COMM_SELF, Nl, &RH2);
		VecCreateSeq(PETSC_COMM_SELF, Nl, &RH3);
		   
		// d w / dy    
		MatMult(DY, UVelW, RH1);
		// d v / dz
		MatMult(DZ, UVelV, RH2);
		// d w / dy - d v / dz
		VecAXPY(RH1,-1,RH2);
		// d rho / dx
		MatMult(DX, URho, RH3);
		// (d w / dy - d v / dz) * d rho / dx
		VecPointwiseMult(RH2, RH1,RH3);
		
		MatMult(DZ, UVelU, RH);
		MatMult(DX, UVelW, RH);
		MatMult(DY, URho, RH);

		MatMult(DX, UVelV, RH);
		MatMult(DY, UVelU, RH);
		MatMult(DZ, URho, RH);		
		*/
		
		
	}
    
    double currTime=0;
    
    int iterations=1;
    
	/*
	IS              isrow,iscol;            
	PetscErrorCode  ierr;
	MatOrderingType rtype = MATORDERINGQMD;//MATORDERINGRCM;
	
	MatGetOrdering(P_,rtype,&isrow,&iscol);
	
	MatPermute(P_,isrow,iscol,&P_);
	MatPermute(Q_,isrow,iscol,&Q_);
	VecPermute(UExact,iscol,PETSC_FALSE);
	*/
	
    KSP ksp;
        // Preconditioner
    PC pc;
   // cout << "Solve"<<endl;
        // Create a solver
    //KSPCreate(PETSC_COMM_SELF, &ksp);
	KSPCreate(PETSC_COMM_WORLD, &ksp);
    
    //cout << "ksp create"<<endl;
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    
        //! Direct solver: only a preconditioner is needed. To change
        //this to an iterative method, you can change the solver to
        //KSPMINRES, KSPGMRES, etc. Similarly, you can choose PCNONE
        //(instead of PCLU) for the preconditioner below if you want to
        //use a nonpreconditioned iterative method.
    
    //KSPSetOperators(ksp, P_, P_, SAME_NONZERO_PATTERN);
    KSPSetOperators(ksp, P_, P_);
    //cout << "Setup operators"<<endl;
    
    //KSPSetType(ksp, KSPGMRES);
	//KSPSetType(ksp, KSPMINRES);
	//KSPSetType(ksp, KSPGCR);
	//KSPSetType(ksp, KSPCG);
	//KSPSetType(ksp, KSPBCGS);
	//KSPSetType(ksp, KSPFGMRES);
    //cout << "Setup solver"<<endl; 
    //cout << "Setup Initial guess"<<endl;
    
    //KSPGetPC(ksp,&pc);
    //PCFactorSetFill(pc,5);
    //cout << "Setup pc"<<endl;
    
	/*if (config->incompressible_ == 0)
	{    
		KSPSetType(ksp, KSPGMRES);
		KSPGetPC(ksp,&pc);
		PCSetType(pc, PCILU);//PCNONE);
		PCFactorSetLevels(pc,1);
		PCFactorSetReuseOrdering(pc, PETSC_TRUE);
	}
	else if (config->incompressible_ == 1)// && config->stratified_ == 1)
	{	*/
	 
		//KSPSetType(ksp, KSPBCGS);
		//KSPGetPC(ksp,&pc);	
		//PCSetType(pc, PCSOR);

		KSPSetType(ksp, KSPGMRES);
		KSPGetPC(ksp,&pc);	
		
		//PCSetType(pc, PCILU);
		PCSetType(pc, PCNONE);
		//PCFactorSetLevels(pc,1);
		//PCFactorSetLevels(pc,3); // doesn't work for WA2D (1,1)
		PCFactorSetReuseOrdering(pc, PETSC_TRUE);	
		// PCILUSetLevels(PC pc,PetscInt levels); test this command
		//PCFactorSetFill(pc,4.40293);//3.02268); // run with -info to find optimal value
		 PCFactorSetReuseFill(pc, PETSC_TRUE); //test this command
		// KSPGMRESSetRestart(KSP ksp, PetscInt restart); can try to play around, default == 30
		
		
		//PCSetType(pc, PCNONE);
		
	//PCSetType(pc,PCCOMPOSITE);
	//PCCompositeSetType(pc,PC_COMPOSITE_ADDITIVE); //[PC COMPOSITE ADDITIVE,	PC COMPOSITE MULTIPLICATIVE]);
	//PCCompositeAddPC(pc,PCJACOBI);
	//PCCompositeAddPC(pc,PCILU);		
	
	

	//}	
		/*
		KSPSetType(ksp, KSPGMRES);
		KSPGetPC(ksp,&pc);
		PCSetType(pc, PCILU);//PCNONE);
		PCFactorSetLevels(pc,1);
		PCFactorSetReuseOrdering(pc, PETSC_TRUE);	
		*/
	//MatOrderingType rtype = MATORDERINGRCM;//MATORDERINGQMD;//
	//PCFactorSetMatOrderingType(pc,rtype);
	//else if (config->incompressible_ == 1 && config->stratified_ == 0)
	//{	
	//	KSPSetType(ksp, KSPBCGS);
	//	KSPGetPC(ksp,&pc);	
	//	PCSetType(pc, PCSOR);
	//}	
	
	

	//PCSetType(pc, PCJACOBI);
	//PCSetType(pc, PCSOR);
		//PCSetType(pc, PCMG);
		//PCSetType(pc, PCICC);
		//PCSetType(pc, PCLU);
        //
        //PCSetType(pc, PCCHOLESKY);
		
		
		// Use Matlab script to estimate the amount of Fill In
		// Important for optimization since most time is spent in iterating
		
        //PCILUSetFill(pc, 4); // outdated command
		//PCFactorSetFill(pc,5);

    KSPSetFromOptions(ksp);
	KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
        	//PCFactorSetFill(pc, 3.83016);
        	//PCFactorSetLevels(pc,4);
/*    
PCFactorSetLevels(PC pc,int levels);
PCFactorSetReuseOrdering(PC pc,PetscBool flag);
PCFactorSetDropTolerance(PC pc,double dt,double dtcol,int dtcount);
PCFactorSetReuseFill(PC pc,PetscBool flag);
PCFactorSetUseInPlace(PC pc);
PCFactorSetAllowDiagonalFill(PC pc);
*/

    //cout << "Setup options"<<endl;
    
    
    double reltol = 1.0e-14;
    double abstol = 1.0e-14;
	int    maxits = 1e8;//1e5
    KSPSetTolerances(ksp, reltol, abstol, PETSC_DEFAULT, maxits);
    //cout << "Setup tolerance"<<endl;
    
    PCSetUp(pc);
    KSPSetUp(ksp);
    //cout << "Setup ksp"<<endl<<endl;

	/*
     MatInfo           matinfo;
	 double fillrationeeded;
	 MatGetInfo(P_,MAT_LOCAL,&matinfo);
     printf("matinfo.nz_used %g\n",matinfo.nz_used);
	 fillrationeeded  = matinfo.fill_ratio_needed;
	 printf("matinfo.fill_ratio_needed %g\n",matinfo.fill_ratio_needed);
	  * */
    
        //! Solving linear system.
	cout << "n_x = "<< config->nx_ << endl;
	cout << "n_b = "<< configData_->numberOfBasisFunctions_ << endl;
	cout << "n_El = " << static_cast<HEulerGlobalVariables*>(globalData_)->nElements_ << endl;
	cout << "Nu = "<< Nu << endl;	
	cout << "dt = " << dt << endl;
    cout << endl << "Solving System" << endl;
	
    while(currTime<=endTime)//
    {     
        //cout<<"construct proper rightHandside"<<endl;
        MatMult(Q_, UExact, RHS);
		/*
		if (config->incompressible_ == 0)
		{
			VecAssemblyBegin(RHS);
			VecAssemblyEnd(RHS);			
			MatMult(R_, RHO, RHS1);
			VecAssemblyBegin(RHS1);
			VecAssemblyEnd(RHS1);
			VecAXPY(RHS,1.0,RHS1);
		}
		 * */
		
        
            ///***********************************************************************************
            ///***************************construct proper rightHandside**************************
        
        
        //cout << "Finalizing vector creation"<<endl;
        
        
            // KSP Solver initializer
            // 	    outputVectorMatlab(RHS, "RHS.dat");
        KSPSolve(ksp, RHS, Lambda);
            //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
		if (iterations%nplotPerPeriod==0)
		{
        if(reason < 0)
        {
            PetscInt its;
            KSPGetIterationNumber(ksp, &its);
            PetscReal rnorm;
            KSPGetResidualNorm(ksp, &rnorm);
            cout << "\t\tPetSc: Solving the linear system has failed. Reason code: "
            << reason << endl //<< "Check KSPConvergedReason for the reason" << endl
            << "\t\tPetSc: Residual after " << int(its) << " iterations : ";
            cout.setf(ios::scientific, ios::floatfield);
            cout.precision(4);
            cout << rnorm << endl;
            cout.setf(ios::fixed, ios::floatfield);
            cout.precision(5);
        }
        
        else
        {
            if(reason > 0)
            {
				
                PetscInt its;
                KSPGetIterationNumber(ksp, &its);
                PetscReal rnorm;
                KSPGetResidualNorm(ksp, &rnorm);
                cout << "\t\tPetSc: Solving the linear system has succeeded. Reason code: "
                << reason << endl //<< "Check KSPConvergedReason for the reason" << endl
                << "\t\tPetsc: Convergence in " << int(its) << " iterations : ";
                cout.setf(ios::scientific, ios::floatfield);
                cout.precision(4);
                cout << rnorm << endl;
                cout.setf(ios::fixed, ios::floatfield);
                cout.precision(5);
				
            }
            else
            {
                cout << "\t\tPetSc: Solving the linear system is still under way" << endl;
            }
        }
		}
        
        //cout << "Solved a timestep" << endl;
        
        
        VecGetArray(Lambda, &XTEMP);
		//VecGetArray(UExact, &XTEMP1);
		//VecGetArray(RHO, &XTEMPRHO);
		
        int pos=0;
        int k=0;
        currTime+=dt;
        //for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
        //{
        //    Base::Element* element= *cit;
        //    int k = element->getID();
        //    
        //    unsigned int pos = k*nb;
        //    
        //    for (int i =0; i<nb;++i)
        //    {
		//		//element->setData(0,3,i, XTEMPRHO[Nu+Nv+pos+i]+0.5*dt*(XTEMP[Nu+Nv+pos+i]+XTEMP1[Nu+Nv+pos+i])); //set rho
		//	}
		//}
		//VecRestoreArray(UExact, &XTEMP1);
		//VecRestoreArray(RHO, &XTEMPRHO);
        pos=0;
        k=0;
		//sum = 0.0;
		//index = 0.0;
        for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
        {
            Base::Element* element= *cit;
            int k = element->getID();
            
            unsigned int pos = k*nb;
            
            for (int i =0; i<nb;++i)
            {	
		
                element->setData(0,3,i, XTEMP[Nu+Nv+Nw+pos+i]); //set rho
				VecSetValue(UExact, Nu+Nv+Nw+pos+i, XTEMP[Nu+Nv+Nw+pos+i], INSERT_VALUES);
				
				element->setData(0,4,i, XTEMP[Nu+Nv+Nw+Nl+pos+i]); //set p
				//sum = sum + XTEMP[Nu+Nv+Nw+Nl+pos+i];
				//index = index + 1;
				if (config->incompressible_ == 0)
				{                
					VecSetValue(UExact, Nu+Nv+Nw+Nl+pos+i, XTEMP[Nu+Nv+Nw+Nl+pos+i], INSERT_VALUES);	
				}
                element->setData(0,0,i, XTEMP[pos+i]);//set U
                VecSetValue(UExact, pos+i,     		XTEMP[pos+i], INSERT_VALUES);
				VecSetValue(UVel, 	pos+i,     		XTEMP[pos+i], INSERT_VALUES);
                
				if (config->DIM_ >= 2)
				{
					element->setData(0,1,i, XTEMP[Nu+pos+i]);//set V
					VecSetValue(UExact, Nu+pos+i,   	XTEMP[Nu+pos+i], INSERT_VALUES);
					VecSetValue(UVel, 	Nu+pos+i,   	XTEMP[Nu+pos+i], INSERT_VALUES);
				}
				if (config->DIM_ == 3)
				{                
					element->setData(0,2,i, XTEMP[Nu+Nv+pos+i]);//set W
					VecSetValue(UExact, Nu+Nv+pos+i, 	XTEMP[Nu+Nv+pos+i], INSERT_VALUES);
					VecSetValue(UVel, 	Nu+Nv+pos+i, 	XTEMP[Nu+Nv+pos+i], INSERT_VALUES);
				}
				
				//rho = element->getData(0,4,i);
				//VecSetValue(RHO, Nu+Nv+pos+i, 	rho, INSERT_VALUES);
				//VecSetValue(RHO, N+pos+i, 		rho, INSERT_VALUES);				
            }
        }
		VecRestoreArray(Lambda, &XTEMP);
		//cout << "Average P = " << sum/index << endl;
        VecAssemblyBegin(UExact);
		VecAssemblyEnd(UExact);	
		VecAssemblyBegin(UVel);
        VecAssemblyEnd(UVel);
        int numOfPeriodsInOnePlotStep  = config->numOfPeriodsInOnePlotStep_;
		int numOfHamInOnePeriod = config->numOfHamInOnePeriod_;
		if ((iterations*numOfHamInOnePeriod)%(numOfPeriodsInOnePlotStep) == 0)
		{
			HEuler::calculateEnergy(Lambda,Ham, Kin, Pot, Int, xMom, yMom, zMom, TotalMass);
			VecSetValue(TotalMassVector, iterations, 	TotalMass, INSERT_VALUES);
			VecSetValue(Hamiltonian,     iterations, 	Ham, INSERT_VALUES);
			VecSetValue(KineticEnergy,   iterations, 	Kin, INSERT_VALUES);
			VecSetValue(InternalEnergy,  iterations, 	Int, INSERT_VALUES);
			VecSetValue(PotentialEnergy, iterations, 	Pot, INSERT_VALUES);
			VecSetValue(xMomentum,     	 iterations, 	xMom, INSERT_VALUES);
			VecSetValue(yMomentum,     	 iterations, 	yMom, INSERT_VALUES);
			VecSetValue(zMomentum,     	 iterations, 	zMom, INSERT_VALUES);	
			if (config->incompressible_ == 1)
			{
				MatMult(globalDIV, UVel, RH);
				PetscInt a(10);
				PetscReal maxDiv;
				VecMax(RH, &a, &maxDiv);
				VecSetValue(Divergence, iterations, 	maxDiv, INSERT_VALUES);
				//PRp(maxDiv);
				
			}			
	        if (iterations%nplotPerPeriod==0)
			{
				output(currTime);
				cout<<"currTime="<<currTime<<endl;
				
				PRp(Ham);
				cout << endl;
			}
		}

        iterations++;
    }
    
    cout << "The end of the time loop"<<endl;
    cout << "iterations =  " << iterations << endl;
	cout << "Nu = "<<Nu<<", nb = "<<nb<<endl;
    cout << "ne = "<<Nu/nb<<endl;		
	HEuler::calculateEnergy(Lambda,Ham, Kin, Pot, Int, xMom, yMom, zMom, TotalMass);
	VecSetValue(TotalMassVector, iterations, 	TotalMass, INSERT_VALUES);
	VecSetValue(Hamiltonian,     iterations, 	Ham, INSERT_VALUES);
	VecSetValue(KineticEnergy,   iterations, 	Kin, INSERT_VALUES);
	VecSetValue(InternalEnergy,  iterations, 	Int, INSERT_VALUES);
	VecSetValue(PotentialEnergy, iterations, 	Pot, INSERT_VALUES);
	VecSetValue(xMomentum,     	 iterations, 	xMom, INSERT_VALUES);
	VecSetValue(yMomentum,     	 iterations, 	yMom, INSERT_VALUES);
	VecSetValue(zMomentum,     	 iterations, 	zMom, INSERT_VALUES);
	if (config->incompressible_ == 1)
	{
		MatMult(globalDIV, UVel, RH);
		PetscInt a(10);
		PetscReal maxDiv;
		VecMax(RH, &a, &maxDiv);
		VecSetValue(Divergence, iterations, 	maxDiv, INSERT_VALUES);
		//cout << "Divergence= " << std::setprecision (20) << maxDiv<<endl;
		//PRp(maxDiv);
	}
	
    KSPDestroy(&ksp);
    VecDestroy(&UExact);
	VecDestroy(&UVel);
    VecDestroy(&RHS);
    VecDestroy(&Lambda);
    VecDestroy(&RH);
	//VecDestroy(&RHS1);
	//VecDestroy(&RHO);
    

	
	PetscScalar* HAMTEMP = new PetscScalar [Num];
	PetscScalar* KINTEMP = new PetscScalar [Num];
	PetscScalar* POTTEMP = new PetscScalar [Num];
	PetscScalar* INTTEMP = new PetscScalar [Num];
    PetscScalar* MASSTEMP = new PetscScalar [Num];
    PetscScalar* DIVTEMP = new PetscScalar [Num];
	PetscScalar* xMomTEMP = new PetscScalar [Num];
	PetscScalar* yMomTEMP = new PetscScalar [Num];
	PetscScalar* zMomTEMP = new PetscScalar [Num];
    VecGetArray(Hamiltonian, &HAMTEMP);
	VecGetArray(KineticEnergy, &KINTEMP);
	VecGetArray(PotentialEnergy, &POTTEMP);
	VecGetArray(InternalEnergy, &INTTEMP);
    VecGetArray(TotalMassVector, &MASSTEMP);
    VecGetArray(Divergence, &DIVTEMP);
	VecGetArray(xMomentum, &xMomTEMP);
	VecGetArray(yMomentum, &yMomTEMP);
	VecGetArray(zMomentum, &zMomTEMP);

    currTime = 0.0;
    ofstream fout("HAM.txt");
	
	fout << "Time" << " " << "Kinetic Energy" << " " << "Potential Energy" << " " << "Internal Energy" << " " << "Total Energy" << " " << "Total Mass" << " " << "Total X-Momentum" << " " << "Total Y-Momentum" << " " << "Total Z-Momentum" << " " << "Maximum Divergence" <<endl; 
        
    for (int i =0; i<Num-1;++i)
    {
        fout << currTime << " " << std::setprecision (20) << KINTEMP[i] << " " << std::setprecision (20) << POTTEMP[i] << " " << std::setprecision (20) << INTTEMP[i] << " " << std::setprecision (20) << HAMTEMP[i] << " " << std::setprecision (20) << MASSTEMP[i] << " " << std::setprecision (20) << xMomTEMP[i] << " " << std::setprecision (20) << yMomTEMP[i] << " " << std::setprecision (20) << zMomTEMP[i] << " " << std::setprecision (20) << DIVTEMP[i] <<endl; 
        currTime += dt; //*config->numOfTimeStepInOnePeriod_/config->numOfHamInOnePeriod_; should then also change the Num, the length of the vectors, otherwise our HAM.txt is mostly zeros.
    }
    fout.close();
    VecRestoreArray(Hamiltonian, &HAMTEMP);
	VecRestoreArray(KineticEnergy, &KINTEMP);
	VecRestoreArray(PotentialEnergy, &POTTEMP);	
	VecRestoreArray(InternalEnergy, &INTTEMP);	
    VecRestoreArray(TotalMassVector, &MASSTEMP);
    VecRestoreArray(Divergence, &DIVTEMP);
	VecRestoreArray(xMomentum, &xMomTEMP);
	VecRestoreArray(yMomentum, &yMomTEMP);
	VecRestoreArray(zMomentum, &zMomTEMP);	

    VecDestroy(&Hamiltonian);
	VecDestroy(&KineticEnergy);
	VecDestroy(&PotentialEnergy);
	VecDestroy(&InternalEnergy);
    VecDestroy(&TotalMassVector);
    VecDestroy(&Divergence);	
	VecDestroy(&xMomentum);	
	VecDestroy(&yMomentum);	
	VecDestroy(&zMomentum);	
    
}

void
HEuler::createCompressibleSystemnew()
{
 	PetscInfoAllow(PETSC_TRUE, "history.txt");
	
	//Mat C;
	Mat BF;
    
	//Mat C1;
	Mat BF1;
	Mat MInv;
	Mat sMInv;
	//Mat A;
	//Mat Ah;
    Mat globalDIV;
	Mat DIV;
	
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);	
	
    unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nw     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
	unsigned int Nv		= 0;
	unsigned int Nu     = 0;
	if (config->DIM_ == 3)
	{
         Nv= Nw;
	}
	if (config->DIM_ >= 2)
	{
         Nu = Nw;
	}	
    
    unsigned int Nl     = Nw;
    unsigned int N      = Nu+Nv+Nw;	
	
	double N2			= config->N2_;
	
    double dt           = static_cast<HEulerGlobalVariables*>(globalData_)->dt_;
    cout << "dt=" <<dt<<endl;	
	
	cout << "**************Starting create Matrices*************"<<endl;
	Mat M1;
	Mat sM1;
	Mat VecNMat;
	Mat NMat;
 	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    1*nb, 				PETSC_NULL, &M1);
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &sM1);
	
 	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    1*nb, 				PETSC_NULL, &MInv);
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &sMInv);	
	
	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nw, 	Nl,          1*nb, 				PETSC_NULL, &VecNMat);
	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nw, 	Nl,          1*nb, 				PETSC_NULL, &NMat);		
	
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nl,         (7)*nb, 			PETSC_NULL, &BF);//number of possible nonzero blocks are 7: element and his 6 neighbours
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 	Nu+Nv+Nw,    (7*3)*nb, 				PETSC_NULL, &globalDIV);//number of possible nonzero blocks are 7: element and his 6 neighbours	
	
	cout << "**************Mass Matrix Calculated*************"<<endl;

//	timer t0;
//	std::cout << t0.elapsed() << " - : Element Integration started\n";
    ElementIntegralData gradMass;
    bool useCache = false;
    ElementIntegralT   elIntegral(useCache);
    typedef void  (HEuler::*Integrand)(const ElementT* , const PointReferenceT&, ElementIntegralData&);
    Integrand gradMassInteg = &HEuler::elementIntegrand;
    
    std::cout << " - : Element Integration started\n";
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {

        int m  = (*cit)->getID();
        unsigned int pos = (*cit)->getID()*nb;
        
        HEulerElementData* elementData = static_cast<HEulerElementData*>((*cit)->getUserData());
        
        ElementT* element =(*cit);
        
        elIntegral.integrate<ElementIntegralData>(element, this, gradMass);
		
		for (unsigned int j = 0; j<nb;++j)
		{
			for (unsigned int i=0; i<nb;++i)//can be optimized later!!!
 			{
                //MatSetValue(NMat, 	Nu+Nv+pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensity_(j,i), 		ADD_VALUES);
				//MatSetValue(VecNMat, 	Nu+Nv+pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensityDeriv_(j,i), 	ADD_VALUES);
                
				//MatSetValue(NMat, 		pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensity_(j,i), 		ADD_VALUES);
				//MatSetValue(VecNMat, 	pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensityDeriv_(j,i), 	ADD_VALUES);	

                //MatSetValue(NMatT, 		pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensity_(i,j), 		ADD_VALUES);
				//MatSetValue(VecNMatT, 	pos+j, 	pos+i, 		elementData->massMatrixTimesBackgroundDensityDeriv_(i,j), 	ADD_VALUES);	
                
                MatSetValue(M1,   		pos+j, 			pos+i,     	elementData->massMatrixOverBackgroundDensity_(j,i), 		ADD_VALUES);
				MatSetValue(sM1,   		pos+j, 			pos+i,      elementData->massMatrixOverBackgroundDensity_(j,i), 		ADD_VALUES);
            
                MatSetValue(MInv, 		pos+j,       	pos+i,    	elementData->invMassMatrix_(j,i) , 							ADD_VALUES);
                MatSetValue(sMInv, 		pos+j, 			pos+i,      elementData->invMassMatrix_(j,i) , 							ADD_VALUES);
                
                MatSetValue(BF, 		pos+j, 		 	pos+i,   	-gradMass.xGrad_(i,j),  									ADD_VALUES);
				MatSetValue(globalDIV, 	pos+j, 			pos+i,   	gradMass.xGrad_(j,i),  										ADD_VALUES);

				if (config->DIM_ >= 2)
				{
					MatSetValue(MInv, 		Nu+pos+j, 		Nu+pos+i, 	elementData->invMassMatrix_(j,i) , 							ADD_VALUES);
					MatSetValue(M1,   		Nu+pos+j, 		Nu+pos+i,   elementData->massMatrixOverBackgroundDensity_(j,i), 		ADD_VALUES);
					MatSetValue(BF, 		Nu+pos+j, 	 	pos+i,   	-gradMass.yGrad_(i,j),  									ADD_VALUES);
					MatSetValue(globalDIV, 	pos+j, 			Nu+pos+i,   gradMass.yGrad_(j,i),  										ADD_VALUES);					
				}				
				if (config->DIM_ == 3)
				{
					MatSetValue(MInv, 		Nu+Nv+pos+j,    Nu+Nv+pos+i,    elementData->invMassMatrix_(j,i), 					ADD_VALUES);	
					MatSetValue(M1,   		Nu+Nv+pos+j, 	Nu+Nv+pos+i,    elementData->massMatrixOverBackgroundDensity_(j,i), ADD_VALUES);					
					MatSetValue(BF, 		Nu+Nv+pos+j, 	pos+i,   		-gradMass.zGrad_(i,j),  							ADD_VALUES);
					MatSetValue(globalDIV, 	pos+j, 			Nu+Nv+pos+i,	gradMass.zGrad_(j,i),  								ADD_VALUES);
				}			
                                                                                          
 			}
		}
    }
	cout << "Mat Assembly for "<< "MINV"<<endl;
    MatAssemblyBegin(MInv,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MInv,MAT_FINAL_ASSEMBLY);
	printFullMatrixInfo(MInv, "MInv");
    cout << "Mat Assembly for "<< "MINVSm"<<endl;
    MatAssemblyBegin(sMInv,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sMInv,MAT_FINAL_ASSEMBLY);
	printFullMatrixInfo(sMInv, "sMInv");
    cout << "Mat Assembly for "<< "M1"<<endl;
    MatAssemblyBegin(M1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for "<< "sM1"<<endl;
    MatAssemblyBegin(sM1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sM1,MAT_FINAL_ASSEMBLY);
    //cout << "Mat Assembly for "<< "NMat"<<endl;
    //MatAssemblyBegin(NMat, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(NMat,MAT_FINAL_ASSEMBLY);
    //cout << "Mat Assembly for "<< "VecNMat"<<endl;
    //MatAssemblyBegin(VecNMat,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(VecNMat,MAT_FINAL_ASSEMBLY);
	
    //cout << "Mat Assembly for "<< "NMatT"<<endl;
    //MatAssemblyBegin(NMatT, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(NMatT,MAT_FINAL_ASSEMBLY);
    //cout << "Mat Assembly for "<< "VecNMatT"<<endl;
    //MatAssemblyBegin(VecNMatT,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(VecNMatT,MAT_FINAL_ASSEMBLY);	
		
    unsigned int posR=0;
    unsigned int posL=0;
    unsigned int eR=0;
    unsigned int eL=0;

    std::cout << " - : Face Integration started\n";
        //
        //
    FluxData fData(nb);
    //typedef void  (HEuler::*FaceIntegrand)(const FaceT*, const PointPhysicalT& normal , const PointReferenceOnTheFaceT&, FluxData&);
    //FaceIntegrand faceInteg = &HEuler::faceIntegrand;
    FaceIntegralT   faceIntegral(useCache);
    
    for (ConstFaceIterator citFe = faceColBegin(); citFe != faceColEnd(); ++citFe)
    {
        
        if ((*citFe)->getPtrElementRight()== NULL) // boundary face
        {
        }
        else
        {
            eR = (*citFe)->getPtrElementRight()->getID();
            eL = (*citFe)->getPtrElementLeft()->getID();
            
            posR  = eR*nb;
            
            posL  = eL*nb;
            
            faceIntegral.integrate((*citFe), this, fData);

            for (unsigned int j=0;j<nb;++j)
            {
                for (unsigned int i=0;i<nb;++i)
                {
                    MatSetValue(BF, posR+j, 		posL+i,     fData.right_[j](0,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posR+j, 		posR+i,    -fData.right_[j](3,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posL+j, 		posL+i,   	fData.left_[j](0,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posL+j, 		posR+i,    -fData.left_[j](3,i),  	ADD_VALUES);//U coefficient	

                    MatSetValue(globalDIV, posR+j, posL+i,    	  fData.right_[j](6,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posL+j, posL+i,   	 -fData.left_[j](6,i),  	ADD_VALUES);//U coefficient		
                    MatSetValue(globalDIV, posL+j, posR+i,   	 -fData.left_[j](9,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posR+j, posR+i,    	  fData.right_[j](9,i),  	ADD_VALUES);//U coefficient		
						
					if (config->DIM_ >= 2)
					{					
                    MatSetValue(BF, Nu+posR+j,      posL+i,     fData.right_[j](1,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+posR+j,      posR+i,    -fData.right_[j](4,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+posL+j,      posL+i,   	fData.left_[j](1,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+posL+j,      posR+i,    -fData.left_[j](4,i),  	ADD_VALUES);//V coefficient
				
                    MatSetValue(globalDIV, posL+j, Nu+posL+i,    -fData.left_[j](7,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posR+j, Nu+posL+i,     fData.right_[j](7,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posR+j, Nu+posR+i,     fData.right_[j](10,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posL+j, Nu+posR+i,    -fData.left_[j](10,i),  	ADD_VALUES);//V coefficient
					}
					
					if (config->DIM_ >= 3) // changed 2 to 3
					{						
                    MatSetValue(BF, Nu+Nv+posL+j, 	posR+i,    -fData.left_[j](5,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posL+j, 	posL+i,   	fData.left_[j](2,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posR+j, 	posR+i,    -fData.right_[j](5,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posR+j, 	posL+i,     fData.right_[j](2,i),  	ADD_VALUES);//W coefficient
					
                    MatSetValue(globalDIV, posL+j, Nu+Nv+posL+i, -fData.left_[j](8,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posR+j, Nu+Nv+posL+i,  fData.right_[j](8,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posL+j, Nu+Nv+posR+i, -fData.left_[j](11,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posR+j, Nu+Nv+posR+i,  fData.right_[j](11,i),  	ADD_VALUES);//W coefficient
					}
                }
            }
                //                cout <<"***********************************"<<endl;
        }
    }
	cout << "Mat Assembly for " << "globalDIV"<<endl;
    MatAssemblyBegin(globalDIV,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalDIV,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for " << "BF"<<endl;
    MatAssemblyBegin(BF,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(BF,MAT_FINAL_ASSEMBLY);
	
    printFullMatrixInfo(globalDIV, "globalDIV");
    printFullMatrixInfo(BF, "BF");	
	
    double fillBF=1;//(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
    //MatMatMult(MInv, BF, MAT_INITIAL_MATRIX, fillBF, &BF1);
	//MatMatMult(M1, BF1, MAT_INITIAL_MATRIX, fillBF, &BF1);
	if (config->stratified_ == 1)
	{	
    //MatMatMult(MInv, BF, MAT_INITIAL_MATRIX, fillBF, &BF1);
	//MatMatMult(M1, BF1, MAT_INITIAL_MATRIX, fillBF, &BF1);		
	
		MatMatMult(BF, sMInv, MAT_INITIAL_MATRIX, fillBF, &BF1);
		MatMatMult(BF1, sM1, MAT_INITIAL_MATRIX, fillBF, &BF1);
		MatMatMult(MInv, BF1, MAT_INITIAL_MATRIX, fillBF, &BF1);
	}
	else
	{
		MatMatMult(MInv, BF, MAT_INITIAL_MATRIX, fillBF, &BF1);	
	}

	double fillDIV=1;//(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
    //MatMatMult(sMInv, globalDIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
	//MatMatMult(sM1, DIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
	
	if (config->stratified_ == 1)
	{
    //MatMatMult(sMInv, globalDIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
	//MatMatMult(sM1, DIV, MAT_INITIAL_MATRIX, fillBF, &DIV);		
	
		MatMatMult(globalDIV, MInv, MAT_INITIAL_MATRIX, fillBF, &DIV);
		MatMatMult(DIV, M1, MAT_INITIAL_MATRIX, fillBF, &DIV);	
		MatMatMult(sMInv, DIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
	}
	else
	{		
		MatMatMult(sMInv, globalDIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
	}
	/*
	if (config->stratified_ == 1)
	{	// 1D Case:
		MatMatMult(sMInv, NMat, MAT_INITIAL_MATRIX, 1., &C);
		MatMatMult(sM1, C, MAT_INITIAL_MATRIX, 1., &C);		
		MatMatMult(sMInv, C, MAT_INITIAL_MATRIX, 1., &C);

		MatMatMult(sMInv, VecNMat, MAT_INITIAL_MATRIX, 1., &C1);
		MatMatMult(sM1, C1, MAT_INITIAL_MATRIX, 1., &C1);		
		MatMatMult(sMInv, C1, MAT_INITIAL_MATRIX, 1., &C1);	
	}
	else	
	{
		MatMatMult(sMInv, NMat, MAT_INITIAL_MATRIX, 1., &C);
		MatMatMult(sMInv, VecNMat, MAT_INITIAL_MATRIX, 1., &C1);
	}
	*/
	
		//MatMatMult(sMInv, NMatT, MAT_INITIAL_MATRIX, 1., &CT);
		//MatMatMult(sM1, CT, MAT_INITIAL_MATRIX, 1., &CT);		
		//MatMatMult(sMInv, CT, MAT_INITIAL_MATRIX, 1., &CT);

		//MatMatMult(sMInv, VecNMatT, MAT_INITIAL_MATRIX, 1., &C1T);
		//MatMatMult(sM1, C1T, MAT_INITIAL_MATRIX, 1., &C1T);		
		//MatMatMult(sMInv, C1T, MAT_INITIAL_MATRIX, 1., &C1T);			
		
		// 2D Case:
		//MatMatMult(NMat, sMInv, MAT_INITIAL_MATRIX, fillBF, &C);
		//MatMatMult(C, sM1, MAT_INITIAL_MATRIX, fillBF, &C);
		//MatMatMult(MInv, C, MAT_INITIAL_MATRIX, fillBF, &C);

		//MatMatMult(VecNMat, sMInv, MAT_INITIAL_MATRIX, fillBF, &C1);
		//MatMatMult(C1, sM1, MAT_INITIAL_MATRIX, fillBF, &C1);
		//MatMatMult(MInv, C1, MAT_INITIAL_MATRIX, fillBF, &C1);		
	
	//outputMatrix(sMInv, "sMInv.dat");
	//outputMatrix(sM1, "sM1.dat");
	
	//outputMatrix(NMat, "NMat.dat");
	//outputMatrix(NMatT, "NMatT.dat");	
	
	//outputMatrix(VecNMat, "VecNMat.dat");
	//outputMatrix(VecNMatT, "VecNMatT.dat");
	
	MatDestroy(&BF);
	//MatDestroy(&NMat);
	//MatDestroy(&VecNMat);
	MatDestroy(&globalDIV);
	//MatDestroy(&C1);
	//MatDestroy(&C);
	MatDestroy(&sMInv);
	MatDestroy(&MInv);
	MatDestroy(&sM1);
	MatDestroy(&M1);
	
    cout << "Mat Assembly for "<< "BF1"<<endl;
    MatAssemblyBegin(BF1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(BF1,MAT_FINAL_ASSEMBLY);
    
    cout << "Mat Assembly for "<< "DIV"<<endl;
    MatAssemblyBegin(DIV,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(DIV,MAT_FINAL_ASSEMBLY);
	
    //cout << "Mat Assembly for "<< "C"<<endl;
    //MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);

	//cout << "Mat Assembly for "<< "C1"<<endl;
    //MatAssemblyBegin(C1,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(C1,MAT_FINAL_ASSEMBLY);	
    
    printFullMatrixInfo(DIV, "DIV");
    printFullMatrixInfo(BF1, "BF1");	
	//printFullMatrixInfo(C, "C");
	//printFullMatrixInfo(C1, "C1");
	
    //cout << "Mat Assembly for "<< "CT"<<endl;
    //MatAssemblyBegin(CT,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(CT,MAT_FINAL_ASSEMBLY);

	//cout << "Mat Assembly for "<< "C1T"<<endl;
    //MatAssemblyBegin(C1T,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(C1T,MAT_FINAL_ASSEMBLY);	
	
	//outputMatrix(C, "C.dat");
	//outputMatrix(C1, "C1.dat");	
	
	//outputMatrix(C, "CT.dat");
	//outputMatrix(C1, "C1T.dat");

			
	
    double dummy=0;
    const PetscInt* cols;
    const PetscScalar* values;
    PetscInt 	numberOfNonZeros;
    std::cout << " - : Started Creating P matrix!\n";

	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nl+Nl,         (7)*nb, 			PETSC_NULL, &B_);//number of possible nonzero blocks are 7: element and his 6 neighbours
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl+Nl, 	Nu+Nv+Nw,    (7*3)*nb, 				PETSC_NULL, &C_);
	//Mat& B = static_cast<HEulerGlobalVariables*>(globalData_)->B_;
	//Mat& C = static_cast<HEulerGlobalVariables*>(globalData_)->C_;	
	
    for (unsigned int i =0;i<N+Nl+Nl;++i)
    {

		// replacing C and C1:
		if (i<Nl)	
		{

			MatSetValue(B_, 	Nu+Nv+i, 		i, -0.5*dt*1.0/N2, 	ADD_VALUES);

			MatSetValue(B_, 	Nu+Nv+i, 		Nl+i, 0.5*dt*(1.0/N2+1.0), 	ADD_VALUES);

			MatSetValue(C_, 	Nl+i, 		Nu+Nv+i, -0.5*dt*1.0/N2, 	ADD_VALUES);

			MatSetValue(B_, 	Nu+Nv+i, 		i, -0.5*dt*1.0/N2*-3.0, 	ADD_VALUES);

			MatSetValue(B_, 	Nu+Nv+i, 		Nl+i, 0.5*dt*(1.0/N2)*-3.0, 	ADD_VALUES);


			MatSetValue(C_, 	i, 		Nu+Nv+i, -0.5*dt*1.0/N2, 	ADD_VALUES);		
						
		}	

    }
    
    cout <<endl<< "N = " << N<<endl;
	
    for (unsigned int i=0;i<N;++i)
    {
        MatGetRow(BF1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {

				//MatSetValue(P_, 	i,     N+cols[j],     -0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
				//MatSetValue(Q_, 	i,     N+cols[j],     0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
				
				//MatSetValue(P_, 	i,     N+cols[j],     0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
				//MatSetValue(Q_, 	i,     N+cols[j],     -0.5*dt*1.0/N2*dummy, 	ADD_VALUES);

				//MatSetValue(P_, 	i,     N+Nl+cols[j],     0.5*dt*(1.0/N2)*dummy, 	ADD_VALUES);
                //MatSetValue(Q_, 	i,     N+Nl+cols[j],     -0.5*dt*(1.0/N2)*dummy, 		ADD_VALUES);	
				
				//MatSetValue(P_, 	i,     N+Nl+cols[j],     -0.5*dt*(1.0/N2+1.0)*dummy, 	ADD_VALUES);
                //MatSetValue(Q_, 	i,     N+Nl+cols[j],     0.5*dt*(1.0/N2+1.0)*dummy, 		ADD_VALUES);	

				MatSetValue(B_, 	i, 		Nl+cols[j], -0.5*dt*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);

        if (i<Nl)
        {
			// C and C1 now added directly: NMat and VecNMat no longer needed
			/*
			MatGetRow(C, i, &numberOfNonZeros, &cols, &values);
			for (int j=0;j<numberOfNonZeros;++j)
			{
				dummy = (values[j]);
				
				if (dummy!=0.0)
				{
					MatSetValue(P_, 	Nu+Nv+i,     N+cols[j],     -0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	Nu+Nv+i,     N+cols[j],     0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
					
					MatSetValue(P_, 	Nu+Nv+i,     N+Nl+cols[j],     0.5*dt*(1.0/N2+1.0)*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	Nu+Nv+i,     N+Nl+cols[j],     -0.5*dt*(1.0/N2+1.0)*dummy, 	ADD_VALUES);
					
					MatSetValue(P_, 	N+Nl+i,    Nu+Nv+cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	N+Nl+i,    Nu+Nv+cols[j],     0.5*dt*dummy, 	ADD_VALUES);	
				}
			}
			MatRestoreRow(C, i, &numberOfNonZeros, &cols, &values);		

			MatGetRow(C1, i, &numberOfNonZeros, &cols, &values);
			for (int j=0;j<numberOfNonZeros;++j)
			{
				cout << dummy << endl;
				dummy = (values[j]);
				if (dummy!=0.0)
				{
					
					MatSetValue(P_, 	Nu+Nv+i,     N+cols[j],     -0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	Nu+Nv+i,     N+cols[j],     0.5*dt*1.0/N2*dummy, 	ADD_VALUES);
					
					MatSetValue(P_, 	Nu+Nv+i,     N+Nl+cols[j],     0.5*dt*(1.0/N2)*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	Nu+Nv+i,     N+Nl+cols[j],     -0.5*dt*(1.0/N2)*dummy, 	ADD_VALUES);
					
					MatSetValue(P_, 	N+i,    Nu+Nv+cols[j],     0.5*dt*dummy, 	ADD_VALUES);
					MatSetValue(Q_, 	N+i,    Nu+Nv+cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
				}
			}
			MatRestoreRow(C1, i, &numberOfNonZeros, &cols, &values);	
			*/ 
            MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
                dummy = (values[j]);
                if (dummy!=0.0)
                {		
					
					MatSetValue(C_, 	i, 		cols[j], -0.5*dt*dummy, 	ADD_VALUES);
					MatSetValue(C_, 	Nl+i, 		cols[j], -0.5*dt*dummy, 	ADD_VALUES);
                }
            }
            MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);
			
        }
    }
	
    cout << "Mat Assembly for "<< "B_"<<endl;
    MatAssemblyBegin(B_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B_,MAT_FINAL_ASSEMBLY);	
	printFullMatrixInfo(B_, "B_");
    cout << "Mat Assembly for "<< "C_"<<endl;
    MatAssemblyBegin(C_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(C_,MAT_FINAL_ASSEMBLY);	
	printFullMatrixInfo(C_, "C_");	
	
	MatDestroy(&BF1);
	MatDestroy(&DIV);
	
}

void
HEuler::solvenew()
{
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
    HEulerGlobalVariables* globalData = static_cast<HEulerGlobalVariables*>(globalData_);
    
	double dt           = static_cast<HEulerGlobalVariables*>(globalData_)->dt_;
    
    unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nw     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
	unsigned int Nv		= 0;
	unsigned int Nu		= 0;
	if (config->DIM_ == 3)
	{
         Nv = Nw;
	}	
	if (config->DIM_ >= 2)
	{
         Nu = Nw;
	}	
    unsigned int Nl     = Nw;
    unsigned int N      = Nu+Nv+Nw;
    
    double endTime      = config->numOfPeriods_*config->onePeriod_ ;
    int nplotPerPeriod  = config->numOfPeriodsInOnePlotStep_;
    
	double Ham, TotalMass, xMom, yMom, zMom;
	double Kin, Pot, Int;
    Mat& globalDIV      = static_cast<HEulerGlobalVariables*>(globalData_)->DivergenceFreeMatrix_;
	
    cout << "endTime="<<endTime<<endl;
    Vec UExact, Lambda, RHS, RH, UVel; //RHS1,RHO,
    
	Vec URP, RHS1, RHS2, Lambda1, Lambda2, Lambda3, Lambda4;
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Destroy these vectors!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	VecCreateSeq(PETSC_COMM_SELF, Nl+Nl, &URP);
	VecCreateSeq(PETSC_COMM_SELF, Nl+Nl, &RHS1);
	VecCreateSeq(PETSC_COMM_SELF, Nl+Nl, &Lambda1);
	VecCreateSeq(PETSC_COMM_SELF, N, &Lambda2);
	VecCreateSeq(PETSC_COMM_SELF, N, &UVel);
	
	VecCreateSeq(PETSC_COMM_SELF, Nl+Nl, &RHS2);
	VecCreateSeq(PETSC_COMM_SELF, Nl+Nl, &Lambda3);
	VecCreateSeq(PETSC_COMM_SELF, N, &Lambda4);
	
	Vec RHSTEMP;
	VecCreateSeq(PETSC_COMM_SELF, N, &RHSTEMP);
	Vec LambdaTEMP;
	VecCreateSeq(PETSC_COMM_SELF, Nl+Nl, &LambdaTEMP);
	//VecCreateSeq(PETSC_COMM_SELF, N+Nl+Nl, &RHO);
    VecCreateSeq(PETSC_COMM_SELF, N+Nl+Nl, &Lambda);
    VecCreateSeq(PETSC_COMM_SELF, N+Nl+Nl, &RHS);
	//VecCreateSeq(PETSC_COMM_SELF, N+Nl+Nl, &RHS1);
    
	 VecCreateSeq(PETSC_COMM_SELF, N+Nl+Nl, &UExact);
	 VecCreateSeq(PETSC_COMM_SELF, N, &UVel);
	 
	 VecCreateSeq(PETSC_COMM_SELF, Nl, &RH);
	 
    int Num =  config->numOfTimeStepInOnePeriod_*config->numOfPeriods_+2;
	
    Vec Hamiltonian, TotalMassVector, Divergence, xMomentum, yMomentum, zMomentum;
	
	Vec KineticEnergy, PotentialEnergy, InternalEnergy;

	
	VecCreateSeq(PETSC_COMM_SELF, Num, &KineticEnergy);
	VecCreateSeq(PETSC_COMM_SELF, Num, &PotentialEnergy);
	VecCreateSeq(PETSC_COMM_SELF, Num, &InternalEnergy);
	
    VecCreateSeq(PETSC_COMM_SELF, Num, &Hamiltonian);
    VecCreateSeq(PETSC_COMM_SELF, Num, &TotalMassVector);
    VecCreateSeq(PETSC_COMM_SELF, Num, &Divergence);
	VecCreateSeq(PETSC_COMM_SELF, Num, &xMomentum);
	VecCreateSeq(PETSC_COMM_SELF, Num, &yMomentum);
	VecCreateSeq(PETSC_COMM_SELF, Num, &zMomentum);

    double u, v, w, rho, p;
    
    //cout<<"Initializing the solver with velocity..."<<endl;
    
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {
        Base::Element* element= *cit;
        int m = element->getID();
        
        unsigned int pos = m*nb;
        
        
        for (int i =0; i<nb;++i)
        {
            u = element->getData(0,0,i);
            VecSetValue(Lambda, pos+i,     u, INSERT_VALUES);
            VecSetValue(UExact, pos+i,     u, INSERT_VALUES);
			VecSetValue(UVel, pos+i,     u, INSERT_VALUES);
            
            
			if (config->DIM_ >= 2)
			{
				v = element->getData(0,1,i);
				VecSetValue(Lambda, Nu+pos+i,   v, INSERT_VALUES);
				VecSetValue(UExact, Nu+pos+i,   v, INSERT_VALUES);
				VecSetValue(UVel, Nu+pos+i,   v, INSERT_VALUES);
			}
            if (config->DIM_ == 3)
			{
				w = element->getData(0,2,i);
				VecSetValue(UExact, Nu+Nv+pos+i, w, INSERT_VALUES);
				VecSetValue(UVel, Nu+Nv+pos+i, w, INSERT_VALUES);
				VecSetValue(Lambda, Nu+Nv+pos+i, w, INSERT_VALUES);
			}
            
            rho = element->getData(0,3,i);
            VecSetValue(UExact, N+pos+i, rho, INSERT_VALUES);
			VecSetValue(URP, pos+i, rho, INSERT_VALUES);
            VecSetValue(Lambda, N+pos+i, rho, INSERT_VALUES);			
			
            p = element->getData(0,4,i);
            VecSetValue(UExact, N+Nl+pos+i, p, INSERT_VALUES);
			VecSetValue(URP, Nl+pos+i, p, INSERT_VALUES);
            VecSetValue(Lambda, N+Nl+pos+i, p, INSERT_VALUES);

		
        }
    }
    VecAssemblyBegin(UExact);
    VecAssemblyEnd(UExact);
    VecAssemblyBegin(UVel);
    VecAssemblyEnd(UVel);	
    VecAssemblyBegin(URP);
    VecAssemblyEnd(URP);		
    //VecAssemblyBegin(RHO);
    //VecAssemblyEnd(RHO);	
    VecAssemblyBegin(Lambda);
    VecAssemblyEnd(Lambda);		
	
    PetscScalar* XTEMP = new PetscScalar [Nu+Nv+Nw+Nl+Nl];
	PetscScalar* XTEMP1 = new PetscScalar [Nu+Nv+Nw+Nl+Nl];
	PetscScalar* XTEMPRHO = new PetscScalar [Nu+Nv+Nw+Nl+Nl];
	
    PetscScalar* XTEMPVel = new PetscScalar [Nu+Nv+Nw];
	PetscScalar* XTEMPRP = new PetscScalar [Nl+Nl];	
	
	VecGetArray(Lambda, &XTEMP);
	//VecGetArray(RHO, &XTEMPRHO);
	TotalMass = 0.0;
	Ham = 0.0;
	Kin = 0.0;
	Pot = 0.0;
	Int = 0.0;
	xMom = 0.0;
	yMom = 0.0;
	zMom = 0.0;
	double yTerm = 0.0;
	double xTerm = 0.0;
	for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
	{

		Base::Element* element = *cit;
		int k = element->getID();
		unsigned int pos = k*nb;
		//LinearAlgebra::NumericalVector sol;
		//element->getSolution(0,p, sol);
		HEulerElementData* elementData = static_cast<HEulerElementData*>((*cit)->getUserData());

		for (int j =0; j<nb;++j)
		{
			for (int i =0; i<nb;++i)
			{
				if (config->DIM_ == 3)
				{
					yTerm = XTEMP[Nu+Nv+pos+i]*XTEMP[Nu+Nv+pos+j];
					yMom 	= yMom + elementData->massMatrix_(j,i)*XTEMP[Nu+Nv+pos+i];
				}
				if (config->DIM_ >= 2)
				{
					xTerm = XTEMP[Nu+pos+i]*XTEMP[Nu+pos+j];
					xMom 	= xMom + elementData->massMatrix_(j,i)*XTEMP[Nu+pos+i];
				}				
				
				if (config->incompressible_ == 0)
				{ //+yTerm+xTerm
					Kin = Kin + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(XTEMP[pos+i]*XTEMP[pos+j]+yTerm+xTerm);
					// Inlcudes Internal Energy
					Pot = Pot + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(1.0/config->N2_*XTEMP[N+pos+i]*XTEMP[N+pos+j]);
					Int = Int + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(-2.0/config->N2_*XTEMP[N+pos+i]*XTEMP[N+Nl+pos+j]+1.0/config->N2_*XTEMP[N+Nl+pos+i]*XTEMP[N+Nl+pos+j]+XTEMP[N+Nl+pos+i]*XTEMP[N+Nl+pos+j]);
					Ham = Ham + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(XTEMP[pos+i]*XTEMP[pos+j]+yTerm+xTerm+1.0/config->N2_*XTEMP[N+pos+i]*XTEMP[N+pos+j]-2.0/config->N2_*XTEMP[N+pos+i]*XTEMP[N+Nl+pos+j]+1.0/config->N2_*XTEMP[N+Nl+pos+i]*XTEMP[N+Nl+pos+j]+XTEMP[N+Nl+pos+i]*XTEMP[N+Nl+pos+j]);
					//Ham = Ham + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(XTEMP[pos+i]*XTEMP[pos+j]+yTerm+xTerm+XTEMP[N+Nl+pos+i]*XTEMP[N+Nl+pos+j]); 
										
					TotalMass = TotalMass + elementData->massMatrix_(j,i)*XTEMP[N+pos+i];//*XTEMP[N+pos+j]
				}
				else 
				{
					Kin = Kin + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(XTEMP[pos+i]*XTEMP[pos+j]+yTerm+xTerm);
					Pot = Pot + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*1.0/config->N2_*XTEMP[N+pos+i]*XTEMP[N+pos+j];
					Ham = Ham + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(XTEMP[pos+i]*XTEMP[pos+j]+yTerm+xTerm+1.0/config->N2_*XTEMP[N+pos+i]*XTEMP[N+pos+j]);
					TotalMass = TotalMass + elementData->massMatrix_(j,i)*XTEMP[N+pos+i];//*XTEMP[N+pos+j]
				}				
				zMom 	= zMom + elementData->massMatrix_(j,i)*XTEMP[pos+i];
			}
		}
	}
	cout << "Ham = " << std::setprecision (20) << Ham << endl;
	//cout << "TotalMass = " << std::setprecision (20) << TotalMass << endl;
    VecSetValue(TotalMassVector, 0, 	TotalMass, INSERT_VALUES);
    VecSetValue(KineticEnergy,   0, 	Kin, INSERT_VALUES);
    VecSetValue(PotentialEnergy, 0, 	Pot, INSERT_VALUES);	
	VecSetValue(InternalEnergy,  0, 	Int, INSERT_VALUES);	
    VecSetValue(Hamiltonian,     0, 	Ham, INSERT_VALUES);
	VecSetValue(xMomentum,     	 0, 	xMom, INSERT_VALUES);
	VecSetValue(yMomentum,     	 0, 	yMom, INSERT_VALUES);
	VecSetValue(zMomentum,     	 0, 	zMom, INSERT_VALUES);
		
	VecRestoreArray(Lambda, &XTEMP);
	//VecRestoreArray(RHO, &XTEMPRHO);	

	if (config->incompressible_ == 1)
	{
		MatMult(globalDIV, UVel, RH);
		PetscInt a(10);
		PetscReal max;
		VecMax(RH, &a, &max);
		cout << "Divergence = " << std::setprecision (20) << max<<endl;	
		VecSetValue(Divergence, 0, 	max, INSERT_VALUES);
	}
    
    double currTime=0;
    
    int iterations=1;
    
	/*
	IS              isrow,iscol;            
	PetscErrorCode  ierr;
	MatOrderingType rtype = MATORDERINGQMD;//MATORDERINGRCM;
	
	MatGetOrdering(P_,rtype,&isrow,&iscol);
	
	MatPermute(P_,isrow,iscol,&P_);
	MatPermute(Q_,isrow,iscol,&Q_);
	VecPermute(UExact,iscol,PETSC_FALSE);
	*/
	
    KSP ksp;
        // Preconditioner
    PC pc;
   // cout << "Solve"<<endl;
        // Create a solver
    //KSPCreate(PETSC_COMM_SELF, &ksp);
	KSPCreate(PETSC_COMM_WORLD, &ksp);
    
    //cout << "ksp create"<<endl;
    
        //! Direct solver: only a preconditioner is needed. To change
        //this to an iterative method, you can change the solver to
        //KSPMINRES, KSPGMRES, etc. Similarly, you can choose PCNONE
        //(instead of PCLU) for the preconditioner below if you want to
        //use a nonpreconditioned iterative method.
    
    //KSPSetOperators(ksp, P_, P_, SAME_NONZERO_PATTERN);
    
	
	Mat I2, I3;
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,        1, 			PETSC_NULL, &I3);
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv, 	Nu+Nv,        1, 			PETSC_NULL, &I2);
	
	for (int i=0;i<Nu;i++)
	{
		MatSetValue(I2, i,     i,     1.0, ADD_VALUES);
		MatSetValue(I2, Nu+i,     Nu+i,     1.0, ADD_VALUES);
		MatSetValue(I3, i,     i,     1.0, ADD_VALUES);
		MatSetValue(I3, Nu+i,     Nu+i,     1.0, ADD_VALUES);
		MatSetValue(I3, Nu+Nv+i,     Nu+Nv+i,     1.0, ADD_VALUES);
	}
	
	
	
	MatAssemblyBegin(I3,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(I3,MAT_FINAL_ASSEMBLY);
	
	MatAssemblyBegin(I2,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(I2,MAT_FINAL_ASSEMBLY);
	
	Mat S;
	MatCreate(PETSC_COMM_WORLD,&S);
	//MatCreateSchurComplement(I3,I3,B_,C_,I2,&S);
	
	MatMatMult(C_, B_, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &S);
	MatScale(S,-1.0);
	////MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv, 	Nu+Nv,        3*7*nb, 			PETSC_NULL, &S);
	MatShift(S,1.0);
	MatAssemblyBegin(S,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(S,MAT_FINAL_ASSEMBLY);
	//printFullMatrixInfo(S, "S")	;
	
	Mat S1;
	MatCreate(PETSC_COMM_WORLD,&S1);
	MatMatMult(C_, B_, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &S1);
	//MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv, 	Nu+Nv,        3*7*nb, 			PETSC_NULL, &S);
	MatShift(S1,1.0);
	MatAssemblyBegin(S,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(S,MAT_FINAL_ASSEMBLY);
	printFullMatrixInfo(S1, "S1")	;	
	
	//KSPSetOperators(ksp, S, S, SAME_NONZERO_PATTERN);
	KSPSetOperators(ksp, S, S);
	//MatAXPY(I2, S, DIFFERENT_NONZERO_PATTERN);
	
    //cout << "Setup operators"<<endl;
    
    //KSPSetType(ksp, KSPGMRES);
	//KSPSetType(ksp, KSPMINRES);
	//KSPSetType(ksp, KSPGCR);
	//KSPSetType(ksp, KSPCG);
	//KSPSetType(ksp, KSPBCGS);
	//KSPSetType(ksp, KSPFGMRES);
    //cout << "Setup solver"<<endl; 
    //cout << "Setup Initial guess"<<endl;
    
    //KSPGetPC(ksp,&pc);
    //PCFactorSetFill(pc,5);
    //cout << "Setup pc"<<endl;
    
	/*if (config->incompressible_ == 0)
	{    
		KSPSetType(ksp, KSPGMRES);
		KSPGetPC(ksp,&pc);
		PCSetType(pc, PCILU);//PCNONE);
		PCFactorSetLevels(pc,1);
		PCFactorSetReuseOrdering(pc, PETSC_TRUE);
	}
	else if (config->incompressible_ == 1)// && config->stratified_ == 1)
	{	*/
	 
		KSPSetType(ksp, KSPBCGS);
		KSPGetPC(ksp,&pc);	
		//PCSetType(pc, PCSOR);

		//KSPSetType(ksp, KSPGMRES);
		//KSPGetPC(ksp,&pc);	
		//PCSetType(pc, PCILU);
		//PCFactorSetLevels(pc,1);
		////PCFactorSetLevels(pc,10); // doesn't work for WA2D (1,1)
		//PCFactorSetReuseOrdering(pc, PETSC_TRUE);		
		PCSetType(pc, PCNONE);

	//}	
		/*
		KSPSetType(ksp, KSPGMRES);
		KSPGetPC(ksp,&pc);
		PCSetType(pc, PCILU);//PCNONE);
		PCFactorSetLevels(pc,1);
		PCFactorSetReuseOrdering(pc, PETSC_TRUE);	
		*/
	//MatOrderingType rtype = MATORDERINGRCM;//MATORDERINGQMD;//
	//PCFactorSetMatOrderingType(pc,rtype);
	//else if (config->incompressible_ == 1 && config->stratified_ == 0)
	//{	
	//	KSPSetType(ksp, KSPBCGS);
	//	KSPGetPC(ksp,&pc);	
	//	PCSetType(pc, PCSOR);
	//}	
	
	

	//PCSetType(pc, PCJACOBI);
	//PCSetType(pc, PCSOR);
		//PCSetType(pc, PCMG);
		//PCSetType(pc, PCICC);
		//PCSetType(pc, PCLU);
        //
        //PCSetType(pc, PCCHOLESKY);
		
		
		// Use Matlab script to estimate the amount of Fill In
		// Important for optimization since most time is spent in iterating
		
        //PCILUSetFill(pc, 4); // outdated command
		//PCFactorSetFill(pc,5);

    KSPSetFromOptions(ksp);
	KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
        	//PCFactorSetFill(pc, 3.83016);
        	//PCFactorSetLevels(pc,4);
/*    
PCFactorSetLevels(PC pc,int levels);
PCFactorSetReuseOrdering(PC pc,PetscBool flag);
PCFactorSetDropTolerance(PC pc,double dt,double dtcol,int dtcount);
PCFactorSetReuseFill(PC pc,PetscBool flag);
PCFactorSetUseInPlace(PC pc);
PCFactorSetAllowDiagonalFill(PC pc);
*/

    //cout << "Setup options"<<endl;
    
    
    double reltol = 1.0e-14;
    double abstol = 1.0e-14;
	int    maxits = 1e8;//1e5
    KSPSetTolerances(ksp, reltol, abstol, PETSC_DEFAULT, maxits);
    //cout << "Setup tolerance"<<endl;
    
    PCSetUp(pc);
    KSPSetUp(ksp);
    //cout << "Setup ksp"<<endl<<endl;

	/*
     MatInfo           matinfo;
	 double fillrationeeded;
	 MatGetInfo(P_,MAT_LOCAL,&matinfo);
     printf("matinfo.nz_used %g\n",matinfo.nz_used);
	 fillrationeeded  = matinfo.fill_ratio_needed;
	 printf("matinfo.fill_ratio_needed %g\n",matinfo.fill_ratio_needed);
	  * */
    
        //! Solving linear system.
	cout << "n_x = "<< config->nx_ << endl;
	cout << "n_b = "<< configData_->numberOfBasisFunctions_ << endl;
	cout << "n_El = " << static_cast<HEulerGlobalVariables*>(globalData_)->nElements_ << endl;
	cout << "Nu = "<< Nu << endl;	
	cout << "dt = " << dt << endl;
    cout << endl << "Solving System" << endl;
	
    while(currTime<=endTime)//
    {     
		// U^n terms
		//cout << "Step 1 "<< endl;
		//Step 1
		MatMult(C_,UVel, RHS1);
		//cout << "Step 2 "<< endl;
		VecScale(RHS1,2.0);
		//cout << "Step 3 "<< endl;
		//Step 2
		KSPSolve(ksp, RHS1, Lambda1);
		//cout << "Step 4 "<< endl;
		//Step 3
		MatMult(B_, Lambda1, Lambda2);
		//cout << "Step 5 "<< endl;
		//Step 4
		VecAXPY(Lambda2,1.0,UVel);
		//cout << "Step 6 "<< endl;
		//Step 5
		//VecScale(Lambda1,-1.0);
		
		// [R^n,P^n]^T terms
		
		// Step 1
		//MatMatMult(C, B, MATTemp);
		//MatShift(MATTemp, 1.0);
		//MatMult(MATTemp, URP, RHS2);
		////MatMult(B_, URP, RHSTEMP);
		//cout << "Step 7 "<< endl;
		// Computes v3 = v2 + A * v1
		////MatMultAdd(C_,RHSTEMP, URP,RHS2);
		
		//MatMult(C_, RHS2, RHS2);
		//cout << "Step 8 "<< endl;
		//VecAXPY(RHS2,1.0,URP);
		//cout << "Step 9 "<< endl;
		
		MatMult(S1, URP, RHS2);
		
		// Step 2
		KSPSolve(ksp, RHS2, Lambda3);
		//cout << "Step 10 "<< endl;
		// Step 3 n
		//Computes w = alpha x + y. 
		VecWAXPY(LambdaTEMP,1.0,Lambda3,URP);
		//VecAXPY(Lambda3,1.0,URP);
		//cout << "Step 11 "<< endl;
		// Step 4n
		MatMult(B_, LambdaTEMP, Lambda4);
		//cout << "Step 12 "<< endl;
		//// Step 3
		//MatMult(B, Lambda3, Lambda4);
		///// Step 4
		//MatMult(B, URP, RHSTEMP);
		//VecAXPY(Lambda4,1.0,RHSTEMP);
		// Step 5 
		//VecScale(Lambda4,-1.0);
		
		// Combine Results
		//U^n+1 = Lambda1 + Lambda3;
		VecWAXPY(UVel,-1.0,Lambda4,Lambda2);
		//cout << "Step 13 "<< endl;
		//[R^n+1,P^n+1]^T = Lambda2 + Lambda4;
		VecWAXPY(URP,-1.0,Lambda1,Lambda3);
		//cout << "Step 14 "<< endl;
		
		
		
		
        //cout<<"construct proper rightHandside"<<endl;
        //MatMult(Q_, UExact, RHS);
		/*
		if (config->incompressible_ == 0)
		{
			VecAssemblyBegin(RHS);
			VecAssemblyEnd(RHS);			
			MatMult(R_, RHO, RHS1);
			VecAssemblyBegin(RHS1);
			VecAssemblyEnd(RHS1);
			VecAXPY(RHS,1.0,RHS1);
		}
		 * */
		
        
            ///***********************************************************************************
            ///***************************construct proper rightHandside**************************
        
        
        //cout << "Finalizing vector creation"<<endl;
        
        
            // KSP Solver initializer
            // 	    outputVectorMatlab(RHS, "RHS.dat");
        //KSPSolve(ksp, RHS, Lambda);
            //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
		if (iterations%nplotPerPeriod==0)
		{
        if(reason < 0)
        {
            PetscInt its;
            KSPGetIterationNumber(ksp, &its);
            PetscReal rnorm;
            KSPGetResidualNorm(ksp, &rnorm);
            cout << "\t\tPetSc: Solving the linear system has failed. Reason code: "
            << reason << endl //<< "Check KSPConvergedReason for the reason" << endl
            << "\t\tPetSc: Residual after " << int(its) << " iterations : ";
            cout.setf(ios::scientific, ios::floatfield);
            cout.precision(4);
            cout << rnorm << endl;
            cout.setf(ios::fixed, ios::floatfield);
            cout.precision(5);
        }
        
        else
        {
            if(reason > 0)
            {
				
                PetscInt its;
                KSPGetIterationNumber(ksp, &its);
                PetscReal rnorm;
                KSPGetResidualNorm(ksp, &rnorm);
                cout << "\t\tPetSc: Solving the linear system has succeeded. Reason code: "
                << reason << endl //<< "Check KSPConvergedReason for the reason" << endl
                << "\t\tPetsc: Convergence in " << int(its) << " iterations : ";
                cout.setf(ios::scientific, ios::floatfield);
                cout.precision(4);
                cout << rnorm << endl;
                cout.setf(ios::fixed, ios::floatfield);
                cout.precision(5);
				
            }
            else
            {
                cout << "\t\tPetSc: Solving the linear system is still under way" << endl;
            }
        }
		}
        
        //cout << "Solved a timestep" << endl;
        
        
        //VecGetArray(Lambda, &XTEMP);
		//VecGetArray(UExact, &XTEMP1);
		//VecGetArray(RHO, &XTEMPRHO);
		
		VecGetArray(UVel, &XTEMPVel);
		VecGetArray(URP, &XTEMPRP);
		
        int pos=0;
        int k=0;
        currTime+=dt;
        //for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
        //{
        //    Base::Element* element= *cit;
        //    int k = element->getID();
        //    
        //    unsigned int pos = k*nb;
        //    
        //    for (int i =0; i<nb;++i)
        //    {
		//		//element->setData(0,3,i, XTEMPRHO[Nu+Nv+pos+i]+0.5*dt*(XTEMP[Nu+Nv+pos+i]+XTEMP1[Nu+Nv+pos+i])); //set rho
		//	}
		//}
		//VecRestoreArray(UExact, &XTEMP1);
		//VecRestoreArray(RHO, &XTEMPRHO);
		/*
        pos=0;
        k=0;
        for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
        {
            Base::Element* element= *cit;
            int k = element->getID();
            
            unsigned int pos = k*nb;
            
            for (int i =0; i<nb;++i)
            {	
		
                element->setData(0,3,i, XTEMP[Nu+Nv+Nw+pos+i]); //set rho
				VecSetValue(UExact, Nu+Nv+Nw+pos+i, XTEMP[Nu+Nv+Nw+pos+i], INSERT_VALUES);
                element->setData(0,4,i, XTEMP[Nu+Nv+Nw+Nl+pos+i]); //set p
				VecSetValue(UExact, Nu+Nv+Nw+Nl+pos+i, XTEMP[Nu+Nv+Nw+Nl+pos+i], INSERT_VALUES);	
	
                element->setData(0,0,i, XTEMP[pos+i]);//set U
                VecSetValue(UExact, pos+i,     		XTEMP[pos+i], INSERT_VALUES);
				VecSetValue(UVel, 	pos+i,     		XTEMP[pos+i], INSERT_VALUES);
                
				if (config->DIM_ >= 2)
				{
					element->setData(0,1,i, XTEMP[Nu+pos+i]);//set V
					VecSetValue(UExact, Nu+pos+i,   	XTEMP[Nu+pos+i], INSERT_VALUES);
					VecSetValue(UVel, 	Nu+pos+i,   	XTEMP[Nu+pos+i], INSERT_VALUES);
				}
				if (config->DIM_ == 3)
				{                
					element->setData(0,2,i, XTEMP[Nu+Nv+pos+i]);//set W
					VecSetValue(UExact, Nu+Nv+pos+i, 	XTEMP[Nu+Nv+pos+i], INSERT_VALUES);
					VecSetValue(UVel, 	Nu+Nv+pos+i, 	XTEMP[Nu+Nv+pos+i], INSERT_VALUES);
				}
				
				//rho = element->getData(0,4,i);
				//VecSetValue(RHO, Nu+Nv+pos+i, 	rho, INSERT_VALUES);
				//VecSetValue(RHO, N+pos+i, 		rho, INSERT_VALUES);				
            }
        }
        */
        //VecRestoreArray(Lambda, &XTEMP);
        //VecAssemblyBegin(UExact);
		//VecAssemblyEnd(UExact);	
		//VecAssemblyBegin(UVel);
        //VecAssemblyEnd(UVel);
        
        	
        //VecAssemblyBegin(RHO);
        //VecAssemblyEnd(RHO);	

        //VecGetArray(Lambda, &XTEMP);
		//VecGetArray(RHO, &XTEMPRHO);

		TotalMass = 0.0;
		Ham = 0.0;
		Kin = 0.0;
		Pot = 0.0;
		Int = 0.0;
		xMom = 0.0;
		yMom = 0.0;
		zMom = 0.0;
		yTerm = 0.0;
		xTerm = 0.0;
        for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
        {

            Base::Element* element = *cit;
            int k = element->getID();
            unsigned int pos = k*nb;
            //LinearAlgebra::NumericalVector sol;
            //element->getSolution(0,p, sol);
            HEulerElementData* elementData = static_cast<HEulerElementData*>((*cit)->getUserData());

            for (int j =0; j<nb;++j)
            {
                for (int i =0; i<nb;++i)
                {
					if (config->DIM_ == 3)
					{
						yTerm = XTEMPVel[Nu+Nv+pos+i]*XTEMPVel[Nu+Nv+pos+j];
						yMom 	= yMom + elementData->massMatrix_(j,i)*XTEMPVel[Nu+Nv+pos+i];
					}
					if (config->DIM_ >= 2)
					{
						xTerm = XTEMPVel[Nu+pos+i]*XTEMPVel[Nu+pos+j];
						xMom 	= xMom + elementData->massMatrix_(j,i)*XTEMPVel[Nu+pos+i];
					}				
					
					if (config->incompressible_ == 0)
					{
						Kin = Kin + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(XTEMPVel[pos+i]*XTEMPVel[pos+j]+yTerm+xTerm);
						Pot = Pot + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(1.0/config->N2_*XTEMPRP[pos+i]*XTEMPRP[pos+j]);
						Int = Int + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(-2.0/config->N2_*XTEMPRP[pos+i]*XTEMPRP[Nl+pos+j]+1.0/config->N2_*XTEMPRP[Nl+pos+i]*XTEMPRP[Nl+pos+j]+XTEMPRP[Nl+pos+i]*XTEMPRP[Nl+pos+j]);						
						Ham = Ham + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(XTEMPVel[pos+i]*XTEMPVel[pos+j]+yTerm+xTerm+1.0/config->N2_*XTEMPRP[pos+i]*XTEMPRP[pos+j]-2.0/config->N2_*XTEMPRP[pos+i]*XTEMPRP[Nl+pos+j]+1.0/config->N2_*XTEMPRP[Nl+pos+i]*XTEMPRP[Nl+pos+j]+XTEMPRP[Nl+pos+i]*XTEMPRP[Nl+pos+j]);
						TotalMass = TotalMass + elementData->massMatrix_(j,i)*XTEMPRP[pos+i];
					}
					else 
					{
						//Kin = Kin + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(XTEMP[pos+i]*XTEMP[pos+j]+yTerm+xTerm);
						//Pot = Pot + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(1.0/config->N2_*XTEMP[N+pos+i]*XTEMP[N+pos+j]);		
						//Ham = Ham + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(XTEMP[pos+i]*XTEMP[pos+j]+yTerm+xTerm+1.0/config->N2_*XTEMP[N+pos+i]*XTEMP[N+pos+j]);
						//TotalMass = TotalMass + elementData->massMatrix_(j,i)*XTEMP[N+pos+i];//*XTEMP[N+pos+j]
					}				
					zMom 	= zMom + elementData->massMatrix_(j,i)*XTEMPVel[pos+i];				
                }
            }
        }
        //cout << "Ham = " << std::setprecision (20) << Ham << endl;
        //cout << "TotalMass = " << std::setprecision (20) << TotalMass << endl;
        //cout << "iterations = " << iterations << endl;
		
        VecSetValue(TotalMassVector, iterations, 	TotalMass, INSERT_VALUES);
        VecSetValue(Hamiltonian,     iterations, 	Ham, INSERT_VALUES);
		VecSetValue(KineticEnergy,   iterations, 	Kin, INSERT_VALUES);
		VecSetValue(InternalEnergy,  iterations, 	Int, INSERT_VALUES);
		VecSetValue(PotentialEnergy, iterations, 	Pot, INSERT_VALUES);
		VecSetValue(xMomentum,     	 iterations, 	xMom, INSERT_VALUES);
		VecSetValue(yMomentum,     	 iterations, 	yMom, INSERT_VALUES);
		VecSetValue(zMomentum,     	 iterations, 	zMom, INSERT_VALUES);		
		
        //VecRestoreArray(Lambda, &XTEMP);
		//VecRestoreArray(RHO, &XTEMPRHO);	
		VecRestoreArray(UVel, &XTEMPVel);
		VecRestoreArray(URP, &XTEMPRP);		


	
	if (config->incompressible_ == 1)
	{
	    MatMult(globalDIV, UVel, RH);
        PetscInt a(10);
        PetscReal max;
        VecMax(RH, &a, &max);
        //cout << "Divergence= " << std::setprecision (20) << max<<endl;
		VecSetValue(Divergence, iterations, 	max, INSERT_VALUES);
	}
	
        if (iterations%nplotPerPeriod==0)
			{
				output(currTime);
				cout<<"currTime="<<currTime<<endl;
				cout << "Ham = " << std::setprecision (20) << Ham << endl; 
				if (config->incompressible_ == 1)
				{
					MatMult(globalDIV, UVel, RH);
					PetscInt a(10);
					PetscReal max;
					VecMax(RH, &a, &max);
					cout << "Divergence= " << std::setprecision (20) << max<<endl;
				}
				cout << endl;
			}
			
        
        iterations++;
    }
	//cout << "Ham = " << std::setprecision (20) << Ham << endl;
    
    KSPDestroy(&ksp);
    VecDestroy(&UExact);
	VecDestroy(&UVel);
    VecDestroy(&RHS);
    VecDestroy(&Lambda);
    VecDestroy(&RH);
	//VecDestroy(&RHS1);
	//VecDestroy(&RHO);
	
    VecDestroy(&URP);
	VecDestroy(&RHS1);
	VecDestroy(&RHS2);
	VecDestroy(&RHSTEMP);
	VecDestroy(&Lambda1);
	VecDestroy(&Lambda2);
	VecDestroy(&Lambda3);
	VecDestroy(&Lambda4);
	
	MatDestroy(&I2);
	MatDestroy(&I3);
	MatDestroy(&S);
	
    cout << "The end of the time loop"<<endl;
    cout << "iterations =  " << iterations << endl;
	cout << "Nu = "<<Nu<<", nb = "<<nb<<endl;
    cout << "ne = "<<Nu/nb<<endl;	
	
	PetscScalar* HAMTEMP = new PetscScalar [Num];
	PetscScalar* KINTEMP = new PetscScalar [Num];
	PetscScalar* POTTEMP = new PetscScalar [Num];
	PetscScalar* INTTEMP = new PetscScalar [Num];
    PetscScalar* MASSTEMP = new PetscScalar [Num];
    PetscScalar* DIVTEMP = new PetscScalar [Num];
	PetscScalar* xMomTEMP = new PetscScalar [Num];
	PetscScalar* yMomTEMP = new PetscScalar [Num];
	PetscScalar* zMomTEMP = new PetscScalar [Num];
    VecGetArray(Hamiltonian, &HAMTEMP);
	VecGetArray(KineticEnergy, &KINTEMP);
	VecGetArray(PotentialEnergy, &POTTEMP);
	VecGetArray(InternalEnergy, &INTTEMP);
    VecGetArray(TotalMassVector, &MASSTEMP);
    VecGetArray(Divergence, &DIVTEMP);
	VecGetArray(xMomentum, &xMomTEMP);
	VecGetArray(yMomentum, &yMomTEMP);
	VecGetArray(zMomentum, &zMomTEMP);
    
	/*for (int i = 0; i<Num-1;++i)
    {
        cout  << std::setprecision (20) << HAMTEMP[i] << endl;
    }*/
    currTime = 0.0;
    ofstream fout("HAM.txt");
	
	fout << "Time" << " " << "Kinetic Energy" << " " << "Potential Energy" << " " << "Internal Energy" << " " << "Total Energy" << " " << "Total Mass" << " " << "Total X-Momentum" << " " << "Total Y-Momentum" << " " << "Total Z-Momentum" << " " << "Maximum Divergence" <<endl; 
        
    for (int i =0; i<Num-1;++i)
    {
        fout << currTime << " " << std::setprecision (20) << KINTEMP[i] << " " << std::setprecision (20) << POTTEMP[i] << " " << std::setprecision (20) << INTTEMP[i] << " " << std::setprecision (20) << HAMTEMP[i] << " " << std::setprecision (20) << MASSTEMP[i] << " " << std::setprecision (20) << xMomTEMP[i] << " " << std::setprecision (20) << yMomTEMP[i] << " " << std::setprecision (20) << zMomTEMP[i] << " " << std::setprecision (20) << DIVTEMP[i] <<endl; 
        currTime +=dt;
    }
    fout.close();
    VecRestoreArray(Hamiltonian, &HAMTEMP);
	VecRestoreArray(KineticEnergy, &KINTEMP);
	VecRestoreArray(PotentialEnergy, &POTTEMP);	
	VecRestoreArray(InternalEnergy, &INTTEMP);	
    VecRestoreArray(TotalMassVector, &MASSTEMP);
    VecRestoreArray(Divergence, &DIVTEMP);
	VecRestoreArray(xMomentum, &xMomTEMP);
	VecRestoreArray(yMomentum, &yMomTEMP);
	VecRestoreArray(zMomentum, &zMomTEMP);	

    VecDestroy(&Hamiltonian);
	VecDestroy(&KineticEnergy);
	VecDestroy(&PotentialEnergy);
	VecDestroy(&InternalEnergy);
    VecDestroy(&TotalMassVector);
    VecDestroy(&Divergence);	
	VecDestroy(&xMomentum);	
	VecDestroy(&yMomentum);	
	VecDestroy(&zMomentum);	
    
}

void
HEuler::output(double time)
{
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
	
    string dxName  = "x.dat";
    string dName   = "u.dat";
    string dExName = "uEx.dat";
    char outF[20];
    ofstream    	l2ErrorFile("l2error.dat");
    ofstream    	energyFile("energy.dat");
    ofstream    	divFile("div.dat");
    ofstream    	energyExfile("energyEx.dat");
    
    string          outputFormat("u,uExact, uError, v, vExact, vError, w, wExact, wError, rho, rhoExact, rhoError, p, pExact, pError, rho0, rho0Exact, rho0Error, N2, N2Exact, N2Error, Energy");
    
	if (config->solutionType_==HEulerConfigurationData::WA2D)
	{
		sprintf(outF,"outputTec%f.dat",time);
	}
	else
	{
		if (time==0)
		{
			sprintf(outF,"outputTec%f.dat",  time);
		}
		else
		{
		   sprintf(outF,"outputTec%f.dat",1.0);
		}
	}
    
    string          outFile(outF);
    
    ofstream    file1Dx0(dxName.c_str());
    ofstream    file1D0(dName.c_str());
    ofstream    file1DEx0(dExName.c_str());
    std::ofstream file3D;
    file3D.open (outFile.c_str());
    
	// Need to interface with mesh generation: not always "RectangularMesh"
    Output::TecplotDiscontinuousSolutionWriter out(file3D,"RectangularMesh","012",outputFormat);
    
    
    TecplotWriteFunction userWriteFunction(&energyFile, &divFile, &energyExfile, exactSolution_, &file1Dx0, &file1D0,&file1DEx0, &l2ErrorFile);
	
	//time = 0.0;
    userWriteFunction.time_=time;//time;
	
	//userWriteFunction.outputErrors();
    
    cout <<"OUTPUTING...."<<endl;
    ostringstream ostr;
    ostr << "t = " << time;
    out.write(meshes_[0],std::string(ostr.str()),false, &userWriteFunction);
    cout <<"FINISH OUTPUTING...."<<endl;
 
    file1Dx0.close();
    file1D0.close();
    file1DEx0.close();
    divFile.close();
    energyFile.close();
    energyExfile.close();
    l2ErrorFile.close();
    
    file3D.close();

    cout <<"FINISH...."<<endl;
}

void HEuler::calculateEnergy(Vec Lambda,double& Ham, double& Kin, double& Pot, double&  Int, double& xMom, double& yMom, double& zMom, double& TotalMass)
{
	// takes the vector of coefficients (Lambda = (U, V, W, R, P)^T 
	// calculates the total energy, kinetic energy, potential energy, internal energy, momentum in x-direction, momentum in y-direction, momentumin z-direction and mass
	Ham = 0.0; TotalMass = 0.0; xMom = 0.0; yMom = 0.0; zMom = 0.0; Kin = 0.0; Pot = 0.0; Int = 0.0;
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
    HEulerGlobalVariables* globalData = static_cast<HEulerGlobalVariables*>(globalData_);
    
    unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nw     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;	
	unsigned int Nl     = Nw;	
	unsigned int Nv		= 0;
	unsigned int Nu		= 0;
	if (config->DIM_ == 3)
	{
         Nv = Nw;
	}	
	if (config->DIM_ >= 2)
	{
         Nu = Nw;
	}	
	unsigned int N = Nu+Nv+Nw;

	//double Ham = 0.0, TotalMass = 0.0, xMom = 0.0, yMom = 0.0, zMom = 0.0, Kin = 0.0, Pot = 0.0, Int = 0.0,
	double xTerm = 0.0, yTerm = 0.0, zTerm =0.0;
	double massji = 0.0;
    //Mat& globalDIV      = static_cast<HEulerGlobalVariables*>(globalData_)->DivergenceFreeMatrix_;	
	
	PetscScalar* XTEMP = new PetscScalar [Nu+Nv+Nw+Nl+Nl];
	VecGetArray(Lambda, &XTEMP);
	
	for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
	{

		Base::Element* element = *cit;
		int k = element->getID();
		unsigned int pos = k*nb;
		HEulerElementData* elementData = static_cast<HEulerElementData*>((*cit)->getUserData());

		for (int j =0; j<nb;++j)
		{
			for (int i =0; i<nb;++i)
			{
				massji = elementData->massMatrix_(j,i);
				
				xMom = xMom + massji*XTEMP[pos+i];
				xTerm = XTEMP[pos+i]*XTEMP[pos+j];
				if (config->DIM_ >= 2)
				{
					yMom = yMom + massji*XTEMP[Nu+pos+i];
					yTerm = XTEMP[Nu+pos+i]*XTEMP[Nu+pos+j];
				}
				if (config->DIM_ == 3)
				{
					zMom 	= zMom + massji*XTEMP[Nu+Nv+pos+i];
					zTerm = XTEMP[Nu+Nv+pos+i]*XTEMP[Nu+Nv+pos+j];
				}

				
				
				if (config->incompressible_ == 0) // compressible case
				{
					Pot = Pot - 0.5*elementData->massMatrixOverBackgroundDensityDeriv_(j,i)*(XTEMP[N+pos+i]*XTEMP[N+pos+j]);
					Kin = Kin + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(xTerm+yTerm+zTerm);
					// ~~~ THIS IS TO BE UPDATED ~~~ \\  
					// Only works for constant N^2
					// see below, the incompressible case for varying N^2
					Pot = Pot + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(1.0/config->N2_*XTEMP[N+pos+i]*XTEMP[N+pos+j]);
					Int = Int + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(-2.0/config->N2_*XTEMP[N+pos+i]*XTEMP[N+Nl+pos+j]+1.0/config->N2_*XTEMP[N+Nl+pos+i]*XTEMP[N+Nl+pos+j]+XTEMP[N+Nl+pos+i]*XTEMP[N+Nl+pos+j]);	
				}
				else							  
				{
					if (config->EB_ == 0)		// incompressible case
					{				
						//Pot = Pot - 0.5*elementData->massMatrixOverN2_(j,i)*(XTEMP[N+pos+i]*XTEMP[N+pos+j]); //massMatrixOverBackgroundDensityDeriv_
						Pot = Pot - 0.5*elementData->massMatrixOverBackgroundDensityDeriv_(j,i)*(XTEMP[N+pos+i]*XTEMP[N+pos+j]); //
						Kin = Kin + 0.5*elementData->massMatrixOverBackgroundDensity_(j,i)*(xTerm+yTerm+zTerm);
					}
						
					else						// EB case
					{
						Pot = Pot + 0.5*elementData->massMatrixOverN2_(j,i)*(XTEMP[N+pos+i]*XTEMP[N+pos+j]);
						Kin = Kin + 0.5*massji*(xTerm+yTerm+zTerm);
					}
					
				}
				
				TotalMass = TotalMass + massji*XTEMP[N+pos+i];
			}
		}
	}			
	Ham = Pot + Kin + Int;
	
	VecRestoreArray(Lambda, &XTEMP);
}
