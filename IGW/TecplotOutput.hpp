//
//  TecplotOutput.h
//  
//
//  Created by Shavarsh Nurijanyan on 7/4/13.
//
//

#ifndef TECPLOTOUTPUT_HH
#define TECPLOTOUTPUT_HH
    // To declare Physical and reference space points
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PointReference.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "Output/TecplotSingleElementWriter.hpp"
#include "InitialConditions.hpp"

using Geometry::PointPhysical;
using Geometry::PointReference;
#include <iomanip>
#include <iostream>
#include <fstream>
using std::ofstream;
using std::ifstream;
using std::istream;
using std::ostream;


class TecplotWriteFunction:public Output::TecplotSingleElementWriter
{
	
    
        // 	double exactSolutionP(const PhysSpacePoint<dim>& pPhys)
        // 	{
        // //std::cout<<"currentTime"<<CurrentTime<<endl;
        // 		return 2*cos(2*pi*(pPhys[0] + pPhys[1] + pPhys[2]) + sqrt(3)*currTime/3)/(4*pi*pi);
        // 	}
        //
        // 	double exactSolutionU(const PhysSpacePoint<dim>& pPhys)
        // 	{
        // 		return (sqrt(3)*cos(2*pi*(pPhys[0] + pPhys[1] + pPhys[2]) + sqrt(3)*currTime/3)+3*sin(2*pi*(pPhys[0] + pPhys[1] + pPhys[2])  + sqrt(3)*currTime/3))/pi/2;
        // 	}
        //
        // 	double exactSolutionV(const PhysSpacePoint<dim>& pPhys)
        // 	{
        // 		return (sqrt(3)*cos(2*pi*(pPhys[0] + pPhys[1] + pPhys[2]) + sqrt(3)*currTime/3)-3*sin(2*pi*(pPhys[0] + pPhys[1] + pPhys[2])  + sqrt(3)*currTime/3))/pi/2;
        // 	}
        //
        // 	double exactSolutionW(const PhysSpacePoint<dim>& pPhys)
        // 	{
        // 		return -(2*sqrt(3)*cos(2*pi*(pPhys[0] + pPhys[1] + pPhys[2]) + sqrt(3)*currTime/3))/pi/2;
        // 	}
	
	
    
public:
    
    TecplotWriteFunction(ofstream* energyFile, ofstream* divFile, ofstream* energyExfile, const ExactSolutionBase* vel, ofstream* file1Dx, ofstream* file1D, ofstream* file1DEx, ofstream* l2File):
    velocity_(vel),
    maxError_(0),
    maxErrorU_(0),
    file1Dx_(file1Dx),
    file1D_(file1D),
    file1DEx_(file1DEx),
    energyFile_(energyFile),
    divFile_(divFile),
    energyExfile_(energyExfile),
    l2errorfile_(l2File),
    maxErrorPoint_(3)
    {
        *energyFile_ 	<< "Time"<<" , "<< "Energy" <<endl;
        *energyExfile_ 	<< "Time"<<" , "<< "EnergyEx" <<endl;
        *divFile_ 		<< "Time"<<" , "<< "DIV" <<endl;
    }
    void writeToTecplotFile(const Base::Element* element, const PointReference& pRef, ostream& os)
    {
		double u, v, w, rho, p, rho0, N2, Energy, uExact, vExact, wExact, rhoExact, pExact, rho0Exact, N2Exact;

        LinearAlgebra::NumericalVector sol;
        element->getSolution(0,pRef, sol);
	
		u = sol[0];
		v = sol[1];
		w = sol[2];
		rho = sol[3];
		p = sol[4];		

		PointPhysical pPhys(1);
		element->referenceToPhysical(pRef, pPhys);
		
		uExact		= velocity_->getU(pPhys, time_);
		vExact		= velocity_->getV(pPhys, time_);
		wExact      = velocity_->getW(pPhys, time_);
		rhoExact 	= velocity_->getR(pPhys, time_);
		pExact		= velocity_->getP(pPhys, time_);
		
		rho0 = 0.0;
		N2 = 0.0;
		rho0Exact 	= velocity_->getrho0(pPhys, time_);
		N2Exact 	= velocity_->getN2(pPhys, time_);	

		Energy = 0.0;//0.5*(u*u+v*v+w*w+rho*rho);
		
		os << u<<"\t"<<uExact<<"\t"<<(u-uExact)<<"\t"<<v<<"\t"<<vExact<<"\t"<<(v-vExact)<<"\t"<<w<<"\t"<<wExact<<"\t"<<(w-wExact)<<"\t"<<  rho << "\t"<< rhoExact << "\t" << (rho-rhoExact)<<"\t"<<  p << "\t"<< pExact << "\t" << (p-pExact)<<"\t"<<  rho0 << "\t"<< rho0Exact << "\t" << (rho0-rho0Exact)<<"\t"<<  N2 << "\t"<< N2Exact << "\t" << (N2-N2Exact)<<"\t"<<Energy;		
		
		//VecSetValue(UError, i,     (u-uExact), INSERT_VALUES);
		//i = i+1;
/*
        double lambda, u, uExact, lambdaExact, v, w, energy, vExact, wExact;
        double energyExact, uError, vError, wError;        
        lambda      = sol[3];
        u           = sol[0];
        v           = sol[1];
        w           = sol[2];
//        uError	    = data.getUError();
//        vError	    = data.getVError();
//        wError	    = data.getWError();
       energy      = 0;
//        energyExact = 0;//data.energyExact_;
//        
        PointPhysical pPhys(3);      // point in physical space
        element->referenceToPhysical(pRef, pPhys);
//        
        uExact      = velocity_->getU(pPhys, time_);//exactSolutionU(pPhys);
        vExact      = velocity_->getV(pPhys, time_);//exactSolutionV(pPhys);
        wExact      = velocity_->getW(pPhys, time_);//exactSolutionW(pPhys);/
        lambdaExact = velocity_->getP(pPhys, time_);//exactSolutionP(pPhys);//
//
//
         os << u<<"\t"<<uExact<<"\t"<<(u-uExact)<<"\t"<<v<<"\t"<<vExact<<"\t"<<(v-vExact)<<"\t"<<w<<"\t"<<wExact<<"\t"<<(w-wExact)<<"\t"<<  lambda << "\t"<< lambdaExact << "\t" << (lambda-lambdaExact)<<"\t"<<energy;
*/		  
//        os << u<<"\t"<<uExact<<"\t"<<uError<<"\t"<<v<<"\t"<<vExact<<"\t"<<vError<<"\t"<<w<<"\t"<<wExact<<"\t"<<wError<<"\t"<<  lambda << "\t"<< lambdaExact << "\t" << (lambda-lambdaExact)<<"\t"<<energy;


		/*PointPhysical pPhys(3);
		element->referenceToPhysical(pRef, pPhys);
	
		if (maxError_< std::abs(p-pExact))
		{
			maxError_ = std::abs(p-pExact);
			maxErrorPoint_= pPhys;
		}*/
        
		//if (data.isInternal())
		//{
		/*
			if (maxErrorU_<std::abs(u-uExact))
				maxErrorU_ = std::abs(u-uExact);
			if (maxErrorV_<std::abs(v-vExact))
				maxErrorV_ = std::abs(v-vExact);
			if (maxErrorW_<std::abs(w-wExact))
				maxErrorW_ = std::abs(w-wExact);
			if (maxErrorR_<std::abs(rho-rhoExact))
				maxErrorR_ = std::abs(rho-rhoExact);	
			if (maxErrorP_<std::abs(p-pExact))
				maxErrorP_ = std::abs(p-pExact);	
		*/			
		//}
    }
    void writeDiv(double max)
    {
        *divFile_<<time_<<" "<<max<<endl;
        cout << "Max error in discrete divergence="<< max <<endl;
    }
    void outputErrors()
    {
        //cout << "Maximum error of Pressure is ="<<maxError_<<" at point ("<<maxErrorPoint_[0]<<", "<<maxErrorPoint_[1]<<", "<<maxErrorPoint_[2]<<endl;
        cout << "Maximum error of U is = "<< maxErrorU_<<endl;
        cout << "Maximum error of V is = "<< maxErrorV_<<endl;
		cout << "Maximum error of W is = "<< maxErrorW_<<endl;
		cout << "Maximum error of R is = "<< maxErrorR_<<endl;
		cout << "Maximum error of P is = "<< maxErrorP_<<endl;
			//maxError_=maxErrorU_=maxErrorV_=maxErrorW_=0;
    }
    
    void writeEnergy()
    {
//        double totalEnergy=0;
//        double totalEnergyEx=0;
//        double totalPressure=0;
//        double totalErrorU=0;
//        double totalErrorV=0;
//        double totalErrorW=0;
//        double totalErrorP=0;
//        for (ElItType it = mesh.pubElCo.begin(); it != mesh.pubElCo.end(); ++it)
//        {
//            ElementData& elData = data[(*it)->id()];
//            
//            totalEnergy 	+= 	elData.getEnergy();
//            totalEnergyEx 	+= 	elData.getExactEnergy();
//            totalPressure 	+= 	elData.getPressureIntegral();
//			
//                //				if (elData.internal_)
//            {
//                totalErrorU		+=	elData.getUError();
//                totalErrorV		+=	elData.getVError();
//                totalErrorW		+=	elData.getWError();
//                totalErrorP     += 	elData.getPError();
//            }
//        }
//        
//        *energyFile_ << std::setprecision (6)<<time_ <<" "<<totalEnergy <<endl;
//        *energyExfile_<< std::setprecision (6) <<time_ <<" "<< totalPressure <<endl;
//        *l2errorfile_<< std::setprecision (4)<<time_<<"\t\t"<< std::setprecision (10)<< sqrt(totalErrorU) <<"\t\t"<< sqrt(totalErrorV)<<"\t\t"<< sqrt(totalErrorW)<<"\t\t"<< sqrt(totalErrorP)<<endl;
    }
    
public:
	double 									time_;
	double 									maxError_;
	double 									maxErrorU_;
	double 									maxErrorV_;
	double 									maxErrorW_;
	double 									maxErrorR_;
	double 									maxErrorP_;
    
	PointPhysical						maxErrorPoint_;
	
	Vec										UError;
    //int										i = 0;
	
	void setOutputFiles(ofstream* file1Dx, ofstream* file1D, ofstream* file1DEx){file1Dx_=file1Dx;file1D_=file1D;file1DEx_=file1DEx;}
    
private:
	const ExactSolutionBase*  		     velocity_;
	ofstream*                                file1Dx_;
	ofstream*                                file1D_;
	ofstream*                                file1DEx_;
	
	ofstream*                                energyFile_;
	ofstream*                                divFile_;
	ofstream*                                energyExfile_;
	ofstream*                                l2errorfile_;
};

#endif

