#ifndef HEULER_HH
#define HEULER_HH

#define PR(x) cout << #x " = " << x << "\n";
#define PRp(x) cout << #x " = " << std::setprecision (20) << x << "\n";


#include <string>
using std::string;

//#include <petscmat.h>
#include <petscksp.h>


#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "Base/GlobalData.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/HpgemUI.hpp"
#include "Base/Norm2.hpp"
#include "Base/PhysGradientOfBasisFunction.hpp"
#include "Base/UserData.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Integration/FaceIntegral.hpp"

#include "InitialConditions.hpp"
#include "TecplotOutput.hpp"
#include "Integration/ElementIntegrandBase.hpp"
#include "Integration/FaceIntegrandBase.hpp"
#include "Base/FaceCacheData.hpp"
#include "Base/ElementCacheData.hpp"

#include "Base/MeshManipulator.hpp"

using Base::RectangularMeshDescriptor;
using Base::HpgemUI;
using Base::GlobalData;
using Base::ConfigurationData;
using namespace Base;

#include "airy.hpp"


typedef std::vector<LinearAlgebra::Matrix>   VectorOfMatrices;
struct ElementIntegralData
{
        //optimize later!
    ElementIntegralData operator*= (const double& scalar){xGrad_*=scalar; yGrad_*=scalar; zGrad_*=scalar; return *this;}
    void axpy(double a, const ElementIntegralData& x){ xGrad_.axpy(a, x.xGrad_); yGrad_.axpy(a, x.yGrad_); zGrad_.axpy(a, x.zGrad_);}
    
    LinearAlgebra::Matrix xGrad_;
    LinearAlgebra::Matrix yGrad_;
    LinearAlgebra::Matrix zGrad_;
};

struct ElementIntegralBackgroundDensity
{
    ElementIntegralBackgroundDensity operator*= (const double& scalar){massMatrixR_*=scalar; massMatrixOverBackgroundDensityR_*=scalar; massMatrixTimesBackgroundDensityR_*=scalar; massMatrixTimesBackgroundDensityDerivR_*=scalar; massMatrixOverBackgroundDensityDerivR_*=scalar; massMatrixTimesN2R_*=scalar;  massMatrixOverN2R_*=scalar; return *this;}
    void axpy(double a, const ElementIntegralBackgroundDensity& x){ massMatrixR_.axpy(a, x.massMatrixR_); massMatrixOverBackgroundDensityR_.axpy(a, x.massMatrixOverBackgroundDensityR_); massMatrixTimesBackgroundDensityR_.axpy(a, x.massMatrixTimesBackgroundDensityR_); massMatrixTimesBackgroundDensityDerivR_.axpy(a, x.massMatrixTimesBackgroundDensityDerivR_); massMatrixOverBackgroundDensityDerivR_.axpy(a, x.massMatrixOverBackgroundDensityDerivR_); massMatrixTimesN2R_.axpy(a,x.massMatrixTimesN2R_); massMatrixOverN2R_.axpy(a,x.massMatrixOverN2R_);}
    LinearAlgebra::Matrix massMatrixR_;
    LinearAlgebra::Matrix massMatrixOverBackgroundDensityR_;
	LinearAlgebra::Matrix massMatrixTimesBackgroundDensityR_;
	LinearAlgebra::Matrix massMatrixTimesBackgroundDensityDerivR_;
	LinearAlgebra::Matrix massMatrixOverBackgroundDensityDerivR_; 
	LinearAlgebra::Matrix massMatrixTimesN2R_;
	LinearAlgebra::Matrix massMatrixOverN2R_;	
};

struct FluxData
{
    FluxData(unsigned int nb)
    {
        left_.resize(nb);
        right_.resize(nb);
        
        for (unsigned int i=0; i<nb; ++i)
        {
            LinearAlgebra::Matrix& left= left_[i];
            LinearAlgebra::Matrix& right= right_[i];
            
            left.resize(12,nb);
            right.resize(12,nb);
        }
    }
    void print()
    {
        cout <<"left="<<endl;
        for (unsigned int n=0; n < left_.size();++n)
        {
            cout <<"n="<<n<<", "<<left_[n]<<endl;
        }
        
        cout <<"right="<<endl;
        for (unsigned int n=0; n < right_.size();++n)
        {
            cout <<"n="<<n<<", "<<right_[n]<<endl;
        }
        
    }
    FluxData operator*= (const double& scalar)
    {
        for (unsigned int n=0; n < left_.size();++n)
        {
            left_[n]*=scalar;
            right_[n]*=scalar;
        }
        return *this;
    }
    void axpy(double a, const FluxData& x)
    {
        for (unsigned int n=0; n < left_.size();++n)
        {
            left_[n].axpy(a, x.left_[n]);
            right_[n].axpy(a, x.right_[n]);
        }
    }
    
    VectorOfMatrices    left_;
    VectorOfMatrices    right_;
};



struct HEulerElementData: public UserElementData
{
    HEulerElementData(unsigned int ndof):
        massMatrix_(ndof, ndof),
        invMassMatrix_(ndof, ndof),
		massMatrixOverBackgroundDensity_(ndof, ndof),
		massMatrixTimesBackgroundDensity_(ndof, ndof),
		massMatrixTimesBackgroundDensityDeriv_(ndof, ndof),
		massMatrixOverBackgroundDensityDeriv_(ndof, ndof),
		massMatrixTimesN2_(ndof, ndof),
		massMatrixOverN2_(ndof, ndof)
    {
    }
    
    LinearAlgebra::Matrix   massMatrix_;
    LinearAlgebra::Matrix   invMassMatrix_;
	LinearAlgebra::Matrix   massMatrixOverBackgroundDensity_;
	LinearAlgebra::Matrix   massMatrixTimesBackgroundDensity_;
	LinearAlgebra::Matrix   massMatrixTimesBackgroundDensityDeriv_;
	LinearAlgebra::Matrix   massMatrixOverBackgroundDensityDeriv_;
	LinearAlgebra::Matrix   massMatrixTimesN2_;
	LinearAlgebra::Matrix   massMatrixOverN2_;	
};

struct HEulerGlobalVariables: public GlobalData
{
    unsigned int    nElements_;
    Mat             DivergenceFreeMatrix_;
    
    double          dt_;
	
	Mat				DX_;
	Mat				DY_;
	Mat				DZ_;
};



struct HEulerConfigurationData: public ConfigurationData
{
    enum  SolutionType 		{EBTrapezoidWAtau3Over2, EBbuckete, EBbucketd, EBbucketb, EBsemiellipse11mode, CS3DSingleMode, CS2DSingleMode, RE0, RES0, RES1, RES2, RESDO0, REMTC0, CS1D, CS1DSingleMode, RES02D, CS2D, CS3D, ICS2D, EBSingleMode, EB2DSingleMode, IGW2D, WA2D, Lambd3DMixed, IGW2Dprolonged, IGWN2LinearAiry}; //{INCOMPRESSIBLE_WALLS, INCOMPRESSIBLE_ONETHIRDPERIODIC, INCOMPRESSIBLE_PERIODIC, COMPRESSIBLE_WALLS, COMPRESSIBLE_PERIODIC};

    HEulerConfigurationData(unsigned int DIM, unsigned int numberOfUnknowns, unsigned int polynomialOrder, unsigned int  numberOfTimeLevels=1, SolutionType type=EBSingleMode):
        ConfigurationData(DIM, numberOfUnknowns, polynomialOrder, numberOfTimeLevels=1),
        solutionType_(type)
    {
            /// reading from a file
        theta_=0.5;
        //numOfPeriods_=5.0-1.0+0.9999999;
        //numOfTimeStepInOnePeriod_=100;
        //numOfPeriodsInOnePlotStep_=100;
        //onePeriod_=7.695298980971185;
		DIM_ = DIM;
		polynomialOrder_ = polynomialOrder;
		numberOfUnknowns_ = numberOfUnknowns;
        //stratified_ = 1;
		//N2_ = 2.0;
		//incompressible_ = 0;
    }
    
public:
    SolutionType    solutionType_;
    
    unsigned int    nx_;
    unsigned int    ny_;
    unsigned int    nz_;
    
    double          lx_;
    double          ly_;
    double          lz_;
    
    double          theta_;
	unsigned int	DIM_;
	unsigned int	polynomialOrder_;
	unsigned int	numberOfUnknowns_;
    
    double          numOfPeriods_;
    double          numOfTimeStepInOnePeriod_;
    double          numOfPeriodsInOnePlotStep_;
    double          onePeriod_;
	double			numOfHamInOnePeriod_;
	
	bool			stratified_;
	bool			incompressible_;
	bool			WA_;
	bool			EB_;
	
	double 			gamma_;
	double 			N2_;
};

class HEuler: public HpgemUI,public Integration::ElementIntegrandBase<ElementIntegralData>,public Integration::ElementIntegrandBase<ElementIntegralBackgroundDensity>,public Integration::FaceIntegrandBase<FluxData>,public Integration::ElementIntegrandBase<LinearAlgebra::Matrix>
{
public:
    typedef HpgemUI                        HpgemUIT;
    typedef Integration::ElementIntegral   ElementIntegralT;
    typedef Integration::FaceIntegral      FaceIntegralT;
    typedef ExactSolutionBase              ExactSolutionT;
    typedef PointReference               PointReferenceOnTheFaceT;
    
public:
    HEuler(HEulerGlobalVariables* global, const HEulerConfigurationData* config);

    ~HEuler();
public:
    
    void printFullMatrixInfo(Mat& matrix, const string& name);
    
    bool initialiseMesh();

    ///calculates mass matrix
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& massMatrix);
	
	void elementIntegrand(const ElementT* element, const PointReferenceT& p, ElementIntegralBackgroundDensity& returnObject);	


    //void calculateLocalEnergy(const ElementT& element, const PointReferenceT& p, double& returnValue);
    
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, ElementIntegralData& returnObject);
    
    void faceIntegrand(const FaceT* face,          const LinearAlgebra::NumericalVector& normal,
                       const PointReferenceOnTheFaceT& p,  FluxData& ret);
    
    void createCompressibleSystem();
    void createIncompressibleSystem();
	
	void calculateEnergy(Vec Lambda, double& Ham, double& Kin, double& Pot, double&  Int, double& xMom, double& yMom, double& zMom, double& TotalMass);

        /// create Mass Matrices, store them as a User Defined Element Data
        /// calculate projection of the every unknown on the FEM spaces.
    void initialConditions();
    
    void solve();
	
	 void createCompressibleSystemnew();
	void solvenew();
    
    void output(double time=0.0);
    
private:///utilities
    void outputMatrix(Mat& matrix, const string& name);
    void outputMatrix(const Mat& matrix, const string& name)const;

    
    void outputVectorMatlab(Vec& vec, const string& name);
    void outputVectorMatlab(const Vec& vec, const string& name)const;
    
    
    void correctInitialProjectionOfVelocity(const Vec& UInit, Vec& UCorrected)const;
    
    //void calculatePressure(const Mat& A, const Mat& Ah,const Vec& UCorrected);
    
    
private:
    ExactSolutionT*          exactSolution_;
    Mat                      P_;
    Mat                      Q_;
	Mat 					 R_;
	Mat						B_;
	Mat						C_;
};
#endif
