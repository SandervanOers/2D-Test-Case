static char help[] = "Solves the 2D EB WA Viscous Forced System. Creates Convergence Tables\n \n\n";

#include <petscksp.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "initial_cond.hpp"
#include "mesh_gen.hpp"
#include "Legendre_Gauss_Lobatto.hpp"
#include "HIGW.hpp"
#include "Elements.hpp"
#include <chrono>
//#include <slepceps.h>
//#include <slepcsys.h>
#include "Cubature2D.hpp"

int main(int argc,char **args)
{
PetscInitialize(&argc,&args,(char*)0,help);

//for (unsigned int I = 1; I <= 6; I++)

auto t1 = std::chrono::high_resolution_clock::now();
// Read in options from command line
PetscInt   Number_Of_Elements_Petsc=2, Number_Of_TimeSteps_In_One_Period=10, Method=1; // Method is meant for time integrator: Midpoint vs Stormer-Verlet vs Third Order
PetscInt   Number_Of_Periods = 64; // 48 periods forced, 16 decay
PetscInt   kmode=1;
PetscScalar N2 = 0.0;//(0.37*5.7558)*(0.37*5.7558);
PetscScalar   theta = 0.5;
PetscInt    N_Petsc = 0, N_Q=0;
PetscScalar nu = 1;
PetscInt    Dimensions = 2;
PetscScalar F0 = 3.4*pow(10.0,-6.0);
PetscScalar omega = 0.5*std::sqrt(2.0);//0.16*5.7558;
PetscScalar Fr = 1;
PetscScalar Re = 1.8*pow(10.0,4.0);
PetscScalar gamma = std::atan2(1,3);// PETSC_PI/20.0;
/*// Read in options from command line
PetscInt   Number_Of_Elements_Petsc=2, Number_Of_TimeSteps_In_One_Period=100, Method=1;
PetscInt   Number_Of_Periods = 64;//30,
PetscInt   kmode=1;
PetscScalar N2 = (0.37*1.79/0.325)*(0.37*1.79/0.325); //(0.37*0.052/0.325)*(0.37*0.052/0.325); //(0.37*1.79/0.325)*(0.37*1.79/0.325); //(0.37*0.56/0.325)*(0.37*0.56/0.325); //1;//0.13*0.13;  //   //0.1625*0.1625; // N2 = beta-1; beta = 1/rho_0 drho_0/dz
PetscScalar   theta = 0.5;
PetscInt    N_Petsc = 0, N_Q=0;
PetscScalar nu = 1; //0;//
PetscInt    Dimensions = 2;
PetscScalar F0 = 3*pow(10.0,-7.0);//0.00001; //0.001;//
PetscScalar omega = 0.16*1.79/0.325;// 0.04*1.79/0.325;//0.16*0.052/0.325;//0.14*1.79/0.325;//0.245*0.56/0.325;//0.7071;//0.42; ////0.05518;//0.17;//0.06139;//        //0.052; //0.06825;//
PetscScalar Fr = 1;//0.56;//0.0291;//
PetscScalar Re = 5.82*pow(10.0,5.0); // 1.82 16900;//
PetscScalar gamma = 0.0;// PETSC_PI/20.0;*/

PetscOptionsGetInt(NULL, NULL, "-n", &Number_Of_Elements_Petsc, NULL);
PetscOptionsGetInt(NULL, NULL, "-k", &kmode, NULL);
PetscOptionsGetInt(NULL, NULL, "-t", &Number_Of_TimeSteps_In_One_Period, NULL);
PetscOptionsGetInt(NULL, NULL, "-P", &Number_Of_Periods, NULL);
PetscOptionsGetInt(NULL, NULL, "-Method", &Method, NULL);
PetscOptionsGetScalar(NULL, NULL, "-N2", &N2, NULL);
PetscOptionsGetScalar(NULL, NULL, "-theta", &theta, NULL);
PetscOptionsGetScalar(NULL, NULL, "-F0", &F0, NULL);
PetscOptionsGetScalar(NULL, NULL, "-omega", &omega, NULL);
PetscOptionsGetInt(NULL, NULL, "-Order", &N_Petsc, NULL);
PetscOptionsGetInt(NULL, NULL, "-QuadratureAdded", &N_Q, NULL);
PetscOptionsGetScalar(NULL, NULL, "-nu", &nu, NULL);
PetscOptionsGetInt(NULL, NULL, "-dim", &Dimensions, NULL);
PetscOptionsGetScalar(NULL, NULL, "-Fr", &Fr, NULL);
PetscOptionsGetScalar(NULL, NULL, "-Re", &Re, NULL);
PetscOptionsGetScalar(NULL, NULL, "-gamma", &gamma, NULL);

PetscLogStage stage;

std::vector<std::vector<double>> L2Errors;
// Viewer Full Digits
PetscViewer viewer;
PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
PetscViewerSetType(viewer, PETSCVIEWERASCII);
PetscViewer viewer_dense;
PetscViewerCreate(PETSC_COMM_WORLD, &viewer_dense);
PetscViewerPushFormat(viewer_dense, PETSC_VIEWER_ASCII_DENSE);
PetscViewerSetType(viewer_dense, PETSCVIEWERASCII);
PetscViewer viewer_info;
PetscViewerASCIIOpen(PETSC_COMM_WORLD, NULL, &viewer_info);
PetscViewerPushFormat(viewer_info, PETSC_VIEWER_ASCII_INFO);

double Eold = 0;
//for (int Number_Of_Polynomial_Steps = 0; Number_Of_Polynomial_Steps < 3; Number_Of_Polynomial_Steps ++ ) //Number_Of_Polynomial_Steps += 2
//for (int Number_Of_Polynomial_Steps = 1; Number_Of_Polynomial_Steps < 2; Number_Of_Polynomial_Steps ++ ) //Number_Of_Polynomial_Steps += 2
//{

    Eold = 0.0;
//int Number_Of_Polynomial_Steps = N_Petsc    ;
//for (int Number_Of_Spatial_Steps = 1; Number_Of_Spatial_Steps < 8-Number_Of_Polynomial_Steps; Number_Of_Spatial_Steps++) //std::max(4,7-Number_Of_Polynomial_Steps)
//{

int Number_Of_Spatial_Steps = 0;
auto t0 = std::chrono::high_resolution_clock::now();
//Number_Of_Elements_Petsc = pow(2.0, (double)Number_Of_Spatial_Steps);
//N_Petsc = Number_Of_Polynomial_Steps;

std::string mesh_name;

//mesh_name = mesh_name_trapezoid(Number_Of_Elements_Petsc);
//mesh_name = "Mesh/Trapezoid_Sensor_5mm.msh";
//mesh_name = "Mesh/Trapezoid_Sensor_2p5mm.msh";
mesh_name = "Mesh/cuboid_"+std::to_string(Number_Of_Elements_Petsc)+".msh";
//mesh_name = "Mesh/square_"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+".msh";
//mesh_name = "Mesh/Bucket_TopRight_"+std::to_string(Number_Of_Elements_Petsc)+".msh";
//mesh_name = "Mesh/Bucket_TopLeft_"+std::to_string(Number_Of_Elements_Petsc)+".msh";
//mesh_name = "Mesh/bucket_"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+".msh";
//mesh_name = "Mesh/Bucket_ConformingNew_"+std::to_string(Number_Of_Elements_Petsc)+".msh";

//std::vector<Squares2D> List_Of_Elements;
std::vector<std::unique_ptr<Element>> List_Of_Elements;
std::vector<std::unique_ptr<Boundary>> List_Of_Boundaries;
std::vector<std::unique_ptr<Vertex>>   List_Of_Vertices;

//Vec VX, VY;
Mat EToV, EToV2;
int element_num, node_num;
//load_msh_mesh2D(mesh_name, VX, VY, EToV, List_Of_Vertices, List_Of_Elements, element_num, node_num);
//load_msh_mesh2D(mesh_name, EToV, List_Of_Vertices, List_Of_Elements, element_num, node_num);

load_msh_mesh(mesh_name, EToV, List_Of_Vertices, List_Of_Elements, element_num,  node_num);



print(List_Of_Vertices);

std::cout << "EToV = " << std::endl;
MatView(EToV, viewer);

Mat EToE, EToF;
Connect_3D(EToV, element_num, node_num, EToE, EToF, List_Of_Boundaries);
std::cout << "EToE = " << std::endl;
MatView(EToE, viewer_dense);
std::cout << "EToF = " << std::endl;
MatView(EToF, viewer_dense);
MatDestroy(&EToE);
MatDestroy(&EToF);
MatDestroy(&EToV);
//VecDestroy(&VX);
//VecDestroy(&VY);

print(List_Of_Boundaries);


set_Order_Polynomials_Uniform(List_Of_Elements, N_Petsc, N_Petsc, N_Petsc);
set_theta_Uniform(List_Of_Boundaries, theta);
set_Node_Coordinates_Uniform(List_Of_Elements, List_Of_Vertices, N_Petsc, N_Petsc, N_Petsc);
unsigned int N_Nodes = get_Number_Of_Nodes(List_Of_Elements);
std::cout << "Total Number of Nodes = " << N_Nodes << std::endl;
unsigned int N_Elements = List_Of_Elements.size();
std::cout << "Total Number of Elements = " << N_Elements << std::endl;


PetscScalar DeltaT = 2.0*PETSC_PI/(double)Number_Of_TimeSteps_In_One_Period;
unsigned int Number_Of_TimeSteps = Number_Of_TimeSteps_In_One_Period*Number_Of_Periods;
std::cout << "Number_Of_TimeSteps_In_One_Period  =  " << Number_Of_TimeSteps_In_One_Period << " => Delta T = " << DeltaT << std::endl;

Mat E, ET, invM, M1, M2, M2_small, NMat, NDerivMat, invM_small, M1_small;
create_Matrices_Cuboids(List_Of_Vertices, List_Of_Boundaries, List_Of_Elements, N_Nodes, N_Petsc, N_Petsc, N_Petsc, N_Q, Fr, E, ET, invM, invM_small, M1, M1_small, M2, M2_small, NMat, NDerivMat);

MatDestroy(&E);
MatDestroy(&ET);
MatDestroy(&invM);
MatDestroy(&invM_small);
MatDestroy(&M2);
MatDestroy(&NMat);
MatDestroy(&NDerivMat);






MatDestroy(&M1);
MatDestroy(&M1_small);
MatDestroy(&M2_small);
auto t3 = std::chrono::high_resolution_clock::now();


std::cout << "Execution took " ;
std::cout << std::chrono::duration_cast<std::chrono::seconds>(t3-t1).count() << " seconds ";
std::cout << std::chrono::duration_cast<std::chrono::minutes>(t3-t1).count() << " minutes ";
std::cout << std::chrono::duration_cast<std::chrono::hours>(t3-t1).count()   << " hours\n";




PetscViewerDestroy(&viewer);
PetscViewerDestroy(&viewer_dense);
PetscViewerDestroy(&viewer_info);

PetscFinalize();
return 1;
}
