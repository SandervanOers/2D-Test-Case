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

    auto t1 = std::chrono::high_resolution_clock::now();
    // Read in options from command line
    // Bucket
    PetscInt   Number_Of_Elements_Petsc=2, Number_Of_TimeSteps_In_One_Period=200, Method=1;
    PetscInt   Number_Of_Periods = 1;
    PetscInt   kmode=1;
    PetscScalar N2 = 1.0;
    PetscScalar   theta = 0.5;
    PetscInt    N_Petsc = 0, N_Q=0;
    PetscScalar nu = 0;
    PetscInt    Dimensions = 2;
    PetscScalar F0 = 0.0;
    PetscScalar omega = 0.5*sqrt(2.0);
    PetscScalar Fr = 1;
    PetscScalar Re = 1.0;
    PetscScalar gamma = 0.0;

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

    std::string mesh_name;
    mesh_name = "Mesh/cuboid_"+std::to_string(Number_Of_Elements_Petsc)+".msh";

    std::vector<Cuboid> List_Of_Elements;
    std::vector<InternalBoundariesCuboid> List_Of_Boundaries;
    std::vector<VertexCoordinates3D> List_Of_Vertices;
    Vec VX, VY, VZ;
    Mat EToV;
    int element_num, node_num;
    load_msh_mesh3D(mesh_name, VX, VY, VZ, EToV, List_Of_Vertices, List_Of_Elements, element_num, node_num);

    // Check det J > 0: counterclockwise ordering
    //Calculate_Jacobian_RectangularCuboid(List_Of_Elements, List_Of_Vertices);
/*
    std::cout << "List of Vertices "  << std::endl;
    std::cout << "ID : x y z"  << std::endl;
    for(auto i = List_Of_Vertices.begin(); i < List_Of_Vertices.end(); i++)
        std::cout << (*i).getID() << ": " << (*i).getxCoordinate() << " " << (*i).getyCoordinate() << " " << (*i).getzCoordinate() << std::endl;

    std::cout << "VX = " << std::endl;
    VecView(VX, viewer);
    std::cout << "VY = " << std::endl;
    VecView(VY, viewer);
    std::cout << "VZ = " << std::endl;
    VecView(VZ, viewer);
    std::cout << "EToV = " << std::endl;
    MatView(EToV, viewer);
*/

    Mat EToE, EToF;
    Connect3D(EToV, element_num, node_num, EToE, EToF, List_Of_Boundaries);
    std::cout << "EToE = " << std::endl;
    MatView(EToE, viewer_dense);
    std::cout << "EToF = " << std::endl;
    MatView(EToF, viewer_dense);
    MatDestroy(&EToE);
    MatDestroy(&EToF);
    MatDestroy(&EToV);
    VecDestroy(&VX);
    VecDestroy(&VY);
    VecDestroy(&VZ);

    std::cout << "List of Elements "  << std::endl;
    std::cout << "ID : V1 V2 V3 V4 V5 V6 V7 V8"  << std::endl;
    for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
    {
        std::cout << (*i).getID() << ": " << (*i).getVertex_V1() << " " << (*i).getVertex_V2() << " " << (*i).getVertex_V3() << " " << (*i).getVertex_V4() << " " << (*i).getVertex_V5() << " " << (*i).getVertex_V6() << " " << (*i).getVertex_V7() << " " << (*i).getVertex_V8() << std::endl;
    }

    std::cout << "List of Boundaries "  << std::endl;
    std::cout << "ID : LeftElement_ID RightElement_ID LeftElement_face RightElement_face "  << std::endl;
    for(auto i = List_Of_Boundaries.begin(); i < List_Of_Boundaries.end(); i++)
    {
        std::cout << (*i).getID() <<  " : " << (*i).getLeftElementID() << " " << (*i).getRightElementID() << " " << (*i).get_Type_Left() << " " << (*i).get_Type_Right() << std::endl;
    }

    Calculate_CuboidFaceNormals(List_Of_Elements, List_Of_Boundaries, List_Of_Vertices);
    set_Order_Polynomials_Uniform(List_Of_Elements, N_Petsc);

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


