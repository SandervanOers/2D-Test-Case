    static char help[] = "Solves the 2D Nonviscous Full System. Creates Convergence Tables\n \n\n";
/*
**********************************************************
Number_Of_TimeSteps_In_One_Period = 100
    2       0             0     -nan        0
    4       0        0.4713     -inf        0
    8       0        0.3898     0.27        0
   16       0        0.1199      1.7        0
   32       0       0.05411      1.1        2
   64       0       0.02501      1.1       10
  128       0        0.0125        1       40
  256       0      0.006911     0.85      170
    2       1        0.5607     -inf        0
    4       1        0.3954      0.5        0
    8       1        0.2336     0.76        0
   16       1        0.1962     0.25        1
   32       1        0.1981    -0.014        5
   64       1        0.2015    -0.024       22
  128       1        0.2009    0.0042      106
    2       2        0.2787     -inf        0
    4       2        0.2244     0.31        0
    8       2       0.04207      2.4        0
   16       2       0.01208      1.8        3
   32       2      0.003305      1.9       13
   64       2      0.001362      1.3       65
    2       3       0.04634     -inf        0
    4       3        0.1066     -1.2        0
    8       3       0.08104      0.4        1
   16       3       0.04919     0.72        7
   32       3        0.0604     -0.3       34
    2       4       0.00442     -inf        0
    4       4        0.1257     -4.8        1
    8       4       0.02534      2.3        4
   16       4      0.007179      1.8       17
    2       5     0.0007959     -inf        0
    4       5       0.09166     -6.8        3
    8       5       0.06724     0.45       13
    2       6     0.0005679     -inf        0
    4       6        0.1036     -7.5        5
    2       7     0.0005509     -inf        1
**********************************************************
    */
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
    PetscInt   Number_Of_Elements_Petsc=2, Number_Of_TimeSteps_In_One_Period=10, Method=1;
    PetscInt   Number_Of_Periods=1, kmode=1;
    PetscScalar N2 = 0.0;//1.0; // N2 = beta-1; beta = 1/rho_0 drho_0/dz
    PetscScalar   theta = 0.5;
    PetscInt    N_Petsc = 0, N_Q=0;
    PetscScalar nu = 0.0;
    PetscInt    Dimensions = 2;

    PetscOptionsGetInt(NULL, NULL, "-n", &Number_Of_Elements_Petsc, NULL);
    PetscOptionsGetInt(NULL, NULL, "-k", &kmode, NULL);
    PetscOptionsGetInt(NULL, NULL, "-t", &Number_Of_TimeSteps_In_One_Period, NULL);
    PetscOptionsGetInt(NULL, NULL, "-P", &Number_Of_Periods, NULL);
    PetscOptionsGetInt(NULL, NULL, "-Method", &Method, NULL);
    PetscOptionsGetScalar(NULL, NULL, "-N2", &N2, NULL);
    PetscOptionsGetScalar(NULL, NULL, "-theta", &theta, NULL);
    PetscOptionsGetInt(NULL, NULL, "-Order", &N_Petsc, NULL);
    PetscOptionsGetInt(NULL, NULL, "-QuadratureAdded", &N_Q, NULL);
    PetscOptionsGetScalar(NULL, NULL, "-nu", &nu, NULL);
    PetscOptionsGetInt(NULL, NULL, "-dim", &Dimensions, NULL);

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
    //for (int Number_Of_Polynomial_Steps = 0; Number_Of_Polynomial_Steps < 10; Number_Of_Polynomial_Steps++)
    {
        Eold = 0.0;
    int Number_Of_Polynomial_Steps = N_Petsc    ;
    for (int Number_Of_Spatial_Steps = 1; Number_Of_Spatial_Steps < std::max(5,9-Number_Of_Polynomial_Steps); Number_Of_Spatial_Steps++) //std::max(5,7-Number_Of_Polynomial_Steps)
    {

    //int Number_Of_Spatial_Steps = 0;
    auto t0 = std::chrono::high_resolution_clock::now();
    Number_Of_Elements_Petsc = pow(2.0, (double)Number_Of_Spatial_Steps);
    N_Petsc = Number_Of_Polynomial_Steps;

    //std::string mesh_name = "Mesh/square_"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(2*Number_Of_Elements_Petsc)+".msh"; //2* switched/
    //std::string mesh_name = "Mesh/square_unstructured_"+std::to_string(Number_Of_Elements_Petsc)+".msh"; //2* switched/
    //std::string mesh_name = "Mesh/square_tilted_quads_"+std::to_string(2*Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+".msh"; //2* switched/
    std::string mesh_name;
    if (Method == 0)
    {
        mesh_name = "Mesh/square_tilted_quads_"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+".msh"; //2* switched/
    }
    /**************************************************************** /
    // square_tilted_quads and manual_square (+ square) have opposite sign for E
    / ****************************************************************/
    else if (Method == 1)
    {
        mesh_name = "Mesh/manual_square_"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+".msh"; //2* switched/
    }
    else if (Method == 2)
    {
        mesh_name = "Mesh/square_unstructured_"+std::to_string(Number_Of_Elements_Petsc)+".msh"; //2* switched/
    }
    else if (Method == 3)
    {
        mesh_name = "Mesh/manual_square_tilted_quads_"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+".msh"; //2* switched/
    }
    else if (Method == 4)
    {
        mesh_name = "Mesh/square_"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+".msh"; //2* switched/
    }
    std::vector<Squares2D> List_Of_Elements;
    std::vector<InternalBoundariesSquares2D> List_Of_Boundaries;
    std::vector<VertexCoordinates2D> List_Of_Vertices;
    Vec VX, VY;
    Mat EToV;
    int element_num, node_num;
    load_msh_mesh2D(mesh_name, VX, VY, EToV, List_Of_Vertices, List_Of_Elements, element_num, node_num);

    /*
    std::cout << "List of Vertices "  << std::endl;
    std::cout << "ID : x y"  << std::endl;
    for(auto i = List_Of_Vertices.begin(); i < List_Of_Vertices.end(); i++)
        std::cout << (*i).getID() << ": " << (*i).getxCoordinate() << " " << (*i).getyCoordinate() << std::endl;
    */
    /*
    std::cout << "VX = " << std::endl;
    VecView(VX, viewer);
    std::cout << "VY = " << std::endl;
    VecView(VY, viewer);
    std::cout << "EToV = " << std::endl;
    MatView(EToV, viewer);
    */

    Mat EToE, EToF;
    Connect2D(EToV, element_num, node_num, EToE, EToF, List_Of_Boundaries);
    //std::cout << "EToE = " << std::endl;
    //MatView(EToE, viewer_dense);
    MatDestroy(&EToE);
    MatDestroy(&EToF);
    MatDestroy(&EToV);
    VecDestroy(&VX);
    VecDestroy(&VY);

    /*
    std::cout << "List of Elements "  << std::endl;
    std::cout << "ID : V1 V2 V3 V4"  << std::endl;
    for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
    {
        std::cout << (*i).getID() << ": " << (*i).getVertex_V1() << " " << (*i).getVertex_V2() << " " << (*i).getVertex_V3() << " " << (*i).getVertex_V4() << std::endl;
    }

    std::cout << "List of Boundaries "  << std::endl;
    std::cout << "ID : LeftElement RightElement"  << std::endl;
    for(auto i = List_Of_Boundaries.begin(); i < List_Of_Boundaries.end(); i++)
    {
        std::cout << (*i).getID() <<  " : " << (*i).getLeftElementID() << " " << (*i).getRightElementID() << std::endl;
    }
    */
    Calculate_Jacobian_Square(List_Of_Elements, List_Of_Vertices);
    Calculate_Jacobian_Boundaries_Square(List_Of_Elements, List_Of_Boundaries, List_Of_Vertices);
    set_Order_Polynomials_Uniform(List_Of_Elements, N_Petsc);
    set_theta_Uniform(List_Of_Boundaries, theta);
    set_Node_Coordinates_Uniform_Square2D(List_Of_Elements, List_Of_Vertices, N_Petsc);

    std::cout << "N = " << N_Petsc << std::endl;// << ", Np = " << Np << std::endl;
    unsigned int N_Nodes = get_Number_Of_Nodes(List_Of_Elements);
    std::cout << "Total Number of Nodes = " << N_Nodes << std::endl;
    unsigned int N_Elements = List_Of_Elements.size();
    std::cout << "Total Number of Elements = " << N_Elements << std::endl;
    PetscInt kxmode, kzmode;
    kxmode = 1;
    kzmode = 1;
    unsigned int rho_0_Deriv = N2;// + 1.0; // = beta
    /// Estimate the required time step
    PetscScalar   sigma;
    sigma = calculate_sigma_2DIC(rho_0_Deriv, kxmode, kzmode);
    PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);

    Number_Of_TimeSteps_In_One_Period = 100;//10*10*pow((double)N_Elements, (N_Petsc+1.0)/2.0);
    PetscScalar DeltaT = 1.0/(double)Number_Of_TimeSteps_In_One_Period/sigma;

    unsigned int Number_Of_TimeSteps = Number_Of_TimeSteps_In_One_Period*Number_Of_Periods;

    std::cout << "Number_Of_TimeSteps_In_One_Period  =  " << Number_Of_TimeSteps_In_One_Period << " => Delta T = " << DeltaT << std::endl;

    //PetscLogStageRegister("Assembly", &stage);
    //PetscLogStagePush(stage);
    Mat E, ET, invM, M1, M2, M2_small, NMat, NDerivMat, invM_small, M1_small;
    create_Matrices_Quads_Full_IC(List_Of_Vertices, List_Of_Boundaries, List_Of_Elements, N_Nodes, N_Petsc, N_Q, rho_0_Deriv, E, ET, invM, invM_small, M1, M1_small, M2, M2_small, NMat, NDerivMat);

    //std::cout <<"E = " << std::endl;
    //MatView(E, viewer_dense);
    //std::cout <<"M1_small = " << std::endl;
    //MatView(M1_small, viewer_dense);
    //std::cout <<"M2_small = " << std::endl;
    //MatView(M2_small, viewer_dense);
    //std::cout <<"invM_small = " << std::endl;
    //MatView(invM_small, viewer_dense);
    //std::cout <<"invM = " << std::endl;
    //MatView(invM, viewer);

    Mat A, B;
    Mat ALaplacian, ALaplacian_h;
    Mat DIV;
    // Send List of Elements -> Get Np per Element for preallocation
    create_Incompressible_System_MidPoint_Full(E, ET, invM, invM_small, M1, M1_small, M2, M2_small, NMat, NDerivMat, N_Nodes, N_Petsc, DeltaT, nu, A, B, ALaplacian, ALaplacian_h, DIV);

    //std::cout << "Store Global Matrices" << std::endl;
    /// TO DO
    //std::cout << "A = " << std::endl;
    //MatView(A, viewer_dense);
    //std::cout << "B = " << std::endl;
    //MatView(B, viewer_dense);
           // PetscViewer    viewerstore;
            //PetscViewerPushFormat(viewerstore,PETSC_VIEWER_ASCII_MATLAB);
           // std::string name = "A_N2_h16x16.dat";
            //PetscViewerBinaryOpen(PETSC_COMM_WORLD,name.c_str(),FILE_MODE_WRITE,&viewerstore);
           // MatView(A,viewerstore);
           // PetscViewerDestroy(&viewerstore);

    //std::cout << "B = " << std::endl;
    //MatView(B, viewer_dense);

    MatDestroy(&E);
    MatDestroy(&ET);
    MatDestroy(&invM);
    MatDestroy(&invM_small);
    MatDestroy(&M2);
    MatDestroy(&NMat);
    MatDestroy(&NDerivMat);

    Vec Initial_Condition;
    compute_InitialCondition_Incompressible(List_Of_Elements, N_Nodes, rho_0_Deriv, kxmode, kzmode, Initial_Condition, Number_Of_Elements_Petsc, Number_Of_TimeSteps_In_One_Period, N_Petsc, ALaplacian, ALaplacian_h, DIV);

    MatDestroy(&ALaplacian);
    MatDestroy(&ALaplacian_h);

    Vec Sol;
    Simulate_IC(A, B, M1_small, M2_small, DIV, Initial_Condition, List_Of_Elements, N_Nodes, Number_Of_TimeSteps, DeltaT, Sol, 4);
    MatDestroy(&DIV);

    double Error = calculate_Error2D(Initial_Condition, Sol, 2, 1.0/Number_Of_Elements_Petsc, 1.0/Number_Of_Elements_Petsc, (N_Petsc+1)*(N_Petsc+1));

        char szFileName[255] = {0};
        std::string store_solution = "Solution/Solutions/Sol_n"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+"N"+std::to_string(N_Petsc)+"Ts"+std::to_string(Number_Of_TimeSteps_In_One_Period)+".txt";
        std::string store_IC       = "Solution/Solutions/IC_n"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+"N"+std::to_string(N_Petsc)+"Ts"+std::to_string(Number_Of_TimeSteps_In_One_Period)+".txt";
        const char *store_solution_char = store_solution.c_str();
        const char *store_IC_char = store_IC.c_str();
        PetscViewer viewer2;
        sprintf(szFileName, store_solution_char);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
        VecView(Sol, viewer2);
        PetscViewerDestroy(&viewer2);
        PetscViewer viewer3;
        sprintf(szFileName, store_IC_char);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer3);
        VecView(Initial_Condition, viewer3);
        PetscViewerDestroy(&viewer3);


    VecDestroy(&Sol);

    VecDestroy(&Initial_Condition);
    MatDestroy(&A);
    MatDestroy(&B);

    MatDestroy(&M1);
    MatDestroy(&M1_small);
    MatDestroy(&M2_small);
    auto t2 = std::chrono::high_resolution_clock::now();

    double Enew = Error;
    std::vector<double> L2Error;
    L2Error.push_back(Number_Of_Elements_Petsc);
    L2Error.push_back(N_Petsc);
    L2Error.push_back(Enew);
    double Order = 0.0;
    if (Number_Of_Spatial_Steps > 0)
    {
        Order = log(Eold/Enew)/log(2);
    }
    L2Error.push_back(Order);
    L2Error.push_back(std::chrono::duration_cast<std::chrono::seconds>(t2-t0).count());
    L2Errors.push_back(L2Error);

    //FILE *g = fopen("L2Errors.txt", "a");
    //fprintf(g, "%d \t %d \t %1.16e \t %1.1e \t %1.1e \n", Number_Of_Elements, N, Enew, Order, std::chrono::duration_cast<std::chrono::milliseconds>(t2-t0).count()/1000.0);
    //fclose(g);

    Eold = Enew;
    std::cout << "**********************************************************"<< std::endl;
    std::cout << "L2 Errors " << std::endl;
    std::cout << "\t h \t N \t Er \t Order \t t " << std::endl;
    for (auto i = L2Errors.begin(); i != L2Errors.end(); ++i)
    {
        std::cout << std::setw(5) << std::setprecision(4) << (*i)[0] << "    " << std::setw(4) << (*i)[1] << "    " << std::setw(10) << std::setprecision(4) << (*i)[2] << "    " << std::setw(5) << std::setprecision(2) << (*i)[3] << "    " << std::setw(5) << std::setprecision(4) << (*i)[4] << std::endl;
    }
    std::cout << "**********************************************************"<< std::endl;
    }
    }

    std::cout << "**********************************************************"<< std::endl;
    std::cout << "L2 Errors " << std::endl;
    std::cout << "\t h \t N \t Er \t Order \t t " << std::endl;
    for (auto i = L2Errors.begin(); i != L2Errors.end(); ++i)
    {
        std::cout << std::setw(5) << std::setprecision(4) << (*i)[0] << "    " << std::setw(4) << (*i)[1] << "    " << std::setw(10) << std::setprecision(4) << (*i)[2] << "    " << std::setw(5) << std::setprecision(2) << (*i)[3] << "    " << std::setw(5) << std::setprecision(4) << (*i)[4] << std::endl;
    }
    std::cout << "**********************************************************"<< std::endl;

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


