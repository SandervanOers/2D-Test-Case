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
    // Bucket
    PetscInt   Number_Of_Elements_Petsc=2, Number_Of_TimeSteps_In_One_Period=100, Method=1;
    PetscInt   Number_Of_Periods = 64; // 48 periods forced, 16 decay
    PetscInt   kmode=1;
    PetscScalar N2 = (0.37*5.57)*(0.37*5.57);
    PetscScalar   theta = 0.5;
    PetscInt    N_Petsc = 0, N_Q=0;
    PetscScalar nu = 1;//
    PetscInt    Dimensions = 2;
    PetscScalar F0 = 3.4*pow(10.0,-6.0);
    PetscScalar omega = 0.16*5.57;//0.16*5.57;
    PetscScalar Fr = 1;//0.56;//0.0291;//
    PetscScalar Re = 1.8*pow(10.0,4.0);
    PetscScalar gamma = 0.0;// PETSC_PI/20.0;
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
    {

        Eold = 0.0;
    int Number_Of_Polynomial_Steps = N_Petsc    ;
    //for (int Number_Of_Spatial_Steps = 1; Number_Of_Spatial_Steps < 8-Number_Of_Polynomial_Steps; Number_Of_Spatial_Steps++) //std::max(4,7-Number_Of_Polynomial_Steps)
    {

    int Number_Of_Spatial_Steps = 0;
    auto t0 = std::chrono::high_resolution_clock::now();
    //Number_Of_Elements_Petsc = pow(2.0, (double)Number_Of_Spatial_Steps);
    N_Petsc = Number_Of_Polynomial_Steps;

    std::string mesh_name;

    //mesh_name = mesh_name_trapezoid(Number_Of_Elements_Petsc);
    mesh_name = "Mesh/Trapezoid_Sensor_5mm.msh";
    //mesh_name = "Mesh/square_"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+".msh";
    //mesh_name = "Mesh/Bucket_TopRight_"+std::to_string(Number_Of_Elements_Petsc)+".msh";
    //mesh_name = "Mesh/Bucket_TopLeft_"+std::to_string(Number_Of_Elements_Petsc)+".msh";
    //mesh_name = "Mesh/bucket_"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+".msh";
    //mesh_name = "Mesh/Bucket_ConformingNew_"+std::to_string(Number_Of_Elements_Petsc)+".msh";

    std::vector<Squares2D> List_Of_Elements;
    s//td::vector<std::unique_ptr<Boundary>> List_Of_Boundaries;
    std::vector<InternalBoundariesSquares2D> List_Of_Boundaries;
    std::vector<std::unique_ptr<Vertex>> List_Of_Vertices;

    Vec VX, VY;
    Mat EToV;
    int element_num, node_num;
    load_msh_mesh2D(mesh_name, VX, VY, EToV, List_Of_Vertices, List_Of_Elements, element_num, node_num);

    //print(List_Of_Vertices);


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

    Calculate_Area_Square(List_Of_Elements, List_Of_Vertices);
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
    double rho_0_Deriv = N2;// + 1.0; // = beta
    /// Estimate the required time step
    PetscScalar   sigma;
    //sigma = calculate_sigma_2DEB(rho_0_Deriv, kxmode, kzmode, Fr);
    //sigma = calculate_sigma_2DEB_Bucket();
    sigma = omega;
    PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);
    PetscScalar Period = period_2DEB_Bucket();
    PetscPrintf(PETSC_COMM_SELF,"Period %6.4e\n", Period);
    PetscPrintf(PETSC_COMM_SELF,"Number of Periods %6.4e\n",(double)Number_Of_Periods);

    //Number_Of_TimeSteps_In_One_Period = 1000;//10*10*pow((double)N_Elements, (N_Petsc+1.0)/2.0);
    //PetscScalar DeltaT = 1.0/(double)Number_Of_TimeSteps_In_One_Period/sigma;


    //sigma = omega;
    PetscScalar DeltaT = 2.0*PETSC_PI/(double)Number_Of_TimeSteps_In_One_Period/sigma;


    unsigned int Number_Of_TimeSteps = Number_Of_TimeSteps_In_One_Period*Number_Of_Periods;

    std::cout << "Number_Of_TimeSteps_In_One_Period  =  " << Number_Of_TimeSteps_In_One_Period << " => Delta T = " << DeltaT << std::endl;

    //PetscLogStageRegister("Assembly", &stage);
    //PetscLogStagePush(stage);
    Mat E, ET, invM, M1, M2, M2_small, NMat, NDerivMat, invM_small, M1_small;
    create_Matrices_Quads_EB(List_Of_Vertices, List_Of_Boundaries, List_Of_Elements, N_Nodes, N_Petsc, N_Q, rho_0_Deriv, Fr, E, ET, invM, invM_small, M1, M1_small, M2, M2_small, NMat, NDerivMat);

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

    Vec Forcing_a;
    ComputeForcing(List_Of_Elements, N_Nodes, Forcing_a);

    Mat A, B;
    Mat DIV;
    //double Re = 100000;//pow(10,5);//3.25*100000;
    // Send List of Elements -> Get Np per Element for preallocation
    create_WA_System_Forced_MidPoint(E, ET, invM, invM_small, M1, M1_small, M2, M2_small, NMat, NDerivMat, Forcing_a, N_Nodes, N_Petsc, DeltaT, nu, A, B, DIV, Re, Fr, gamma);
    //create_EB_System_MidPoint(E, ET, invM, invM_small, M1, M1_small, M2, M2_small, NMat, NDerivMat, N_Nodes, N_Petsc, DeltaT, 0.0, A, B, DIV, Re, Fr);

    VecDestroy(&Forcing_a);

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

    Vec VecNodes;
    compute_InitialCondition_WA(List_Of_Elements, N_Nodes, rho_0_Deriv, Fr, kxmode, kzmode, Initial_Condition, VecNodes, Number_Of_Elements_Petsc, Number_Of_TimeSteps_In_One_Period, N_Petsc, DIV);
    //compute_InitialCondition_EB(List_Of_Elements, N_Nodes, rho_0_Deriv, Fr, kxmode, kzmode, Initial_Condition, VecNodes, Number_Of_Elements_Petsc, Number_Of_TimeSteps_In_One_Period, N_Petsc, DIV);
    //compute_InitialCondition_EB_Bucket(List_Of_Elements, N_Nodes, rho_0_Deriv, Fr, kxmode, kzmode, Initial_Condition, Number_Of_Elements_Petsc, Number_Of_TimeSteps_In_One_Period, N_Petsc, DIV);


    Vec Sol;
    //Number_Of_TimeSteps = 10*Number_Of_TimeSteps_In_One_Period;
    //Simulate_WA_Forced_Continuous(A, B, M1_small, M2_small, DIV, Initial_Condition, List_Of_Elements, N_Nodes, Number_Of_TimeSteps, DeltaT, Sol, 4, F0, omega);
    //Number_Of_TimeSteps = 1;
    Simulate_WA_Forced(A, B, M1_small, M2_small, DIV, Initial_Condition, List_Of_Elements, N_Nodes, Number_Of_TimeSteps, DeltaT, Sol, 4, F0, omega);


    //Simulate_EB(A, B, M1_small, M2_small, DIV, Initial_Condition, VecNodes, List_Of_Elements, N_Nodes, Number_Of_TimeSteps, DeltaT, Sol, 4, F0, omega);
    MatDestroy(&DIV);
    //double Error = calculate_Error2D(Initial_Condition, Sol, 2, 1.0/Number_Of_Elements_Petsc, 1.0/Number_Of_Elements_Petsc, (N_Petsc+1)*(N_Petsc+1));
    //double Error = calculate_Error2D(Initial_Condition, Sol, 2, List_Of_Elements, N_Nodes);
    //double Error = calculate_Error2D(Initial_Condition, Sol, 2, 1.0/Number_Of_Elements_Petsc, 1.0/Number_Of_Elements_Petsc, (N_Petsc+1)*(N_Petsc+1));
    double Error = calculate_Error2D(Initial_Condition, Sol, 2, 1.0/sqrt(Number_Of_Elements_Petsc), 1.0/sqrt(Number_Of_Elements_Petsc), (N_Petsc+1)*(N_Petsc+1));

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
    VecDestroy(&VecNodes);
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
    std::cout << "   h \t N \t \t Er \t Order \t t " << std::endl;
    for (auto i = L2Errors.begin(); i != L2Errors.end(); ++i)
    {
        std::cout << std::setw(5) << std::setprecision(4) << (*i)[0] << "    " << std::setw(4) << (*i)[1] << "    " << std::setw(10) << std::setprecision(4) << (*i)[2] << "    " << std::setw(5) << std::setprecision(2) << (*i)[3] << "    " << std::setw(5) << std::setprecision(4) << (*i)[4] << std::endl;
    }
    std::cout << "**********************************************************"<< std::endl;
    }
    }

    std::cout << "**********************************************************"<< std::endl;
    std::cout << "L2 Errors " << std::endl;
    std::cout << "   h \t N \t \t Er \t Order \t t " << std::endl;
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


