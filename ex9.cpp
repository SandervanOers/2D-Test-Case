    static char help[] = "Solves the 1D Compressible Nonviscous System. Creates Convergence Tables\n \n\n";

#include <petscksp.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "initial_cond.hpp"
#include "mesh_gen.hpp"
#include "Legendre_Gauss_Lobatto.hpp"
#include "HIGW.hpp"
#include "Elements.hpp"
//#include "meshgen2D.hpp"
#include <chrono>
//#include <slepceps.h>
//#include <slepcsys.h>
#include "Cubature2D.hpp"

int main(int argc,char **args)
{
    PetscInitialize(&argc,&args,(char*)0,help);

    //for (unsigned int I = 1; I <= 6; I++)
    {
    auto t0 = std::chrono::high_resolution_clock::now();
    // Read in options from command line
    PetscInt   Number_Of_Elements_Petsc=2, Number_Of_TimeSteps_In_One_Period=10, Method=1;
    PetscInt   Number_Of_Periods=1, kmode=1;
    PetscScalar N2 = -1.0;//1.0; // N2 = beta-1; beta = 1/rho_0 drho_0/dz
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

    // Viewer Full Digits
    PetscViewer viewer;
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    PetscViewerSetType(viewer, PETSCVIEWERASCII);
    PetscViewer viewer_dense;
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer_dense);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
    PetscViewerSetType(viewer_dense, PETSCVIEWERASCII);
    PetscViewer viewer_info;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, NULL, &viewer_info);
    PetscViewerPushFormat(viewer_info, PETSC_VIEWER_ASCII_INFO);

    /*--------------------------------------------------------------------------*/
    /* Compute and Store Nodes and Vandermonde Matrices                         */
    //store_Nodes_Reference_Triangle();
    /* Compute and Store Lagrange Polynomials on Cubature Nodes                 */
    //store_LagrangePolynomial_Cubature();
    //store_DerivativeLagrangePolynomial_Cubature();
    /*--------------------------------------------------------------------------*/

    //for (unsigned int Nel = 1; Nel < 10; Nel++ )
    {
    //Number_Of_Elements_Petsc = pow(2,Nel);
    std::vector<Elements2D> List_Of_Elements2D;
    std::vector<Boundaries2D> List_Of_Boundaries2D;
    std::vector<VertexCoordinates2D> List_Of_Vertices;

    // Only Run when Recomputing the Nodes on Reference Triangles
    //store_Nodes_Reference_Triangle();

    //Number_Of_Elements_Petsc = pow(2.0,(double)I);
    unsigned int Nel_x = Number_Of_Elements_Petsc;//2;
    unsigned int Nel_y = Number_Of_Elements_Petsc;//2;
    Compute_Vertex_Coordinates_Uniform_Rectangle_2D(0,1, 0,1, Nel_x, Nel_y, List_Of_Vertices, List_Of_Boundaries2D, List_Of_Elements2D);
    Calculate_Jacobian(List_Of_Elements2D, List_Of_Vertices);
    Calculate_Jacobian_boundaries(List_Of_Boundaries2D, List_Of_Vertices);
    set_Order_Polynomials_Uniform(List_Of_Elements2D, N_Petsc);
    set_theta_Uniform(List_Of_Boundaries2D, theta);

    //set_Node_Coordinates_Uniform(List_Of_Elements2D, List_Of_Vertices, N_Petsc);

    /*
    auto tn1 = std::chrono::high_resolution_clock::now();
    set_Node_Coordinates_Uniform(List_Of_Elements2D, List_Of_Vertices, N_Petsc);
    auto tn2 = std::chrono::high_resolution_clock::now();
    set_Node_Coordinates_NonUniform(List_Of_Elements2D, List_Of_Vertices);
    auto tn3 = std::chrono::high_resolution_clock::now();
    */
    set_Node_Coordinates_ReadNonUniform(List_Of_Elements2D, List_Of_Boundaries2D, List_Of_Vertices);
    /*
    auto tn4 = std::chrono::high_resolution_clock::now();


    std::cout << "N = " << N_Petsc << ": Uniform Execution took "
          << std::chrono::duration_cast<std::chrono::milliseconds>(tn2-tn1).count()
          << " milliseconds\n";
              std::cout << "N = " << N_Petsc << ": Non Uniform Execution took "
          << std::chrono::duration_cast<std::chrono::milliseconds>(tn3-tn2).count()
          << " milliseconds\n";
              std::cout << "N = " << N_Petsc << ": Read Non Uniform Execution took "
          << std::chrono::duration_cast<std::chrono::milliseconds>(tn4-tn3).count()
          << " milliseconds\n";
    */
    std::cout << "N = " << N_Petsc << std::endl;// << ", Np = " << Np << std::endl;
    unsigned int N_Nodes = get_Number_Of_Nodes(List_Of_Elements2D);
    std::cout << "Total Number of Nodes = " << N_Nodes << std::endl;
    unsigned int N_Elements = List_Of_Elements2D.size();
    std::cout << "Total Number of Elements = " << N_Elements << std::endl;
/*
    std::cout << "List of Vertices "  << std::endl;
    std::cout << "ID : x y isInternal "  << std::endl;
    for(auto i = List_Of_Vertices.begin(); i < List_Of_Vertices.end(); i++)
        std::cout << (*i).getID() << ": " << (*i).getxCoordinate() << " " << (*i).getyCoordinate() << " " << (*i).isInternal() << std::endl;


    std::cout << "List of Boundaries "  << std::endl;
    std::cout << "ID : isInternal : V1 V2: LeftElement RightElement"  << std::endl;
    for(auto i = List_Of_Boundaries2D.begin(); i < List_Of_Boundaries2D.end(); i++)
        std::cout << (*i).getID() << ": " << (*i).isInternal() << " : " << (*i).getVertex_V1() << " " << (*i).getVertex_V2() << " : " << (*i).getLeftElementID() << " " << (*i).getRightElementID() << std::endl;



    std::cout << "List of Elements "  << std::endl;
    std::cout << "ID : V1 V2 V3 : B1 B2 B3 J N"  << std::endl;
    for(auto i = List_Of_Elements2D.begin(); i < List_Of_Elements2D.end(); i++)
    {
        std::cout << (*i).getID() << ": " << (*i).getVertex_V1() << "  " << (*i).getVertex_V2() << " " << (*i).getVertex_V3() << ": " << (*i).getBoundary_B1() << " " <<  (*i).getBoundary_B2() << " " << (*i).getBoundary_B3() << " " << (*i).getJacobian() << " " << (*i).get_Order_Of_Polynomials()<< std::endl;
    }
*/
/*



    std::cout << "List of Elements "  << std::endl;
    std::cout << "ID : V1 V2 V3 : B1 B2 B3 J N"  << std::endl;
    for(auto i = List_Of_Elements2D.begin(); i < List_Of_Elements2D.end(); i++)
    {
        std::cout << (*i).getID() << std::endl;
        std::cout << "X, Y = " << std::endl;
        std::vector<double> xx = (*i).get_node_coordinates_x();
        std::vector<double> yy = (*i).get_node_coordinates_y();
        for(unsigned int j = 0; j < xx.size() ; j++)
        {
            std::cout << xx[j] << " " << yy[j] << std::endl;
        }
        std::cout << std::endl;
    }
    */
/*
    {
        // Compute the nodes for a equilateral element
        Vec X, Y;
        Nodes2D(N_Petsc, X, Y);
        // Convert nodes to a reference element
        Vec R, S;
        XYtoRS(X, Y, R, S);

        Mat V2Dr, V2Ds;
        GradVandermonde2D(N_Petsc, R, S, V2Dr, V2Ds);
        std::cout << "V2Dr = " << std::endl;
        MatView(V2Dr, viewer);
        std::cout << "V2Ds = " << std::endl;
        MatView(V2Ds, viewer);

        Mat V2D;
        V2D = Vandermonde2D(N_Petsc, R, S);
        std::cout << "V2D = " << std::endl;
        MatView(V2D, viewer);
        Mat Dr, Ds;
        DMatrices2D(N_Petsc, R, S, V2D, Dr, Ds);
        std::cout << "Dr = " << std::endl;
        MatView(Dr, viewer);
        std::cout << "Ds = " << std::endl;
        MatView(Ds, viewer);

        unsigned int Order_Gaussian_Quadrature = 10;
        // Get Quadrature Rules Reference Triangle
        Vec cubR, cubS, cubW;
        unsigned int Ncub;
        Cubature2D(Order_Gaussian_Quadrature, cubR, cubS, cubW, Ncub);
        Mat cubDr, cubDs;
        DMatrices2D(N_Petsc, cubR, cubS, V2D, cubDr, cubDs);

        std::cout << "cubDr = " << std::endl;
        MatView(cubDr, viewer);
        std::cout << "cubDs = " << std::endl;
        MatView(cubDs, viewer);

        Mat cV = InterpMatrix2D(N_Petsc, cubR, cubS);
        Mat cVT ;
        MatTranspose(cV, MAT_INITIAL_MATRIX, &cVT);

        std::cout << "cV = " << std::endl;
        MatView(cV,  viewer);
        unsigned int Np = (N_Petsc+1)*(N_Petsc+2)/2;

        // build mass matrix
        //cMM = cV’*diag(cJ.*cW)*cV; cinfo(c).MM = cMM;
        // build physical derivative matrices
        //cinfo(c).Dx = cMM\(cV’*diag(cW.*cJ)*(diag(crx)*cDr+diag(csx)*cDs));
        //cinfo(c).Dy = cMM\(cV’*diag(cW.*cJ)*(diag(cry)*cDr+diag(csy)*cDs));

            Mat W;

            Mat D;
            MatCreate(PETSC_COMM_WORLD,&D);
            MatSetType(D,MATSEQAIJ);
            MatSetSizes(D, Ncub, Ncub, PETSC_DECIDE, PETSC_DECIDE);
            MatSeqAIJSetPreallocation(D, 1, NULL);
            std::cout << "cubW = " << std::endl;
            VecView(cubW, viewer);

            MatDiagonalSet(D, cubW, INSERT_VALUES);
            std::cout << "D = " << std::endl;
            MatView(D, viewer);

            Mat Intermediate;
            MatMatMult(cVT, D, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Intermediate);
            MatMatMult(Intermediate, cV, MAT_INITIAL_MATRIX, PETSC_DEFAULT,  &W);

            std::cout << "W = " << std::endl;
            MatView(W, viewer);



        MatDestroy(&cubDr);
        MatDestroy(&cubDs);
        VecDestroy(&cubR);
        VecDestroy(&cubS);
        VecDestroy(&cubW);
        MatDestroy(&V2Dr);
        MatDestroy(&V2Ds);
        MatDestroy(&V2D);
        MatDestroy(&Dr);
        MatDestroy(&Ds);
        VecDestroy(&X);
        VecDestroy(&Y);
        VecDestroy(&R);
        VecDestroy(&S);
    }
*/
    /*--------------------------------------------------------------------------*/
    /* Initial Condition / Exact Solution                                       */
    /*--------------------------------------------------------------------------*/
    /*
    std::cout << "Start Elemental Calculations " << std::endl;
    for (auto e = List_Of_Elements2D.begin(); e < List_Of_Elements2D.end(); e++)
    {
        //unsigned int ID = (*e).getID()-1;
        unsigned int Np = (*e).get_Number_Of_Nodes();
        double J = (*e).getJacobian();
        double drdx = (*e).get_rx();
        double drdy = (*e).get_ry();
        double dsdx = (*e).get_sx();
        double dsdy = (*e).get_sy();
        //std::cout << "Np = " << Np << std::endl;
        unsigned int pos = (*e).getPosition();

        unsigned int Order_Polynomials = (*e).get_Order_Of_Polynomials();

        Vec X, Y;
        Nodes2D(Order_Polynomials, X, Y);
        // Convert nodes to a reference element
        Vec R, S;
        XYtoRS(X, Y, R, S);

        Mat V2D;
        V2D = Vandermonde2D(N_Petsc, R, S);
        std::cout << "V2D = " << std::endl;
        MatView(V2D, viewer);
        Mat Dr, Ds;
        DMatrices2D(N_Petsc, R, S, V2D, Dr, Ds);

        std::vector<double> x = (*e).get_node_coordinates_x();
        std::vector<double> y = (*e).get_node_coordinates_y();
        double x_arr[Np];
        std::copy(x.begin(), x.end(), x_arr);
        double y_arr[Np];
        std::copy(y.begin(), y.end(), y_arr);
        Vec xx, yy;
        VecCreateSeqWithArray(PETSC_COMM_SELF, Np, Np, x_arr, &xx);
        VecCreateSeqWithArray(PETSC_COMM_SELF, Np, Np, y_arr, &yy);

        Vec rx, ry, sx, sy;
        double JJ;
        GeometricFactors2D(xx, yy, Dr, Ds, Np, rx, ry, sx, sy, JJ);

        std::cout << "rx = " << std::endl;
        VecView(rx, viewer);
        std::cout << "ry = " << std::endl;
        VecView(ry, viewer);
        std::cout << "sx = " << std::endl;
        VecView(sx, viewer);
        std::cout << "sy = " << std::endl;
        VecView(sy, viewer);
        std::cout << " old: " << std::endl;
        std::cout << J << " " << drdx <<  " " << drdy << " " << dsdx << " " << dsdy << std::endl;
    }*/

    PetscInt kxmode, kzmode;
    kxmode = 1;
    kzmode = 1;
    unsigned int rho_0_Deriv = N2 + 1.0; // = beta
    PetscScalar   sigma;
    sigma = calculate_sigma_2D_system1(rho_0_Deriv, kxmode, kzmode);
    PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);

    PetscScalar DeltaX = 1.0/((double)Nel_x/2.0);
    PetscScalar DeltaY = 1.0/((double)Nel_y/2.0);
    //Number_Of_TimeSteps_In_One_Period = 10*pow((Nel_x+Nel_y)/2, (N_Petsc+1)/2);//Number_Of_Elements*Number_Of_Elements;
    //Number_Of_TimeSteps_In_One_Period = 10*pow(N_Elements, (N_Petsc+1.0)/2.0);//Number_Of_Elements*Number_Of_Elements;
    Number_Of_TimeSteps_In_One_Period = 10.0 * pow(std::ceil(sqrt(sqrt(N_Elements))), N_Petsc+1.0);
    std::cout << "Number_Of_TimeSteps_In_One_Period = " <<  Number_Of_TimeSteps_In_One_Period << std::endl;

    PetscScalar DeltaT = 0.05;//1.0/(double)Number_Of_TimeSteps_In_One_Period/sigma;

    std::cout << "DeltaX = " << DeltaX << std::endl;
    std::cout << "DeltaY = " << DeltaY << std::endl;
    std::cout << Number_Of_TimeSteps_In_One_Period << " => " << DeltaT << std::endl;

    /// Check for CFL condition (when explicit)

    Mat E, ET, invM, M1, M2, NMat, NDerivMat;
    Mat Ex, ExT, Ey, EyT;
    Mat invM_small, M1_small;
     // Np can be variable
    double Np = (N_Petsc+1)*(N_Petsc+2)/2; // N_Nodes/N_Elements // Array of N_Nodes per Element
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, 5*Np, NULL, &E);  // 2*N_Nodes x N_Nodes //number of possible nonzero blocks are 5: element and his 4 neighbours (2D)
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, 2*N_Nodes, 2*5*Np, NULL, &ET); // N_Nodes x 2*N_Nodes
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 5*Np, NULL, &Ex);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 5*Np, NULL, &ExT);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 5*Np, NULL, &Ey);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 5*Np, NULL, &EyT);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &invM);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &invM_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &M1);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M1_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NMat);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M2);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NDerivMat);


    std::cout << "Start Elemental Calculations " << std::endl;
    unsigned int Order_Polynomials_old = 1000;
    for (auto e = List_Of_Elements2D.begin(); e < List_Of_Elements2D.end(); e++)
    {
        //unsigned int ID = (*e).getID()-1;
        unsigned int Np = (*e).get_Number_Of_Nodes();
        double J = (*e).getJacobian();
        double drdx = (*e).get_rx();
        double drdy = (*e).get_ry();
        double dsdx = (*e).get_sx();
        double dsdy = (*e).get_sy();
        //std::cout << "Np = " << Np << std::endl;
        unsigned int pos = (*e).getPosition();

        unsigned int Order_Polynomials = (*e).get_Order_Of_Polynomials();

        unsigned int Order_Gaussian_Quadrature = 2*Order_Polynomials+3+N_Q;//ceil(Order_Polynomials+3+N_Q); // + higher order for rho_0 term
        Order_Gaussian_Quadrature = 10  ;// 28;//std::max((uint)10, Order_Gaussian_Quadrature);


        // Elemental Matrices
        Mat Ex_el, ExT_el, Ey_el, EyT_el,  M1_el, NMat_el, M2_el, NDerivMat_el; //invM_el,
        MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &Ex_el);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &ExT_el);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &Ey_el);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &EyT_el);
        //MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &invM_el);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &M1_el);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &NMat_el);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &M2_el);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &NDerivMat_el);

        //if (Order_Polynomials != Order_Polynomials_old)
        //{

            PetscInt in[Np];
            for (unsigned int n=0;n<Np; n++)
            {
                in[n] = n+pos;
            }

            // Compute New Quadrature Points and Weights
            // Get Nodal Coordinates Reference Triangle
            Vec R, S;
            Read_RS_Coordinates_Reference_Triangle(Order_Polynomials, R, S);
                //VecView(R, viewer);
                //VecView(S, viewer);

            // Get Quadrature Rules Reference Triangle
            Vec cubR, cubS, cubW;
            unsigned int Ncub;
            Cubature2D(Order_Gaussian_Quadrature, cubR, cubS, cubW, Ncub);


            //drdx = 1;
            //drdy = 0;
            //dsdx = 0;
            //dsdy = 1;
            //J = 1;
      /*  {
        Mat V2D;
        V2D = Vandermonde2D(N_Petsc, R, S);
        Mat cubDr, cubDs;
        DMatrices2D(N_Petsc, cubR, cubS, V2D, cubDr, cubDs);

        std::cout << "cubDr = " << std::endl;
        MatView(cubDr, viewer);
        std::cout << "cubDs = " << std::endl;
        MatView(cubDs, viewer);

        Mat cV = InterpMatrix2D(N_Petsc, cubR, cubS);
        Mat cVT ;
        MatTranspose(cV, MAT_INITIAL_MATRIX, &cVT);

        std::cout << "cV = " << std::endl;
        MatView(cV,  viewer);
        unsigned int Np = (N_Petsc+1)*(N_Petsc+2)/2;

            Mat W;

            Mat D;
            MatCreate(PETSC_COMM_WORLD,&D);
            MatSetType(D,MATSEQAIJ);
            MatSetSizes(D, Ncub, Ncub, PETSC_DECIDE, PETSC_DECIDE);
            MatSeqAIJSetPreallocation(D, 1, NULL);
            std::cout << "cubW = " << std::endl;
            VecView(cubW, viewer);

            MatDiagonalSet(D, cubW, INSERT_VALUES);
            //std::cout << "D = " << std::endl;
            //MatView(D, viewer);

            Mat Intermediate;
            MatMatMult(cVT, D, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Intermediate);
            MatMatMult(Intermediate, cV, MAT_INITIAL_MATRIX, PETSC_DEFAULT,  &W);

            std::cout << "W = " << std::endl;
            MatView(W, viewer);

            Mat W1, W2;
            Mat Dx, Dy;
            MatDuplicate(cubDr, MAT_COPY_VALUES, &Dx);
            MatScale(Dx, drdx);
            MatAXPY(Dx,dsdx,cubDs,SAME_NONZERO_PATTERN);
            MatMatMult(Intermediate, Dx, MAT_INITIAL_MATRIX, PETSC_DEFAULT,  &W1);
            MatDuplicate(cubDr, MAT_COPY_VALUES, &Dy);
            MatScale(Dy, drdy);
            MatAXPY(Dy,dsdy,cubDs,SAME_NONZERO_PATTERN);
            MatMatMult(Intermediate, Dy, MAT_INITIAL_MATRIX, PETSC_DEFAULT,  &W2);
            std::cout << "W1 = " << std::endl;
            MatView(W1, viewer);
            std::cout << "W2 = " << std::endl;
            MatView(W2, viewer);


            Mat invMM = InverseMassMatrix2D(N_Petsc);
            Mat Dx1, Dy1;
            MatMatMult(invMM, W1, MAT_INITIAL_MATRIX, PETSC_DEFAULT,  &Dx1);
            MatMatMult(invMM, W2, MAT_INITIAL_MATRIX, PETSC_DEFAULT,  &Dy1);

            std::cout << "Dx = " << std::endl;
            MatView(Dx1, viewer);
            std::cout << "Dy = " << std::endl;
            MatView(Dy1, viewer);

            MatDestroy(&cV);
            MatDestroy(&W);
            MatDestroy(&cVT);
            MatDestroy(&D);
            MatDestroy(&Dx);
            MatDestroy(&Dy);
            MatDestroy(&Dx1);
            MatDestroy(&Dy1);
            MatDestroy(&W1);
            MatDestroy(&W2);
            MatDestroy(&Intermediate);
            MatDestroy(&invMM);
            MatDestroy(&V2D);

        MatDestroy(&cubDr);
        MatDestroy(&cubDs);
        }*/


            Mat L = load_LagrangePolynomial_Cubature(Order_Polynomials, Order_Gaussian_Quadrature); // Watch Transpose // GetRow faster than GetColumn in Petsc
            Mat dLdr = load_DerivativeLagrangePolynomial_Cubature(Order_Polynomials, Order_Gaussian_Quadrature, 1);
            Mat dLds = load_DerivativeLagrangePolynomial_Cubature(Order_Polynomials, Order_Gaussian_Quadrature, 0);

            //std::cout << "drdx = " << drdx << std::endl;
            //std::cout << "drdy = " << drdy << std::endl;
            //std::cout << "dsdx = " << dsdx << std::endl;
            //std::cout << "dsdy = " << dsdy << std::endl;

            std::cout << "J = " <<  J << std::endl;

            //std::cout << "L = " <<   std::endl;
            //MatView(L, viewer);
            //std::cout << "dLdr = " << std::endl;
            //MatView(dLdr, viewer);
            //std::cout << "dLds = " << std::endl;
            //MatView(dLds, viewer);


            /// Not Lagrange Polynomials, but Legendre Polynomials

            // Transform from (r,s) derivatives to (x,y) derivatives
            // dL/dx = dr/dx dL/dr + ds/dx dL/ds
            Mat dLdx, dLdy;
            MatDuplicate(dLdr, MAT_COPY_VALUES, &dLdx);
            MatScale(dLdx, drdx);
            MatAXPY(dLdx,dsdx,dLds,SAME_NONZERO_PATTERN);

            MatDuplicate(dLdr, MAT_COPY_VALUES, &dLdy);
            ///MatDuplicate(dLds, MAT_COPY_VALUES, &dLdy);
            MatScale(dLdy, drdy);
            MatAXPY(dLdy,dsdy,dLds,SAME_NONZERO_PATTERN);

            MatDestroy(&dLdr);
            MatDestroy(&dLds);

            //std::cout << "dLdx = " << std::endl;
            //MatView(dLdx, viewer);
            //std::cout << "dLdy = " << std::endl;
            //MatView(dLdy, viewer);

            PetscScalar *r_a, *s_a;
            VecGetArray(cubR, &r_a);
            VecGetArray(cubS, &s_a);
            double x_v1 = List_Of_Vertices[(*e).getVertex_V1()-1].getxCoordinate();
            double y_v1 = List_Of_Vertices[(*e).getVertex_V1()-1].getyCoordinate();
            double x_v2 = List_Of_Vertices[(*e).getVertex_V2()-1].getxCoordinate();
            double y_v2 = List_Of_Vertices[(*e).getVertex_V2()-1].getyCoordinate();
            double x_v3 = List_Of_Vertices[(*e).getVertex_V3()-1].getxCoordinate();
            double y_v3 = List_Of_Vertices[(*e).getVertex_V3()-1].getyCoordinate();

            std::vector<double> X, Y;
            for(unsigned int k = 0; k < Ncub; k++)
            {
                double x = 0.5*(-(r_a[k]+s_a[k])*x_v1+(1.0+r_a[k])*x_v2+(1.0+s_a[k])*x_v3);
                double y = 0.5*(-(r_a[k]+s_a[k])*y_v1+(1.0+r_a[k])*y_v2+(1.0+s_a[k])*y_v3);
                //std::cout << k << " " << x << " " << y << " " << rho_0_2D_system1(y, 1.0) << std::endl;
                X.push_back(x);
                Y.push_back(y);
            }
            // Get Varying Coefficients on Cubature Nodes
            std::vector<double> Rho0 = rho_0_2D_system1(Y, rho_0_Deriv);
            std::vector<double> Rho0Deriv = rho_0_deriv_2D_system1(Y, rho_0_Deriv);
            std::vector<double> N2Val = N_2_2D_system1(Y, rho_0_Deriv);

            //std::cout << "Y El = ";
            //for (const double& i : Y)
            //{
            //    std::cout << i << " " ;
            //}
            //std::cout << std::endl;
            //std::cout << "Rho0Deriv  El = ";
            //for (const double& i : Rho0Deriv )
            //{
            //    std::cout << i << " " ;
            //}
            //std::cout << std::endl;
            VecRestoreArray(cubR, &r_a);
            VecRestoreArray(cubS, &s_a);
            VecDestroy(&cubR);
            VecDestroy(&cubS);

            Mat W;
            MatCreate(PETSC_COMM_WORLD,&W);
            MatSetType(W,MATSEQAIJ);
            MatSetSizes(W, Np, Np, PETSC_DECIDE, PETSC_DECIDE);
            MatSeqAIJSetPreallocation(W, Np, NULL);

            PetscScalar *cubW_a;
            VecGetArray(cubW, &cubW_a);
            for (unsigned int k = 0; k < Np; k++)
            {
                PetscInt    nck; // nck == Ncub
                const PetscInt    *ak;
                const PetscScalar *Lk;
                MatGetRow(L, k, &nck, &ak, &Lk);
                for (unsigned int l = 0; l < Np; l++)
                {
                    PetscInt    ncl; // ncl == Ncub
                    const PetscInt    *al;
                    const PetscScalar *Ll;
                    const PetscScalar *dLdxl, *dLdyl;
                    MatGetRow(dLdx, l, &ncl, &al, &dLdxl);
                    MatGetRow(dLdy, l, &ncl, &al, &dLdyl);
                    MatGetRow(L, l, &ncl, &al, &Ll);

                    double value = 0;
                    double value_ex = 0.0;
                    double value_ey = 0.0;
                    double value_m = 0.0;
                    double value_n = 0.0;
                    double value_m2 = 0.0;
                    double value_n_deriv = 0.0;

                    for (unsigned int i = 0; i < nck; i++)
                    {
                        value += cubW_a[i]*Lk[i]*Ll[i];
                        value_m += cubW_a[i]*Lk[i]*Ll[i]/Rho0[i]*J;
                        value_n += cubW_a[i]*Lk[i]*Ll[i]*Rho0[i]*J;
                        value_n_deriv += cubW_a[i]*Lk[i]*Ll[i]*Rho0Deriv[i]*J;
                        value_ey += cubW_a[i]*Lk[i]*Ll[i]*Rho0Deriv[i]*J;
                        value_m2 += cubW_a[i]*Lk[i]*Ll[i]/Rho0[i]/N2Val[i]*J;
                        // Lk * d(rho_0 * Ll)/dx = Lk * rho_0 * dLl/dx
                        // Lk * d(rho_0 * Ll)/dy = Lk * rho_0 * dLl/dy + Lk * Ll *drho_0/dy
                        if (Np > 1)
                        {
                            value_ex += cubW_a[i]*Lk[i]*dLdxl[i]*Rho0[i]*J;
                            value_ey += cubW_a[i]*Lk[i]*dLdyl[i]*Rho0[i]*J;
                        }
                        //std::cout << "Rho0[i] = " << Rho0[i] << std::endl;
                        //std::cout << "Rho0Deriv[i] = " << Rho0Deriv[i] << std::endl;
                        //std::cout << "cubW_a[i] = " << cubW_a[i] << std::endl;
                        //std::cout << "Lk[i] = " << Lk[i] << std::endl;
                        //std::cout << "dLdxl[i] = " << dLdxl[i] << std::endl;

                    }

                    MatSetValue(W, k, l, value, INSERT_VALUES);
                    MatRestoreRow(dLdx, l, &ncl, &al, &dLdxl);
                    MatRestoreRow(dLdy, l, &ncl, &al, &dLdyl);
                    MatRestoreRow(L, l, &ncl, &al, &Ll);

                    double factor = -1.0;
                    MatSetValue(Ex_el, k, l, factor*value_ex, ADD_VALUES);
                    MatSetValue(Ey_el, k, l, factor*value_ey, ADD_VALUES);
                    value_ex = - value_ex;
                    value_ey = - value_ey;
                    MatSetValue(ExT_el, l, k, factor*value_ex, ADD_VALUES);
                    MatSetValue(EyT_el, l, k, factor*value_ey, ADD_VALUES);

                    MatSetValue(M1_el, k, l, value_m, ADD_VALUES);
                    MatSetValue(NMat_el, k, l, value_n, ADD_VALUES);
                    MatSetValue(NDerivMat_el, k, l, value_n_deriv, ADD_VALUES);
                    MatSetValue(M2_el, k, l, value_m2, ADD_VALUES);

                    /// ==> Elemental Matrices ==> Multiply V2DInv^T * W * V2DInv for Lagrange Values ==> Collect into global Matrices
                    //std::cout << k << " " << l << " " << value_ex << std::endl; // row k // col l
                }
                MatRestoreRow(L, k, &nck, &ak, &Lk);
            }
            VecRestoreArray(cubW, &cubW_a);
            VecDestroy(&cubW);
            MatDestroy(&L);
            MatDestroy(&dLdx);
            MatDestroy(&dLdy);

            MatAssemblyBegin(Ex_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(Ex_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyBegin(Ey_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(Ey_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyBegin(ExT_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(ExT_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyBegin(EyT_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(EyT_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyBegin(M1_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(M1_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyBegin(NMat_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(NMat_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyBegin(NDerivMat_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(NDerivMat_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyBegin(M2_el, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(M2_el, MAT_FINAL_ASSEMBLY);


            MatAssemblyBegin(W, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(W, MAT_FINAL_ASSEMBLY);
            //std::cout << "W = " << std::endl;
            //MatView(W, viewer);
            ///std::cout << "Ex_el = " << std::endl;
            ///MatView(Ex_el, viewer);
            //std::cout << "Ey_el = " << std::endl;
            //MatView(Ey_el, viewer);

            Mat Ex_el_L, Ey_el_L, ExT_el_L, EyT_el_L, M1_el_L, NMat_el_L, NDerivMat_el_L, M2_el_L;
            Mat V2DInv = load_InverseVandermondeMatrix(Order_Polynomials);
            M1_el_L = VandermondeMultiply(V2DInv, M1_el);
            M2_el_L = VandermondeMultiply(V2DInv, M2_el);
            NMat_el_L = VandermondeMultiply(V2DInv, NMat_el);
            NDerivMat_el_L = VandermondeMultiply(V2DInv, NDerivMat_el);
            Mat W_L;
            W_L = VandermondeMultiply(V2DInv, W);

            Mat V2D = load_VandermondeMatrix(Order_Polynomials);
            Ex_el_L = VandermondeMultiply(V2DInv, Ex_el);
            Ey_el_L = VandermondeMultiply(V2DInv, Ey_el);
            ExT_el_L = VandermondeMultiply(V2DInv, ExT_el);
            EyT_el_L = VandermondeMultiply(V2DInv, EyT_el);
            /*
            Mat TestIntermediateExel, TestIntermediateEyel, TestIntermediateExTel, TestIntermediateEyTel;

            MatMatMult(V2D, Ex_el, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &TestIntermediateExel);
            MatMatMult(TestIntermediateExel, V2DInv, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Ex_el_L);
            MatMatMult(V2D, Ey_el, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &TestIntermediateEyel);
            MatMatMult(TestIntermediateEyel, V2DInv, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Ey_el_L);
            MatMatMult(V2D, ExT_el, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &TestIntermediateExTel);
            MatMatMult(TestIntermediateExTel, V2DInv, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &ExT_el_L);
            MatMatMult(V2D, EyT_el, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &TestIntermediateEyTel);
            MatMatMult(TestIntermediateEyTel, V2DInv, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &EyT_el_L);
            MatDestroy(&TestIntermediateExel);
            MatDestroy(&TestIntermediateExTel);
            MatDestroy(&TestIntermediateEyel);
            MatDestroy(&TestIntermediateEyTel);

            std::cout << "Ex_el_L = " << std::endl;
            MatView(Ex_el_L, viewer);
            */

            /*
            {
            Mat TestDx,
            std::cout << "Dx = " << std::endl;
            MatView(TestDx, viewer);
            MatDestroy(&TestDx);
            MatDestroy(&TestIntermediate);
            }
            {
            Mat TestDx, TestIntermediate;
            MatMatMult(V2D, Ey_el, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &TestIntermediate);
            MatMatMult(TestIntermediate, V2DInv, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &TestDx);
            std::cout << "Dy = " << std::endl;
            MatView(TestDx, viewer);
            MatDestroy(&TestDx);
            MatDestroy(&TestIntermediate);
            }*/
            MatDestroy(&V2D);
            //Mat Test;
            //Mat Vinv;
            //Vinv = load_InverseVandermondeMatrix(N_Petsc);
            //MatMatMult(Ex_el, Vinv, MAT_INITIAL_MATRIX, 1, &Test);
            //std::cout << "Test = " << std::endl;
            //MatView(Test, viewer);
            //MatDestroy(&Test);

            //Mat Test2;
            //Mat Vinv;
            //Vinv = load_InverseVandermondeMatrix(N_Petsc);
            //MatMatMult(Ex_el_L, Vinv, MAT_INITIAL_MATRIX, 1, &Test2);
            //std::cout << "Test = " << std::endl;
            //MatView(Test2, viewer);
            //MatDestroy(&Test2);






            //MatDestroy(&Vinv);
            //std::cout << "Ex_el = " << std::endl;
            //MatView(Ex_el, viewer);
            //std::cout << "Ex_el_L = " << std::endl;
            //MatView(Ex_el_L, viewer);
            //std::cout << "V2DInv = " << std::endl;
            //MatView(V2DInv, viewer);

            //std::cout << "W = " << std::endl;
            //MatView(W, viewer);
            //std::cout << "W_L = " << std::endl;
            //MatView(W_L, viewer);
            //std::cout << "V2DInv = " << std::endl;
            //MatView(V2DInv, viewer);

            MatDestroy(&W_L);

            Mat invM_el = InverseMassMatrix2D(N_Petsc);


/*
            Mat Dx1, Dy1;
            MatMatMult(invM_el, Ex_el_L, MAT_INITIAL_MATRIX, PETSC_DEFAULT,  &Dx1);
            MatMatMult(invM_el, Ey_el_L, MAT_INITIAL_MATRIX, PETSC_DEFAULT,  &Dy1);

            std::cout << "Dx = " << std::endl;
            MatView(Dx1, viewer);
            std::cout << "Dy = " << std::endl;
            MatView(Dy1, viewer);



            std::vector<double> x = (*e).get_node_coordinates_x();
            std::vector<double> y = (*e).get_node_coordinates_y();

            std::vector<double> FunctionVal = rho_0_2D_system1(y, 1.0);
            std::vector<double> FunctionValDeriv = rho_0_deriv_2D_system1(y, 1.0);
            //for (auto i = FunctionVal.begin(); i < FunctionVal.end(); i++)
            //    std::cout << (*i) << " " ;
            //std::cout << std::endl;
            Vec FunctionValVec;
            Vec FunctionValVecDeriv;
            stdVectorToPetscVec(FunctionVal,FunctionValVec);
            stdVectorToPetscVec(FunctionValDeriv,FunctionValVecDeriv);
            std::cout << "Function Val = " << std::endl;
            VecView(FunctionValVec, viewer);
            std::cout << "FunctionDeriv Val = " << std::endl;
            VecView(FunctionValVecDeriv, viewer);
            std::cout << "x derivative of function = " << std::endl;
            Vec Result;
            VecDuplicate(FunctionValVec , &Result);
            MatMult(Dx1, FunctionValVec , Result);
            VecView(Result, viewer);

            std::cout << "y derivative of function = " << std::endl;
            Vec Result2;
            VecDuplicate(FunctionValVec , &Result2);
            MatMult(Dy1, FunctionValVec , Result2);
            VecView(Result2, viewer);


            Vec Difference;
            VecDuplicate(Result2, &Difference);

            VecWAXPY(Difference, 1.0, FunctionValVecDeriv, Result2);
            //VecWAXPY(Difference, -1.0, FunctionValVecDeriv, Result2);
            std::cout << "Difference = " << std::endl;
            VecView(Difference, viewer);
        PetscReal error;
        VecNorm(Difference, NORM_2, &error);

        std::cout << "L2 Error = " << error << std::endl;
        std::cout << "L2 Error / sqrt(Np)= " << error/sqrt(Np) << std::endl;


            MatDestroy(&Dx1);
            MatDestroy(&Dy1);
            VecDestroy(&FunctionValVec);
            VecDestroy(&FunctionValVecDeriv);
            VecDestroy(&Result);
            VecDestroy(&Result2);
            VecDestroy(&Difference);

*/

            //std::cout << "invM_el = " << std::endl;
            //MatView(invM_el, viewer);

            //PetscInt nc;
            //const PetscInt    *aj;
            //const PetscScalar *values;
            for (PetscInt i=0; i<Np; i++)
            {
                PetscInt ncols;
                const PetscInt *cols;
                const PetscScalar *vals_Ex, *vals_Ey, *vals_M1, *vals_invM, *vals_ExT, *vals_EyT;
                MatGetRow(Ex_el_L, i, &ncols, &cols, &vals_Ex);
                MatGetRow(ExT_el_L, i, &ncols, &cols, &vals_ExT);
                MatGetRow(Ey_el_L, i, &ncols, &cols, &vals_Ey);
                MatGetRow(EyT_el_L, i, &ncols, &cols, &vals_EyT);
                MatGetRow(M1_el_L, i, &ncols, &cols, &vals_M1);
                MatGetRow(invM_el, i, &ncols, &cols, &vals_invM);
                const PetscInt IndexI = pos+i;
                PetscInt GlobalIndexCol[ncols];
                PetscInt GlobalIndexCol2[ncols];
                for (unsigned int j=0; j < ncols; j++)
                {
                    GlobalIndexCol[j] = cols[j]+pos;
                    GlobalIndexCol2[j] = N_Nodes+cols[j]+pos;
                }
                PetscInt GlobalIndex[1] = {i + pos};
                PetscInt GlobalIndex2[1] = {N_Nodes+i + pos};
                //std::cout << "1" << std::endl;
                MatSetValues(Ex,1, GlobalIndex, ncols, GlobalIndexCol,vals_Ex,ADD_VALUES); // When changing to E and Et, indices should be changed
                //std::cout << "2" << std::endl;
                MatSetValues(E,1, GlobalIndex, ncols, GlobalIndexCol,vals_Ex,ADD_VALUES);
                //std::cout << "3" << std::endl;
                MatSetValues(ExT,1, GlobalIndex, ncols, GlobalIndexCol,vals_ExT,ADD_VALUES);
                //std::cout << "4" << std::endl;
                MatSetValues(ET,1, GlobalIndex, ncols, GlobalIndexCol,vals_ExT,ADD_VALUES);
                //std::cout << "5" << std::endl;
                MatSetValues(Ey,1, GlobalIndex, ncols, GlobalIndexCol,vals_Ey,ADD_VALUES);
                //std::cout << "6" << std::endl;
                MatSetValues(E,1, GlobalIndex2, ncols, GlobalIndexCol,vals_Ey,ADD_VALUES);
                //std::cout << "7" << std::endl;
                MatSetValues(EyT,1, GlobalIndex, ncols, GlobalIndexCol,vals_EyT,ADD_VALUES);
                //std::cout << "8" << std::endl;
                MatSetValues(ET,1, GlobalIndex, ncols, GlobalIndexCol2,vals_EyT,ADD_VALUES);
                //std::cout << "9" << std::endl;
                MatSetValues(M1_small,1, GlobalIndex, ncols, GlobalIndexCol,vals_M1,ADD_VALUES);
                //std::cout << "10" << std::endl;
                MatSetValues(M1,1, GlobalIndex, ncols, GlobalIndexCol,vals_M1,ADD_VALUES);
                //std::cout << "11" << std::endl;
                MatSetValues(M1,1, GlobalIndex2, ncols, GlobalIndexCol2,vals_M1,ADD_VALUES);
                //std::cout << "12" << std::endl;
                MatSetValues(invM_small,1, GlobalIndex, ncols, GlobalIndexCol,vals_invM,ADD_VALUES);
                //std::cout << "13" << std::endl;
                MatSetValues(invM,1, GlobalIndex, ncols, GlobalIndexCol,vals_invM,ADD_VALUES);
                //std::cout << "14" << std::endl;
                MatSetValues(invM,1, GlobalIndex2, ncols, GlobalIndexCol2,vals_invM,ADD_VALUES);
                //std::cout << "15" << std::endl;
                MatRestoreRow(Ex_el_L, i, &ncols, &cols, &vals_Ex);
                MatRestoreRow(Ey_el_L, i, &ncols, &cols, &vals_Ey);
                MatRestoreRow(M1_el_L, i, &ncols, &cols, &vals_M1);
                MatRestoreRow(invM_el, i, &ncols, &cols, &vals_invM);
                MatRestoreRow(ExT_el_L, i, &ncols, &cols, &vals_ExT);
                MatRestoreRow(EyT_el_L, i, &ncols, &cols, &vals_EyT);
            }

            //std::cout << "Np = " << Np << std::endl;
            //std::cout << "pos = " << pos << std::endl;
            //std::cout << "N_Nodes = " << N_Nodes << std::endl;


/*

            //std::cout << "W = " << std::endl;
            //MatView(W, viewer);


            Mat V2D;
            V2D = Vandermonde2D(N_Petsc, R, S);

            std::cout << "V2D = " << std::endl;
            MatView(V2D, viewer);
            std::cout << "V2D Inv= " << std::endl;
            MatView(V2DInv, viewer);

            //Mat Dr, Ds;
            //GradVandermonde2D(N_Petsc, R, S, Dr, Ds);
            //DMatrices2D(N_Petsc, R, S, V2D, Dr, Ds);

            //std::cout << "Dr = " << std::endl;
            //MatView(Dr, viewer);

            Mat Intermediate, cubDr;
            MatTransposeMatMult(V2DInv, W, MAT_INITIAL_MATRIX, 1.0, &Intermediate);
            MatMatMult(Intermediate,V2DInv,MAT_INITIAL_MATRIX,1.0, &cubDr);
            std::cout << "MM = " << std::endl;
            MatView(cubDr, viewer);

            Mat MM = MassMatrix2D(N_Petsc);
            std::cout << "MM = " << std::endl;
            MatView(MM, viewer);

            Mat  MMInv = Inverse_Matrix(MM);
            std::cout << "MMInv = " << std::endl;
            MatView(MMInv, viewer);
            MatDestroy(&MMInv);


            //MatDestroy(&V2DInv);
            MatDestroy(&cubDr);
            MatDestroy(&Intermediate);

            MatDestroy(&V2D);
            //MatDestroy(&Dr);
            //MatDestroy(&Ds);
            ///VecDestroy(&X);
            ///VecDestroy(&Y);
            ///VecDestroy(&R);
            ///VecDestroy(&S);

            MatDestroy(&MM);
            //MatDestroy(&Sr);
            //MatDestroy(&Ss);



*/






            MatDestroy(&W);
            MatDestroy(&V2DInv);










            /*
            Mat cubV = InterpMatrix2D(N_Petsc, cubR, cubS);

            Mat MM =  MassMatrix2D(N_Petsc);
            Mat MMcub = MassMatrix2D_Cubature(N_Petsc, cubV, cubW, Ncub);


            //custom mass matrix per element
            std::cout << "MM = " << std::endl;
            MatView(MM, viewer);

            std::cout << "MMcub = " << std::endl;
            MatView(MMcub, viewer);

            //std::cout << "cubW = " << std::endl;
            //VecView(cubW, viewer);

                    PetscInt size_V1, size_V2;
                    MatGetSize(cubV, &size_V1, &size_V2);
                    std::cout << "Size cub V = " << size_V1 << " x " << size_V2 << std::endl;
           // //MatDestroy(&cubV);
            //MatDestroy(&MM);
            //MatDestroy(&MMcub);



*/

        VecDestroy(&R);
        VecDestroy(&S);

        MatDestroy(&Ex_el);
        MatDestroy(&ExT_el);
        MatDestroy(&Ey_el);
        MatDestroy(&EyT_el);
        MatDestroy(&invM_el);
        MatDestroy(&M1_el);
        MatDestroy(&NMat_el);
        MatDestroy(&M2_el);
        MatDestroy(&NDerivMat_el);

        MatDestroy(&Ex_el_L);
        MatDestroy(&ExT_el_L);
        MatDestroy(&Ey_el_L);
        MatDestroy(&EyT_el_L);
        MatDestroy(&M1_el_L);
        MatDestroy(&NMat_el_L);
        MatDestroy(&M2_el_L);
        MatDestroy(&NDerivMat_el_L);
    }
    MatAssemblyBegin(M1, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(M1_small, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1_small, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(invM, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(invM, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(invM_small, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(invM_small, MAT_FINAL_ASSEMBLY);
/*
    MatAssemblyBegin(E, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd(E, MAT_FINAL_ASSEMBLY );
    MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY );

    std::cout << "M1 = " << std::endl;
    MatView(M1, viewer);
    std::cout << "M1_small = " << std::endl;
    MatView(M1_small, viewer);
    std::cout << "E = " << std::endl;
    MatView(E, viewer);
    std::cout << "ET = " << std::endl;
    MatView(ET, viewer);
    */
    /*
    MatAssemblyBegin(M2, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M2, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(NMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(NMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(NDerivMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(NDerivMat, MAT_FINAL_ASSEMBLY);

    MatDestroy(&invM_Elemental);
    */



    std::cout << "Start Faces Calculations " << std::endl;

    /////////////////////////////////////////////////////////////////////
    // What happens when N = 0 => Np = 1, one node per element         //
    /////////////////////////////////////////////////////////////////////
    // Is the order of the node numbers correct                        //
    /////////////////////////////////////////////////////////////////////
    // one node => on entire element, # of nodes on each face is 1     //
    /////////////////////////////////////////////////////////////////////

    for (auto f = List_Of_Boundaries2D.begin(); f < List_Of_Boundaries2D.end(); f++)
    {
        if ((*f).isInternal())
        {
            //std::cout << (*f).getID() << " is Internal" << std::endl;

            double Jacobian = (*f).getJacobian();
            unsigned int Type_Boundary = (*f).getType();
            double theta = (*f).get_theta();
            double nx = (*f).get_nx();
            double ny = (*f).get_ny();
            ///std::cout << " nx = " << nx << ", ny = " << ny << std::endl;

            int left = (*f).getLeftElementID()-1;
            int right = (*f).getRightElementID()-1;

            std::vector<unsigned int> Node_Numbers_On_Boundary_Left = List_Of_Elements2D[left].get_nodes_on_boundary(Type_Boundary);
            std::vector<unsigned int> Node_Numbers_On_Boundary_Right = List_Of_Elements2D[right].get_nodes_on_boundary(Type_Boundary);


            //std::cout << std::endl << "Internal Boundary " << (*f).getID() << " Type = " << Type_Boundary << std::endl;
                /*
            std::cout << "Left: " ;
            for (auto i = Node_Numbers_On_Boundary_Left.begin(); i < Node_Numbers_On_Boundary_Left.end(); i++)
            {
                std::cout << (*i) << " ";
            }
            std::cout << "Right: " ;
            for (auto i = Node_Numbers_On_Boundary_Right.begin(); i < Node_Numbers_On_Boundary_Right.end(); i++)
            {
                std::cout << (*i) << " ";
            }
            std::cout << std::endl << std::endl;
            */

            if (Type_Boundary == 3)
            {
                std::reverse(Node_Numbers_On_Boundary_Left.begin(),Node_Numbers_On_Boundary_Left.end());
                //std::reverse(Node_Numbers_On_Boundary_Right.begin(),Node_Numbers_On_Boundary_Right.end());
            }
            if (Type_Boundary == 1)
            {
                std::reverse(Node_Numbers_On_Boundary_Right.begin(),Node_Numbers_On_Boundary_Right.end());
            }
            if (Type_Boundary == 2)
            {
                std::reverse(Node_Numbers_On_Boundary_Right.begin(),Node_Numbers_On_Boundary_Right.end());
                //std::reverse(Node_Numbers_On_Boundary_Left.begin(),Node_Numbers_On_Boundary_Left.end());
            }


            /*
            std::cout << "Internal Boundary " << (*f).getID() << std::endl;
            std::cout << "Left: " ;
            for (auto i = Node_Numbers_On_Boundary_Left.begin(); i < Node_Numbers_On_Boundary_Left.end(); i++)
            {
                std::cout << (*i) << " ";
            }
            std::cout << "Right: " ;
            for (auto i = Node_Numbers_On_Boundary_Right.begin(); i < Node_Numbers_On_Boundary_Right.end(); i++)
            {
                std::cout << (*i) << " ";
            }
            std::cout << std::endl;
            */

            unsigned int Np_left = List_Of_Elements2D[left].get_Number_Of_Nodes(); // Assumes Ordering
            unsigned int Np_right = List_Of_Elements2D[right].get_Number_Of_Nodes(); // Assumes Ordering

            unsigned int posL = List_Of_Elements2D[left].getPosition();
            unsigned int posR = List_Of_Elements2D[right].getPosition();

            //std::cout << "posL = " << posL << std::endl;
            //std::cout << "posR = " << posR << std::endl;

            unsigned int Order_Polynomials_left = (List_Of_Elements2D[left]).get_Order_Of_Polynomials();
            unsigned int Order_Polynomials_right = (List_Of_Elements2D[right]).get_Order_Of_Polynomials();

            /// Or use two different gaussian quadratures
            // unsigned int Order_Gaussian_Quadrature_L
            // unsigned int Order_Gaussian_Quadrature_R
            unsigned int Order_Gaussian_Quadrature  = ceil(std::max(2*Order_Polynomials_left, 2*Order_Polynomials_right)+3+N_Q);
            Order_Gaussian_Quadrature = 10;//std::max((uint)10, Order_Gaussian_Quadrature);
            // Order_Gaussian_Quadrature+1 = Number of Points
            Vec Weights;
            Vec QuadraturePoints;
            QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);

            unsigned int Vertex_ID_1 = (*f).getVertex_V1()-1;
            unsigned int Vertex_ID_2 = (*f).getVertex_V2()-1;

            double x_v1 = List_Of_Vertices[Vertex_ID_1].getxCoordinate(); // Physical Coordinates
            double y_v1 = List_Of_Vertices[Vertex_ID_1].getyCoordinate();
            double x_v2 = List_Of_Vertices[Vertex_ID_2].getxCoordinate();
            double y_v2 = List_Of_Vertices[Vertex_ID_2].getyCoordinate();

            PetscScalar *w_a, *r_a;
            VecGetArray(QuadraturePoints, &r_a);
            VecGetArray(Weights, &w_a);
            std::vector<double> X, Y;
            for(unsigned int k = 0; k <= Order_Gaussian_Quadrature; k++) // left and right different quadrature order?
            {
                double x = x_v1 + (x_v2-x_v1)*(1.0+r_a[k])/2.0;
                double y = y_v1 + (y_v2-y_v1)*(1.0+r_a[k])/2.0;
                //std::cout << k << " " << x << " " << y << " " << rho_0_2D_system1(y, 1.0) << std::endl;
                X.push_back(x);
                Y.push_back(y);
            }
            // Get Varying Coefficients on Cubature Nodes
            std::vector<double> Rho0_L = rho_0_2D_system1(Y, rho_0_Deriv);
            std::vector<double> Rho0_R = rho_0_2D_system1(Y, rho_0_Deriv);

            Vec ri_left, ri_right;
            ri_left = JacobiGL(0, 0, Order_Polynomials_left);
            ri_right = JacobiGL(0, 0, Order_Polynomials_right);

            //std::cout << "ri_left = " << std::endl;
            //VecView(ri_left, viewer);

            //std::cout << "Order_Polynomials_left = " << Order_Polynomials_left << std::endl;
            std::cout << "Jacobian = " << Jacobian << std::endl;

            // GLL
            Mat GLL;
            MatCreate(PETSC_COMM_WORLD,&GLL);
            MatSetType(GLL,MATSEQAIJ);
            MatSetSizes(GLL, Np_left, Np_left, PETSC_DECIDE, PETSC_DECIDE);
            MatSeqAIJSetPreallocation(GLL, Order_Polynomials_left+1, NULL);
            //std::cout << "MM 1D = " << std::endl;
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++) // Should be N_Left + 1, etc.
            {
                for (unsigned int j = 0; j <= Order_Polynomials_left; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_L
                    {
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += (1.0-theta)*w_a[q]*Rho0_L[q]*Li*Lj*Jacobian;

                    }
                    MatSetValue(GLL, Node_Numbers_On_Boundary_Left[i], Node_Numbers_On_Boundary_Left[j], value_e, ADD_VALUES);
                    MatSetValue(Ex,  posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ExT, posL+Node_Numbers_On_Boundary_Left[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(Ey,  posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(EyT, posL+Node_Numbers_On_Boundary_Left[j], posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    //std::cout << posL+Node_Numbers_On_Boundary_Left[i] << " " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    //std::cout << value_e << " " ;
                }
                //std::cout << std::endl;
            }
                //std::cout << std::endl;
            MatAssemblyBegin(GLL, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(GLL, MAT_FINAL_ASSEMBLY);
            //std::cout << "GLL = " << std::endl;
            //MatView(GLL, viewer);
            MatDestroy(&GLL);
            // GLR
            Mat GLR;
            MatCreate(PETSC_COMM_WORLD,&GLR);
            MatSetType(GLR,MATSEQAIJ);
            MatSetSizes(GLR, Np_left, Np_right, PETSC_DECIDE, PETSC_DECIDE);
            MatSeqAIJSetPreallocation(GLR, Order_Polynomials_left+1, NULL);
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -(1.0-theta)*w_a[q]*Rho0_R[q]*Li*Lj*Jacobian;
                    }
                    MatSetValue(GLR, Node_Numbers_On_Boundary_Left[i], Node_Numbers_On_Boundary_Right[j], value_e, ADD_VALUES);
                    MatSetValue(Ex,  posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ExT, posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(Ey,  posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(EyT, posR+Node_Numbers_On_Boundary_Right[j],  posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    //std::cout << posL+Node_Numbers_On_Boundary_Left[i] << " " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                }
            }
            MatAssemblyBegin(GLR, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(GLR, MAT_FINAL_ASSEMBLY);
            //std::cout << "GLR = " << std::endl;
            //MatView(GLR, viewer);
            MatDestroy(&GLR);
            // GRL
            Mat GRL;
            MatCreate(PETSC_COMM_WORLD,&GRL);
            MatSetType(GRL,MATSEQAIJ);
            MatSetSizes(GRL, Np_right, Np_left, PETSC_DECIDE, PETSC_DECIDE);
            MatSeqAIJSetPreallocation(GRL, Order_Polynomials_right+1, NULL);
            for (unsigned int i = 0; i <= Order_Polynomials_right; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_left; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += theta*w_a[q]*Rho0_L[q]*Li*Lj*Jacobian;
                    }
                    MatSetValue(GRL, Node_Numbers_On_Boundary_Right[i], Node_Numbers_On_Boundary_Left[j], value_e, ADD_VALUES);
                    MatSetValue(Ex,  posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ExT, posL+Node_Numbers_On_Boundary_Left[j],  posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(Ey,  posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(EyT, posL+Node_Numbers_On_Boundary_Left[j],  posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);

                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);

                    //std::cout << posR+Node_Numbers_On_Boundary_Right[i] << " " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                }
            }
            MatAssemblyBegin(GRL, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(GRL, MAT_FINAL_ASSEMBLY);
            //std::cout << "GRL = " << std::endl;
            //MatView(GRL, viewer);
            MatDestroy(&GRL);
            // GRR
            Mat GRR;
            MatCreate(PETSC_COMM_WORLD,&GRR);
            MatSetType(GRR,MATSEQAIJ);
            MatSetSizes(GRR, Np_right, Np_right, PETSC_DECIDE, PETSC_DECIDE);
            MatSeqAIJSetPreallocation(GRR, Order_Polynomials_right+1, NULL);
            for (unsigned int i = 0; i <= Order_Polynomials_right; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_R
                    {
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -theta*w_a[q]*Rho0_R[q]*Li*Lj*Jacobian;
                    }
                    MatSetValue(GRR, Node_Numbers_On_Boundary_Right[i], Node_Numbers_On_Boundary_Right[j], value_e, ADD_VALUES);
                    MatSetValue(Ex,  posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ExT, posR+Node_Numbers_On_Boundary_Right[j], posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(Ey,  posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(EyT, posR+Node_Numbers_On_Boundary_Right[j], posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);

                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);

                    //std::cout << posR+Node_Numbers_On_Boundary_Right[i] << " " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                }
            }
            MatAssemblyBegin(GRR, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(GRR, MAT_FINAL_ASSEMBLY);
            //std::cout << "GRR = " << std::endl;
            //MatView(GRR, viewer);
            MatDestroy(&GRR);

            VecRestoreArray(QuadraturePoints, &r_a);
            VecRestoreArray(Weights, &w_a);
            VecDestroy(&QuadraturePoints);
            VecDestroy(&Weights);
            VecDestroy(&ri_right);
            VecDestroy(&ri_left);


            /////////////////////////////////////////////
            // vec n = (nx, ny) -> still split into Ex, Ey?
            // if nx = 0 -> where in E do the boundary terms show up -> in Ey.
            // if nx, ny /= 0, then in both Ex and Ey
            /////////////////////////////////////////////

        }
        else
        {
            //std::cout << (*f).getID() << " is External" << std::endl;
        }
    }

    MatAssemblyBegin(Ex, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Ex, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Ey, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Ey, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(ExT, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(ExT, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(EyT, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EyT, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(E, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(E, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY);

    //std::cout << "List of Boundaries "  << std::endl;
    //std::cout << "ID : isInternal LeftElement RightElement"  << std::endl;
    //for(auto i = List_Of_Boundaries2D.begin(); i < List_Of_Boundaries2D.end(); i++)
    //    std::cout << (*i).getID() << ": " << (*i).isInternal()  << " " << (*i).getLeftElementID() << " " << (*i).getRightElementID() << std::endl;

    MatAssemblyBegin(E, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd(E, MAT_FINAL_ASSEMBLY );
    MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY );

    //std::cout << "M1 = " << std::endl;
    //MatView(M1, viewer);
    //std::cout << "M1_small = " << std::endl;
    //MatView(M1_small, viewer);
    ///std::cout << "E = " << std::endl;
    ///MatView(E, viewer);
    //std::cout << "ET = " << std::endl;
    //MatView(ET, viewer);

    std::cout << "Start Global Matrices Construction" << std::endl;

    double fillBF = 1;

    Mat BF1, BF1_TEMP1, BF1_TEMP2;
    MatMatMult(E, invM_small, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP1);
    MatMatMult(BF1_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP2);
    MatMatMult(invM, BF1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF1);
    MatDestroy(&BF1_TEMP1);
    MatDestroy(&BF1_TEMP2);
/*
    Mat BF2, BF2_TEMP1, BF2_TEMP2;
    MatMatMult(E, invM, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP1);
    MatMatMult(BF2_TEMP1, M2, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP2);
    MatMatMult(invM, BF2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF2);
    MatDestroy(&BF2_TEMP1);
    MatDestroy(&BF2_TEMP2);
*/
    Mat DIV, DIV_TEMP1, DIV_TEMP2;
    MatMatMult(ET, invM, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP1);
    MatMatMult(DIV_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP2);
    MatMatMult(invM_small, DIV_TEMP2, MAT_INITIAL_MATRIX, fillBF, &DIV);
    MatDestroy(&DIV_TEMP1);
    MatDestroy(&DIV_TEMP2);
/*
    Mat C, C_TEMP1, C_TEMP2;
    MatMatMult(NMat, invM, MAT_INITIAL_MATRIX, fillBF, &C_TEMP1);
    MatMatMult(C_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &C_TEMP2);
    MatMatMult(invM, C_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C);
    MatDestroy(&C_TEMP1);
    MatDestroy(&C_TEMP2);

    Mat C2, C2_TEMP1, C2_TEMP2;
    MatMatMult(NMat, invM, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP1);
    MatMatMult(C2_TEMP1, M2, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP2);
    MatMatMult(invM, C2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C2);
    MatDestroy(&C2_TEMP1);
    MatDestroy(&C2_TEMP2);

    Mat D1, D1_TEMP1, D1_TEMP2;
    MatMatMult(NDerivMat, invM, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP1);
    MatMatMult(D1_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP2);
    MatMatMult(invM, D1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D1);
    MatDestroy(&D1_TEMP1);
    MatDestroy(&D1_TEMP2);

    Mat D2, D2_TEMP1, D2_TEMP2;
    MatMatMult(NDerivMat, invM, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP1);
    MatMatMult(D2_TEMP1, M2, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP2);
    MatMatMult(invM, D2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D2);
    MatDestroy(&D2_TEMP1);
    MatDestroy(&D2_TEMP2);

    MatDestroy(&E);
    MatDestroy(&ET);
    MatDestroy(&invM);
    MatDestroy(&NMat);
    MatDestroy(&NDerivMat);

    Mat Laplacian;
	MatMatMult(BF1, DIV, MAT_INITIAL_MATRIX, 1, &Laplacian);
	MatScale(Laplacian, nu);

        //MatView(Laplacian, viewer_info);
    */
    Mat A, B;    // factor 4 = Number of Variables
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 6*Np+1+3*Np*Np,  NULL, &A);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 6*Np+1+3*Np*Np,  NULL, &B);
    std::cout << " Global Matrices Preallocated" << std::endl;
    for (unsigned int i = 0; i < 2*N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        MatGetRow(BF1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            //if (dummy!=0.0)
            if (abs(dummy)>1e-12)
            {
                MatSetValue(A, 	i,     3*N_Nodes+cols[j],    -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     3*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	    ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);
    }
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;

        // Fill Diagonals
        MatSetValue(A, i, i, 1.0, ADD_VALUES);
        MatSetValue(A, N_Nodes+i, N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(A, 2*N_Nodes+i, 2*N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(A, 3*N_Nodes+i, 3*N_Nodes+i, 1.0, ADD_VALUES);

        MatSetValue(B, i, i, 1.0, ADD_VALUES);
        MatSetValue(B, N_Nodes+i, N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(B, 2*N_Nodes+i, 2*N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(B, 3*N_Nodes+i, 3*N_Nodes+i, 1.0, ADD_VALUES);

        /*
        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        */

        /*
        MatGetRow(BF2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     2*Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     2*Np*Number_Of_Elements+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	i,     Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     Np*Number_Of_Elements+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF2, i, &numberOfNonZeros, &cols, &values);
        * /
        MatGetRow(C, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     2*Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     2*Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	2*Np*Number_Of_Elements+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	2*Np*Number_Of_Elements+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(C2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	i,     2*Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     2*Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C2, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(D1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	Np*Number_Of_Elements+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	Np*Number_Of_Elements+i,     cols[j],    -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D1, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	i,     2*Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     2*Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);
        */
        MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            //if (dummy!=0.0)
            if (abs(dummy)>1e-12)
            {
                MatSetValue(A, 	3*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                //MatSetValue(A, 	Np*Number_Of_Elements+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	Np*Number_Of_Elements+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);

    }

    std::cout << " Global Matrices Constructed" << std::endl;
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
    /*
    MatDestroy(&BF2);
    MatDestroy(&C);
    MatDestroy(&C2);
    MatDestroy(&D1);
    MatDestroy(&D2);
    MatDestroy(&Laplacian);
    */

    //MatView(BF1, viewer);
    MatDestroy(&BF1);
    MatDestroy(&DIV);

    MatDestroy(&Ex);
    MatDestroy(&Ey);
    MatDestroy(&E);
    MatDestroy(&ET);
    MatDestroy(&M1);
    MatDestroy(&invM);
    MatDestroy(&invM_small);
    std::cout << "Store Global Matrices" << std::endl;


    std::cout << "Computing Initial Condition " << std::endl;
    // Initial Condition
    Vec Initial_Condition, VecU, VecW, VecR, VecP;
    // Size = sum_i Np_i, i = 1 .. Nel
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecU);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecW);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecR);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecP);

    for (auto k = List_Of_Elements2D.begin(); k < List_Of_Elements2D.end(); k++)
    {

        //unsigned int ID = (*k).getID()-1;
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int pos = (*k).getPosition();
        std::vector<double> xCoor, zCoor;
        xCoor = (*k).get_node_coordinates_x();
        zCoor = (*k).get_node_coordinates_y();
        int i = 0;
        double t = 0;
        for (unsigned int n = 0; n < Np; n++)
        {
            i = n;
            //std::cout << "i = " << i << ", pos = " << pos << ", xCoor[n] = " << xCoor[n] << ", zCoor[n] = " << zCoor[n] << std::endl;
            double value = Exact_Solution_mx_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecU, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, pos + i, value, INSERT_VALUES);

            value = Exact_Solution_mz_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecW, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, N_Nodes+ pos + i, value, INSERT_VALUES);

            value = Exact_Solution_r_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecR, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, 2*N_Nodes+ pos + i, value, INSERT_VALUES);

            value = Exact_Solution_p_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecP, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, 3*N_Nodes+ pos + i, value, INSERT_VALUES);
        }

    }
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    VecAssemblyBegin(VecU);
    VecAssemblyEnd(VecU);
    VecAssemblyBegin(VecW);
    VecAssemblyEnd(VecW);
    VecAssemblyBegin(VecR);
    VecAssemblyEnd(VecR);
    VecAssemblyBegin(VecP);
    VecAssemblyEnd(VecP);


    double H0 = calculate_Hamiltonian2D(M1_small, Initial_Condition, List_Of_Elements2D, N_Nodes);
    std::cout << "Initial Energy = " << std::setprecision(16) << H0 << std::endl;

    std::cout << "Start Simulations " << std::endl;




    //MatDestroy(&E);
    //MatDestroy(&ET);
    //MatDestroy(&Ex);
    MatDestroy(&ExT);
    //MatDestroy(&Ey);
    MatDestroy(&EyT);
    MatDestroy(&NMat);
    MatDestroy(&NDerivMat);
    //MatDestroy(&M1);
    MatDestroy(&M2);

    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    KSPSetUp(ksp);
    //KSPSetTolerances(ksp, 1e-12, 1e-12, 1e-12, PETSC_DEFAULT);
    KSPSetTolerances(ksp, 1e-12, 1e-12, 1e30, PETSC_DEFAULT);

    KSPSetType(ksp,KSPPREONLY);
    //KSPSetType(ksp,KSPCG);
    //KSPSetType(ksp,KSPGMRES);
    //KSPSetType(ksp,KSPBCGS);
    //KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    //PCSetType(pc,PCILU);
    //PCSetType(pc,PCNONE);
    //PCSetType(pc,PCSOR);

    KSPSetFromOptions(ksp);

    Vec Sol, QX;
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes, &Sol);
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes, &QX);
    VecCopy(Initial_Condition, Sol);
    double H1 = 0.0;

    //PetscPrintf(PETSC_COMM_SELF,"Size Global Matrices %6.4e\n",(double)sigma);
    // MatView(A, viewer_info);

    char szFileName[255] = {0};
    FILE *f = fopen("Energy.txt", "w");
    // Solve Linear System
    std::cout << "Start Time Stepping" << std::endl;
    double time = 0.0;


    //std::cout << "Initial Condition = " << std::endl;
    //VecView(Initial_Condition, viewer);

    ///std::cout << "P = " << std::endl;
    //MatView(A, viewer);
    ///MatView(A, viewer_dense);
    //std::cout << "Q = " << std::endl;
    //MatView(B, viewer);

    double Hold = H0;
    for (unsigned int t = 0; t < Number_Of_Periods*Number_Of_TimeSteps_In_One_Period; t++) //
    {

        fprintf(f, "%1.16e \t %1.16e\n", time, H1);
        time = (t+1)*DeltaT;
            MatMult(B, Sol, QX);
            KSPSolve(ksp, QX, Sol);

            H1 = calculate_Hamiltonian2D(M1_small, Sol, List_Of_Elements2D, N_Nodes);

            //std::cout << "Energy Diff= " << std::setprecision(16) << H1-Hold <<std::endl;
            Hold = H1;

            //std::cout << "Solution = " << std::endl;
            //VecView(Sol, viewer);
        /*
        PetscViewer viewer2;
        sprintf(szFileName, "solution%d.txt", t);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
        VecView(Sol, viewer2);
        PetscViewerDestroy(&viewer2);
        */

    }
    H1 = calculate_Hamiltonian2D(M1_small, Sol, List_Of_Elements2D, N_Nodes);

     ///std::cout << "Initial Condition = " << std::endl;
    ///VecView(Initial_Condition, viewer);
    ///std::cout << "Final Solution = " << std::endl;
    ///VecView(Sol, viewer);

    fprintf(f, "%1.16e \t %1.16e\n", time, H1);
    fclose(f);

    std::cout << "End Time Stepping" << std::endl;
    KSPDestroy(&ksp);
    VecDestroy(&QX);
    MatDestroy(&M1);
    MatDestroy(&M2);
    std::cout << "Initial Energy    = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Final Energy      = " << std::setprecision(16) << H1 << std::endl;
    std::cout << "Difference Energy = " << std::setprecision(16) << H1-H0 << std::endl;

    //MatView(A, PETSC_VIEWER_STDOUT_SELF);

    /*
    MatConvert(A, MATSEQDENSE,  MAT_INPLACE_MATRIX, &A);
    MatConvert(B, MATSEQDENSE,  MAT_INPLACE_MATRIX, &B);

    PetscViewer Aviewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrixA.txt", &Aviewer);
    MatView(A,Aviewer);
    PetscViewerDestroy(&Aviewer);
    PetscViewer Bviewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrixB.txt", &Bviewer);
    MatView(B,Bviewer);
    PetscViewerDestroy(&Bviewer);
*/
    double Enew = calculate_Error2D(Initial_Condition, Sol, 2, DeltaX, DeltaY, Np);
    PetscPrintf(PETSC_COMM_WORLD,"L2-Norm of error new %1.3e\n",(double)Enew);

/*
    PetscScalar Solution_Centroid_a[4*N_Elements];

    PetscScalar *Sol_a;
    VecGetArray(Sol, &Sol_a);

    PetscInt size_r;
    VecGetSize(Sol, &size_r);
    std::cout << "Size Sol = " << size_r << std::endl;
    for (auto k = List_Of_Elements2D.begin(); k < List_Of_Elements2D.end(); k++)
    {

        unsigned int ID = (*k).getID()-1;
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int pos = (*k).getPosition();
        unsigned int Order_Polynomials = (*k).get_Order_Of_Polynomials();
        Mat V2DInv = load_InverseVandermondeMatrix(Order_Polynomials);

        // Coefficients Lagrange Polynomial // Nodal
        PetscScalar SolElem_u[Np];
        PetscScalar SolElem_w[Np];
        PetscScalar SolElem_r[Np];
        PetscScalar SolElem_p[Np];
        for (unsigned int l = 0; l < Np; l++)
        {
            SolElem_u[l] = Sol_a[pos+l]; //u
            SolElem_w[l] = Sol_a[N_Nodes+pos+l]; //w
            SolElem_r[l] = Sol_a[2*N_Nodes+pos+l]; //r
            SolElem_p[l] = Sol_a[3*N_Nodes+pos+l]; //p
        }
        Vec CoeffLagU;
        Vec CoeffLagW;
        Vec CoeffLagR;
        Vec CoeffLagP;
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, Np, SolElem_u, &CoeffLagU);
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, Np, SolElem_w, &CoeffLagW);
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, Np, SolElem_r, &CoeffLagR);
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, Np, SolElem_p, &CoeffLagP);
        // Coefficients Legendre Polynomial // Modal
        Vec CoeffLegU;
        Vec CoeffLegW;
        Vec CoeffLegR;
        Vec CoeffLegP;
        VecDuplicate(CoeffLagU, &CoeffLegU);
        VecDuplicate(CoeffLagW, &CoeffLegW);
        VecDuplicate(CoeffLagR, &CoeffLegR);
        VecDuplicate(CoeffLagP, &CoeffLegP);
        MatMult(V2DInv,CoeffLagU,CoeffLegU);
        MatMult(V2DInv,CoeffLagW,CoeffLegW);
        MatMult(V2DInv,CoeffLagR,CoeffLegR);
        MatMult(V2DInv,CoeffLagP,CoeffLegP);


        PetscScalar *CoeffLegU_a;
        PetscScalar *CoeffLegW_a;
        PetscScalar *CoeffLegR_a;
        PetscScalar *CoeffLegP_a;
        VecGetArray(CoeffLegU, &CoeffLegU_a);
        VecGetArray(CoeffLegW, &CoeffLegW_a);
        VecGetArray(CoeffLegR, &CoeffLegR_a);
        VecGetArray(CoeffLegP, &CoeffLegP_a);

        double xc = (*k).get_centroid_x();
        double yc = (*k).get_centroid_y();

        Vec X, Y;
        PetscScalar X_a[1];
        PetscScalar Y_a[1];
        X_a[0] = xc;
        Y_a[0] = yc;
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, 1, X_a, &X);
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, 1, Y_a, &Y);
        Vec R, S;
        XYtoRS(X, Y, R, S);
        Vec A, B;
        RStoAB(R, S, A, B);

        std::cout << "B = " << std::endl;
        VecView(B, viewer);
        //Vec LengendreinCentroid;
        PetscScalar LengendreinCentroid_a[Np];
        unsigned int sk = 0;
        for (unsigned int i = 0; i <= Order_Polynomials; i++)
        {
            for (unsigned int j = 0; j <= Order_Polynomials-i; j++)
            {
                Vec P;
                P = Simplex2DP(A, B, i, j);
                PetscScalar *P_a;
                VecGetArray(P, &P_a);
                std::cout << "P = " << std::endl;
                VecView(P, viewer);
                LengendreinCentroid_a[sk] = P_a[0];
                VecRestoreArray(P, &P_a);
                VecDestroy(&P);
                sk += 1;
            }
        }
        double ValueInCentroidu = 0;
        double ValueInCentroidw = 0;
        double ValueInCentroidr = 0;
        double ValueInCentroidp = 0;
        for (unsigned int i = 0; i < Np; i++)
        {
            ValueInCentroidu += LengendreinCentroid_a[i]*CoeffLegU_a[i];
            std::cout << "LengendreinCentroid_a[i] = " << LengendreinCentroid_a[i] << std::endl;
            ValueInCentroidw += LengendreinCentroid_a[i]*CoeffLegW_a[i];
            ValueInCentroidr += LengendreinCentroid_a[i]*CoeffLegR_a[i];
            ValueInCentroidp += LengendreinCentroid_a[i]*CoeffLegP_a[i];
        }
        Solution_Centroid_a[ID] = ValueInCentroidu;
        Solution_Centroid_a[ID+N_Elements] = ValueInCentroidw;
        Solution_Centroid_a[ID+2*N_Elements] = ValueInCentroidr;
        Solution_Centroid_a[ID+3*N_Elements] = ValueInCentroidp;
        VecRestoreArray(CoeffLegU, &CoeffLegU_a);
        VecRestoreArray(CoeffLegW, &CoeffLegW_a);
        VecRestoreArray(CoeffLegR, &CoeffLegR_a);
        VecRestoreArray(CoeffLegP, &CoeffLegP_a);
        MatDestroy(&V2DInv);
        VecDestroy(&CoeffLagU);
        VecDestroy(&CoeffLagW);
        VecDestroy(&CoeffLagR);
        VecDestroy(&CoeffLagP);
        VecDestroy(&CoeffLegU);
        VecDestroy(&CoeffLegW);
        VecDestroy(&CoeffLegR);
        VecDestroy(&CoeffLegP);
        VecDestroy(&X);
        VecDestroy(&Y);
        VecDestroy(&R);
        VecDestroy(&S);
        VecDestroy(&A);
        VecDestroy(&B);
    }
    VecRestoreArray(Sol, &Sol_a);

    Vec Solution_Centroid;
    VecCreateSeqWithArray(PETSC_COMM_SELF,1, 4*N_Elements,Solution_Centroid_a, &Solution_Centroid);

    std::cout << "Computing Initial Condition " << std::endl;
    // Initial Condition
    Vec Initial_Condition_Centroid;
    // Size = sum_i Np_i, i = 1 .. Nel
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Elements,&Initial_Condition_Centroid);

    std::cout << "4*N_Elements = " << 4*N_Elements << std::endl;
    for (auto k = List_Of_Elements2D.begin(); k < List_Of_Elements2D.end(); k++)
    {

        unsigned int ID = (*k).getID()-1;
        double xc = (*k).get_centroid_x();
        double zc = (*k).get_centroid_y();
        double t = 0;

            double value = Exact_Solution_mx_2D_system1(xc, zc, t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(Initial_Condition_Centroid, ID, value, INSERT_VALUES);

            value = Exact_Solution_mz_2D_system1(xc, zc, t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(Initial_Condition_Centroid, N_Elements+ID, value, INSERT_VALUES);

            value = Exact_Solution_r_2D_system1(xc, zc, t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(Initial_Condition_Centroid, 2*N_Elements+ID, value, INSERT_VALUES);

            value = Exact_Solution_p_2D_system1(xc, zc, t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(Initial_Condition_Centroid, 3*N_Elements+ID, value, INSERT_VALUES);
    }
    VecAssemblyBegin(Initial_Condition_Centroid);
    VecAssemblyEnd(Initial_Condition_Centroid);


    double EnewCentroid = calculate_Error2D(Initial_Condition_Centroid, Solution_Centroid, 2, DeltaX, DeltaY, 1.0);
    PetscPrintf(PETSC_COMM_WORLD,"L2-Norm of error Centroid %1.3e\n",(double)Enew);

    VecDestroy(&Initial_Condition_Centroid);
    VecDestroy(&Solution_Centroid);



    std::cout << "Initial Condition Centroid = " << std::endl;
    VecView(Initial_Condition_Centroid, viewer);
    std::cout << "Final Solution Centroid = " << std::endl;
    VecView(Solution_Centroid, viewer);
    /*
    // Exact Solution
    Vec Exact_Solution;
    VecCreateSeq(PETSC_COMM_WORLD, 2*Number_Of_Elements*Np,&Exact_Solution);
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int ID = (*k).getID();
        unsigned int pos = ID*Np;
        std::vector<double> xCoor;
        xCoor = (*k).get_VertexCoordinates();
        int i = 0;
        for (auto c = xCoor.begin(); c < xCoor.end(); c++)
        {
            double value = Exact_Solution_m(*c, time, N2, sigma, kmode);
            VecSetValue(Exact_Solution, pos + i, value, INSERT_VALUES);
            value = Exact_Solution_p(*c, time, N2, sigma, kmode);
            VecSetValue(Exact_Solution, Number_Of_Elements*Np+pos + i, value, INSERT_VALUES);
            i++;
        }
    }
    VecAssemblyBegin(Exact_Solution);
    VecAssemblyEnd(Exact_Solution);

    PetscReal      norm2;
    VecAXPY(Sol2,-1.0,Exact_Solution);
    VecNorm(Sol2,NORM_2,&norm2);
    norm2 *= sqrt(DeltaX);
    PetscPrintf(PETSC_COMM_WORLD,"L2-Norm of error %1.9e\n",(double)norm2);
    VecDestroy(&Exact_Solution);
    */


    /*
    std::cout << "Computing Initial Condition " << std::endl;
    // Initial Condition
    Vec Initial_Condition, VecU, VecW, VecR, VecP;
    // Size = sum_i Np_i, i = 1 .. Nel
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecU);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecW);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecR);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecP);

    Mat Values, Mat HowOftenAdded;

    for (auto k = List_Of_Elements2D.begin(); k < List_Of_Elements2D.end(); k++)
    {

        //unsigned int ID = (*k).getID()-1;
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int pos = (*k).getPosition();
        std::vector<double> xCoor, zCoor;
        xCoor = (*k).get_node_coordinates_x();
        zCoor = (*k).get_node_coordinates_y();
        int i = 0;
        double t = 0;
        for (unsigned int n = 0; n < Np; n++)
        {
            i = n;
            std::cout << "i = " << i << ", pos = " << pos << ", xCoor[n] = " << xCoor[n] << ", zCoor[n] = " << zCoor[n] << std::endl;
            double value = Exact_Solution_mx_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecU, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, pos + i, value, INSERT_VALUES);

            value = Exact_Solution_mz_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecW, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, N_Nodes+ pos + i, value, INSERT_VALUES);

            value = Exact_Solution_r_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecR, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, 2*N_Nodes+ pos + i, value, INSERT_VALUES);

            value = Exact_Solution_p_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecP, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, 3*N_Nodes+ pos + i, value, INSERT_VALUES);
        }

    }
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    VecAssemblyBegin(VecU);
    VecAssemblyEnd(VecU);
    VecAssemblyBegin(VecW);
    VecAssemblyEnd(VecW);
    VecAssemblyBegin(VecR);
    VecAssemblyEnd(VecR);
    VecAssemblyBegin(VecP);
    VecAssemblyEnd(VecP);
    */



    VecDestroy(&Sol);
    //MatDestroy(&V);
    VecDestroy(&Initial_Condition);
    VecDestroy(&VecU);
    VecDestroy(&VecW);
    VecDestroy(&VecR);
    VecDestroy(&VecP);

    MatDestroy(&A);
    MatDestroy(&B);
    MatDestroy(&M1_small);


    auto t2 = std::chrono::high_resolution_clock::now();

    //std::cout << "Time Loop took "
    std::cout << "Execution took "
              << std::chrono::duration_cast<std::chrono::seconds>(t2-t0).count()
              << " seconds\n";
              }

/*
    std::vector<double> L2Error;
    L2Error.push_back(Number_Of_Elements);
    L2Error.push_back(N);
    L2Error.push_back(Enew);
    double Order = 0.0;
    if (Number_Of_Spatial_Steps > 0)
    {
        Order = log(Eold/Enew)/log(2);
    }
    L2Error.push_back(Order);
    L2Error.push_back(std::chrono::duration_cast<std::chrono::seconds>(t2-t0).count());
    L2Errors.push_back(L2Error);

    FILE *g = fopen("L2Errors.txt", "a");
    fprintf(g, "%d \t %d \t %1.16e \t %1.1e \t %1.1e \n", Number_Of_Elements, N, Enew, Order, std::chrono::duration_cast<std::chrono::milliseconds>(t2-t0).count()/1000.0);
    fclose(g);

    Eold = Enew;
    }
    }

    std::cout << "**********************************************************"<< std::endl;
    std::cout << "L2 Errors " << std::endl;
    std::cout << "h N Er O t " << std::endl;
    for (auto i = L2Errors.begin(); i != L2Errors.end(); ++i)
    {
        for (auto j = (*i).begin(); j != (*i).end(); ++j)
        {
            std::cout << *j <<  " ";
        }
        std::cout << std::endl;
    }
    //    std::cout << *i << ' ' << std::endl;
    std::cout << "**********************************************************"<< std::endl;

    auto t2 = std::chrono::high_resolution_clock::now();

    //std::cout << "Total took "
    std::cout << "Execution took "
              << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count()
              << " seconds\n";
*/
    PetscViewerDestroy(&viewer);
    PetscViewerDestroy(&viewer_dense);
    PetscViewerDestroy(&viewer_info);
    }
    PetscFinalize();
    return 1;
}

