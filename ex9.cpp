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



int main(int argc,char **args)
{
    PetscInitialize(&argc,&args,(char*)0,help);

    auto t1 = std::chrono::high_resolution_clock::now();
    // Read in options from command line
    PetscInt   Number_Of_Elements_Petsc=10, Number_Of_TimeSteps_In_One_Period=10, Method=1;
    PetscInt   Number_Of_Periods=1, kmode=1;
    PetscScalar N2 = 0.0;//1.0;
    PetscScalar   theta = 0.5;
    PetscInt    N_Petsc = 1, N_Q=0;
    PetscScalar nu = 0.0;

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



    std::vector<Elements2D> List_Of_Elements2D;
    std::vector<Boundaries2D> List_Of_Boundaries2D;
    std::vector<VertexCoordinates2D> List_Of_Vertices;

    Compute_Vertex_Coordinates_Uniform_Rectangle_2D(0,1, 0,1, 12, 6, List_Of_Vertices, List_Of_Boundaries2D, List_Of_Elements2D);

    std::cout << "List of Vertices "  << std::endl;
    std::cout << "ID : x y isInternal "  << std::endl;
    for(auto i = List_Of_Vertices.begin(); i < List_Of_Vertices.end(); i++)
        std::cout << (*i).getID() << ": " << (*i).getxCoordinate() << "  " << (*i).getyCoordinate() << " " << (*i).isInternal() << std::endl;

    std::cout << "List of Boundaries "  << std::endl;
    std::cout << "ID : isInternal LeftElement RightElement"  << std::endl;
    for(auto i = List_Of_Boundaries2D.begin(); i < List_Of_Boundaries2D.end(); i++)
        std::cout << (*i).getID() << ": " << (*i).isInternal()  << " " << (*i).getLeftElementID() << " " << (*i).getRightElementID() << std::endl;


    std::cout << "List of Elements "  << std::endl;
    std::cout << "ID : V1 V2 V3 : B1 B2 B3"  << std::endl;
    for(auto i = List_Of_Elements2D.begin(); i < List_Of_Elements2D.end(); i++)
        std::cout << (*i).getID() << ": " << (*i).getVertex_V1() << "  " << (*i).getVertex_V2() << " " << (*i).getVertex_V3() << ": " << (*i).getBoundary_B1() << " " <<  (*i).getBoundary_B2() << " " << (*i).getBoundary_B3() << std::endl;


    Vec X, Y;
    Nodes2D(10, X, Y);
    std::cout << "X = " << std::endl;
    VecView(X,viewer);
    std::cout << "Y = " << std::endl;
    VecView(Y, viewer);

    Vec R, S;
    XYtoRS(X, Y, R, S);
    std::cout << "R = " << std::endl;
    VecView(R,viewer);
    std::cout << "S = " << std::endl;
    VecView(S, viewer);

    VecDestroy(&X);
    VecDestroy(&Y);
    VecDestroy(&R);
    VecDestroy(&S);
    std::cout << "**********************************************************"<< std::endl;

    std::vector<std::vector<double>> L2Errors;
    FILE *g = fopen("L2Errors.txt", "w");
    fprintf(g, "1/h \t N \t L2 Error \t \t \t Order \t \t time [s] \n");
    fclose(g);
    for (unsigned int Number_Of_Polynomial_Steps = 0; Number_Of_Polynomial_Steps < 1; Number_Of_Polynomial_Steps++)
    {
    double Eold = 0;
    FILE *g = fopen("L2Errors.txt", "a");
    fprintf(g, "\n");
    fclose(g);
    for (unsigned int Number_Of_Spatial_Steps = 0; Number_Of_Spatial_Steps < 1; Number_Of_Spatial_Steps++) //11-Number_Of_Polynomial_Steps
    {
    auto t0 = std::chrono::high_resolution_clock::now();
    unsigned int Number_Of_Elements = pow(2, Number_Of_Spatial_Steps); //Number_Of_Elements_Petsc;
    unsigned int N = Number_Of_Polynomial_Steps; //N_Petsc;


    std::cout << "**********************************************************"<< std::endl;
    std::cout << "N = " << N << ", 1/h = " << Number_Of_Elements << std::endl;
    std::cout << "**********************************************************"<< std::endl;

    unsigned int Np = N + 1;
    unsigned int Nfp = 1;
    unsigned int Nfaces = 2;


    double xmin = 0;
    double xmax = 1;

    std::vector<Elements> List_Of_Elements;
    std::vector<Boundaries> List_Of_Boundaries;
    for (unsigned int i=0; i<Number_Of_Elements+1; i++)
    {
        double xCoordinate = (xmax-xmin)*(double)(i)/Number_Of_Elements+xmin;
        bool isInternal = 1;
        int leftID = i-1;
        int rightID = i;

        if (i==0)
        {
            isInternal = 0;
            leftID = -1;
        }
        if (i==Number_Of_Elements)
        {
            isInternal = 0;
            rightID = -1;
        }

        Boundaries B(i, leftID, rightID, isInternal, xCoordinate);
        List_Of_Boundaries.push_back(B);
    }

    for (unsigned int i=0; i<Number_Of_Elements; i++)
    {
        // Evenly-Spaced
        double left_xCoordinate = (xmax-xmin)*(double)(i)/Number_Of_Elements+xmin;
        double right_xCoordinate = (xmax-xmin)*(double)(i+1)/Number_Of_Elements+xmin;
        double Jacobian = (right_xCoordinate-left_xCoordinate)/2.0;
        Elements E(i, Jacobian, i, i+1, left_xCoordinate, right_xCoordinate, N);
        List_Of_Elements.push_back(E);
    }

    /// Different Order N for different Elements
    /// Call JacobiGL for different N values
    Vec r;
    r = JacobiGL(0, 0, N);

    Mat V;
    V = Vandermonde1D(r, N);

    PetscInt MS, NS;
    MatGetSize(V,&MS,&NS);
    std::cout << "V Size = " << MS << " x " << NS << std::endl;

    Mat Dr;
    Dr = DMatrix1D(r, N, V);

    // Affine mapping
    // (physical) coordinates of all nodes
    Mat x;
    MatCreate(PETSC_COMM_WORLD, &x);
    MatSetSizes(x, PETSC_DECIDE, PETSC_DECIDE, N+1, Number_Of_Elements);
    MatSetType(x,MATSEQAIJ);
    MatSeqAIJSetPreallocation(x, Number_Of_Elements, NULL);

    PetscScalar *raa;
    VecGetArray(r, &raa);
    std::cout << "x Coordinates Nodes " << std::endl;
    int i = 0;
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k ++)
    {
        for (unsigned int n = 0; n <= N; n ++)
        {
            double val = (*k).get_xCoordinateLeft()+0.5*(raa[n]+1.0)*((*k).get_xCoordinateRight()-(*k).get_xCoordinateLeft());
            MatSetValue(x, n, i, val, INSERT_VALUES);
            (*k).set_VertexCoordinates(val);
        }
        i++;
    }
    VecRestoreArray(r, &raa);

    MatAssemblyBegin(x, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(x, MAT_FINAL_ASSEMBLY);
    //MatView(x, PETSC_VIEWER_STDOUT_SELF);
    MatDestroy(&Dr);
    VecDestroy(&r);

    MatConvert(x, MATSEQDENSE,  MAT_INPLACE_MATRIX, &x);
    PetscViewer xviewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "coordinates.txt", &xviewer);
    MatView(x,xviewer);
    PetscViewerDestroy(&xviewer);
    MatDestroy(&x);
    /*--------------------------------------------------------------------------*/
    /* Problem Specific */
    /*--------------------------------------------------------------------------*/

    PetscScalar   sigma;
    sigma = calculate_sigma_system1dcom(N2, kmode);
    PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);
    PetscScalar   DeltaX = 1.0/(double)Number_Of_Elements;
    Number_Of_TimeSteps_In_One_Period = 10*pow(Number_Of_Elements, (Np+1)/2);//Number_Of_Elements*Number_Of_Elements;
    PetscScalar DeltaT=1.0/(double)Number_Of_TimeSteps_In_One_Period/sigma;
    std::cout << Number_Of_Elements << " => " << DeltaX << std::endl;
    std::cout << Number_Of_TimeSteps_In_One_Period << " => " << DeltaT << std::endl;
    /// Check for CFL condition (when explicit)

    // Initial Condition
    Vec Initial_Condition, VecU, VecR, VecP;
    VecCreateSeq(PETSC_COMM_WORLD, 3*Number_Of_Elements*Np,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements*Np, &VecU);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements*Np, &VecR);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements*Np, &VecP);

    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int ID = (*k).getID();
        unsigned int pos = ID*Np;
        std::vector<double> xCoor;
        xCoor = (*k).get_VertexCoordinates();
        int i = 0;
        double t = 0;
        for (auto c = xCoor.begin(); c < xCoor.end(); c++)
        {
            double value = Exact_Solution_m_system1dcom(*c, t, N2, sigma, kmode);
            VecSetValue(VecU, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, pos + i, value, INSERT_VALUES);
            value = Exact_Solution_p_system1dcom(*c, t, N2, sigma, kmode);
            VecSetValue(VecP, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, 2*Number_Of_Elements*Np+pos + i, value, INSERT_VALUES);
            value = Exact_Solution_r_system1dcom(*c, t, N2, sigma, kmode);
            VecSetValue(VecR, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, Number_Of_Elements*Np+pos + i, value, INSERT_VALUES);
            i++;
        }
    }
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    VecAssemblyBegin(VecU);
    VecAssemblyEnd(VecU);
    VecAssemblyBegin(VecR);
    VecAssemblyEnd(VecR);
    VecAssemblyBegin(VecP);
    VecAssemblyEnd(VecP);

    /// In the future:
    /// Get Order of Polynomials from Element
    /// Calculate Mass Matrix for each element
    //  Local Matrices
    Mat M;
    M = MassMatrix_local(V);
    // Physical Mass Matrix
    Mat Mk;
    MatDuplicate(M, MAT_COPY_VALUES, &Mk);
    MatScale(Mk, DeltaX/2.0);
    MatDestroy(&M);

    MatDestroy(&Mk);
    // Inverse Mass Matrix (local): blocked
    Mat invM_Elemental;
    invM_Elemental = MassMatrix_inverse_local(V);
    MatScale(invM_Elemental, 2.0/DeltaX);
    PetscScalar invM_values[Np*Np];
    PetscInt in[Np];
    for (unsigned int n=0;n<Np; n++)
    {
        in[n] = n;
    }
    MatAssemblyBegin(invM_Elemental, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(invM_Elemental, MAT_FINAL_ASSEMBLY);
    MatGetValues(invM_Elemental, Np, in, Np, in, invM_values);

    Mat E, ET, invM, M1, M2, NMat, NDerivMat;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,3*Np, NULL, &E);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,3*Np, NULL, &ET);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &invM);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &M1);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &NMat);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &M2);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &NDerivMat);

    std::cout << " Start Elemental Calculations " << std::endl;
    for (auto e = List_Of_Elements.begin(); e < List_Of_Elements.end(); e++)
    {
        unsigned int ID = (*e).getID();
        unsigned int pos = ID*Np;   // Assumes Ordering
        unsigned int Order_Polynomials = (*e).getOrderOfPolynomials();
        unsigned int Order_Gaussian_Quadrature = ceil(Order_Polynomials+3+N_Q); // + higher order for rho_0 term

        for (unsigned int n=0;n<=Np; n++)
        {
            in[n] = n+pos;
        }
        // Quadrature Rules
        Vec Weights;
        Vec QuadraturePoints;
        QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
        //QuadraturePoints = JacobiGL_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);

        //VecView(Weights, PETSC_VIEWER_STDOUT_SELF);
        //VecView(QuadraturePoints, PETSC_VIEWER_STDOUT_SELF);
        PetscScalar *w, *qp;
        VecGetArray(Weights, &w);
        VecGetArray(QuadraturePoints, &qp);

        Vec ri;
        ri = JacobiGL(0, 0, Order_Polynomials);
        //std::cout << "VecView " << std::endl;
        //VecView(ri, PETSC_VIEWER_STDOUT_SELF);
        //VecView(QuadraturePoints, PETSC_VIEWER_STDOUT_SELF);
        for (unsigned int i = 0; i <= Order_Polynomials; i++)
        {
            for (unsigned int j = 0; j <= Order_Polynomials; j++)
            {
                // E Matrix
                double value_e = 0.0;
                double value_m = 0.0;
                double value_n = 0.0;
                double value_m2 = 0.0;
                double value_n_deriv = 0.0;
                for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                {
                    //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dr
                    double Li = LagrangePolynomial(ri, qp[q], i);
                    double physical_x = (*e).get_xCoordinateLeft()+0.5*(1.0+qp[q])*((*e).get_xCoordinateRight()-(*e).get_xCoordinateLeft());
                    double rho0 = rho_0_compressible(physical_x, N2);
                    if (Order_Polynomials > 0)
                    {
                        double Ljderiv = LagrangePolynomialDeriv(ri, qp[q], j);
                        value_e += w[q]*rho0*Li*Ljderiv;
                    }
                    // w_q drho_0(r_q)/dr l_i(r_q) l_j(r_q)
                    double Lj = LagrangePolynomial(ri, qp[q], j);
                    double rho0deriv = rho_0_deriv_compressible(physical_x, N2);
                    value_e += w[q]*rho0deriv*Li*Lj*DeltaX/2.0;

                    value_m += w[q]*Li*Lj/rho0*DeltaX/2.0;

                    value_n += w[q]*rho0*Li*Lj*DeltaX/2.0;

                    double N2_val = N_2_compressible(physical_x, N2); // N2 is actually rate: rho(-beta*x) => N2 = beta-1

                    value_m2 += w[q]*Li*Lj/rho0/N2_val*DeltaX/2.0;

                    value_n_deriv += w[q]*Li*Lj*rho0deriv*DeltaX/2.0; /// Should we rewrite this: second term value_e => add and subtract same term
                }
                double factor = -1.0;
                MatSetValue(E, pos+i, pos+j, factor*value_e, ADD_VALUES);
                value_e = - value_e;
                MatSetValue(ET, pos+j, pos+i, factor*value_e, ADD_VALUES);

                MatSetValue(M1, pos+i, pos+j, value_m, ADD_VALUES);
                MatSetValue(NMat, pos+i, pos+j, value_n, ADD_VALUES);
                MatSetValue(NDerivMat, pos+i, pos+j, value_n_deriv, ADD_VALUES);
                MatSetValue(M2, pos+i, pos+j, value_m2, ADD_VALUES);
            }
        }
        VecDestroy(&ri);
        VecDestroy(&QuadraturePoints);
        VecDestroy(&Weights);
        MatSetValues(invM,Np, in,Np, in, invM_values, INSERT_VALUES);
    }
    MatAssemblyBegin(M1, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(M2, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M2, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(NMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(NMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(NDerivMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(NDerivMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(invM, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(invM, MAT_FINAL_ASSEMBLY);

    MatDestroy(&invM_Elemental);

    for (auto f = List_Of_Boundaries.begin(); f < List_Of_Boundaries.end(); f++)
    {
        if ((*f).isInternal())
        {
            //std::cout << "Internal" << std::endl;
            int left = (*f).getLeftElementID();
            int right = (*f).getRightElementID();

            unsigned int posL = left*Np;    // Assumes Ordering
            unsigned int posR = right*Np;

            unsigned int Order_Polynomials_left = (List_Of_Elements[left]).getOrderOfPolynomials();
            unsigned int Order_Polynomials_right = (List_Of_Elements[right]).getOrderOfPolynomials();

            double physical_x = (*f).getxCoordinate();
            double rho0 = rho_0_compressible(physical_x, N2);
            unsigned int i = Order_Polynomials_left;
            unsigned int j = Order_Polynomials_left;

            double factor = 1.0;
            // GLL
            MatSetValue(E, posL+i, posL+j, factor*rho0*(1-theta), ADD_VALUES);
            MatSetValue(ET, posL+j, posL+i, factor*-rho0*(1-theta), ADD_VALUES);
            // GLR
            MatSetValue(E, posL+i, posR+0, factor*-rho0*(1-theta), ADD_VALUES);
            MatSetValue(ET, posR+0, posL+i, factor*rho0*(1-theta), ADD_VALUES);
            // GRL
            MatSetValue(E, posR+0, posL+j, factor*rho0*theta, ADD_VALUES);
            MatSetValue(ET, posL+j, posR+0, factor*-rho0*theta, ADD_VALUES);
            // GRR
            MatSetValue(E, posR+0, posR+0, factor*-rho0*theta, ADD_VALUES);
            MatSetValue(ET, posR+0, posR+0, factor*rho0*theta, ADD_VALUES);
        }
        //else
            //std::cout << "External" << std::endl;
    }

    std::cout << "Started Assembly DIV Matrices" << std::endl;
    MatAssemblyBegin(E, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(E, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY);
    //MatView(E, viewer_dense);
    //MatView(ET, viewer_dense);
    std::cout << "Finished Assembly DIV Matrices" << std::endl;
    double fillBF = 1;

    Mat BF1, BF1_TEMP1, BF1_TEMP2;
    MatMatMult(E, invM, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP1);
    MatMatMult(BF1_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP2);
    MatMatMult(invM, BF1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF1);
    MatDestroy(&BF1_TEMP1);
    MatDestroy(&BF1_TEMP2);

    Mat BF2, BF2_TEMP1, BF2_TEMP2;
    MatMatMult(E, invM, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP1);
    MatMatMult(BF2_TEMP1, M2, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP2);
    MatMatMult(invM, BF2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF2);
    MatDestroy(&BF2_TEMP1);
    MatDestroy(&BF2_TEMP2);

    Mat DIV, DIV_TEMP1, DIV_TEMP2;
    MatMatMult(ET, invM, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP1);
    MatMatMult(DIV_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP2);
    MatMatMult(invM, DIV_TEMP2, MAT_INITIAL_MATRIX, fillBF, &DIV);
    MatDestroy(&DIV_TEMP1);
    MatDestroy(&DIV_TEMP2);

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

    Mat A, B;    // factor 3 = Number of Variables
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Np*Number_Of_Elements*3, Np*Number_Of_Elements*3, 6*Np+1+3*Np*Np,  NULL, &A);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Np*Number_Of_Elements*3, Np*Number_Of_Elements*3, 6*Np+1+3*Np*Np,  NULL, &B);
    std::cout << " Global Matrices Preallocated" << std::endl;

    for (unsigned int i = 0; i < Np*Number_Of_Elements; i++)
    {
        MatSetValue(A, i, i, 1.0, ADD_VALUES);
        MatSetValue(A, Np*Number_Of_Elements+i, Np*Number_Of_Elements+i, 1.0, ADD_VALUES);
        MatSetValue(A, 2*Np*Number_Of_Elements+i, 2*Np*Number_Of_Elements+i, 1.0, ADD_VALUES);
        MatSetValue(B, i, i, 1.0, ADD_VALUES);
        MatSetValue(B, Np*Number_Of_Elements+i, Np*Number_Of_Elements+i, 1.0, ADD_VALUES);
        MatSetValue(B, 2*Np*Number_Of_Elements+i, 2*Np*Number_Of_Elements+i, 1.0, ADD_VALUES);

        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;

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

        MatGetRow(BF1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     2*Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     2*Np*Number_Of_Elements+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);
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
        */
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

        MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	2*Np*Number_Of_Elements+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	2*Np*Number_Of_Elements+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	Np*Number_Of_Elements+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	Np*Number_Of_Elements+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);


    }
    std::cout << " Global Matrices Constructed" << std::endl;
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
    MatDestroy(&BF1);
    MatDestroy(&BF2);
    MatDestroy(&DIV);
    MatDestroy(&C);
    MatDestroy(&C2);
    MatDestroy(&D1);
    MatDestroy(&D2);
    MatDestroy(&Laplacian);

    double H0 = calculate_Hamiltonian_comp(M1, M2, Initial_Condition, Number_Of_Elements, Np);
    std::cout << "Initial Energy      = " << std::setprecision(16) << H0 << std::endl;

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
    VecCreateSeq(PETSC_COMM_WORLD, 3*Number_Of_Elements*Np, &Sol);
    VecCreateSeq(PETSC_COMM_WORLD, 3*Number_Of_Elements*Np, &QX);
    VecCopy(Initial_Condition, Sol);
    double H1 = 0.0;

    //PetscPrintf(PETSC_COMM_SELF,"Size Global Matrices %6.4e\n",(double)sigma);
    // MatView(A, viewer_info);

    char szFileName[255] = {0};
    FILE *f = fopen("Energy.txt", "w");
    // Solve Linear System
    std::cout << "Start Time Stepping" << std::endl;
    double time = 0.0;
    for (unsigned int t = 0; t < Number_Of_Periods*Number_Of_TimeSteps_In_One_Period; t++)
    {

        H1 = calculate_Hamiltonian_comp(M1, M2, Sol, Number_Of_Elements, Np);
        fprintf(f, "%1.16e \t %1.16e\n", time, H1);
        time = (t+1)*DeltaT;
            MatMult(B, Sol, QX);
            KSPSolve(ksp, QX, Sol);

        //std::cout << "Energy Diff= " << std::setprecision(16) << H1-calculate_Hamiltonian(M1, Sol, Number_Of_Elements, Np) <<std::endl;

        /*
        PetscViewer viewer2;
        sprintf(szFileName, "solution%d.txt", t);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
        VecView(Sol, viewer2);
        PetscViewerDestroy(&viewer2);
        */

    }
    H1 = calculate_Hamiltonian_comp(M1, M2, Sol, Number_Of_Elements, Np);

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
    double Enew = calculate_Error(Initial_Condition, Sol, Number_Of_Elements, Np, DeltaX);
    PetscPrintf(PETSC_COMM_WORLD,"L2-Norm of error new %1.3e\n",(double)Enew);

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

    VecDestroy(&Sol);
    MatDestroy(&V);
    VecDestroy(&Initial_Condition);
    VecDestroy(&VecU);
    VecDestroy(&VecR);
    VecDestroy(&VecP);

    MatDestroy(&A);
    MatDestroy(&B);


    auto t2 = std::chrono::high_resolution_clock::now();

    //std::cout << "Time Loop took "
    std::cout << "Execution took "
              << std::chrono::duration_cast<std::chrono::seconds>(t2-t0).count()
              << " seconds\n";

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

    PetscViewerDestroy(&viewer);
    PetscViewerDestroy(&viewer_dense);
    PetscViewerDestroy(&viewer_info);
    PetscFinalize();
    return 1;
}

