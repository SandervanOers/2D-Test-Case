static char help[] = "Solves a 1D test case for a viscous homogeneous fluid.\n\n";

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



int main(int argc,char **args)
{
    PetscInitialize(&argc,&args,(char*)0,help);
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

    unsigned int Number_Of_Elements = Number_Of_Elements_Petsc;
    unsigned int N = N_Petsc;

    unsigned int Np = N + 1;
    unsigned int Nfp = 1;
    unsigned int Nfaces = 2;

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
    MatDestroy(&x);
    MatDestroy(&Dr);
    VecDestroy(&r);
    /*--------------------------------------------------------------------------*/
    /* Problem Specific */
    /*--------------------------------------------------------------------------*/

    PetscScalar   sigma;
    //sigma = calculate_sigma(N2, kmode);
    sigma = viscous_calculate_sigma(nu, kmode);
    PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);
    PetscScalar   DeltaX = 1.0/(double)Number_Of_Elements;
    //Number_Of_TimeSteps_In_One_Period = 10*Number_Of_Elements*Number_Of_Elements;//*Number_Of_Elements;

    Number_Of_TimeSteps_In_One_Period = 100*pow(Number_Of_Elements, (Np+1)/2);
    PetscScalar DeltaT=1.0/(double)Number_Of_TimeSteps_In_One_Period/sigma;
    std::cout << Number_Of_Elements << " => " << DeltaX << std::endl;
    std::cout << Number_Of_TimeSteps_In_One_Period << " => " << DeltaT << std::endl;
    /// Check for CFL condition (when explicit)

    // Initial Condition
    Vec Initial_Condition, VecU, VecP;
    VecCreateSeq(PETSC_COMM_WORLD, 2*Number_Of_Elements*Np,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements*Np, &VecU);
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
            //double value = Exact_Solution_m(*c, t, N2, sigma, kmode);
            double value = viscous_Exact_Solution_m(*c, t, nu, sigma, kmode);
            VecSetValue(VecU, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, pos + i, value, INSERT_VALUES);
            //value = Exact_Solution_p(*c, t, N2, sigma, kmode);
            value = viscous_Exact_Solution_p(*c, t, nu, sigma, kmode);
            VecSetValue(VecP, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, Number_Of_Elements*Np+pos + i, value, INSERT_VALUES);
            i++;
        }
    }
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    VecAssemblyBegin(VecU);
    VecAssemblyEnd(VecU);
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

    Mat E, ET, invM, M1;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,3*Np, NULL, &E);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,3*Np, NULL, &ET);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &invM);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &M1);

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
                for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                {
                    //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dr
                    double Li = LagrangePolynomial(ri, qp[q], i);
                    double physical_x = (*e).get_xCoordinateLeft()+0.5*(1.0+qp[q])*((*e).get_xCoordinateRight()-(*e).get_xCoordinateLeft());
                    double rho0 = rho_0(physical_x, N2);
                    if (Order_Polynomials > 0)
                    {
                        double Ljderiv = LagrangePolynomialDeriv(ri, qp[q], j);
                        value_e += w[q]*rho0*Li*Ljderiv;
                    }
                    // w_q drho_0(r_q)/dr l_i(r_q) l_j(r_q)
                    double Lj = LagrangePolynomial(ri, qp[q], j);
                    double rho0deriv = rho_0_deriv(physical_x, N2);
                    value_e += w[q]*rho0deriv*Li*Lj*DeltaX/2.0;

                    value_m += w[q]*Li*Lj/rho0*DeltaX/2.0;
                }
                double factor = -1.0;
                MatSetValue(E, pos+i, pos+j, factor*value_e, ADD_VALUES);
                value_e = - value_e;
                MatSetValue(ET, pos+j, pos+i, factor*value_e, ADD_VALUES);

                MatSetValue(M1, pos+i, pos+j, value_m, ADD_VALUES);
            }
        }
        VecDestroy(&ri);
        VecDestroy(&QuadraturePoints);
        VecDestroy(&Weights);
        MatSetValues(invM,Np, in,Np, in, invM_values, INSERT_VALUES);
    }
    MatAssemblyBegin(M1, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(invM, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(invM, MAT_FINAL_ASSEMBLY);

    //MatView(M1, PETSC_VIEWER_STDOUT_SELF);
    //MatView(invM, PETSC_VIEWER_STDOUT_SELF);

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
            double rho0 = rho_0(physical_x, N2);
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

    Mat BF1, DIV;
    double fillBF = 1;
    Mat BF1_TEMP1, BF1_TEMP2;
    MatMatMult(E, invM, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP1);
    MatMatMult(BF1_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP2);
    MatMatMult(invM, BF1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF1);
    MatDestroy(&BF1_TEMP1);
    MatDestroy(&BF1_TEMP2);

    Mat DIV_TEMP1, DIV_TEMP2;
    MatMatMult(ET, invM, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP1);
    MatMatMult(DIV_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP2);
    MatMatMult(invM, DIV_TEMP2, MAT_INITIAL_MATRIX, fillBF, &DIV);
    MatDestroy(&DIV_TEMP1);
    MatDestroy(&DIV_TEMP2);

    Mat Laplacian;
	//MatMatMult(DIV, BF1, MAT_INITIAL_MATRIX, 1, &Laplacian);
	MatMatMult(BF1, DIV, MAT_INITIAL_MATRIX, 1, &Laplacian);
	MatScale(Laplacian, nu);

    MatDestroy(&E);
    MatDestroy(&ET);
    MatDestroy(&invM);

    //MatView(BF1, viewer_info);
    //MatView(DIV, viewer_info);
    //MatView(Laplacian, viewer_info);

    std::cout << 3*Np+1 << std::endl;
    std::cout << 3*Np+1+Np*Np << std::endl;
    Mat A, B;    // factor 2 = Number of Variables
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Np*Number_Of_Elements*2, Np*Number_Of_Elements*2, 3*Np+1+3*Np*Np,  NULL, &A);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Np*Number_Of_Elements*2, Np*Number_Of_Elements*2, 3*Np+1+3*Np*Np,  NULL, &B);
    std::cout << " Global Matrices Preallocated" << std::endl;

    for (unsigned int i = 0; i < Np*Number_Of_Elements; i++)
    {
        MatSetValue(A, i, i, 1.0, ADD_VALUES);
        MatSetValue(A, Np*Number_Of_Elements+i, Np*Number_Of_Elements+i, 1.0, ADD_VALUES);
        MatSetValue(B, i, i, 1.0, ADD_VALUES);
        MatSetValue(B, Np*Number_Of_Elements+i, Np*Number_Of_Elements+i, 1.0, ADD_VALUES);

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
                MatSetValue(A, 	i,     Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
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
    MatDestroy(&DIV);

    MatDestroy(&Laplacian);


    //MatView(A, viewer_info);

    double H0 = calculate_Hamiltonian(M1, Initial_Condition, Number_Of_Elements, Np);
    std::cout << "Initial Energy      = " << std::setprecision(16) << H0 << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
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
    VecCreateSeq(PETSC_COMM_WORLD, 2*Number_Of_Elements*Np, &Sol);
    VecCreateSeq(PETSC_COMM_WORLD, 2*Number_Of_Elements*Np, &QX);
    VecCopy(Initial_Condition, Sol);
    double H1 = 0.0;

    char szFileName[255] = {0};
    //PetscPrintf(PETSC_COMM_SELF,"Size Global Matrices %6.4e\n",(double)sigma);
    // MatView(A, viewer_info);

    // Solve Linear System
    std::cout << "Start Time Stepping" << std::endl;
    double time = 0.0;
    for (unsigned int t = 0; t < Number_Of_Periods*Number_Of_TimeSteps_In_One_Period; t++)
    {
        time = (t+1)*DeltaT;
            MatMult(B, Sol, QX);
            KSPSolve(ksp, QX, Sol);

        //H1 = calculate_Hamiltonian(M1, Sol, Number_Of_Elements, Np);
        //std::cout << "Energy = " << std::setprecision(16) << H1 ;
        //std::cout << "  Time = " << time << " " << std::endl;
        /*
        PetscViewer viewer2;
        sprintf(szFileName, "solution%d.txt", t);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
        VecView(Sol, viewer2);
        PetscViewerDestroy(&viewer2);
        */

    }
    H1 = calculate_Hamiltonian(M1, Sol, Number_Of_Elements, Np);
    std::cout << "End Time Stepping" << std::endl;
    KSPDestroy(&ksp);
    VecDestroy(&QX);
    MatDestroy(&M1);
    MatDestroy(&Mk);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Initial Energy    = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Final Energy      = " << std::setprecision(16) << H1 << std::endl;
    std::cout << "Difference Energy = " << std::setprecision(16) << H1-H0 << std::endl;

    std::cout << "Time Loop took "
              << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count()
              << " seconds\n";
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
            //double value = Exact_Solution_m(*c, time, N2, sigma, kmode);
            double value = viscous_Exact_Solution_m(*c, time, nu, sigma, kmode);
            VecSetValue(Exact_Solution, pos + i, value, INSERT_VALUES);
            //value = Exact_Solution_p(*c, time, N2, sigma, kmode);
            value = viscous_Exact_Solution_p(*c, time, nu, sigma, kmode);
            VecSetValue(Exact_Solution, Number_Of_Elements*Np+pos + i, value, INSERT_VALUES);
            i++;
        }
    }
    VecAssemblyBegin(Exact_Solution);
    VecAssemblyEnd(Exact_Solution);
    double Enew = calculate_Error(Exact_Solution, Sol, Number_Of_Elements, Np, DeltaX);
    Vec Sol2;
    VecDuplicate(Sol, &Sol2);
    VecCopy(Sol, Sol2);
    PetscReal      norm;
    VecAXPY(Sol,-1.0,Initial_Condition);
    VecNorm(Sol,NORM_2,&norm);
    norm *= sqrt(DeltaX);
    //PetscPrintf(PETSC_COMM_WORLD,"L2-Norm of error %1.9e\n",(double)norm);

    PetscReal      norm2;
    VecAXPY(Sol2,-1.0,Exact_Solution);
    VecNorm(Sol2,NORM_2,&norm2);
    norm2 *= sqrt(DeltaX);
    PetscPrintf(PETSC_COMM_WORLD,"L2-Norm of error %1.9e\n",(double)norm2);

    PetscPrintf(PETSC_COMM_WORLD,"L2-Norm of error new %1.3e\n",(double)Enew);
    VecDestroy(&Exact_Solution);

    VecDestroy(&Sol2);
    VecDestroy(&Sol);
    MatDestroy(&V);
    VecDestroy(&Initial_Condition);
    VecDestroy(&VecU);
    VecDestroy(&VecP);

    MatDestroy(&A);
    MatDestroy(&B);

    PetscViewerDestroy(&viewer);
    PetscViewerDestroy(&viewer_dense);
    PetscViewerDestroy(&viewer_info);

    PetscFinalize();
    return 1;
}

