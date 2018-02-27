static char help[] = "Solves a 1D test case for a stratified fluid.\n\n";

#include <petscksp.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "initial_cond.hpp"
#include "mesh_gen.hpp"
#include "Legendre_Gauss_Lobatto.hpp"
#include "HIGW.hpp"
#include "Elements.hpp"
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

    PetscOptionsGetInt(NULL, NULL, "-n", &Number_Of_Elements_Petsc, NULL);
    PetscOptionsGetInt(NULL, NULL, "-k", &kmode, NULL);
    PetscOptionsGetInt(NULL, NULL, "-t", &Number_Of_TimeSteps_In_One_Period, NULL);
    PetscOptionsGetInt(NULL, NULL, "-P", &Number_Of_Periods, NULL);
    PetscOptionsGetInt(NULL, NULL, "-Method", &Method, NULL);
    PetscOptionsGetScalar(NULL, NULL, "-N2", &N2, NULL);
    PetscOptionsGetScalar(NULL, NULL, "-theta", &theta, NULL);
    PetscOptionsGetInt(NULL, NULL, "-Order", &N_Petsc, NULL);
    PetscOptionsGetInt(NULL, NULL, "-QuadratureAdded", &N_Q, NULL);

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

    //Vec VX;
    //VX = mesh_generation_1D_VX(xmin,xmax, Number_Of_Elements);
    //VecView(VX, PETSC_VIEWER_STDOUT_SELF);

    Mat EtoV;
    EtoV = mesh_generation_1D_EtoV(xmin, xmax, Number_Of_Elements);
    //MatView(EtoV, PETSC_VIEWER_STDOUT_SELF);

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
    //VecView(r, PETSC_VIEWER_STDOUT_SELF);

    Mat V;
    V = Vandermonde1D(r, N);
    //MatView(V, PETSC_VIEWER_STDOUT_SELF);

    Mat Dr;
    Dr = DMatrix1D(r, N, V);
    //MatView(Dr, PETSC_VIEWER_STDOUT_SELF);


    //Mat Lift;
    //Lift = Lift1D(N, V);
    //MatView(Lift, PETSC_VIEWER_STDOUT_SELF);

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
    MatView(x, PETSC_VIEWER_STDOUT_SELF);

    Mat J;
    J = GeometricFactors1D(x, Dr);

    //std::cout << "Jacobian Matrix Start" << std::endl;
    //MatView(J, PETSC_VIEWER_STDOUT_SELF);
    //std::cout << "Jacobian Matrix End" << std::endl;

    // Fx = 2 x Number_Of_Elements: physical coordinates of edge nodes
    // r does not need to be ordered for this
    PetscInt    ind_min, ind_max;
    PetscScalar val_min, val_max;

    VecMin(r, &ind_min, &val_min);
    VecMax(r, &ind_max, &val_max);

    PetscInt nc;
    const PetscScalar *row_min;
    const PetscScalar *row_max;
    MatGetRow(x,ind_min, &nc, NULL, &row_min);
    MatGetRow(x,ind_max, &nc, NULL, &row_max);
    Mat Fx;
    MatCreateSeqAIJ(PETSC_COMM_WORLD,2, nc, nc, NULL, &Fx);
    PetscInt ix[Number_Of_Elements];
    for (unsigned int k=0;k<=Number_Of_Elements; k++)
    {
        ix[k] = k;
    }
    PetscInt jx[1];
    jx[0]=0;
    MatSetValues(Fx, 1, jx, Number_Of_Elements, ix, row_min, INSERT_VALUES);
    jx[0]=1;
    MatSetValues(Fx, 1, jx, Number_Of_Elements, ix, row_max, INSERT_VALUES);
    MatRestoreRow(x,ind_min, &nc, NULL, &row_min);
    MatRestoreRow(x,ind_max, &nc, NULL, &row_max);
    MatAssemblyBegin(Fx, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Fx, MAT_FINAL_ASSEMBLY);
    //MatView(Fx, PETSC_VIEWER_STDOUT_SELF);
    MatDestroy(&Fx);
    VecDestroy(&r);

    MatGetRow(J,ind_min, &nc, NULL, &row_min);
    MatGetRow(J,ind_max, &nc, NULL, &row_max);
    Mat Fscale;
    MatCreateSeqAIJ(PETSC_COMM_WORLD,2, nc, nc, NULL, &Fscale);
    for (int k=0; k<nc; k++)
    {
        MatSetValue(Fscale, 0, k, 1/row_min[k], INSERT_VALUES);
        MatSetValue(Fscale, 1, k, 1/row_max[k], INSERT_VALUES);
    }
    MatRestoreRow(J,ind_min, &nc, NULL, &row_min);
    MatRestoreRow(J,ind_max, &nc, NULL, &row_max);
    MatAssemblyBegin(Fscale, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Fscale, MAT_FINAL_ASSEMBLY);
    //MatView(Fscale, PETSC_VIEWER_STDOUT_SELF);

    MatDestroy(&Fscale);

    Mat nx;
    nx = normals1D(N, Number_Of_Elements);
    /*
    //MatView(nx, PETSC_VIEWER_STDOUT_SELF);
    Mat EtoEF;
    EtoEF = FaceToFace_1D(Number_Of_Elements, EtoV);
    */
    MatDestroy(&EtoV);
    MatDestroy(&nx);


    MatDestroy(&Dr);
    MatDestroy(&x);
    MatDestroy(&J);
    MatDestroy(&nx);

    /*--------------------------------------------------------------------------*/
    /* Problem Specific */
    /*--------------------------------------------------------------------------*/

    PetscScalar   sigma;
    sigma = calculate_sigma(N2, kmode);
    PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);
    PetscScalar   DeltaX = 1.0/(double)Number_Of_Elements;
    std::cout << Number_Of_Elements << " => " << DeltaX << std::endl;
    PetscScalar DeltaT=1.0/(double)Number_Of_TimeSteps_In_One_Period/sigma;
    /// Check for CFL condition (when explicit)

    // Initial Condition
    Vec Initial_Condition, VecU, VecP;
    VecCreateSeq(PETSC_COMM_WORLD, 2*Number_Of_Elements*Np,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements*Np, &VecU);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements*Np, &VecP);
    /*
    for (unsigned int i=0; i<Number_Of_Elements*Np; i++)
    {
        PetscScalar coordinate;
        const PetscInt idxm[1]={(i%Np)}, idxn[1]={floor(i/Np)};
        MatGetValues(x, 1, idxm, 1, idxn, &coordinate);
        double value = Exact_Solution_m(PetscRealPart(coordinate), 0, N2, sigma, kmode);
        VecSetValue(Initial_Condition,i,value, INSERT_VALUES);
        VecSetValue(VecU,i,value, INSERT_VALUES);
        value = Exact_Solution_p(PetscRealPart(coordinate), 0, N2, sigma, kmode);
        VecSetValue(Initial_Condition,Np*Number_Of_Elements+i,value, INSERT_VALUES);
        VecSetValue(VecP,i,value, INSERT_VALUES);
    }
    */
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int ID = (*k).getID();
        unsigned int pos = ID*Np;
        std::vector<double> xCoor;
        xCoor = (*k).get_VertexCoordinates();
        int i = 0;
        for (auto c = xCoor.begin(); c < xCoor.end(); c++)
        {
            double value = Exact_Solution_m(*c, 0, N2, sigma, kmode);
            VecSetValue(VecU, pos + i, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, pos + i, value, INSERT_VALUES);
            value = Exact_Solution_p(*c, 0, N2, sigma, kmode);
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
    /*
    //MatView(M, viewer_dense);
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        (*k).set_MassMatrix(M);
    }


    Mat invM_Elemental;
    invM_Elemental = MassMatrix_inverse_local(V);
    //MatView(invM_Elemental, viewer_dense);
*/
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
        unsigned int pos = ID*Np;
        unsigned int Order_Polynomials = (*e).getOrderOfPolynomials();
        //unsigned int Order_Gaussian_Quadrature = ceil(Order_Polynomials+1+N_Q); // + higher order for rho_0 term
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
                    double rho0deriv = 0;
                    value_e += w[q]*rho0deriv*Li*Lj*DeltaX/2.0;

                    value_m += w[q]*Li*Lj/rho0*DeltaX/2.0;   /// Scale: det Jac = DeltaX/2.0
                }
                double factor = 1.0;
                //MatSetValue(E, pos+i, pos+j, factor*value_e, ADD_VALUES);
                value_e = - value_e;
                //MatSetValue(ET, pos+j, pos+i, factor*value_e, ADD_VALUES);

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

    MatView(M1, PETSC_VIEWER_STDOUT_SELF);

    MatDestroy(&invM_Elemental);

    //Mat GLL, GLR, GRL, GRR; // Global Flux Matrices
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &GLL);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &GLR);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &GRL);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements*Np, Number_Of_Elements*Np,2*Np, NULL, &GRR);
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

            //MatSetValue(GLL, posL+i, posL+j, rho0*(1-theta), ADD_VALUES);
            double factor = -1.0;
            MatSetValue(E, posL+i, posL+j, factor*rho0*(1-theta), ADD_VALUES);
            MatSetValue(ET, posL+j, posL+i, factor*-rho0*(1-theta), ADD_VALUES); /// Exchange i and j?

            //MatSetValue(GLR, posL+i, posR+0, -rho0*(1-theta), ADD_VALUES);
            MatSetValue(E, posL+i, posR+0, factor*-rho0*(1-theta), ADD_VALUES);
            MatSetValue(ET, posR+0, posL+i, factor*rho0*(1-theta), ADD_VALUES);

            //MatSetValue(GRL, posR+0, posL+j, rho0*theta, ADD_VALUES);
            MatSetValue(E, posR+0, posL+j, factor*rho0*theta, ADD_VALUES);
            MatSetValue(ET, posL+j, posR+0, factor*-rho0*theta, ADD_VALUES);

            //MatSetValue(GRR, posR+0, posR+0, -rho0*theta, ADD_VALUES);
            MatSetValue(E, posR+0, posR+0, factor*-rho0*theta, ADD_VALUES);
            MatSetValue(ET, posR+0, posR+0, factor*rho0*theta, ADD_VALUES);


        }
        //else
            //std::cout << "External" << std::endl;
    }

    //MatAssemblyBegin(GLL, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(GLL, MAT_FINAL_ASSEMBLY);
    //MatView(GLL, viewer_dense);
    //MatAssemblyBegin(GLR, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(GLR, MAT_FINAL_ASSEMBLY);
    //MatView(GLR, viewer_dense);
    //MatAssemblyBegin(GRL, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(GRL, MAT_FINAL_ASSEMBLY);
    //MatView(GRL, viewer_dense);
    //MatAssemblyBegin(GRR, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(GRR, MAT_FINAL_ASSEMBLY);
    //MatView(GRR, viewer_dense);

    std::cout << "Started Assembly DIV Matrices" << std::endl;
    MatAssemblyBegin(E, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(E, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY);

    MatView(E, viewer_dense);
    MatView(ET, viewer_dense);
    std::cout << "Finished Assembly DIV Matrices" << std::endl;

    //MatAXPY(E, 1.0, GLL, DIFFERENT_NONZERO_PATTERN);
    //MatAXPY(ET, -1.0, GLL, DIFFERENT_NONZERO_PATTERN);

    //MatAXPY(E, 1.0, GRR, DIFFERENT_NONZERO_PATTERN);
    //MatAXPY(ET, -1.0, GRR, DIFFERENT_NONZERO_PATTERN);

    /*
    MatAXPY(E, 1.0, GLR, DIFFERENT_NONZERO_PATTERN);
    Mat GLR_T, GRL_T;
    MatTranspose(GLR, MAT_INITIAL_MATRIX, &GLR_T);
    MatAXPY(E, -1.0, GLR_T, DIFFERENT_NONZERO_PATTERN);

    MatAXPY(E, 1.0, GRL, DIFFERENT_NONZERO_PATTERN);
    MatTranspose(GRL, MAT_INITIAL_MATRIX, &GRL_T);
    MatAXPY(E, -1.0, GRL_T, DIFFERENT_NONZERO_PATTERN);
    */
    //MatDestroy(&GLL);
    //MatDestroy(&GLR);
    //MatDestroy(&GRL);
    //MatDestroy(&GRR);

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

    MatDestroy(&E);
    MatDestroy(&ET);
    MatDestroy(&invM);

    Mat A, B;    // factor 2 = Number of Variables
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Np*Number_Of_Elements*2, Np*Number_Of_Elements*2, 3*Np+1,  NULL, &A);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Np*Number_Of_Elements*2, Np*Number_Of_Elements*2, 3*Np+1,  NULL, &B);
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

    double H0 = calculate_Hamiltonian(M1, Initial_Condition, Number_Of_Elements, Np);
    std::cout << "Initial Energy      = " << std::setprecision(16) << H0 << std::endl;

    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    //KSPSetOperators(ksp,AA,AA);
    KSPGetPC(ksp,&pc);
    KSPSetUp(ksp);
    KSPSetTolerances(ksp, 1e-12, 1e-12, 1e-12, PETSC_DEFAULT);

    //KSPSetType(ksp,KSPCG);
    //KSPSetType(ksp,KSPGMRES);
    //KSPSetType(ksp,KSPBCGS);
    KSPSetType(ksp,KSPPREONLY);
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
    //MatDuplicate(Initial_Condition,MAT_COPY_VALUES,&Sol);
    double H1 = 0.0;

    // Solve Linear System
    std::cout << "Start Time Stepping" << std::endl;
    for (unsigned int t = 0; t< Number_Of_Periods*Number_Of_TimeSteps_In_One_Period; t++)
    {
        MatMult(B, Sol, QX);
        KSPSolve(ksp, QX, Sol);
        H1 = calculate_Hamiltonian(M1, Sol, Number_Of_Elements, Np);
        std::cout << "Energy      = " << std::setprecision(16) << H1 << std::endl;
    }
    std::cout << "End Time Stepping" << std::endl;
    KSPDestroy(&ksp);
    //PCDestroy(&pc);
    VecDestroy(&QX);
    MatDestroy(&M1);
    MatDestroy(&Mk);


    std::cout << "Initial Energy    = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Final Energy      = " << std::setprecision(16) << H1 << std::endl;
    std::cout << "Difference Energy = " << std::setprecision(16) << H1-H0 << std::endl;
    MatView(A, PETSC_VIEWER_STDOUT_SELF);
    MatConvert(A, MATSEQDENSE,  MAT_INPLACE_MATRIX, &A);
    MatView(A, viewer_dense);

    PetscReal      norm;
    VecAXPY(Sol,-1.0,Initial_Condition);
    VecNorm(Sol,NORM_2,&norm);
    norm *= sqrt(DeltaX);
    PetscPrintf(PETSC_COMM_WORLD,"L2-Norm of error %1.9e\n",(double)norm);


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

