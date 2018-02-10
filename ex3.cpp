static char help[] = "Solves a 1D test case for a stratified fluid (DG(1)).\n\n";

#include <petscksp.h>
#include <iostream>
#include <iomanip>
#include "initial_cond.hpp"
#include "mesh_gen.hpp"
#include "Legendre_Gauss_Lobatto.hpp"
#include "HIGW.hpp"
//#include <slepceps.h>
//#include <slepcsys.h>



int main(int argc,char **args)
{
    PetscInitialize(&argc,&args,(char*)0,help);
    // Read in options from command line
    PetscInt   Number_Of_Elements=10, Number_Of_TimeSteps_In_One_Period=10, Method=1;
    PetscInt   Number_Of_Periods=1, k=1;
    PetscScalar N2 = 0.0;//1.0;
    PetscScalar   theta = 0.5;
    PetscInt    N = 1;

    PetscOptionsGetInt(NULL, NULL, "-n", &Number_Of_Elements, NULL);
    PetscOptionsGetInt(NULL, NULL, "-k", &k, NULL);
    PetscOptionsGetInt(NULL, NULL, "-t", &Number_Of_TimeSteps_In_One_Period, NULL);
    PetscOptionsGetInt(NULL, NULL, "-P", &Number_Of_Periods, NULL);
    PetscOptionsGetInt(NULL, NULL, "-Method", &Method, NULL);
    PetscOptionsGetScalar(NULL, NULL, "-N2", &N2, NULL);
    PetscOptionsGetScalar(NULL, NULL, "-theta", &theta, NULL);
    PetscOptionsGetInt(NULL, NULL, "-Order", &N, NULL);

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

    Vec VX;
    VX = mesh_generation_1D_VX(xmin,xmax, Number_Of_Elements);
    //VecView(VX, PETSC_VIEWER_STDOUT_SELF);

    Mat EtoV;
    EtoV = mesh_generation_1D_EtoV(xmin, xmax, Number_Of_Elements);
    //MatView(EtoV, PETSC_VIEWER_STDOUT_SELF);

    Vec r;
    r = JacobiGL(0, 0, N);
    //VecView(r, PETSC_VIEWER_STDOUT_SELF);

    Mat V;
    V = Vandermonde1D(r, N);
    //MatView(V, PETSC_VIEWER_STDOUT_SELF);

    Mat Dr;
    Dr = DMatrix1D(r, N, V);
    //MatView(Dr, PETSC_VIEWER_STDOUT_SELF);


    Mat Lift;
    Lift = Lift1D(N, V);
    //MatView(Lift, PETSC_VIEWER_STDOUT_SELF);

    // build coordinates of all the nodes
    //va = EToV(:,1)'; vb = EToV(:,2)';
    Vec va, vb;
    VecCreate(PETSC_COMM_WORLD,&va);
    VecSetSizes(va,PETSC_DECIDE,Number_Of_Elements);
    VecSetFromOptions(va);
    VecCreate(PETSC_COMM_WORLD,&vb);
    VecSetSizes(vb,PETSC_DECIDE,Number_Of_Elements);
    VecSetFromOptions(vb);
    MatGetColumnVector(EtoV, va, 0);
    MatGetColumnVector(EtoV, vb, 1);
    //VecView(vb, PETSC_VIEWER_STDOUT_SELF);

    // Affine mapping
    // (physical) coordinates of all nodes
    Mat x;
    MatCreate(PETSC_COMM_WORLD, &x);
    MatSetSizes(x, PETSC_DECIDE, PETSC_DECIDE, N+1, Number_Of_Elements);
    MatSetType(x,MATSEQAIJ);
    MatSeqAIJSetPreallocation(x, Number_Of_Elements, NULL);

    PetscScalar *VXa;
    PetscScalar *ra;
    PetscScalar *vaa;
    PetscScalar *vba;
    VecGetArray(VX, &VXa);
    VecGetArray(r, &ra);
    VecGetArray(va, &vaa);
    VecGetArray(vb, &vba);
    for (unsigned int k = 0; k < Number_Of_Elements; k++)
    {
        for (unsigned int n = 0; n <= N; n ++)
        {
            double val = VXa[(int)vaa[k]]+0.5*(ra[n]+1.0)*(VXa[(int)vba[k]]-VXa[(int)vaa[k]]);
            MatSetValue(x, n, k, val, INSERT_VALUES);
        }
    }
    VecRestoreArray(VX, &VXa);
    VecRestoreArray(r, &ra);
    VecRestoreArray(va, &vaa);
    VecRestoreArray(vb, &vba);

    MatAssemblyBegin(x, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(x, MAT_FINAL_ASSEMBLY);
    //MatView(x, PETSC_VIEWER_STDOUT_SELF);

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

    MatGetRow(J,ind_min, &nc, NULL, &row_min);
    MatGetRow(J,ind_max, &nc, NULL, &row_max);
    Mat Fscale;
    MatCreateSeqAIJ(PETSC_COMM_WORLD,2, nc, nc, NULL, &Fscale);
    for (unsigned int k=0; k<nc; k++)
    {
        MatSetValue(Fscale, 0, k, 1/row_min[k], INSERT_VALUES);
        MatSetValue(Fscale, 1, k, 1/row_max[k], INSERT_VALUES);
    }
    MatRestoreRow(J,ind_min, &nc, NULL, &row_min);
    MatRestoreRow(J,ind_max, &nc, NULL, &row_max);
    MatAssemblyBegin(Fscale, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Fscale, MAT_FINAL_ASSEMBLY);
    //MatView(Fscale, PETSC_VIEWER_STDOUT_SELF);

    Mat nx;
    nx = normals1D(N, Number_Of_Elements);
    //MatView(nx, PETSC_VIEWER_STDOUT_SELF);
    Mat EtoEF;
    EtoEF = FaceToFace_1D(Number_Of_Elements, EtoV);
    //MatView(EtoEF, PETSC_VIEWER_STDOUT_SELF);
    Mat EtoE, EtoF;
    IS isrow1, isrow2, iscol;
    ISCreateStride(PETSC_COMM_WORLD, Number_Of_Elements, 0,1,&isrow1);
    ISCreateStride(PETSC_COMM_WORLD, Nfaces,0, 1,&iscol);
    MatGetSubMatrix(EtoEF,isrow1,iscol,MAT_INITIAL_MATRIX,&EtoE);
    //MatView(EtoE, PETSC_VIEWER_STDOUT_SELF);
    ISCreateStride(PETSC_COMM_WORLD, Number_Of_Elements, Number_Of_Elements,1,&isrow2);
    MatGetSubMatrix(EtoEF,isrow2,iscol,MAT_INITIAL_MATRIX,&EtoF);
    //MatView(EtoF, PETSC_VIEWER_STDOUT_SELF);

    // BuildMaps
    PetscScalar *v;
    IS nodeids;
    ISCreateStride(PETSC_COMM_WORLD, Number_Of_Elements*(N+1), 1,1,&nodeids);
    const PetscInt *indices;
    ISGetIndices(nodeids,&indices);
    IS vmapM, vmapP;

    PetscInt aint[Nfp*Nfaces*Number_Of_Elements], bint[Nfp*Nfaces*Number_Of_Elements];
    for (unsigned int k=0; k<Number_Of_Elements; k++)
    {
        for (unsigned int f=0; f<Nfaces; f++)
        {
            PetscInt inde;
            if (f==0)
                inde = ind_min;
            else
                inde = ind_max;

            aint[f+k*Nfaces] = indices[inde+k*(N+1)];
        }
    }
    ISRestoreIndices(nodeids,&indices);
    ISDestroy(&nodeids);
    ISCreateGeneral(PETSC_COMM_WORLD, Nfp*Nfaces*Number_Of_Elements,aint, PETSC_COPY_VALUES, &vmapM);
    //ISView(vmapM, PETSC_VIEWER_STDOUT_SELF);
    for (PetscInt k=0; k<Number_Of_Elements; k++)
    {
        PetscInt ncEtoE, ncEtoF;
        const PetscScalar *rowEtoE;
        const PetscScalar *rowEtoF;
        MatGetRow(EtoE,k, &ncEtoE, NULL, &rowEtoE);
        MatGetRow(EtoF,k, &ncEtoF, NULL, &rowEtoF);
        for (PetscInt f=0; f<Nfaces; f++)
        {
            PetscInt k2, f2, vidM, vidP;
            /// Added minus one
            k2 = rowEtoE[f]-1;
            f2 = rowEtoF[f]-1;
            vidM = aint[f+k*Nfaces];
            vidP = aint[f2+k2*Nfaces];
            PetscScalar x1, x2;
            const PetscInt idxm[1]={f*N}, idxn[1]={k};
            MatGetValues(x, 1, idxm, 1, idxn, &x1);
            const PetscInt idxm2[1]={f*N}, idxn2[1]={k};
            MatGetValues(x, 1, idxm2, 1, idxn2, &x2);
            // Compute the distance
            PetscScalar D = (x1-x2)*(x1-x2);
            if (D<10.0*std::numeric_limits<double>::epsilon())
            {
                bint[f+k*Nfaces] = vidP;
            }
        }
        MatRestoreRow(EtoE,k, &ncEtoE, NULL, &rowEtoE);
        MatRestoreRow(EtoF,k, &ncEtoF, NULL, &rowEtoF);
    }
    ISCreateGeneral(PETSC_COMM_WORLD, Nfp*Nfaces*Number_Of_Elements,bint, PETSC_COPY_VALUES, &vmapP);
    //ISView(vmapP, PETSC_VIEWER_STDOUT_SELF);

    // Create list of boundary nodes
    PetscInt mapB[2];
    int a = 0;
    for (unsigned int i=0; i<Nfp*Nfaces*Number_Of_Elements; i ++)
    {
        if (abs(aint[i] - bint[i])<10.0*std::numeric_limits<double>::epsilon())
        {
            mapB[a] = aint[i];
            a++;
        }
    }
    // Create specific left (inflow) and right (outflow) maps
    PetscInt mapI = 1, mapO = Number_Of_Elements*Nfaces, vmapI = 1, vmapO = Number_Of_Elements*Np;

    /*--------------------------------------------------------------------------*/
    /* Problem Specific */
    /*--------------------------------------------------------------------------*/

    PetscScalar   sigma;
    sigma = calculate_sigma(N2, k);
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
    for (unsigned int i=0; i<Number_Of_Elements*Np; i++)
    {
        PetscScalar coordinate;
        const PetscInt idxm[1]={(i%Np)}, idxn[1]={floor(i/Np)};
        MatGetValues(x, 1, idxm, 1, idxn, &coordinate);
        double value = Exact_Solution_m(PetscRealPart(coordinate), 0, N2, sigma, k);
        VecSetValue(Initial_Condition,i,value, INSERT_VALUES);
        VecSetValue(VecU,i,value, INSERT_VALUES);
        value = Exact_Solution_p(PetscRealPart(coordinate), 0, N2, sigma, k);
        VecSetValue(Initial_Condition,Np*Number_Of_Elements+i,value, INSERT_VALUES);
        VecSetValue(VecP,i,value, INSERT_VALUES);
    }
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    VecAssemblyBegin(VecU);
    VecAssemblyEnd(VecU);
    VecAssemblyBegin(VecP);
    VecAssemblyEnd(VecP);
    //VecView(Initial_Condition,viewer);

    //  Local Matrices
    Mat M;
    M = MassMatrix_local(V, DeltaX, N);
    //MatView(M, viewer_dense);

    Mat invM;
    invM = MassMatrix_inverse_local(V, DeltaX);
    //MatView(invM, viewer_dense);
    //MatView(Dr, viewer_dense);




    //Mat Mk;
    //MatDuplicate(M, MAT_COPY_VALUES, &Mk);
    //MatScale(Mk, DeltaX/2.0);
    double H0 =0.0;
    H0 = calculate_Hamiltonian(M, Initial_Condition, Number_Of_Elements, Np);

    std::cout << "Initial Energy    = " << std::setprecision(16) << H0 << std::endl;
    //MatDestroy(&Mk);

    // Need inverse Jacobian: rx=1./J
    // to transform from the reference element r in [-,1,1] to the physical element x in[x_l, x_r]
    Mat RHS;
    Mat S;
    MatMatMult(M, Dr, MAT_INITIAL_MATRIX, 1.0, &S);
    //MatScale(S, 2.0/DeltaX);
    //MatView(S, viewer_dense);
    MatMatMult(invM, S, MAT_INITIAL_MATRIX, 1.0, &RHS);
    {
        //MatDuplicate(Dr,MAT_COPY_VALUES,&RHS);
        MatScale(RHS, 2.0/DeltaX); // 2.0/DeltaX = rx = 1/J
    }
    /*{
    Mat RHS_intermediate1;
    MatMatMult(M, S, MAT_INITIAL_MATRIX, 1.0, &RHS_intermediate1);
    Mat RHS_intermediate2;
    MatMatMult(RHS_intermediate1, invM, MAT_INITIAL_MATRIX, 1.0, &RHS_intermediate2);
    MatMatMult(RHS_intermediate2, invM, MAT_INITIAL_MATRIX, 1.0, &RHS);
    MatDestroy(&RHS_intermediate1);
    MatDestroy(&RHS_intermediate2);
    MatDestroy(&S);
    std::cout << "~~~~~~~~~~~~ RHS ~~~~~~~~~~~~" << std::endl;
    MatView(RHS, viewer_dense);
    }*/


    // Implicit Midpoint
    // Construct Global Matrices A and B
    //MatView(RHS, viewer_info);
    Mat RHS_T;
    Mat S_T;
    MatTranspose(S, MAT_INITIAL_MATRIX, &S_T);
    //MatTranspose(S, MAT_INITIAL_MATRIX, &RHS_T);
    MatMatMult(invM, S_T, MAT_INITIAL_MATRIX, 1.0, &RHS_T);
    MatScale(RHS_T, 2.0/DeltaX);
    MatScale(RHS_T, -1.0);



    std::cout << std::endl << std::endl << "Start Viewing Matrices " << std::endl;
    //MatView(V, viewer_dense);
    //MatView(M, viewer_dense);
    /*MatView(invM, viewer_dense);/*
    MatView(Dr, viewer_dense);
    MatView(S, viewer_dense);
    MatView(S_T, viewer_dense);
    MatView(RHS, viewer_dense);
    MatView(RHS_T, viewer_dense);
    std::cout << "End Viewing Matrices " << std::endl << std::endl<< std::endl;
*/

    MatDestroy(&S);
    MatDestroy(&S_T);

    Mat A, B;
    // factor 2 = Number of Variables
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Np*Number_Of_Elements*2, Np*Number_Of_Elements*2, Np+1,  NULL, &A);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, Np*Number_Of_Elements*2, Np*Number_Of_Elements*2, Np+1,  NULL, &B);
    std::cout << " Global Matrices Preallocated" << std::endl;

    PetscScalar rhs[Np*Np];
    PetscScalar rhs_t[Np*Np];
    PetscInt id[N+1];
    for (unsigned int n=0; n<Np; n++)
        id[n] = n;

    MatGetValues(RHS, Np, id, Np, id, rhs);
    MatGetValues(RHS_T, Np, id, Np, id, rhs_t);

    for (unsigned int k=0; k<Number_Of_Elements; k++)
    {
        PetscInt ie[Np], ief[N];
        for (unsigned int n=0; n<Np; n++)
        {
            ie[n] = n+k*Np;
            ief[n] = n+k*Np+Np*Number_Of_Elements;
            // We want to keep the nonzero pattern after assembly
            MatSetValue(A, n+k*Np, n+k*Np, 0.0, ADD_VALUES);
            MatSetValue(A, Np*Number_Of_Elements+n+k*Np, Np*Number_Of_Elements+n+k*Np, 0.0, ADD_VALUES);
            //MatSetValue(B, n+k*Np, n+k*Np, 1.0, ADD_VALUES);
            //MatSetValue(B, Np*Number_Of_Elements+n+k*Np, Np*Number_Of_Elements+n+k*Np, 1.0, ADD_VALUES);
        }
        MatSetValues(A, Np, ie, Np, ief, rhs, ADD_VALUES);
        MatSetValues(A, Np, ief, Np, ie, rhs_t, ADD_VALUES);
        //MatSetValues(B, Np, ie, Np, ief, rhs, ADD_VALUES);
        //MatSetValues(B, Np, ief, Np, ie, rhs_t, ADD_VALUES);

    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    //std::cout << " Matrix A before Scaling " << std::endl;
    //MatView(A, viewer_dense);
    MatScale(A, DeltaT/2.0);
    //std::cout << " Matrix A after Scaling " << std::endl;
    //MatView(A, viewer_dense);


    MatDuplicate(A, MAT_COPY_VALUES, &B);
    MatScale(B, -1.0);
    MatShift(A,1.0);
    MatShift(B,1.0);
    //MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

    //std::cout << "A" << std::endl;
    //MatView(A, viewer_dense);
    //std::cout << "B" << std::endl;
    //MatView(B, viewer_dense);

    //VecView(Initial_Condition, viewer_dense);
    //MatView(x, PETSC_VIEWER_STDOUT_SELF);
    //VecView(VecP, viewer_dense);


    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    KSPSetUp(ksp);

    //KSPSetType(ksp,KSPCG);
    KSPSetType(ksp,KSPGMRES);
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

    KSPGetPC(ksp,&pc);
    //PCSetType(pc,PCILU);
    PCSetType(pc,PCNONE);

    KSPSetFromOptions(ksp);

    Vec Sol, QX;
    VecCreateSeq(PETSC_COMM_WORLD, 2*Number_Of_Elements*Np, &Sol);
    VecCreateSeq(PETSC_COMM_WORLD, 2*Number_Of_Elements*Np, &QX);
    VecCopy(Initial_Condition, Sol);
    //MatDuplicate(Initial_Condition,MAT_COPY_VALUES,&Sol);
    double H1 = 0.0;
    // Solve Linear System
    for (unsigned int t = 0; t< Number_Of_Periods*Number_Of_TimeSteps_In_One_Period; t++)
    {


    //VecView(Initial_Condition, viewer_dense);
    MatMult(B, Sol, QX);

    //MatConvert(B, MATSEQDENSE, MAT_INITIAL_MATRIX, &B);
    //MatConvert(A, MATSEQDENSE, MAT_INITIAL_MATRIX, &A);
    //MatView(B, viewer_dense);
    //MatView(A, viewer_dense);
    //VecView(QX, viewer_dense);
    KSPSolve(ksp, QX, Sol);

    //VecView(Sol, viewer_dense);
    H1 = calculate_Hamiltonian(M, Sol, Number_Of_Elements, Np);
    std::cout << "Energy      = " << std::setprecision(16) << H1 << std::endl;

    //VecCopy()


    }
    //VecView(Initial_Condition, viewer_dense);
    //VecView(Sol, viewer_dense);
    KSPDestroy(&ksp);
    //PCDestroy(&pc);
    VecDestroy(&Sol);
    VecDestroy(&QX);


    std::cout << "Initial Energy    = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Final Energy      = " << std::setprecision(16) << H1 << std::endl;
    std::cout << "Difference Energy = " << std::setprecision(16) << H1-H0 << std::endl;

    VecDestroy(&VX);
    VecDestroy(&r);
    VecDestroy(&va);
    VecDestroy(&vb);
    MatDestroy(&EtoV);
    MatDestroy(&V);
    MatDestroy(&Dr);
    MatDestroy(&Lift);
    MatDestroy(&x);
    MatDestroy(&J);
    MatDestroy(&Fx);
    MatDestroy(&Fscale);
    MatDestroy(&nx);
    MatDestroy(&EtoEF);
    MatDestroy(&EtoF);
    MatDestroy(&EtoE);
    ISDestroy(&isrow1);
    ISDestroy(&isrow2);
    ISDestroy(&iscol);
    ISDestroy(&vmapM);
    ISDestroy(&vmapP);
    VecDestroy(&Initial_Condition);
    VecDestroy(&VecU);
    VecDestroy(&VecP);
    MatDestroy(&M);
    MatDestroy(&invM);

    MatDestroy(&A);
    MatDestroy(&B);

    //MatDestroy(&S);
    MatDestroy(&RHS);
    MatDestroy(&RHS_T);

    PetscViewerDestroy(&viewer);
    PetscViewerDestroy(&viewer_dense);
    PetscViewerDestroy(&viewer_info);

    PetscFinalize();
    return 1;
}

