#include "HIGW.hpp"
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix_local(const Mat &V)
{
    Mat Product;
    // Calculate inverse Vandermonde Matrix
    Mat A, B, X;
    MatDuplicate(V,MAT_COPY_VALUES,&A);
    MatDuplicate(V,MAT_DO_NOT_COPY_VALUES,&X);
    MatDuplicate(V,MAT_DO_NOT_COPY_VALUES,&B);
    MatConvert(B, MATSEQDENSE, MAT_INPLACE_MATRIX, &B);
    MatConvert(X, MATSEQDENSE, MAT_INPLACE_MATRIX, &X);

    MatShift(B, 1.0);
    MatOrderingType rtype = MATORDERINGNATURAL;
    IS row, col;
    MatGetOrdering(A, rtype, &row, &col);
    MatFactorInfo info;
    MatFactorInfoInitialize(&info);
    info.fill=1.0;
    MatLUFactor(A, row, col, &info);
    MatMatSolve(A, B, X);
    //MatView(X, PETSC_VIEWER_STDOUT_SELF);
    // X is the inverse

    //M = inv(V^T*V)
    MatTransposeMatMult(X, X, MAT_INITIAL_MATRIX, 1.0, &Product);
    MatDestroy(&A);
    MatDestroy(&B);
    MatDestroy(&X);
    ISDestroy(&row);
    ISDestroy(&col);

    return Product;
}
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix_inverse_local(const Mat &V)
{
    Mat M;
    MatMatTransposeMult(V, V, MAT_INITIAL_MATRIX, 1, &M);
    //MatScale(M, 2.0/deltaX);
    MatConvert(M, MATSEQDENSE, MAT_INPLACE_MATRIX, &M);
    return M;
}
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian(const Mat &M1, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np)
{

    double H = 0.0;
    //PetscInt size_s = Number_Of_Elements*Np;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (unsigned int k = 0; k < Number_Of_Elements; k++)
    {
        // Get Submatrix from M1
        Mat M1_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, k*Np, 1, &isrow);
        MatGetSubMatrix(M1,isrow, isrow,MAT_INITIAL_MATRIX, &M1_Elemental);
        for (unsigned int i = 0; i < Np; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                double massij;
                const PetscInt idxm[1] = {i};
                const PetscInt idxn[1] = {j};

                MatGetValues(M1_Elemental, 1, idxm, 1,  idxn, &massij);
                //massij = 1;
                H += 0.5*massij*(XTemp[k*Np+i]*XTemp[k*Np+j]+XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[Np*Number_Of_Elements+k*Np+j]);
            }
        }
        MatDestroy(&M1_Elemental);
        ISDestroy(&isrow);
    }

    VecRestoreArray(Solution, &XTemp);
    return H;

}
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian_vec(const Mat &M1, const Vec &Solution)
{
    //double H = 0.0;

    Vec SolSquared;
    VecDuplicate(Solution, &SolSquared);
    VecPointwiseMult(SolSquared, Solution, Solution);
    //Vec Energy_Elemental;
    //VecDuplicate(Solution, &Energy_Elemental);

    // M1: Np*Number_Of_Elements x Np*Number_Of_Elements
    // Solution: 2*Np*Numer_Of_Elements
    PetscInt size_r;
    VecGetSize(Solution, &size_r);
    size_r /= 2;
    IS stride;
    ISCreateStride(PETSC_COMM_SELF, size_r, 0, 1, &stride);
    Vec U;
    VecGetSubVector(SolSquared, stride, &U);
    Vec Energy_Elemental_U;
    VecDuplicate(U, &Energy_Elemental_U);
    MatMult(M1, U, Energy_Elemental_U);
    ISCreateStride(PETSC_COMM_SELF, size_r, size_r+1, 1, &stride);
    Vec P;
    VecGetSubVector(SolSquared, stride, &P);
    Vec Energy_Elemental_P;
    VecDuplicate(P, &Energy_Elemental_P);
    MatMult(M1, P, Energy_Elemental_P);
    VecDestroy(&SolSquared);

    VecAXPY(Energy_Elemental_U, 1.0, Energy_Elemental_P);
    PetscScalar H;
    VecSum(Energy_Elemental_U, &H);
    VecDestroy(&Energy_Elemental_U);
    VecDestroy(&Energy_Elemental_P);
    ISDestroy(&stride);

    H *= 0.5;
    return PetscRealPart(H);
}
/*--------------------------------------------------------------------------*/
