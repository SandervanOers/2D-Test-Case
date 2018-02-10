#include "HIGW.hpp"
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix_local(const Mat &V, const double &deltaX, const unsigned int &N)
{
    //Mat M;
    //unsigned int Np = N + 1;

    Mat Product;
    Mat VV;
    //MatConvert(V, MATSEQDENSE, MAT_INITIAL_MATRIX, &VV);
    //MatMatTransposeMult(V, V, MAT_INITIAL_MATRIX, 1, &Product);
    //MatMatTransposeMult(VV, VV, MAT_INITIAL_MATRIX, 1, &Product);

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
    //MatMatMult(Vr, X, MAT_INITIAL_MATRIX, 1, &Dr);

    MatTransposeMatMult(X, X, MAT_INITIAL_MATRIX, 1.0, &Product);

    //MatView(V, PETSC_VIEWER_STDOUT_SELF);
    //MatView(X, PETSC_VIEWER_STDOUT_SELF);
    // X is the inverse
    //MatScale(X, 0.5*deltaX);
    //MatDestroy(&Product);
    MatDestroy(&A);
    MatDestroy(&B);
    MatDestroy(&X);
    ISDestroy(&row);
    ISDestroy(&col);

    return Product;
}
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix_inverse_local(const Mat &V, const double &deltaX)
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
    PetscInt size_s = Number_Of_Elements*Np;
    //VecGetSize(Solution, size_s);
    //PetscScalar XTemp[2*size_s];
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    //MatView(M1, PETSC_VIEWER_STDOUT_SELF);
    for (unsigned int k = 0; k < Number_Of_Elements; k++)
    {
        for (unsigned int i = 0; i < Np; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                double massij;
                const PetscInt idxm[1] = {i};
                const PetscInt idxn[1] = {j};

                MatGetValues(M1, 1, idxm, 1,  idxn, &massij);
                H += 0.5*massij*(XTemp[k*Np+i]*XTemp[k*Np+j]+XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[Np*Number_Of_Elements+k*Np+j]);
            }
        }
    }
    VecRestoreArray(Solution, &XTemp);
    return H;
}
/*--------------------------------------------------------------------------*/
