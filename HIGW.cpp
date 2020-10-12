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
extern double calculate_Hamiltonian2D(const Mat &M1, const Vec &Solution, const std::vector<Elements2D> &List_Of_Elements2D, const unsigned int &N_Nodes)
{
    double H = 0.0;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (auto k = List_Of_Elements2D.begin(); k < List_Of_Elements2D.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int pos = (*k).getPosition();

        // Get Submatrix from M1
        Mat M1_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, pos, 1, &isrow);
        MatGetSubMatrix(M1,isrow, isrow,MAT_INITIAL_MATRIX, &M1_Elemental);

        for (unsigned int i = 0; i < Np; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                double massij;
                const PetscInt idxm[1] = {i};
                const PetscInt idxn[1] = {j};

                MatGetValues(M1_Elemental, 1, idxm, 1,  idxn, &massij);
                H += 0.5*massij*(XTemp[pos+i]*XTemp[pos+j]+XTemp[N_Nodes+pos+i]*XTemp[N_Nodes+pos+j]+XTemp[3*N_Nodes+pos+i]*XTemp[3*N_Nodes+pos+j]);

            }
        }

        MatDestroy(&M1_Elemental);
        ISDestroy(&isrow);
    }
    VecRestoreArray(Solution, &XTemp);
    return H;
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
extern double calculate_Hamiltonian_pres(const Mat &M1, const Mat &M2, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np)
{

    double H = 0.0;
    //PetscInt size_s = Number_Of_Elements*Np;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (unsigned int k = 0; k < Number_Of_Elements; k++)
    {
        // Get Submatrix from M1 and M2
        Mat M1_Elemental, M2_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, k*Np, 1, &isrow);
        MatGetSubMatrix(M1,isrow, isrow,MAT_INITIAL_MATRIX, &M1_Elemental);
        MatGetSubMatrix(M2,isrow, isrow,MAT_INITIAL_MATRIX, &M2_Elemental);

        for (unsigned int i = 0; i < Np; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                double massij1, massij2;
                const PetscInt idxm[1] = {i};
                const PetscInt idxn[1] = {j};

                MatGetValues(M1_Elemental, 1, idxm, 1,  idxn, &massij1);
                MatGetValues(M2_Elemental, 1, idxm, 1,  idxn, &massij2);
                H += 0.5*massij1*(XTemp[k*Np+i]*XTemp[k*Np+j]+XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[Np*Number_Of_Elements+k*Np+j])+0.5*massij2*(XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[Np*Number_Of_Elements+k*Np+j]);
            }
        }
        MatDestroy(&M1_Elemental);
        MatDestroy(&M2_Elemental);
        ISDestroy(&isrow);
    }

    VecRestoreArray(Solution, &XTemp);
    return H;

}
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian_comp(const Mat &M1, const Mat &M2, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np)
{

    double H = 0.0;
    //PetscInt size_s = Number_Of_Elements*Np;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (unsigned int k = 0; k < Number_Of_Elements; k++)
    {
        // Get Submatrix from M1 and M2
        Mat M1_Elemental, M2_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, k*Np, 1, &isrow);
        MatGetSubMatrix(M1,isrow, isrow,MAT_INITIAL_MATRIX, &M1_Elemental);
        MatGetSubMatrix(M2,isrow, isrow,MAT_INITIAL_MATRIX, &M2_Elemental);

        for (unsigned int i = 0; i < Np; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                double massij1, massij2;
                const PetscInt idxm[1] = {i};
                const PetscInt idxn[1] = {j};

                MatGetValues(M1_Elemental, 1, idxm, 1,  idxn, &massij1);
                MatGetValues(M2_Elemental, 1, idxm, 1,  idxn, &massij2);
                H += 0.5*massij1*(XTemp[k*Np+i]*XTemp[k*Np+j]+XTemp[2*Np*Number_Of_Elements+k*Np+i]*XTemp[2*Np*Number_Of_Elements+k*Np+j])
                    +0.5*massij2*(XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[Np*Number_Of_Elements+k*Np+j]-XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[2*Np*Number_Of_Elements+k*Np+j]
                        -XTemp[Np*Number_Of_Elements+k*Np+j]*XTemp[2*Np*Number_Of_Elements+k*Np+i]+XTemp[2*Np*Number_Of_Elements+k*Np+i]*XTemp[2*Np*Number_Of_Elements+k*Np+j]);
            }
        }
        MatDestroy(&M1_Elemental);
        MatDestroy(&M2_Elemental);
        ISDestroy(&isrow);
    }

    VecRestoreArray(Solution, &XTemp);
    return H;

}
/*--------------------------------------------------------------------------*/
extern double calculate_Error(const Vec &Exact, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np, const double &DeltaX)
{
    PetscInt size_r;
    VecGetSize(Solution, &size_r);
    PetscScalar *ExactTemp, *SolutionTemp;
    VecGetArray(Exact, &ExactTemp);
    VecGetArray(Solution, &SolutionTemp);

    double error = 0.0;
    for (unsigned int i=0; i < size_r; i++)
    {
        if (Np > 1 && i > 0 && i < size_r-1 && (i+1)%(Np*Number_Of_Elements)!=0   && (i+1)%Np == 0 && i%Np==Np-1)
        {
            // Average over Interface
            SolutionTemp[i+1] = 0.5*(SolutionTemp[i]+SolutionTemp[i+1]);
        }
        else
        {
            error += (SolutionTemp[i]-ExactTemp[i])*(SolutionTemp[i]-ExactTemp[i]);
        }
    }
    error = sqrt(error);
    error *= sqrt(DeltaX);


    VecRestoreArray(Exact, &ExactTemp);
    VecRestoreArray(Solution, &SolutionTemp);

    return error;
}
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const double &DeltaX, const double &DeltaY, const unsigned int &Np)
{
    Vec Difference;
    VecDuplicate(Exact, &Difference);

    VecWAXPY(Difference, -1.0, Exact, Solution);

    PetscReal error;
    if (Norm_Type == 1)
    {
        VecNorm(Difference, NORM_1, &error);
    }
    else if (Norm_Type == 2)
    {
        VecNorm(Difference, NORM_2, &error);
    }
    else
    {
        VecNorm(Difference, NORM_INFINITY, &error);
    }

    VecDestroy(&Difference);
    error /= sqrt(Np);
    error *= sqrt(0.5*DeltaX*DeltaY);
    return error;

}
/*--------------------------------------------------------------------------*/
