#include "mesh_gen.hpp"
/*--------------------------------------------------------------------------*/
extern Vec mesh_generation_1D_VX(const double &xmin, const double &xmax, const unsigned int &Number_Of_Elements)
{
    Vec VX; // Vertex coordinates (physical)
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements+1,&VX);
    for (unsigned int i = 0; i <= Number_Of_Elements; i++)
    {
        double value = (xmax-xmin)*(double)(i)/Number_Of_Elements+xmin;
        VecSetValue(VX, i, value, INSERT_VALUES);
    }
    VecAssemblyBegin(VX);
    VecAssemblyEnd(VX);
    //VecView(VX, PETSC_VIEWER_STDOUT_SELF);

    return VX;
}
/*--------------------------------------------------------------------------*/
extern Mat mesh_generation_1D_EtoV(const double &xmin, const double &xmax, const unsigned int &Number_Of_Elements)
{
    Mat EtoV; // Element to Vector connectivity
    MatCreate(PETSC_COMM_WORLD,&EtoV);
    MatSetType(EtoV,MATSEQAIJ);
    MatSetSizes(EtoV, Number_Of_Elements, 2, PETSC_DECIDE, PETSC_DECIDE);
    MatSeqAIJSetPreallocation(EtoV, 2, NULL);
    for (unsigned int i = 0; i < Number_Of_Elements; i++)
    {
        MatSetValue(EtoV, i, 0, i, INSERT_VALUES);
        MatSetValue(EtoV, i, 1, i+1, INSERT_VALUES);
    }
    MatAssemblyBegin(EtoV, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EtoV, MAT_FINAL_ASSEMBLY);
    //MatView(EtoV, PETSC_VIEWER_STDOUT_SELF);

    return EtoV;
}
/*--------------------------------------------------------------------------*/
Mat FaceToFace_1D(const unsigned int &Number_Of_Elements, const Mat &EtoV)
{
    int Nfaces = 2;
    int Total_Faces = Nfaces*Number_Of_Elements;
    int Nv = Number_Of_Elements+1;

    // List of local face to local vertex connections
    ///int nv[2]={1,2};
    // Build global face to node
    Mat FtoV;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Total_Faces, Nv, 1, NULL, &FtoV);
    PetscInt sk = 0;

    PetscInt ir[1]={0};
    PetscInt ic[1]={0};
    for (unsigned int k=0; k<Number_Of_Elements; k++)
    {
        for (unsigned int f=0; f<Nfaces; f++)
        {
            ir[0] = k;
            ic[0] = f;
            PetscScalar val;
            MatGetValues(EtoV, 1, ir, 1, ic, &val);
            MatSetValue(FtoV, sk, val, 1.0, INSERT_VALUES);
            sk++;
        }
    }
    MatAssemblyBegin(FtoV, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(FtoV, MAT_FINAL_ASSEMBLY);

    Mat FtoF;
    MatMatTransposeMult(FtoV, FtoV, MAT_INITIAL_MATRIX, 1, &FtoF);
    MatShift(FtoF, -1.0);
    MatDestroy(&FtoV);

    PetscInt r_start, r_end, ncols;
    const PetscInt *cols;
    const PetscScalar *vals;
    PetscInt index=0;
    MatGetOwnershipRange(FtoF, &r_start, &r_end);
    PetscInt faces1[r_end-r_start], faces2[r_end-r_start];
    for (PetscInt row=r_start; row<r_end; row++)
    {
        MatGetRow(FtoF,row,&ncols,&cols,&vals);
        for (PetscInt i=0; i<ncols; i++)
        {
            // Removing the diagnoal (MatShift -1), yielded entries which are zero, but are present in the nonzero structure of the matrix
            if (abs(vals[i]) > 10.0*std::numeric_limits<double>::epsilon())
            {
                /// Added plus 1 here
                faces1[index] = cols[i]+1;
                faces2[index] = row+1;
                index++;
            }
        }
        MatRestoreRow(FtoF,row,&ncols,&cols,&vals);
    }
    MatDestroy(&FtoF);
    PetscInt element1[index], element2[index], face1[index], face2[index], ind[index];
    for (PetscInt i=0; i<index;i++)
    {
        element1[i]  = floor((faces1[i]-1)/Nfaces)+1;
        face1[i]     = (  (faces1[i]-1)% Nfaces)+1;
        element2[i]  = floor((faces2[i]-1)/Nfaces)+1;
        face2[i]     = ( ( faces2[i]-1)% Nfaces)+1;
        ind[i] = element1[i]+Number_Of_Elements*(face1[i]-1);
    }
    Mat EtoE, EtoF;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements, Nfaces, Nfaces, NULL, &EtoE);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements, Nfaces, Nfaces, NULL, &EtoF);

    for (PetscInt i=0; i<Number_Of_Elements; i++)
    {
        MatSetValue(EtoE, i, 0, i+1, INSERT_VALUES);
        MatSetValue(EtoE, i, 1, i+1, INSERT_VALUES);
        MatSetValue(EtoF, i, 0, 1, INSERT_VALUES);
        MatSetValue(EtoF, i, 1, 2, INSERT_VALUES);
    }
    MatAssemblyBegin(EtoE, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(EtoE,   MAT_FLUSH_ASSEMBLY);
    MatAssemblyBegin(EtoF, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(EtoF,   MAT_FLUSH_ASSEMBLY);

    for (PetscInt i=0; i<index;i++)
    {
        PetscInt row, col;
        /// Added minus 1 here
        row = (ind[i]-1)%Number_Of_Elements;
        col = floor((ind[i]-1)/Number_Of_Elements);
        MatSetValue(EtoE, row, col, element2[i], INSERT_VALUES);
        MatSetValue(EtoF, row, col, face2[i], INSERT_VALUES);
    }
    MatAssemblyBegin(EtoE,  MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EtoE,    MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(EtoF,  MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EtoF,    MAT_FINAL_ASSEMBLY);

    //MatView(EtoE, PETSC_VIEWER_STDOUT_SELF);
    //MatView(EtoF, PETSC_VIEWER_STDOUT_SELF);

    Mat EtoEF;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*Number_Of_Elements, Nfaces, Nfaces, NULL, &EtoEF);
    for (unsigned int row=0; row<Number_Of_Elements; row++)
    {
        MatGetRow(EtoE,row,&ncols,&cols,&vals);
        for (unsigned int col=0; col<ncols; col++)
        {
                MatSetValue(EtoEF, row, col, vals[col], INSERT_VALUES);
        }
        MatRestoreRow(EtoE,row,&ncols,&cols,&vals);
        MatGetRow(EtoF,row,&ncols,&cols,&vals);
        for (PetscInt col=0; col<ncols; col++)
        {
                MatSetValue(EtoEF, Number_Of_Elements+row, col, vals[col], INSERT_VALUES);
        }
        MatRestoreRow(EtoF,row,&ncols,&cols,&vals);
    }

    MatDestroy(&EtoE);
    MatDestroy(&EtoF);
    MatAssemblyBegin(EtoEF,  MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EtoEF,  MAT_FINAL_ASSEMBLY);

    return EtoEF;
}
