#include "mesh_gen.hpp"
/*--------------------------------------------------------------------------*/
void Compute_Vertex_Coordinates_Uniform_Rectangle_2D(const double &xmin, const double &xmax, const double &ymin, const double &ymax, const unsigned int &Number_Of_Elements_X, const unsigned int &Number_Of_Elements_Y, std::vector<VertexCoordinates2D> &List_Of_Vertices, std::vector<Boundaries2D> &List_Of_Boundaries, std::vector<Elements2D> &List_Of_Elements)
{
    // First Create Squares, then halve the squares to get triangles
    unsigned int ID_Vertices = 1;
    unsigned int ID_Boundaries = 1;
    unsigned int ID_Elements = 1;


    bool isInternal_Vertix = 0;
    bool isInternal_Boundary = 0;

    unsigned int Nx = Number_Of_Elements_X/2+1;
    unsigned int Ny = Number_Of_Elements_Y/2+1;

    //std::cout << "Nx = " << Nx << std::endl;
    //std::cout << "Ny = " << Ny << std::endl;
    unsigned int Total_Number_Of_Boundary_Elements = (Nx-1)*Ny+Nx*(Ny-1)+(Nx-1)*(Ny-1);
    unsigned int Total_Number_Of_Square_Elements = (Nx-1)*(Ny-1);
    //std::cout << "Total_Number_Of_Square_Elements = " << Total_Number_Of_Square_Elements << std::endl;
    unsigned int Total_Number_Of_Triangular_Elements = Number_Of_Elements_X/2*Number_Of_Elements_Y/2*2;
    //std::cout << "Total_Number_Of_Triangular_Elements = " << Total_Number_Of_Triangular_Elements << std::endl;

    //std::cout << "Total_Number_Of_Boundary_Elements = " << Total_Number_Of_Boundary_Elements << std::endl;
    for (unsigned int j = 1; j <= Ny; j++)
    {
        double yvalue1 = (ymax-ymin)*(double)(j-1)/(Ny-1)+ymin;
        for (unsigned int i = 1; i <= Nx; i++)
        {
            double xvalue1 = (xmax-xmin)*(double)(i-1)/(Nx-1)+xmin;

            if (yvalue1 == ymin || yvalue1 == ymax || xvalue1 == xmin || xvalue1 == xmax )
            {
                isInternal_Vertix = 0;
            }
            else
            {
                isInternal_Vertix = 1;
            }

            VertexCoordinates2D V(ID_Vertices, xvalue1, yvalue1, isInternal_Vertix);
            List_Of_Vertices.push_back(V);
            ID_Vertices++;
        }
    }

    for (unsigned int I = 1; I <= Total_Number_Of_Square_Elements; I++)
    {
        unsigned int Y = floor((I-1)/(Nx-1));
        unsigned int Is = (I-1)%(Nx-1) + 1;
        unsigned int S[4] = {I+Y, I+1+Y, I+Nx+Y, I+1+Nx+Y};

        unsigned int B1[3] = {Is+Y*(3*Nx-2), Is+Y*(3*Nx-2)+2*Nx-1, Is+Y*(3*Nx-2)+Nx-1};
        unsigned int B2[3] = {Is+(Y+1)*(3*Nx-2), Is+Y*(3*Nx-2)+2*Nx-1, Is+Y*(3*Nx-2)+Nx};

        //std::cout << "I = " << I << ", Y = " << Y << ", Is = " << Is << std::endl;
        //std::cout << "S = " << S[0] << " " << S[1] << " " << S[2] << " " << S[3] << std::endl;

        // top boundaries
        if (List_Of_Vertices[S[2]-1].getyCoordinate() == ymax)
        {
            // External
            Boundaries2D Boundary5(B2[0], 0, S[3], S[2], ID_Elements+1, -1);
            List_Of_Boundaries.push_back(Boundary5);
        }
        else
        {
            // Internal
            Boundaries2D Boundary5(B2[0], 1, S[3], S[2], ID_Elements+1, ID_Elements+2*(Nx-1));
            List_Of_Boundaries.push_back(Boundary5);
        }

        // right boundaries
        if (List_Of_Vertices[S[3]-1].getxCoordinate() == xmax)
        {
            // External
            Boundaries2D Boundary5(B2[2], 0, S[1], S[3], ID_Elements+1, -1);
            List_Of_Boundaries.push_back(Boundary5);
        }
        else
        {
            // Internal
            Boundaries2D Boundary5(B2[2], 1, S[1], S[3], ID_Elements+1, ID_Elements+2);
            List_Of_Boundaries.push_back(Boundary5);
        }
        // bottom boundaries
        if (List_Of_Vertices[S[1]-1].getyCoordinate() == ymin)
        {
            // External
            Boundaries2D Boundary5(B1[0], 0, S[0], S[1], ID_Elements, -1);
            List_Of_Boundaries.push_back(Boundary5);
        }
        else
        {
            // Internal
            //Boundaries2D Boundary5(B1[0], 1, S[0], S[1], ID_Elements, -1);
            //List_Of_Boundaries.push_back(Boundary5);
        }
        // left boundaries
        if (List_Of_Vertices[S[0]-1].getxCoordinate() == xmin)
        {
            Boundaries2D Boundary5(B1[2], 0, S[2], S[0], ID_Elements, -1);
            List_Of_Boundaries.push_back(Boundary5);
        }

        // Internal Boundaries
        // diagonal boundaries
        {
            Boundaries2D Boundary5(B1[1], 1, S[1], S[2], ID_Elements, ID_Elements+1);
            List_Of_Boundaries.push_back(Boundary5);
        }





        Elements2D T1(ID_Elements, B1[0], B1[1], B1[2], S[0], S[1], S[2], 3);
        ID_Elements++;
        List_Of_Elements.push_back(T1);
        Elements2D T2(ID_Elements, B2[0], B2[1], B2[2], S[3], S[2], S[1], 3);
        ID_Elements++;
        List_Of_Elements.push_back(T2);

    }

    std::sort(List_Of_Boundaries.begin(), List_Of_Boundaries.end());

}

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
    Mat EtoV; // Element to Vertex connectivity
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
    unsigned int Nfaces = 2;
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
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian(std::vector<Elements2D> &List_Of_Elements2D, const std::vector<VertexCoordinates2D> &List_Of_Vertices)
{
    for(auto i = List_Of_Elements2D.begin(); i < List_Of_Elements2D.end(); i++)
    {
        double x1 = List_Of_Vertices[(*i).getVertex_V1()-1].getxCoordinate();
        double y1 = List_Of_Vertices[(*i).getVertex_V1()-1].getyCoordinate();
        double x2 = List_Of_Vertices[(*i).getVertex_V2()-1].getxCoordinate();
        double y2 = List_Of_Vertices[(*i).getVertex_V2()-1].getyCoordinate();
        double x3 = List_Of_Vertices[(*i).getVertex_V3()-1].getxCoordinate();
        double y3 = List_Of_Vertices[(*i).getVertex_V3()-1].getyCoordinate();

        double dxdr = (x2-x1)/2.0;
        double dydr = (y2-y1)/2.0;
        double dxds = (x3-x1)/2.0;
        double dyds = (y3-y1)/2.0;

        double Jacobian = dxdr*dyds-dxds*dydr;
        (*i).setJacobian(Jacobian);
    }
}
/*--------------------------------------------------------------------------*/
void set_Order_Polynomials_Uniform(std::vector<Elements2D> &List_Of_Elements2D, const unsigned int &N)
{
    for(auto i = List_Of_Elements2D.begin(); i < List_Of_Elements2D.end(); i++)
    {
        (*i).set_Order_Of_Polynomials(N); //(rand() % 3) + 1
    }
}
/*--------------------------------------------------------------------------*/
extern unsigned int get_Number_Of_Nodes(const std::vector<Elements2D> &List_Of_Elements2D)
{
    unsigned int Number_Of_Nodes = 0;
    for(auto i = List_Of_Elements2D.begin(); i < List_Of_Elements2D.end(); i++)
    {
        Number_Of_Nodes += (*i).get_Number_Of_Nodes();
    }
    return Number_Of_Nodes;
}
/*--------------------------------------------------------------------------*/

