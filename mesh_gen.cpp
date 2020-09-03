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

    std::cout << "Nx = " << Nx << std::endl;
    std::cout << "Ny = " << Ny << std::endl;
    unsigned int Total_Number_Of_Boundary_Elements = (Nx-1)*Ny+Nx*(Ny-1)+(Nx-1)*(Ny-1);
    unsigned int Total_Number_Of_Square_Elements = (Nx-1)*(Ny-1);
    unsigned int Total_Number_Of_Triangular_Elements = Number_Of_Elements_X*Number_Of_Elements_Y;

    std::cout << "Total_Number_Of_Boundary_Elements = " << Total_Number_Of_Boundary_Elements << std::endl;
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

    for (unsigned int B = 0; B < Total_Number_Of_Boundary_Elements; B++)
    {
        //Boundaries2D(unsigned int IDg, int Left_Element, int Right_Element, bool Internal, unsigned int ID_V1, unsigned int ID_V2)
    }

    for (unsigned int I = 1; I <= Total_Number_Of_Square_Elements; I++)
    {
        unsigned int Y = floor(I/Nx);
        unsigned int Is = (I-1)%(Nx-1) + 1;
        //std::cout << "I-1 = " << I-1 << ", Nx-1 = " << Nx-1 << ", rem = " << remainder(I-1, Nx-1) << ", % = " << (I-1)%(Nx-1) << std::endl;
        //std::cout << "I = " << I << ", Y = " << Y << ", Is = " << Is << std::endl;
        unsigned int S[4] = {I+Y, I+1+Y, I+Nx+Y, I+1+Nx+Y};
        // //         int S[4] = {I+Y, I+1+Y, I+1+Nx+Y, I+Nx+Y}; Counterclockwise
        //std::cout << "V" << ": " << S[0] << ", " << S[1] << ", " << S[2] << ", " << S[3] << std::endl;

        unsigned int B1[3] = {Is+Y*(3*Nx-2), Is+Y*(3*Nx-2)+2*Nx-1, Is+Y*(3*Nx-2)+Nx-1};
        unsigned int B2[3] = {Is+(Y+1)*(3*Nx-2), Is+Y*(3*Nx-2)+2*Nx-1, Is+Y*(3*Nx-2)+Nx};
        //std::cout << "B1" << ": " << B1[0] << ", " << B1[1] << ", " << B1[2]  << std::endl;
        //std::cout << "B2" << ": " << B2[0] << ", " << B2[1] << ", " << B2[2]  << std::endl;

        //Boundaries2D(ID_Boundaries, int Left_Element, int Right_Element, bool Internal, unsigned int ID_V1, unsigned int ID_V2);
        //ID_Boundaries++;
        //Boundaries2D(ID_Boundaries, int Left_Element, int Right_Element, bool Internal, unsigned int ID_V1, unsigned int ID_V2);
        //ID_Boundaries++;
        //Boundaries2D(ID_Boundaries, int Left_Element, int Right_Element, bool Internal, unsigned int ID_V1, unsigned int ID_V2);
        //ID_Boundaries++;
        //Boundaries2D(ID_Boundaries, int Left_Element, int Right_Element, bool Internal, unsigned int ID_V1, unsigned int ID_V2);
        //ID_Boundaries++;
        //Boundaries2D(ID_Boundaries, int Left_Element, int Right_Element, bool Internal, unsigned int ID_V1, unsigned int ID_V2);
        //ID_Boundaries++;

        if (List_Of_Vertices[S[0]-1].isInternal() + List_Of_Vertices[S[1]-1].isInternal() == 0)
        {
            isInternal_Boundary = 0;
        }
        else
        {
            isInternal_Boundary = 1;
        }
        Boundaries2D Boundary1(B1[0], isInternal_Boundary, S[0], S[1]);
        if (List_Of_Vertices[S[1]-1].isInternal() + List_Of_Vertices[S[2]-1].isInternal() == 0)
        {
            isInternal_Boundary = 0;
        }
        else
        {
            isInternal_Boundary = 1;
        }
        Boundaries2D Boundary2(B1[1], isInternal_Boundary, S[1], S[2]);
        if (List_Of_Vertices[S[2]-1].isInternal() + List_Of_Vertices[S[0]-1].isInternal() == 0)
        {
            isInternal_Boundary = 0;
        }
        else
        {
            isInternal_Boundary = 1;
        }
        Boundaries2D Boundary3(B1[2], isInternal_Boundary, S[2], S[0]);
        List_Of_Boundaries.push_back(Boundary1);
        List_Of_Boundaries.push_back(Boundary2);
        List_Of_Boundaries.push_back(Boundary3);

        // top boundaries

        // right boundaries

        // Boundary 14 should be internal, even though it is inbetween two boundary vertices

        // Rewrite for three types of boundary elements: horizontal, vertical and diagonal?


        Elements2D T1(ID_Elements, B1[0], B1[1], B1[2], S[0], S[1], S[2]);
        // Set Boundary : Element Data here
        ID_Elements++;
        List_Of_Elements.push_back(T1);
        Elements2D T2(ID_Elements, B2[0], B2[1], B2[2], S[3], S[2], S[1]);
        // Set Boundary : Element Data here
        ID_Elements++;
        List_Of_Elements.push_back(T2);

    }


    /*
    // Define first three vertices
    double xvalue0 = xmin;
    double yvalue0 = ymin;
    unsigned int ID_V1 = ID_Vertices;
    VertexCoordinates2D V1(ID_Vertices, xvalue0, yvalue0, isInternal_Vertix);
    List_Of_Vertices.push_back(V1);
    ID_Vertices++;

    double xvalue1 = (xmax-xmin)*1.0/Number_Of_Elements_X+xmin;
    unsigned int ID_V2 = ID_Vertices;
    VertexCoordinates2D V2(ID_Vertices, xvalue1, yvalue0, isInternal_Vertix);
    List_Of_Vertices.push_back(V2);
    ID_Vertices++;

    unsigned int ID_B1 = ID_Boundaries;
    Boundaries2D(ID_Boundaries, ID_Elements, -1.0, 0.0, ID_V1, ID_V2);
    ID_Boundaries++;

    double yvalue1 = (ymax-ymin)*1/Number_Of_Elements_Y+ymin;
    unsigned int ID_V3 = ID_Vertices;
    VertexCoordinates2D V3(ID_Vertices, xvalue0, yvalue1, isInternal_Vertix);
    List_Of_Vertices.push_back(V3);
    ID_Vertices++;

    unsigned int ID_B2 = ID_Boundaries;
    Boundaries2D(ID_Boundaries, ID_Elements, ID_E1+1, 1.0, ID_V2, ID_V3);
    ID_Boundaries++;

    unsigned int ID_B3 = ID_Boundaries;
    Boundaries2D(ID_Boundaries, ID_Elements, -1.0, 0.0, ID_V3, ID_V1);
    ID_Boundaries++;

    unsigned int ID_E1 = ID_Elements;
    Elements2D(ID_Elements, ID_B1, ID_B2, ID_B3, ID_V1, ID_V2, ID_V3);
    ID_Elements++;

    unsigned int ID_V4;

    for (unsigned int i = 1; i <= Number_Of_Elements_X; i++)
    {
        double xvalue1 = (xmax-xmin)*(double)(i)/Number_Of_Elements_X+xmin;


        for (unsigned int j = 1; j <= Number_Of_Elements_Y; j++)
        {
            double yvalue1 = (ymax-ymin)*(double)(j)/Number_Of_Elements_Y+ymin;

            if (yvalue1 == ymin || yvalue1 == ymax || xvalue1 == xmin || xvalue1 == xmax )
            {
                isInternal_Vertix = 0;
            }
            else
            {
                isInternal_Vertix = 1;
            }

            ID_V4 = ID_Vertices;
            VertexCoordinates2D V4(ID_Vertices, xvalue1, yvalue1, isInternal_Vertix);
            List_Of_Vertices.push_back(V4);
            ID_Vertices++;

            // Vertex 4 -- Boundary 3 -- Vertex 3
            // Boundary 4 Element 1 Boundary 5 Element 2 Boundary 2
            // Vertex 1 -- Boundary 1 -- Vertex 2



            // 1: Check Iterator Steps
            // 2: Define 2D Elements and Boundaries
            // Vertices V1 V2 V3 and V1 V2 V4 form two elements

        }
    }
    */
}
/*--------------------------------------------------------------------------*/

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
