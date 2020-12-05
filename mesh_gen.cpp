#include "mesh_gen.hpp"
/*--------------------------------------------------------------------------*/
void load_msh_mesh2D(const std::string &mesh_name, Vec &VX, Vec &VY, Mat &EToV, std::vector<VertexCoordinates2D> &List_Of_Vertices, std::vector<Squares2D> &List_Of_Elements, int &element_num,  int &node_num)
{
      int *element_node;
      //int element_num;
      int element_order;
      std::string gmsh_filename = mesh_name;
      int m;
      //int node_num;
      double *node_x;

    //
    //  Get the data size.
    //
      gmsh_size_read ( gmsh_filename, node_num, m, element_num,
        element_order );
    //
    //  Print the sizes.
    //
      std::cout << "\n";
      std::cout << "  Node data read from file \"" << gmsh_filename << "\"\n";
      std::cout << "  Number of vertices = " << node_num << "\n";
      std::cout << "  Spatial dimension = " << m << "\n";
      std::cout << "  Number of elements = " << element_num << "\n";
      //std::cout << "  Element order = " << element_order << "\n";
      std::cout << "\n";
    //
    //  Allocate memory.
    //
      //node_x = ( double * ) malloc ( m * node_num * sizeof ( double ) );
      //element_node = ( int * ) malloc ( element_order * element_num * sizeof ( int ) );

      node_x =  new  double[( m * node_num * sizeof ( double ) )];
      element_node =  new  int[( element_order * element_num * sizeof ( int ) )];
    //
    //  Get the data.
    //
      gmsh_data_read ( gmsh_filename, m, node_num, node_x, element_order, element_num, element_node );
    //
    //  Print some of the data.
    //
     // r8mat_transpose_print_some ( m, node_num, node_x, 1, 1, m, 10, "  Coordinates for first 10 nodes:" );

    unsigned int ID_Vertices = 0;
    unsigned int ID_Elements = 0;
    PetscScalar vx[node_num];
    PetscScalar vy[node_num];
    Vec VXT, VYT;
    for (int i = 0; i < node_num; i++)
    {
        vx[i] = node_x[2*i];
        vy[i] = node_x[2*i+1];
        VertexCoordinates2D V(ID_Vertices, node_x[2*i], node_x[2*i+1], 1);
        List_Of_Vertices.push_back(V);
        ID_Vertices++;
    }
    VecCreateSeqWithArray(PETSC_COMM_SELF,1,node_num,vx,&VXT);
    VecCreateSeqWithArray(PETSC_COMM_SELF,1,node_num,vy,&VYT);
    VecAssemblyBegin(VXT);
    VecAssemblyEnd(VXT);
    VecAssemblyBegin(VYT);
    VecAssemblyEnd(VYT);
    VecDuplicate(VXT, &VX);
    VecDuplicate(VYT, &VY);
    VecCopy(VXT, VX);
    VecCopy(VYT, VY);
    VecDestroy(&VXT);
    VecDestroy(&VYT);

    Mat EToVT;
    PetscScalar etov[element_num*4]; //4 = element_order
    for (unsigned int i = 0; i < element_num; i++)
    {
        //std::cout << element_node[4*i] << " " << element_node[4*i+1] << " " << element_node[4*i+2] << " "<< element_node[4*i+3] << " " << std::endl;
        //std::cout << element_node[4*i+3] << " " << element_node[4*i] << " " << element_node[4*i+1] << " "<< element_node[4*i+2] << " " << std::endl;

        double x1 = List_Of_Vertices[element_node[4*i]-1].getxCoordinate();
        double y1 = List_Of_Vertices[element_node[4*i]-1].getyCoordinate();
        double x2 = List_Of_Vertices[element_node[4*i+1]-1].getxCoordinate();
        double y2 = List_Of_Vertices[element_node[4*i+1]-1].getyCoordinate();
        double x3 = List_Of_Vertices[element_node[4*i+2]-1].getxCoordinate();
        double y3 = List_Of_Vertices[element_node[4*i+2]-1].getyCoordinate();
        double x4 = List_Of_Vertices[element_node[4*i+3]-1].getxCoordinate();
        double y4 = List_Of_Vertices[element_node[4*i+3]-1].getyCoordinate();

        // Bottom left is origin (0,0).
        double dist_1 = sqrt(x1*x1+y1*y1);
        double dist_2 = sqrt(x2*x2+y2*y2);
        double dist_3 = sqrt(x3*x3+y3*y3);
        double dist_4 = sqrt(x4*x4+y4*y4);

        double min = dist_1;
        int flag = 0;

        if(dist_2 < min)
        {
            min = dist_2;
            flag = 1;
        }
        if(dist_3 < min)
        {
            min = dist_3;
            flag = 2;
        }
        if(dist_4 < min)
        {
            min = dist_4;
            flag = 3;
        }

       //std::cout << element_node[4*i+flag]-1 << " " << element_node[4*i+(flag+1)%4]-1 << " " << element_node[4*i+(flag+2)%4]-1 << " " << element_node[4*i+(flag+3)%4]-1 << std::endl;
        Squares2D S(ID_Elements, element_node[4*i+flag]-1, element_node[4*i+(flag+1)%4]-1, element_node[4*i+(flag+2)%4]-1, element_node[4*i+(flag+3)%4]-1);
        ID_Elements++;
        List_Of_Elements.push_back(S);

        etov[4*i+1] = (PetscScalar)element_node[4*i+(flag+1)%4]-1;
        etov[4*i+2] = (PetscScalar)element_node[4*i+(flag+2)%4]-1;
        etov[4*i+3] = (PetscScalar)element_node[4*i+(flag+3)%4]-1;
        etov[4*i] = (PetscScalar)element_node[4*i+flag]-1;


       /* Assumes Ordering
        Squares2D S(ID_Elements, element_node[4*i+3]-1, element_node[4*i]-1, element_node[4*i+1]-1, element_node[4*i+2]-1);
        ID_Elements++;
        List_Of_Elements.push_back(S);

        etov[4*i+1] = (PetscScalar)element_node[4*i]-1;
        etov[4*i+2] = (PetscScalar)element_node[4*i+1]-1;
        etov[4*i+3] = (PetscScalar)element_node[4*i+2]-1;
        etov[4*i] = (PetscScalar)element_node[4*i+3]-1;
        */

        // Msh native Format
        //  3 -- 2
        //  |    |
        //  4 -- 1
        //
        // We want a counterclockwise ordering
        //
        // 3 -- 2
        // |    |
        // 0 -- 1
    }
    //MatCreateDense(PETSC_COMM_SELF, element_num, 4, element_num, 4, etov, &EToVT);
    MatCreateSeqDense(PETSC_COMM_SELF, 4, element_num, etov, &EToVT);
    //MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, element_num, 4, rows, columns, etov, &EToVT); //(PetscScalar*)element_node
    MatAssemblyBegin(EToVT, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EToVT, MAT_FINAL_ASSEMBLY);
    MatTranspose(EToVT, MAT_INPLACE_MATRIX, &EToVT);
    //std::cout << "EToVT = " << std::endl;
    //MatView(EToVT, PETSC_VIEWER_STDOUT_SELF);
    MatDuplicate(EToVT, MAT_COPY_VALUES, &EToV);
    MatDestroy(&EToVT);
    //  i4mat_transpose_print_some ( element_order, element_num, element_node, 1, 1, element_order, 10, "  Connectivity for first 10 elements:" );
    //
    //  Clean up.
    //
      delete[] element_node;
      delete[] node_x;
      //free [] element_node;
      //free [] node_x;


  ///
    }


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
        if (List_Of_Vertices[S[2]-1].getyCoordinate() == ymax) // Type == 1
        {
            // External
            Boundaries2D Boundary5(B2[0], 0, S[2], S[3], -1, ID_Elements+1, 1); // Changed
            List_Of_Boundaries.push_back(Boundary5);
        }
        else
        {
            // Internal
            Boundaries2D Boundary5(B2[0], 1, S[2], S[3], ID_Elements+2*(Nx-1), ID_Elements+1, 1); // Changed
            List_Of_Boundaries.push_back(Boundary5);
        }

        // right boundaries
        if (List_Of_Vertices[S[3]-1].getxCoordinate() == xmax) // Type == 3
        {
            // External
            Boundaries2D Boundary5(B2[2], 0, S[3], S[1], -1, ID_Elements+1, 3); // Changed
            List_Of_Boundaries.push_back(Boundary5);
        }
        else
        {
            // Internal
            Boundaries2D Boundary5(B2[2], 1, S[3], S[1], ID_Elements+2, ID_Elements+1, 3); // Changed
            List_Of_Boundaries.push_back(Boundary5);
        }
        // bottom boundaries
        if (List_Of_Vertices[S[1]-1].getyCoordinate() == ymin)
        {
            // External
            Boundaries2D Boundary5(B1[0], 0, S[0], S[1], ID_Elements, -1, 1); // Type == 1
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
            // External
            Boundaries2D Boundary5(B1[2], 0, S[2], S[0], ID_Elements, -1, 3); // Type == 3
            List_Of_Boundaries.push_back(Boundary5);
        }

        // Internal Boundaries
        // diagonal boundaries
        {
            Boundaries2D Boundary5(B1[1], 1, S[1], S[2], ID_Elements, ID_Elements+1, 2); // Type == 2
            List_Of_Boundaries.push_back(Boundary5);
        }





        Elements2D T1(ID_Elements, B1[0], B1[1], B1[2], S[0], S[1], S[2], 3);
        //List_Of_Boundaries[B1[0]-1].setType(1);
        //List_Of_Boundaries[B1[1]-1].setType(2);
        //List_Of_Boundaries[B1[2]-1].setType(3);
        ID_Elements++;
        List_Of_Elements.push_back(T1);
        Elements2D T2(ID_Elements, B2[0], B2[1], B2[2], S[3], S[2], S[1], 3);
        //List_Of_Boundaries[B2[0]-1].setType(1);
        //List_Of_Boundaries[B2[1]-1].setType(2);
        //List_Of_Boundaries[B2[2]-1].setType(3);
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
void Connect2D(const Mat &EToV, const unsigned int &Number_Of_Elements, const unsigned int &Number_Of_Vertices, Mat &EToE, Mat &EToF, std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries)
{
    unsigned int Nfaces = 4;
    unsigned int Total_Faces = Nfaces*Number_Of_Elements;

    // List of local face to local vertex connections
    // int nv[4]=[[1,2], [2,3], [3,4], [1,4]]
    // Build global face to node
    Mat FToV;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Total_Faces, Number_Of_Vertices, 2, NULL, &FToV);
    PetscInt sk = 0;

    PetscInt ir[1]={0};
    PetscInt ic[2]={0};
    for (unsigned int k=0; k<Number_Of_Elements; k++)
    {
        ir[0] = k;
        for (unsigned int f=0; f<Nfaces; f++)
        {
            if (f == 0)
            {
                ic[0] = 0;
                ic[1] = 1;
            }
            else if (f == 1)
            {
                ic[0] = 1;
                ic[1] = 2;
            }
            else if (f == 2)
            {
                ic[0] = 2;
                ic[1] = 3;
            }
            else if (f == 3)
            {
                ic[0] = 0;
                ic[1] = 3;
            }
            const PetscScalar *vals;
            MatGetRow(EToV, k, NULL, NULL, &vals);
            //std::cout << "val = " << vals[0] << " " << vals[1] << " " << vals[2] <<" " << vals[3] << std::endl;
            MatSetValue(FToV, sk, vals[ic[0]], 1.0, INSERT_VALUES);
            MatSetValue(FToV, sk, vals[ic[1]], 1.0, INSERT_VALUES);
            MatRestoreRow(EToV, k, NULL, NULL, &vals);
            sk++;
        }
    }
    MatAssemblyBegin(FToV, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(FToV, MAT_FINAL_ASSEMBLY);

    //std::cout << "FToV = " << std::endl;
    //MatView(FToV, PETSC_VIEWER_STDOUT_SELF);


    // Build global face to global face sparse array
    Mat FToF;
    MatMatTransposeMult(FToV, FToV, MAT_INITIAL_MATRIX, 1, &FToF);
    MatShift(FToF, -2.0);

    //std::cout << "FToF = " << std::endl;
    //MatView(FToF, PETSC_VIEWER_STDOUT_SELF);


    PetscInt ncols;
    const PetscInt *cols;
    const PetscScalar *vals;
    PetscInt InternalFaces = 0;
    PetscInt faces1[Total_Faces], faces2[Total_Faces];
    //std::cout << "faces1 faces2" << std::endl;
    for (PetscInt row=0; row<Total_Faces; row++)
    {
        MatGetRow(FToF,row,&ncols,&cols,&vals);
        for (PetscInt i=0; i<ncols; i++)
        {
            // std::cout << "vals[i] = " << vals[i] << std::endl;
            if (abs(vals[i]-2.0) < NODETOL)
            {
                faces2[InternalFaces] = row;
                faces1[InternalFaces] = cols[i];
                //std::cout << faces1[InternalFaces] << " " << faces2[InternalFaces] << std::endl;
                InternalFaces++;
            }
        }
        MatRestoreRow(FToF,row,&ncols,&cols,&vals);
    }
    std::cout << "Total_Faces = " << Total_Faces << std::endl;
    std::cout << "InternalFaces = " << InternalFaces << std::endl;

    //std::set<pairs, decltype(compare)> Set(compare);
    std::set<BoundaryInfo, decltype(compare_BoundaryInfo)> SetType(compare_BoundaryInfo);
    PetscInt element1[InternalFaces], element2[InternalFaces], face1[InternalFaces], face2[InternalFaces], ind[InternalFaces];
    //std::cout << "element1 element2 face1 face2 ind " << std::endl;
    for (unsigned int i=0; i<InternalFaces;i++)
    {
        element1[i]  = floor((faces1[i])/Nfaces);
        face1[i]     = (  (faces1[i])% Nfaces);
        element2[i]  = floor((faces2[i])/Nfaces);
        face2[i]     = ( ( faces2[i])% Nfaces);
        ind[i] = element1[i]+Number_Of_Elements*(face1[i]);

        //std::cout << element1[i] << " " << element2[i] << " " << face1[i] << " " << face2[i] << " " << ind[i] << std::endl;
        //Set.emplace(pairs{element1[i], element2[i]});
        BoundaryInfo BI(element1[i], element2[i], face1[i]+1, face2[i]+1);
        SetType.emplace(BI);
    }
    unsigned int ID_Boundary = 0;
    //std::cout << "Boundary: Size = " << Set.size() << std::endl;
    // std::cout << "Internal Boundaries from pair: " << std::endl;
    //for(const auto& it: Set)
    //{
        //InternalBoundariesSquares2D B(ID_Boundary, it.first, it.second);
        //ID_Boundary++;
        //List_Of_Boundaries.push_back(B);
        //std::cout << it.first << " " << it.second << std::endl;
    //}
    //std::cout << "Boundary Info: Size = " << SetType.size() << std::endl;
    for(const auto& it: SetType)
    {
        InternalBoundariesSquares2D B(ID_Boundary, it.left, it.right, it.faceleft, it.faceright);
        ID_Boundary++;
        List_Of_Boundaries.push_back(B);
        //std::cout << (it).left << " " << (it).right << " "<< (it).faceleft << " "<< (it).faceright <<  std::endl;
    }

    //MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements, Nfaces, Nfaces, NULL, &EToE);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Nfaces, Number_Of_Elements, Number_Of_Elements, NULL, &EToE);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Number_Of_Elements, Nfaces, Nfaces, NULL, &EToF);


    PetscInt ix[Nfaces];
    PetscScalar iv[Nfaces];
    for (PetscInt k=0;k<Nfaces; k++)
    {
        ix[k] = k;
        iv[k] = k;
    }
    for (PetscInt i=0; i<Number_Of_Elements; i++)
    {
        ir[0] = i;
        MatSetValues(EToF, 1, ir, Nfaces, ix, iv, INSERT_VALUES);
        MatSetValues(EToE, Nfaces, ix, 1, ir, iv, INSERT_VALUES);
    }

    for (PetscInt i=0; i<InternalFaces;i++)
    {
        PetscInt row, col;
        row = (ind[i])%Number_Of_Elements;
        col = floor((ind[i])/Number_Of_Elements);
        //MatSetValue(EToE, row, col, element2[i], INSERT_VALUES);
        MatSetValue(EToE, col, row, element2[i], INSERT_VALUES);
        MatSetValue(EToF, row, col, face2[i], INSERT_VALUES);
    }
    MatAssemblyBegin(EToF, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EToF,   MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(EToE, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EToE,   MAT_FINAL_ASSEMBLY);

    //MatView(EToE, PETSC_VIEWER_STDOUT_SELF);
    //MatView(EToF, PETSC_VIEWER_STDOUT_SELF);

    std::cout << "Number of Internal Boundaries = " << InternalFaces/2.0 << std::endl;

    MatDestroy(&FToV);
    MatDestroy(&FToF);


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

        double Area = 0.5*abs(x1*y2+x2*y3+x3*y1-x2*y1-x3*y2-x1*y3);

        double dxdr = (x2-x1)/2.0;
        double dydr = (y2-y1)/2.0;
        double dxds = (x3-x1)/2.0;
        double dyds = (y3-y1)/2.0;

        double centroid_x = (x1+x2+x3)/3;
        double centroid_y = (y1+y2+y3)/3;

        double Jacobian = dxdr*dyds-dxds*dydr;
        (*i).setJacobian(Jacobian);
        (*i).setArea(Area);
        (*i).setCentroidX(centroid_x);
        (*i).setCentroidY(centroid_y);

        double drdx = dyds/Jacobian;
        double drdy = -dxds/Jacobian;
        double dsdx = -dydr/Jacobian;
        double dsdy = dxdr/Jacobian;
        (*i).set_rx(drdx);
        (*i).set_ry(drdy);
        (*i).set_sx(dsdx);
        (*i).set_sy(dsdy);
    }
}
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Square(std::vector<Squares2D> &List_Of_Elements, const std::vector<VertexCoordinates2D> &List_Of_Vertices)
{

    // Assumes Square Quadrilaterals => J is constant
    for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
    {
        double x1 = List_Of_Vertices[(*i).getVertex_V1()].getxCoordinate();
        double y1 = List_Of_Vertices[(*i).getVertex_V1()].getyCoordinate();
        double x2 = List_Of_Vertices[(*i).getVertex_V2()].getxCoordinate();
        double y2 = List_Of_Vertices[(*i).getVertex_V2()].getyCoordinate();
        double x3 = List_Of_Vertices[(*i).getVertex_V3()].getxCoordinate();
        double y3 = List_Of_Vertices[(*i).getVertex_V3()].getyCoordinate();
        double x4 = List_Of_Vertices[(*i).getVertex_V4()].getxCoordinate();
        double y4 = List_Of_Vertices[(*i).getVertex_V4()].getyCoordinate();

        //double Area1 = 0.5*abs(x1*y2+x2*y3+x3*y1-x2*y1-x3*y2-x1*y3);
        //double Area2 = 0.5*abs(x1*y3+x3*y4+x4*y1-x3*y1-x4*y3-x1*y4);
        double Area = 0.5*abs(x1*y2+x2*y3+x3*y4-x2*y1-x3*y2-x4*y3-x1*y4);

        double dxdr = (x2+x3-x1-x4)/4.0;
        double dxds = (x3+x4-x2-x1)/4.0;
        double dydr = (y2+y3-y1-y4)/4.0;
        double dyds = (y3+y4-y2-y1)/4.0;

        double Jacobian = dxdr*dyds-dxds*dydr;
        (*i).setJacobian(Jacobian);

        double drdx = dyds/Jacobian;
        double drdy = -dxds/Jacobian;
        double dsdx = -dydr/Jacobian;
        double dsdy = dxdr/Jacobian;
        (*i).set_rx(drdx);
        (*i).set_ry(drdy);
        (*i).set_sx(dsdx);
        (*i).set_sy(dsdy);
    }
}
/*--------------------------------------------------------------------------*/
/*void Calculate_Jacobian_Quadrilateral(const Squares2D &Quad, const std::vector<VertexCoordinates2D> &List_Of_Vertices, const double &r_p, const double &s_p, double &Jacobian, double &drdx, double &drdy, double &dsdx, double &dsdy)
{
    //unsigned int Element_ID = Quad.getID();

    double x1 = List_Of_Vertices[Quad.getVertex_V1()].getxCoordinate();
    double y1 = List_Of_Vertices[Quad.getVertex_V1()].getyCoordinate();
    double x2 = List_Of_Vertices[Quad.getVertex_V2()].getxCoordinate();
    double y2 = List_Of_Vertices[Quad.getVertex_V2()].getyCoordinate();
    double x3 = List_Of_Vertices[Quad.getVertex_V3()].getxCoordinate();
    double y3 = List_Of_Vertices[Quad.getVertex_V3()].getyCoordinate();
    double x4 = List_Of_Vertices[Quad.getVertex_V4()].getxCoordinate();
    double y4 = List_Of_Vertices[Quad.getVertex_V4()].getyCoordinate();

    //double Area = 0.5*abs(x1*y2+x2*y3+x3*y4-x2*y1-x3*y2-x4*y3-x1*y4);

    double dxdr = (x2+x3-x1-x4)/4.0 + (x1-x2+x3-x4)*s_p;
    double dxds = (x3+x4-x2-x1)/4.0 + (x1-x2+x3-x4)*r_p;
    double dydr = (y2+y3-y1-y4)/4.0 + (y1-y2+y3-y4)*s_p;
    double dyds = (y3+y4-y2-y1)/4.0 + (y1-y2+y3-y4)*r_p;

    Jacobian = dxdr*dyds-dxds*dydr;

    drdx = dyds/Jacobian;
    drdy = -dxds/Jacobian;
    dsdx = -dydr/Jacobian;
    dsdy = dxdr/Jacobian;

}*/
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Boundaries_Square(const std::vector<Squares2D> &List_Of_Elements, std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<VertexCoordinates2D> &List_Of_Vertices)
{
    for(auto f = List_Of_Boundaries.begin(); f < List_Of_Boundaries.end(); f++)
    {
        int left = (*f).getLeftElementID();
        int right = (*f).getRightElementID();

        double left_x1 = List_Of_Vertices[List_Of_Elements[left].getVertex_V1()].getxCoordinate();
        double left_y1 = List_Of_Vertices[List_Of_Elements[left].getVertex_V1()].getyCoordinate();
        double left_x2 = List_Of_Vertices[List_Of_Elements[left].getVertex_V2()].getxCoordinate();
        double left_y2 = List_Of_Vertices[List_Of_Elements[left].getVertex_V2()].getyCoordinate();
        double left_x3 = List_Of_Vertices[List_Of_Elements[left].getVertex_V3()].getxCoordinate();
        double left_y3 = List_Of_Vertices[List_Of_Elements[left].getVertex_V3()].getyCoordinate();
        double left_x4 = List_Of_Vertices[List_Of_Elements[left].getVertex_V4()].getxCoordinate();
        double left_y4 = List_Of_Vertices[List_Of_Elements[left].getVertex_V4()].getyCoordinate();

        double right_x1 = List_Of_Vertices[List_Of_Elements[right].getVertex_V1()].getxCoordinate();
        double right_y1 = List_Of_Vertices[List_Of_Elements[right].getVertex_V1()].getyCoordinate();
        double right_x2 = List_Of_Vertices[List_Of_Elements[right].getVertex_V2()].getxCoordinate();
        double right_y2 = List_Of_Vertices[List_Of_Elements[right].getVertex_V2()].getyCoordinate();
        double right_x3 = List_Of_Vertices[List_Of_Elements[right].getVertex_V3()].getxCoordinate();
        double right_y3 = List_Of_Vertices[List_Of_Elements[right].getVertex_V3()].getyCoordinate();
        double right_x4 = List_Of_Vertices[List_Of_Elements[right].getVertex_V4()].getxCoordinate();
        double right_y4 = List_Of_Vertices[List_Of_Elements[right].getVertex_V4()].getyCoordinate();
        //std::cout << "Left Element:  (" << left_x1 << ", " << left_y1 << "); (" << left_x2 << ", " << left_y2 << "); (" << left_x3 << ", " << left_y3 << "); (" << left_x4 << ", " << left_y4 << ")" << std::endl;
        //std::cout << "Right Element: (" << right_x1 << ", " << right_y1 << "); (" << right_x2 << ", " << right_y2 << "); (" << right_x3 << ", " << right_y3 << "); (" << right_x4 << ", " << right_y4 << ")" << std::endl;

        double x[2];
        double y[2];
        unsigned int ID_Vertex[2];
        unsigned int index =0;
        if ((left_x1 == right_x1 && left_y1 == right_y1) || (left_x1 == right_x2 && left_y1 == right_y2) || (left_x1 == right_x3 && left_y1 == right_y3) || (left_x1 == right_x4 && left_y1 == right_y4))
        {
            x[index] = left_x1;
            y[index] = left_y1;
            ID_Vertex[index] = List_Of_Elements[left].getVertex_V1();
            index++;
        }
        if ((left_x2 == right_x1 && left_y2 == right_y1) || (left_x2 == right_x2 && left_y2 == right_y2) || (left_x2 == right_x3 && left_y2 == right_y3) || (left_x2 == right_x4 && left_y2 == right_y4))
        {
            x[index] = left_x2;
            y[index] = left_y2;
            ID_Vertex[index] = List_Of_Elements[left].getVertex_V2();
            index++;
        }
        if ((left_x3 == right_x1 && left_y3 == right_y1) || (left_x3 == right_x2 && left_y3 == right_y2) || (left_x3 == right_x3 && left_y3 == right_y3) || (left_x3 == right_x4 && left_y3 == right_y4))
        {
            x[index] = left_x3;
            y[index] = left_y3;
            ID_Vertex[index] = List_Of_Elements[left].getVertex_V3();
            index++;
        }
        if ((left_x4 == right_x1 && left_y4 == right_y1) || (left_x4 == right_x2 && left_y4 == right_y2) || (left_x4 == right_x3 && left_y4 == right_y3) || (left_x4 == right_x4 && left_y4 == right_y4))
        {
             x[index] = left_x4;
             y[index] = left_y4;
             ID_Vertex[index] = List_Of_Elements[left].getVertex_V4();
            // We pass vertices from V4 to V1, but we check V1 first
            if (x[0] == left_x1 && y[0] == left_y1)
            {
                x[index] = x[index-1];
                y[index] = y[index-1];
                x[index-1] = left_x4;
                y[index-1] = left_y4;
                /// Check this reversion
                 ID_Vertex[index] = ID_Vertex[index-1];
                 ID_Vertex[index-1] = List_Of_Elements[left].getVertex_V4();
                /// Check this reversion

            }

            index++;
        }
        //std::cout << "index = " << index-1 << std::endl;
        //std::cout << x[0] << " " << x[1] << std::endl;
        //std::cout << y[0] << " " << y[1] << std::endl;
        //std::cout << ID_Vertex[0] << " " << ID_Vertex[1] << std::endl;

        (*f).set_Vertex_V1(ID_Vertex[0]);
        (*f).set_Vertex_V2(ID_Vertex[1]);

        double dx = x[1]-x[0];
        double dy = y[1]-y[0];
        double Length = sqrt(dx*dx+dy*dy);
        double ReferenceLength = 2.0;
        double Jacobian = Length/ReferenceLength; // Sign Boundary Jacobian. Can this be negative?
        (*f).setJacobian(Jacobian);
        //(*f).setJacobian(Length);
        (*f).set_nx(dy/Length);
        (*f).set_ny(-dx/Length);

        //std::cout << "Boundary ID = " << (*f).getID() << ". Left El = " << left << ", Right El = " << right << std::endl;
        //std::cout << "x[1] = " << x[1] << ", x[0] = " << x[0] << "y[1] = " << y[1] << ", y[0] = " << y[0] << std::endl;
        //std::cout << "J = " << Jacobian << ". nx = " << dy/Length << ", ny = " << -dx/Length << std::endl;
        //std::cout << "Length = " << Length << std::endl;
    }
}
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_boundaries(std::vector<Boundaries2D> &List_Of_Boundaries2D, const std::vector<VertexCoordinates2D> &List_Of_Vertices)
{
    for(auto i = List_Of_Boundaries2D.begin(); i < List_Of_Boundaries2D.end(); i++)
    {
        double x1 = List_Of_Vertices[(*i).getVertex_V1()-1].getxCoordinate();
        double y1 = List_Of_Vertices[(*i).getVertex_V1()-1].getyCoordinate();
        double x2 = List_Of_Vertices[(*i).getVertex_V2()-1].getxCoordinate();
        double y2 = List_Of_Vertices[(*i).getVertex_V2()-1].getyCoordinate();

        double dx = (x2-x1);
        double dy = (y2-y1);
        double Length = sqrt(dx*dx+dy*dy);

        unsigned int Type = (*i).getType();
        double ReferenceLength = 2.0;
        if (Type == 2)
        {
            ReferenceLength = 2.0*sqrt(2.0);
        }

        double Jacobian = Length/ReferenceLength;
        (*i).setJacobian(Jacobian);
        //std::cout << (*i).getID() << " " << Type << " (" << (dy)/Length << ", " << -(dx)/Length << ")" << std::endl;
        (*i).setNormalX(dy/Length);
        (*i).setNormalY(-dx/Length);

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
void set_Order_Polynomials_Uniform(std::vector<Squares2D> &List_Of_Elements2D, const unsigned int &N)
{
    for(auto i = List_Of_Elements2D.begin(); i < List_Of_Elements2D.end(); i++)
    {
        (*i).set_Order_Of_Polynomials(N); //(rand() % 3) + 1
    }
}
/*--------------------------------------------------------------------------*/
extern unsigned int get_Number_Of_Nodes(std::vector<Elements2D> &List_Of_Elements2D)
{
    unsigned int Number_Of_Nodes = 0;
    for(auto i = List_Of_Elements2D.begin(); i < List_Of_Elements2D.end(); i++)
    {
        (*i).set_pos(Number_Of_Nodes);
        Number_Of_Nodes += (*i).get_Number_Of_Nodes();
    }
    return Number_Of_Nodes;
}
/*--------------------------------------------------------------------------*/
extern unsigned int get_Number_Of_Nodes(std::vector<Squares2D> &List_Of_Elements2D)
{
    unsigned int Number_Of_Nodes = 0;
    for(auto i = List_Of_Elements2D.begin(); i < List_Of_Elements2D.end(); i++)
    {
        (*i).set_pos(Number_Of_Nodes);
        Number_Of_Nodes += (*i).get_Number_Of_Nodes();
    }
    return Number_Of_Nodes;
}
/*--------------------------------------------------------------------------*/
void set_theta_Uniform(std::vector<Boundaries2D> &List_Of_Boundaries2D, const double &theta)
{
    for(auto f = List_Of_Boundaries2D.begin(); f < List_Of_Boundaries2D.end(); f++)
    {
        (*f).set_theta(theta); //(rand() % 3) + 1
    }
}
/*--------------------------------------------------------------------------*/
void set_theta_Uniform(std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries2D, const double &theta)
{
    for(auto f = List_Of_Boundaries2D.begin(); f < List_Of_Boundaries2D.end(); f++)
    {
        (*f).set_theta(theta); //(rand() % 3) + 1
    }
}
/*--------------------------------------------------------------------------*/
extern std::string mesh_name_trapezoid(const unsigned int n)
{
    std::string mesh_name;
    switch(n) {
    case 2: mesh_name = "Mesh/trapezoid_2x6.msh"; break;
    case 3: mesh_name = "Mesh/trapezoid_3x8.msh"; break;
    case 4: mesh_name = "Mesh/trapezoid_4x12.msh"; break;
    case 5: mesh_name = "Mesh/trapezoid_5x14.msh"; break;
    case 6: mesh_name = "Mesh/trapezoid_6x16.msh"; break;
    case 7: mesh_name = "Mesh/trapezoid_7x20.msh"; break;
    case 8: mesh_name = "Mesh/trapezoid_8x22.msh"; break;
    case 9: mesh_name = "Mesh/trapezoid_9x26.msh"; break;
    case 10: mesh_name = "Mesh/trapezoid_10x28.msh"; break;
    case 11: mesh_name = "Mesh/trapezoid_11x30.msh"; break;
    case 12: mesh_name = "Mesh/trapezoid_12x34.msh"; break;
    case 15: mesh_name = "Mesh/trapezoid_15x40.msh"; break;
    case 18: mesh_name = "Mesh/trapezoid_18x50.msh"; break;
    case 27: mesh_name = "Mesh/trapezoid_27x72.msh"; break;
    case 33: mesh_name = "Mesh/trapezoid_33x88.msh"; break;
    case 39: mesh_name = "Mesh/trapezoid_39x104.msh"; break;
    case 54: mesh_name = "Mesh/trapezoid_54x140.msh"; break;
    case 105: mesh_name = "Mesh/trapezoid_105x280.msh"; break;
    default: std::cout << "***********************\n* Mesh does not exist *\n* Defaulting to 2x6   *\n***********************"; mesh_name = "Mesh/trapezoid_2x6.msh"; break;
    }

    return mesh_name;
}
/*--------------------------------------------------------------------------*/
