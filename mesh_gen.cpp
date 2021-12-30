#include "mesh_gen.hpp"
/*--------------------------------------------------------------------------*/
void load_msh_mesh(const std::string &mesh_name, Mat &EToV, std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, std::vector<std::unique_ptr<Element>> &List_Of_Elements, int &element_num,  int &node_num)
{
      int *element_node;
      int element_order;
      std::string gmsh_filename = mesh_name;
      int m;
      double *node_x;

    //
    //  Get the data size.
    //
      gmsh_size_read ( gmsh_filename, node_num, m, element_num, element_order );
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
      node_x =  new  double[( m * node_num * sizeof ( double ) )];
      element_node =  new  int[( element_order * element_num * sizeof ( int ) )];
    //
    //  Get the data.
    //
      gmsh_data_read ( gmsh_filename, m, node_num, node_x, element_order, element_num, element_node );
    //
    //  Print some of the data.
    //
    //r8mat_transpose_print_some ( m, node_num, node_x, 1, 1, m, 10, "  Coordinates for first 10 nodes:" );

    unsigned int ID_Vertices = 0;
    unsigned int ID_Elements = 0;
    unsigned int Number_Of_Vertices;
    switch(m)
    {
      case 1:
      {
        Number_Of_Vertices = 2;
        PetscScalar etov[element_num*Number_Of_Vertices];

        for (int i = 0; i < node_num; i++)
          {
            List_Of_Vertices.push_back(std::make_unique<Vertex1D>(ID_Vertices, node_x[m*i]));
            ID_Vertices++;
          }
          for (unsigned int i = 0; i < element_num; i++)
            {
              List_Of_Elements.push_back(std::make_unique<Element1D>(ID_Elements, element_node[Number_Of_Vertices*i]-1, element_node[Number_Of_Vertices*i+1]-1));
              ID_Elements++;
              etov[Number_Of_Vertices*i] = (PetscScalar)element_node[Number_Of_Vertices*i]-1;
              etov[Number_Of_Vertices*i+1] = (PetscScalar)element_node[Number_Of_Vertices*i+1]-1;
            }
            Mat EToVT;
            MatCreateSeqDense(PETSC_COMM_SELF, Number_Of_Vertices, element_num, etov, &EToVT);
            MatAssemblyBegin(EToVT, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(EToVT, MAT_FINAL_ASSEMBLY);
            MatTranspose(EToVT, MAT_INPLACE_MATRIX, &EToVT);
            MatDuplicate(EToVT, MAT_COPY_VALUES, &EToV);
            MatDestroy(&EToVT);
        break;
      }
      case 2:
      {
        Number_Of_Vertices = 4;
        PetscScalar etov[element_num*Number_Of_Vertices];

        for (int i = 0; i < node_num; i++)
          {
            List_Of_Vertices.push_back(std::make_unique<Vertex2D>(ID_Vertices, node_x[m*i], node_x[m*i+1]));
            ID_Vertices++;
          }
          for (unsigned int i = 0; i < element_num; i++)
            {
              List_Of_Elements.push_back(std::make_unique<Element2D>(ID_Elements, element_node[Number_Of_Vertices*i]-1, element_node[Number_Of_Vertices*i+1%Number_Of_Vertices]-1, element_node[Number_Of_Vertices*i+2%Number_Of_Vertices]-1, element_node[Number_Of_Vertices*i+3%Number_Of_Vertices]-1));
              ID_Elements++;
              etov[Number_Of_Vertices*i] = (PetscScalar)element_node[Number_Of_Vertices*i]-1;
              etov[Number_Of_Vertices*i+1] = (PetscScalar)element_node[Number_Of_Vertices*i+1]-1;
              etov[Number_Of_Vertices*i+2] = (PetscScalar)element_node[Number_Of_Vertices*i+2]-1;
              etov[Number_Of_Vertices*i+3] = (PetscScalar)element_node[Number_Of_Vertices*i+3]-1;
            }
            Mat EToVT;
            MatCreateSeqDense(PETSC_COMM_SELF, Number_Of_Vertices, element_num, etov, &EToVT);
            MatAssemblyBegin(EToVT, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(EToVT, MAT_FINAL_ASSEMBLY);
            MatTranspose(EToVT, MAT_INPLACE_MATRIX, &EToVT);
            MatDuplicate(EToVT, MAT_COPY_VALUES, &EToV);
            MatDestroy(&EToVT);
        break;
      }
      case 3:
      {
        Number_Of_Vertices = 8;
        PetscScalar etov[element_num*Number_Of_Vertices];

        for (int i = 0; i < node_num; i++)
          {
            List_Of_Vertices.push_back(std::make_unique<Vertex3D>(ID_Vertices, node_x[m*i], node_x[m*i+1], node_x[m*i+2]));
            ID_Vertices++;
          }
          for (unsigned int i = 0; i < element_num; i++)
            {
              List_Of_Elements.push_back(std::make_unique<Element3D>(ID_Elements, element_node[Number_Of_Vertices*i]-1, element_node[Number_Of_Vertices*i+1]-1, element_node[Number_Of_Vertices*i+2]-1, element_node[Number_Of_Vertices*i+3]-1, element_node[Number_Of_Vertices*i+4]-1, element_node[Number_Of_Vertices*i+5]-1, element_node[Number_Of_Vertices*i+6]-1, element_node[Number_Of_Vertices*i+7]-1));
              ID_Elements++;
              etov[Number_Of_Vertices*i] = (PetscScalar)element_node[Number_Of_Vertices*i]-1;
              etov[Number_Of_Vertices*i+1] = (PetscScalar)element_node[Number_Of_Vertices*i+1]-1;
              etov[Number_Of_Vertices*i+2] = (PetscScalar)element_node[Number_Of_Vertices*i+2]-1;
              etov[Number_Of_Vertices*i+3] = (PetscScalar)element_node[Number_Of_Vertices*i+3]-1;
              etov[Number_Of_Vertices*i+4] = (PetscScalar)element_node[Number_Of_Vertices*i+4]-1;
              etov[Number_Of_Vertices*i+5] = (PetscScalar)element_node[Number_Of_Vertices*i+5]-1;
              etov[Number_Of_Vertices*i+6] = (PetscScalar)element_node[Number_Of_Vertices*i+6]-1;
              etov[Number_Of_Vertices*i+7] = (PetscScalar)element_node[Number_Of_Vertices*i+7]-1;
            }
            Mat EToVT;
            MatCreateSeqDense(PETSC_COMM_SELF, Number_Of_Vertices, element_num, etov, &EToVT);
            MatAssemblyBegin(EToVT, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(EToVT, MAT_FINAL_ASSEMBLY);
            MatTranspose(EToVT, MAT_INPLACE_MATRIX, &EToVT);
            MatDuplicate(EToVT, MAT_COPY_VALUES, &EToV);
            MatDestroy(&EToVT);
        break;
      }
      default:
          std::cout << "Something went wrong in the determination of the vertices" << std::endl;
    }

       /* Gmsh Ordering
        Line:                 Line3:          Line4:

              v
              ^
              |
              |
        0-----+-----1 --> u   0----2----1     0---2---3---1

        Quadrangle:            Quadrangle8:            Quadrangle9:

              v
              ^
              |
        3-----------2          3-----6-----2           3-----6-----2
        |     |     |          |           |           |           |
        |     |     |          |           |           |           |
        |     +---- | --> u    7           5           7     8     5
        |           |          |           |           |           |
        |           |          |           |           |           |
        0-----------1          0-----4-----1           0-----4-----1

        Hexahedron:             Hexahedron20:          Hexahedron27:

               v
        3----------2            3----13----2           3----13----2
        |\     ^   |\           |\         |\          |\         |\
        | \    |   | \          | 15       | 14        |15    24  | 14
        |  \   |   |  \         9  \       11 \        9  \ 20    11 \
        |   7------+---6        |   7----19+---6       |   7----19+---6
        |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
        0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
         \  |    \  \  |         \  17      \  18       \ 17    25 \  18
          \ |     \  \ |         10 |        12|        10 |  21    12|
           \|      w  \|           \|         \|          \|         \|
            4----------5            4----16----5           4----16----5


        * /

*/
    //  i4mat_transpose_print_some ( element_order, element_num, element_node, 1, 1, element_order, 10, "  Connectivity for first 10 elements:" );
    //
    //  Clean up.
    //
      delete[] element_node;
      delete[] node_x;

      // Check the ordering of the vertices. If right-handed -> jacobian is always positive
      //Calculate_Jacobian_RectangularElements(List_Of_Elements, List_Of_Vertices);


}
/*--------------------------------------------------------------------------*/
/*void load_msh_mesh3D(const std::string &mesh_name, Vec &VX, Vec &VY, Vec &VZ, Mat &EToV, std::vector<VertexCoordinates3D> &List_Of_Vertices, std::vector<Cuboid> &List_Of_Elements, int &element_num,  int &node_num)
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
    //r8mat_transpose_print_some ( m, node_num, node_x, 1, 1, m, 10, "  Coordinates for first 10 nodes:" );

    std::cout << "node_num = " << node_num << std::endl;
    unsigned int ID_Vertices = 0;
    unsigned int ID_Elements = 0;
    PetscScalar vx[node_num];
    PetscScalar vy[node_num];
    PetscScalar vz[node_num];
    Vec VXT, VYT, VZT;
    for (int i = 0; i < node_num; i++)
    {
        vx[i] = node_x[3*i];
        vy[i] = node_x[3*i+1];
        vz[i] = node_x[3*i+2];
        std::cout << i <<  " " << vx[i] << " " << vy[i] << " " << vz[i] << std::endl;
        VertexCoordinates3D V(ID_Vertices, node_x[3*i], node_x[3*i+1], node_x[3*i+2]);
        List_Of_Vertices.push_back(V);
        ID_Vertices++;
    }

    VecCreateSeqWithArray(PETSC_COMM_SELF,1,node_num,vx,&VXT);
    VecCreateSeqWithArray(PETSC_COMM_SELF,1,node_num,vy,&VYT);
    VecCreateSeqWithArray(PETSC_COMM_SELF,1,node_num,vz,&VZT);
    VecAssemblyBegin(VXT);
    VecAssemblyEnd(VXT);
    VecAssemblyBegin(VYT);
    VecAssemblyEnd(VYT);
    VecAssemblyBegin(VZT);
    VecAssemblyEnd(VZT);
    VecDuplicate(VXT, &VX);
    VecDuplicate(VYT, &VY);
    VecDuplicate(VZT, &VZ);
    VecCopy(VXT, VX);
    VecCopy(VYT, VY);
    VecCopy(VZT, VZ);
    VecDestroy(&VXT);
    VecDestroy(&VYT);
    VecDestroy(&VZT);

    Mat EToVT;
    PetscScalar etov[element_num*8]; //8 = number of vertices
    for (unsigned int i = 0; i < element_num; i++)
    {
    /*
        // We want all elements to have the same ordering of the vertices//
        double x1 = List_Of_Vertices[element_node[8*i]-1].getxCoordinate();
        double y1 = List_Of_Vertices[element_node[8*i]-1].getyCoordinate();
        double x2 = List_Of_Vertices[element_node[8*i+1]-1].getxCoordinate();
        double y2 = List_Of_Vertices[element_node[8*i+1]-1].getyCoordinate();
        double x3 = List_Of_Vertices[element_node[8*i+2]-1].getxCoordinate();
        double y3 = List_Of_Vertices[element_node[8*i+2]-1].getyCoordinate();
        double x4 = List_Of_Vertices[element_node[8*i+3]-1].getxCoordinate();
        double y4 = List_Of_Vertices[element_node[8*i+3]-1].getyCoordinate();
        double x5 = List_Of_Vertices[element_node[8*i+4]-1].getxCoordinate();
        double y5 = List_Of_Vertices[element_node[8*i+4]-1].getyCoordinate();
        double x6 = List_Of_Vertices[element_node[8*i+5]-1].getxCoordinate();
        double y6 = List_Of_Vertices[element_node[8*i+5]-1].getyCoordinate();
        double x7 = List_Of_Vertices[element_node[8*i+6]-1].getxCoordinate();
        double y7 = List_Of_Vertices[element_node[8*i+6]-1].getyCoordinate();
        double x8 = List_Of_Vertices[element_node[8*i+7]-1].getxCoordinate();
        double y8 = List_Of_Vertices[element_node[8*i+7]-1].getyCoordinate();


        double xorigin = std::min({ x1, x2, x3, x4 });
        double yorigin = std::min({ y1, y2, y3, y4 });
        // Bottom left is origin (0,0).
        double dist_1 = sqrt((x1-xorigin)*(x1-xorigin)+(y1-yorigin)*(y1-yorigin));
        double dist_2 = sqrt((x2-xorigin)*(x2-xorigin)+(y2-yorigin)*(y2-yorigin));
        double dist_3 = sqrt((x3-xorigin)*(x3-xorigin)+(y3-yorigin)*(y3-yorigin));
        double dist_4 = sqrt((x4-xorigin)*(x4-xorigin)+(y4-yorigin)*(y4-yorigin));
        double min = dist_1;* /
        int flag = 0;

        /*
        if(dist_2 <= min)
        {
            min = dist_2;
            flag = 1;
        }
        if(dist_3 <= min)
        {
            min = dist_3;
            flag = 2;
        }
        if(dist_4 <= min)
        {
            min = dist_4;
            flag = 3;
        }*/

        /*
        std::cout << "ID: " << ID_Elements << std::endl;
        std::cout << element_node[4*i]-1 << " " << element_node[4*i+1]-1 << " " << element_node[4*i+2]-1 << " "<< element_node[4*i+3]-1 << std::endl;
        std::cout << "(x0, y0) = " << xorigin << ", " << yorigin << std::endl;
        std::cout << dist_1 << " " << dist_2 << " " << dist_3 << " " << dist_4 << std::endl;
        std::cout << "min = " << min << " => flag = " << flag << std::endl;
        std::cout << std::endl;
        * /



       //std::cout << element_node[4*i+flag]-1 << " " << element_node[4*i+(flag+1)%4]-1 << " " << element_node[4*i+(flag+2)%4]-1 << " " << element_node[4*i+(flag+3)%4]-1 << std::endl;
        Cuboid S(ID_Elements, element_node[8*i+flag]-1, element_node[8*i+(flag+1)%8]-1, element_node[8*i+(flag+2)%8]-1, element_node[8*i+(flag+3)%8]-1, element_node[8*i+(flag+4)%8]-1, element_node[8*i+(flag+5)%8]-1, element_node[8*i+(flag+6)%8]-1, element_node[8*i+(flag+7)%8]-1);
        ID_Elements++;
        List_Of_Elements.push_back(S);

        etov[8*i+1] = (PetscScalar)element_node[8*i+(flag+1)%8]-1;
        etov[8*i+2] = (PetscScalar)element_node[8*i+(flag+2)%8]-1;
        etov[8*i+3] = (PetscScalar)element_node[8*i+(flag+3)%8]-1;
        etov[8*i+4] = (PetscScalar)element_node[8*i+(flag+4)%8]-1;
        etov[8*i+5] = (PetscScalar)element_node[8*i+(flag+5)%8]-1;
        etov[8*i+6] = (PetscScalar)element_node[8*i+(flag+6)%8]-1;
        etov[8*i+7] = (PetscScalar)element_node[8*i+(flag+7)%8]-1;
        etov[8*i] = (PetscScalar)element_node[8*i+flag]-1;
    }

       /* Gmsh Ordering
        Line:                 Line3:          Line4:

              v
              ^
              |
              |
        0-----+-----1 --> u   0----2----1     0---2---3---1

        Quadrangle:            Quadrangle8:            Quadrangle9:

              v
              ^
              |
        3-----------2          3-----6-----2           3-----6-----2
        |     |     |          |           |           |           |
        |     |     |          |           |           |           |
        |     +---- | --> u    7           5           7     8     5
        |           |          |           |           |           |
        |           |          |           |           |           |
        0-----------1          0-----4-----1           0-----4-----1

        Hexahedron:             Hexahedron20:          Hexahedron27:

               v
        3----------2            3----13----2           3----13----2
        |\     ^   |\           |\         |\          |\         |\
        | \    |   | \          | 15       | 14        |15    24  | 14
        |  \   |   |  \         9  \       11 \        9  \ 20    11 \
        |   7------+---6        |   7----19+---6       |   7----19+---6
        |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
        0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
         \  |    \  \  |         \  17      \  18       \ 17    25 \  18
          \ |     \  \ |         10 |        12|        10 |  21    12|
           \|      w  \|           \|         \|          \|         \|
            4----------5            4----16----5           4----16----5


        * /



    MatCreateSeqDense(PETSC_COMM_SELF, 8, element_num, etov, &EToVT);
    MatAssemblyBegin(EToVT, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EToVT, MAT_FINAL_ASSEMBLY);
    MatTranspose(EToVT, MAT_INPLACE_MATRIX, &EToVT);
    MatDuplicate(EToVT, MAT_COPY_VALUES, &EToV);
    MatDestroy(&EToVT);

    //  i4mat_transpose_print_some ( element_order, element_num, element_node, 1, 1, element_order, 10, "  Connectivity for first 10 elements:" );
    //
    //  Clean up.
    //
      delete[] element_node;
      delete[] node_x;


    }

*/
/*--------------------------------------------------------------------------*/
//void load_msh_mesh2D(const std::string &mesh_name, Vec &VX, Vec &VY, Mat &EToV, std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, std::vector<Squares2D> &List_Of_Elements, int &element_num,  int &node_num)
void load_msh_mesh2D(const std::string &mesh_name, Mat &EToV, std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, std::vector<Squares2D> &List_Of_Elements, int &element_num,  int &node_num)
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
    //PetscScalar vx[node_num];
    //PetscScalar vy[node_num];
    Vec VXT, VYT;
    for (int i = 0; i < node_num; i++)
    {
        //vx[i] = node_x[2*i];
        //vy[i] = node_x[2*i+1];
        List_Of_Vertices.push_back(std::make_unique<Vertex2D>(ID_Vertices, node_x[2*i], node_x[2*i+1]));
        ID_Vertices++;
    }
    //VecCreateSeqWithArray(PETSC_COMM_SELF,1,node_num,vx,&VXT);
    //VecCreateSeqWithArray(PETSC_COMM_SELF,1,node_num,vy,&VYT);
    //VecAssemblyBegin(VXT);
    //VecAssemblyEnd(VXT);
    //VecAssemblyBegin(VYT);
    //VecAssemblyEnd(VYT);
    //VecDuplicate(VXT, &VX);
    //VecDuplicate(VYT, &VY);
    //VecCopy(VXT, VX);
    //VecCopy(VYT, VY);
    //VecDestroy(&VXT);
    //VecDestroy(&VYT);

    Mat EToVT;
    PetscScalar etov[element_num*4]; //4 = element_order
    for (unsigned int i = 0; i < element_num; i++)
    {
        //std::cout << element_node[4*i] << " " << element_node[4*i+1] << " " << element_node[4*i+2] << " "<< element_node[4*i+3] << " " << std::endl;

        double x1 = List_Of_Vertices[element_node[4*i]-1]->getxCoordinate();
        double y1 = List_Of_Vertices[element_node[4*i]-1]->getyCoordinate();
        double x2 = List_Of_Vertices[element_node[4*i+1]-1]->getxCoordinate();
        double y2 = List_Of_Vertices[element_node[4*i+1]-1]->getyCoordinate();
        double x3 = List_Of_Vertices[element_node[4*i+2]-1]->getxCoordinate();
        double y3 = List_Of_Vertices[element_node[4*i+2]-1]->getyCoordinate();
        double x4 = List_Of_Vertices[element_node[4*i+3]-1]->getxCoordinate();
        double y4 = List_Of_Vertices[element_node[4*i+3]-1]->getyCoordinate();


        double xorigin = std::min({ x1, x2, x3, x4 });
        double yorigin = std::min({ y1, y2, y3, y4 });
        // Bottom left is origin (0,0).
        double dist_1 = sqrt((x1-xorigin)*(x1-xorigin)+(y1-yorigin)*(y1-yorigin));
        double dist_2 = sqrt((x2-xorigin)*(x2-xorigin)+(y2-yorigin)*(y2-yorigin));
        double dist_3 = sqrt((x3-xorigin)*(x3-xorigin)+(y3-yorigin)*(y3-yorigin));
        double dist_4 = sqrt((x4-xorigin)*(x4-xorigin)+(y4-yorigin)*(y4-yorigin));
        double min = dist_1;
        int flag = 0;

        if(dist_2 <= min)
        {
            min = dist_2;
            flag = 1;
        }
        if(dist_3 <= min)
        {
            min = dist_3;
            flag = 2;
        }
        if(dist_4 <= min)
        {
            min = dist_4;
            flag = 3;
        }

        /*
        std::cout << "ID: " << ID_Elements << std::endl;
        std::cout << element_node[4*i]-1 << " " << element_node[4*i+1]-1 << " " << element_node[4*i+2]-1 << " "<< element_node[4*i+3]-1 << std::endl;
        std::cout << "(x0, y0) = " << xorigin << ", " << yorigin << std::endl;
        std::cout << dist_1 << " " << dist_2 << " " << dist_3 << " " << dist_4 << std::endl;
        std::cout << "min = " << min << " => flag = " << flag << std::endl;
        std::cout << std::endl;
        */



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
void Connect(const Mat &EToV, const unsigned int &Number_Of_Elements, const unsigned int &Number_Of_Vertices, Mat &EToE, Mat &EToF, std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries)
{

    unsigned int Nfaces;
    unsigned int Total_Faces;
    unsigned int m =1;
    switch(m)
    {
          case 1:
          {
            Nfaces = 2;
            Total_Faces = Nfaces*Number_Of_Elements;

            break;
          }
          case 2:
          {
            Nfaces = 4;
            Total_Faces = Nfaces*Number_Of_Elements;

            break;
          }
          case 3:
          {
            Nfaces = 6;
            Total_Faces = Nfaces*Number_Of_Elements;

            break;
          }
          default:
              std::cout << "Something went wrong in the connection of the boundaries" << std::endl;
        }
}
/*--------------------------------------------------------------------------*/
void Connect_3D(const Mat &EToV, const unsigned int &Number_Of_Elements, const unsigned int &Number_Of_Vertices, Mat &EToE, Mat &EToF, std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries)
{
  unsigned int Nfaces = 6;
  unsigned int Total_Faces = Nfaces*Number_Of_Elements;
  MatCreateSeqDense(PETSC_COMM_WORLD, Number_Of_Elements, Nfaces, NULL, &EToE);
  MatCreateSeqDense(PETSC_COMM_WORLD, Number_Of_Elements, Nfaces, NULL, &EToF);

  PetscInt ir[1]={0};
   // Store 1 .. 6 (=Nfaces) in each row of EtoF
  PetscInt array_if[Nfaces];
  PetscScalar iv_f[Nfaces];
  for (PetscInt k=0;k<Nfaces; k++)
  {
      array_if[k] = k;
      iv_f[k] = k;
  }
  for (PetscInt i=0; i<Number_Of_Elements; i++)
  {
      ir[0] = i;
      MatSetValues(EToF, 1, ir, Nfaces, array_if, iv_f, INSERT_VALUES);
  }

   // Store 1 .. K in each column of EtoE
  PetscInt array_ie[Number_Of_Elements];
  PetscScalar iv[Number_Of_Elements];
  for (PetscInt k=0;k<Number_Of_Elements; k++)
  {
      array_ie[k] = k;
      iv[k] = k;
  }
  for (PetscInt i=0; i<Nfaces; i++)
  {
      ir[0] = i;
      MatSetValues(EToE, Number_Of_Elements, array_ie, 1, ir, iv, INSERT_VALUES);
  }

  MatAssemblyBegin(EToF, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(EToF,   MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(EToE, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(EToE,   MAT_FINAL_ASSEMBLY);

  PetscInt sk = 0;
  PetscInt ic[4]={0};
  PetscInt i_insert[4]={0, 1, 2, 3};
  // uniquely number each set of four faces by their node numbers
  std::vector<std::vector<unsigned int>> spNodetoNode;

  for (unsigned int f=0; f<Nfaces; f++)
  {
      if (f == 0)
      {
          ic[0] = 0;
          ic[1] = 1;
          ic[2] = 2;
          ic[3] = 3;
      }
      else if (f == 1)
      {
          ic[0] = 4;
          ic[1] = 5;
          ic[2] = 6;
          ic[3] = 7;
      }
      else if (f == 2)
      {
          ic[0] = 0;
          ic[1] = 4;
          ic[2] = 5;
          ic[3] = 1;
      }
      else if (f == 3)
      {
          ic[0] = 3;
          ic[1] = 7;
          ic[2] = 6;
          ic[3] = 2;
      }
      else if (f == 4)
      {
          //ic[0] = 0;
          //ic[1] = 4;
          //ic[2] = 7;
          //ic[3] = 3;
          ic[0] = 0;
          ic[1] = 3;
          ic[2] = 7;
          ic[3] = 4;
      }
      else if (f == 5)
      {
          //ic[0] = 1;
          //ic[1] = 5;
          //ic[2] = 6;
          //ic[3] = 2;
          ic[0] = 1;
          ic[1] = 2;
          ic[2] = 6;
          ic[3] = 5;
      }
      for (unsigned int k=0; k<Number_Of_Elements; k++)
      {
          //ir[0] = sk;
          const PetscScalar *vals;
          MatGetRow(EToV, k, NULL, NULL, &vals);
          // We do not subtract one
          PetscScalar v[4] = {vals[ic[0]], vals[ic[1]], vals[ic[2]], vals[ic[3]]};
          std::sort(std::begin(v), std::end(v) );
          unsigned int id = v[0]*Number_Of_Vertices*Number_Of_Vertices*Number_Of_Vertices+v[1]*Number_Of_Vertices*Number_Of_Vertices+v[2]*Number_Of_Vertices+v[3];
          //std::cout << id << std::endl;
          MatRestoreRow(EToV, k, NULL, NULL, &vals);
          std::vector<unsigned int> Row {id, sk, array_ie[sk%Number_Of_Elements], array_if[sk/Number_Of_Elements]};
          spNodetoNode.push_back(Row);
          // [ id sk EtoE_a[sk] EtoF_a[sk]]
          sk++;
      }
  }
  sortrows(spNodetoNode,0);

  /*
  std::cout << "spNodetoNode sorted = " << std::endl;
  for ( const std::vector<unsigned int> &v : spNodetoNode )
  {
      for (unsigned int x : v )
      {
          std::cout << x << " " ;
      }
      std::cout << std::endl;
  }*/

  std::vector<unsigned int> indices;
  //std::cout << "indices = " << std::endl;
  int Xold = -1;
  unsigned int Xnew, k = 0;
  for ( const std::vector<unsigned int> &v : spNodetoNode )
  {
     Xnew = v[0];
     if (k > 0)
     {
          if (Xold == Xnew)
          {
              indices.push_back(k-1);
          }
          Xold = Xnew;
     }
     k++;
  }
  //for (unsigned int x : indices ) std::cout << x << ' ';
  //std::cout << std::endl;

  std::vector<unsigned int> matchL2;
  std::vector<unsigned int> matchR3;
  std::vector<unsigned int> matchR4;
  for (unsigned int kk = 0 ; kk < indices.size(); kk++)
  {
      matchL2.push_back(spNodetoNode[indices[kk]][1]);
      matchR3.push_back(spNodetoNode[indices[kk]+1][2]);
      matchR4.push_back(spNodetoNode[indices[kk]+1][3]);
  }
  for (unsigned int kk = 0 ; kk < indices.size(); kk++)
  {

      matchL2.push_back(spNodetoNode[indices[kk]+1][1]);
      matchR3.push_back(spNodetoNode[indices[kk]][2]);
      matchR4.push_back(spNodetoNode[indices[kk]][3]);
  }
  /*
  std::cout << "matchL2 = " << std::endl;
  for(const auto& value: matchL2) {
  std::cout << value << "\n";
  }
  std::cout << std::endl;
  std::cout << "matchR3 = " << std::endl;
  for(const auto& value: matchR3) {
  std::cout << value << "\n";
  }
  std::cout << std::endl;
  std::cout << "matchR4 = " << std::endl;
  for(const auto& value: matchR4) {
  std::cout << value << "\n";
  }
  std::cout << std::endl;
  */
  for (unsigned int kk = 0 ; kk < matchL2.size(); kk++)
  {
      PetscInt row, col;
      row = ((int)matchL2[kk])%((int)Number_Of_Elements);
      col = floor(((int)matchL2[kk])/(int)Number_Of_Elements);
      MatSetValue(EToE, row, col, matchR3[kk], INSERT_VALUES);
      MatSetValue(EToF, row, col, matchR4[kk], INSERT_VALUES);
  }
  MatAssemblyBegin(EToF, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(EToF,   MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(EToE, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(EToE,   MAT_FINAL_ASSEMBLY);

  create_ListInternalBoundaries(Number_Of_Elements, EToE, EToF, List_Of_Boundaries);


}
/*--------------------------------------------------------------------------*/
void create_ListInternalBoundaries(const unsigned int &Number_Of_Elements, Mat &EToE, Mat &EToF, std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries)
{
    unsigned int Nfaces = 6;
    unsigned int TotalFaces = Nfaces * Number_Of_Elements;

    std::set<BoundaryInfo, decltype(compare_BoundaryInfo)> SetType(compare_BoundaryInfo);
    for (unsigned int i = 0; i < Number_Of_Elements; i++)
    {
        const PetscScalar *EToE_row, *EToF_row;
        MatGetRow(EToE, i, NULL, NULL, &EToE_row);
        MatGetRow(EToF, i, NULL, NULL, &EToF_row);
        for (unsigned int j = 0; j < Nfaces; j++)
        {
            if (i==EToE_row[j])
            {
                // External Face
            }
            else
            {
                // Internal Face
                BoundaryInfo BI(i, EToE_row[j], j, EToF_row[j]);
                SetType.emplace(BI);
            }
        }
        MatRestoreRow(EToE, i, NULL, NULL, &EToE_row);
        MatRestoreRow(EToF, i, NULL, NULL, &EToF_row);
    }
    unsigned int ID_Boundary = 0;
    for(const auto& it: SetType)
    {
        List_Of_Boundaries.push_back(std::make_unique<Boundary3D>(ID_Boundary, it.left, it.right, it.faceleft, it.faceright));
        ID_Boundary++;
    }

}
/*--------------------------------------------------------------------------*/


/*void Connect3D(const Mat &EToV, const unsigned int &Number_Of_Elements, const unsigned int &Number_Of_Vertices, Mat &EToE, Mat &EToF, std::vector<InternalBoundariesCuboid> &List_Of_Boundaries)
{
    unsigned int Nfaces = 6;
    unsigned int Total_Faces = Nfaces*Number_Of_Elements;
    MatCreateSeqDense(PETSC_COMM_WORLD, Number_Of_Elements, Nfaces, NULL, &EToE);
    MatCreateSeqDense(PETSC_COMM_WORLD, Number_Of_Elements, Nfaces, NULL, &EToF);

    PetscInt ir[1]={0};
     // Store 1 .. 6 (=Nfaces) in each row of EtoF
    PetscInt array_if[Nfaces];
    PetscScalar iv_f[Nfaces];
    for (PetscInt k=0;k<Nfaces; k++)
    {
        array_if[k] = k;
        iv_f[k] = k;
    }
    for (PetscInt i=0; i<Number_Of_Elements; i++)
    {
        ir[0] = i;
        MatSetValues(EToF, 1, ir, Nfaces, array_if, iv_f, INSERT_VALUES);
    }

     // Store 1 .. K in each column of EtoE
    PetscInt array_ie[Number_Of_Elements];
    PetscScalar iv[Number_Of_Elements];
    for (PetscInt k=0;k<Number_Of_Elements; k++)
    {
        array_ie[k] = k;
        iv[k] = k;
    }
    for (PetscInt i=0; i<Nfaces; i++)
    {
        ir[0] = i;
        MatSetValues(EToE, Number_Of_Elements, array_ie, 1, ir, iv, INSERT_VALUES);
    }

    MatAssemblyBegin(EToF, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EToF,   MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(EToE, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EToE,   MAT_FINAL_ASSEMBLY);

    PetscInt sk = 0;
    PetscInt ic[4]={0};
    PetscInt i_insert[4]={0, 1, 2, 3};
    std::vector<std::vector<unsigned int>> spNodetoNode;

    for (unsigned int f=0; f<Nfaces; f++)
    {
        if (f == 0)
        {
            ic[0] = 0;
            ic[1] = 1;
            ic[2] = 2;
            ic[3] = 3;
        }
        else if (f == 1)
        {
            ic[0] = 4;
            ic[1] = 5;
            ic[2] = 6;
            ic[3] = 7;
        }
        else if (f == 2)
        {
            ic[0] = 0;
            ic[1] = 4;
            ic[2] = 5;
            ic[3] = 1;
        }
        else if (f == 3)
        {
            ic[0] = 3;
            ic[1] = 7;
            ic[2] = 6;
            ic[3] = 2;
        }
        else if (f == 4)
        {
            //ic[0] = 0;
            //ic[1] = 4;
            //ic[2] = 7;
            //ic[3] = 3;
            ic[0] = 0;
            ic[1] = 3;
            ic[2] = 7;
            ic[3] = 4;
        }
        else if (f == 5)
        {
            //ic[0] = 1;
            //ic[1] = 5;
            //ic[2] = 6;
            //ic[3] = 2;
            ic[0] = 1;
            ic[1] = 2;
            ic[2] = 6;
            ic[3] = 5;
        }
        for (unsigned int k=0; k<Number_Of_Elements; k++)
        {
            //ir[0] = sk;
            const PetscScalar *vals;
            MatGetRow(EToV, k, NULL, NULL, &vals);
            // We do not subtract one
            PetscScalar v[4] = {vals[ic[0]], vals[ic[1]], vals[ic[2]], vals[ic[3]]};
            std::sort(std::begin(v), std::end(v) );
            unsigned int id = v[0]*Number_Of_Vertices*Number_Of_Vertices*Number_Of_Vertices+v[1]*Number_Of_Vertices*Number_Of_Vertices+v[2]*Number_Of_Vertices+v[3];
            //std::cout << id << std::endl;
            MatRestoreRow(EToV, k, NULL, NULL, &vals);
            std::vector<unsigned int> Row {id, sk, array_ie[sk%Number_Of_Elements], array_if[sk/Number_Of_Elements]};
            spNodetoNode.push_back(Row);
            // [ id sk EtoE_a[sk] EtoF_a[sk]]
            sk++;
        }
    }
    sortrows(spNodetoNode,0);
    /*
    std::cout << "spNodetoNode sorted = " << std::endl;
    for ( const std::vector<unsigned int> &v : spNodetoNode )
    {
        for (unsigned int x : v )
        {
            std::cout << x << " " ;
        }
        std::cout << std::endl;
    }
    * /
    std::vector<unsigned int> indices;
    //std::cout << "indices = " << std::endl;
    int Xold = -1;
    unsigned int Xnew, k = 0;
    for ( const std::vector<unsigned int> &v : spNodetoNode )
    {
       Xnew = v[0];
       if (k > 0)
       {
            if (Xold == Xnew)
            {
                indices.push_back(k-1);
            }
            Xold = Xnew;
       }
       k++;
    }
    //for (unsigned int x : indices ) std::cout << x << ' ';
    //std::cout << std::endl;

    std::vector<unsigned int> matchL2;
    std::vector<unsigned int> matchR3;
    std::vector<unsigned int> matchR4;
    for (unsigned int kk = 0 ; kk < indices.size(); kk++)
    {
        matchL2.push_back(spNodetoNode[indices[kk]][1]);
        matchR3.push_back(spNodetoNode[indices[kk]+1][2]);
        matchR4.push_back(spNodetoNode[indices[kk]+1][3]);
    }
    for (unsigned int kk = 0 ; kk < indices.size(); kk++)
    {

        matchL2.push_back(spNodetoNode[indices[kk]+1][1]);
        matchR3.push_back(spNodetoNode[indices[kk]][2]);
        matchR4.push_back(spNodetoNode[indices[kk]][3]);
    }
    /*
    std::cout << "matchL2 = " << std::endl;
    for(const auto& value: matchL2) {
    std::cout << value << "\n";
    }
    std::cout << std::endl;
    std::cout << "matchR3 = " << std::endl;
    for(const auto& value: matchR3) {
    std::cout << value << "\n";
    }
    std::cout << std::endl;
    std::cout << "matchR4 = " << std::endl;
    for(const auto& value: matchR4) {
    std::cout << value << "\n";
    }
    std::cout << std::endl;
    * /
    for (unsigned int kk = 0 ; kk < matchL2.size(); kk++)
    {
        PetscInt row, col;
        row = ((int)matchL2[kk])%((int)Number_Of_Elements);
        col = floor(((int)matchL2[kk])/(int)Number_Of_Elements);
        MatSetValue(EToE, row, col, matchR3[kk], INSERT_VALUES);
        MatSetValue(EToF, row, col, matchR4[kk], INSERT_VALUES);
    }
    MatAssemblyBegin(EToF, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EToF,   MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(EToE, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(EToE,   MAT_FINAL_ASSEMBLY);

    create_ListInternalBoundaries(Number_Of_Elements, EToE, EToF, List_Of_Boundaries);

}
/*--------------------------------------------------------------------------*/
/*void create_ListInternalBoundaries(const unsigned int &Number_Of_Elements, Mat &EToE, Mat &EToF, std::vector<InternalBoundariesCuboid> &List_Of_Boundaries)
{
    unsigned int Nfaces = 6;
    unsigned int TotalFaces = Nfaces * Number_Of_Elements;

    std::set<BoundaryInfo, decltype(compare_BoundaryInfo)> SetType(compare_BoundaryInfo);
    for (unsigned int i = 0; i < Number_Of_Elements; i++)
    {
        const PetscScalar *EToE_row, *EToF_row;
        MatGetRow(EToE, i, NULL, NULL, &EToE_row);
        MatGetRow(EToF, i, NULL, NULL, &EToF_row);
        for (unsigned int j = 0; j < Nfaces; j++)
        {
            if (i==EToE_row[j])
            {
                // External Face
            }
            else
            {
                // Internal Face
                BoundaryInfo BI(i, EToE_row[j], j, EToF_row[j]);
                SetType.emplace(BI);
            }
        }
        MatRestoreRow(EToE, i, NULL, NULL, &EToE_row);
        MatRestoreRow(EToF, i, NULL, NULL, &EToF_row);
    }
    unsigned int ID_Boundary = 0;
    //std::cout << "left right left-face right-face" <<  std::endl;
    for(const auto& it: SetType)
    {
        InternalBoundariesCuboid B(ID_Boundary, it.left, it.right, it.faceleft, it.faceright);
        ID_Boundary++;
        List_Of_Boundaries.push_back(B);
        //std::cout << (it).left << " " << (it).right << " "<< (it).faceleft << " "<< (it).faceright <<  std::endl;
    }

} */
/*--------------------------------------------------------------------------*/
void Connect2D(const Mat &EToV, const unsigned int &Number_Of_Elements, const unsigned int &Number_Of_Vertices, Mat &EToE, Mat &EToF, std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries)
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
        List_Of_Boundaries.push_back(std::make_unique<Boundary2D>(ID_Boundary, it.left, it.right, it.faceleft, it.faceright));
        //InternalBoundariesSquares2D B(ID_Boundary, it.left, it.right, it.faceleft, it.faceright);
        ID_Boundary++;
        //List_Of_Boundaries.push_back(B);


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
void Calculate_Jacobian_RectangularElements(std::vector<std::unique_ptr<Element>> &List_Of_Elements, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices)
{

    // only for 3D Cuboids

    // Assumes Rectangular Cuboids => J is constant
    for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
    {
        double x0 = List_Of_Vertices[(*i)->getVertex_V1()]->getxCoordinate();
        double y0 = List_Of_Vertices[(*i)->getVertex_V1()]->getyCoordinate();
        double z0 = List_Of_Vertices[(*i)->getVertex_V1()]->getzCoordinate();
        double x1 = List_Of_Vertices[(*i)->getVertex_V2()]->getxCoordinate();
        double y1 = List_Of_Vertices[(*i)->getVertex_V2()]->getyCoordinate();
        double z1 = List_Of_Vertices[(*i)->getVertex_V2()]->getzCoordinate();
        double x2 = List_Of_Vertices[(*i)->getVertex_V3()]->getxCoordinate();
        double y2 = List_Of_Vertices[(*i)->getVertex_V3()]->getyCoordinate();
        double z2 = List_Of_Vertices[(*i)->getVertex_V3()]->getzCoordinate();
        double x3 = List_Of_Vertices[(*i)->getVertex_V4()]->getxCoordinate();
        double y3 = List_Of_Vertices[(*i)->getVertex_V4()]->getyCoordinate();
        double z3 = List_Of_Vertices[(*i)->getVertex_V4()]->getzCoordinate();
        double x4 = List_Of_Vertices[(*i)->getVertex_V5()]->getxCoordinate();
        double y4 = List_Of_Vertices[(*i)->getVertex_V5()]->getyCoordinate();
        double z4 = List_Of_Vertices[(*i)->getVertex_V5()]->getzCoordinate();
        double x5 = List_Of_Vertices[(*i)->getVertex_V6()]->getxCoordinate();
        double y5 = List_Of_Vertices[(*i)->getVertex_V6()]->getyCoordinate();
        double z5 = List_Of_Vertices[(*i)->getVertex_V6()]->getzCoordinate();
        double x6 = List_Of_Vertices[(*i)->getVertex_V7()]->getxCoordinate();
        double y6 = List_Of_Vertices[(*i)->getVertex_V7()]->getyCoordinate();
        double z6 = List_Of_Vertices[(*i)->getVertex_V7()]->getzCoordinate();
        double x7 = List_Of_Vertices[(*i)->getVertex_V8()]->getxCoordinate();
        double y7 = List_Of_Vertices[(*i)->getVertex_V8()]->getyCoordinate();
        double z7 = List_Of_Vertices[(*i)->getVertex_V8()]->getzCoordinate();


        double dxdr = (-x0+x1+x2-x3-x4+x5+x6-x7)/8.0;
        double dxds = (-x0-x1+x2+x3-x4-x5+x6+x7)/8.0;
        double dxdt = (-x0-x1-x2-x3+x4+x5+x6+x7)/8.0;
        double dydr = (-y0+y1+y2-y3-y4+y5+y6-y7)/8.0;
        double dyds = (-y0-y1+y2+y3-y4-y5+y6+y7)/8.0;
        double dydt = (-y0-y1-y2-y3+y4+y5+y6+y7)/8.0;
        double dzdr = (-z0+z1+z2-z3-z4+z5+z6-z7)/8.0;
        double dzds = (-z0-z1+z2+z3-z4-z5+z6+z7)/8.0;
        double dzdt = (-z0-z1-z2-z3+z4+z5+z6+z7)/8.0;


        double det_Jacobian = dzdt*(dxdr*dyds-dxds*dydr)-dzds*(dxdt*dydr-dxdr*dydt)+dzdr*(dxds*dydt-dxdt*dyds);

        std::cout << "Element " << (*i)->getID() << ", det J = " << det_Jacobian << std::endl;
    }
}
/*--------------------------------------------------------------------------*/
void Calculate_Area_Square(std::vector<Squares2D> &List_Of_Elements, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices)
{
    // Assumes Square Quadrilaterals => J is constant
    for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
    {
        double x1 = List_Of_Vertices[(*i).getVertex_V1()]->getxCoordinate();
        double y1 = List_Of_Vertices[(*i).getVertex_V1()]->getyCoordinate();
        double x2 = List_Of_Vertices[(*i).getVertex_V2()]->getxCoordinate();
        double y2 = List_Of_Vertices[(*i).getVertex_V2()]->getyCoordinate();
        double x3 = List_Of_Vertices[(*i).getVertex_V3()]->getxCoordinate();
        double y3 = List_Of_Vertices[(*i).getVertex_V3()]->getyCoordinate();
        double x4 = List_Of_Vertices[(*i).getVertex_V4()]->getxCoordinate();
        double y4 = List_Of_Vertices[(*i).getVertex_V4()]->getyCoordinate();

        double Area = 0.5*abs(x1*y2+x2*y3+x3*y4-x2*y1-x3*y2-x4*y3-x1*y4);
        //std::cout << "Area = " << Area << std::endl;

        (*i).setArea(Area);
    }
}
/*--------------------------------------------------------------------------*/
void Calculate_CuboidFaceNormals(const std::vector<std::unique_ptr<Element>> &List_Of_Elements, std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices)
{
    for(auto f = List_Of_Boundaries.begin(); f < List_Of_Boundaries.end(); f++)
    {
        int left = (*f)->getLeftElementID();
        int left_face = (*f)->getTypeLeft();
        int right = (*f)->getRightElementID();

        double x1, x2, x3, y1, y2, y3, z1, z2, z3;
        x1 = x2 = x3 = y1 = y2 = y3 = z1 = z2 = z3 = 0.0;
        // We read the vertices counterclockwise so that we obtain the outward normal
        switch(left_face)
        {
            case 0:
                x1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V1()]->getxCoordinate();
                y1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V1()]->getyCoordinate();
                z1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V1()]->getzCoordinate();
                x2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V2()]->getxCoordinate();
                y2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V2()]->getyCoordinate();
                z2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V2()]->getzCoordinate();
                x3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V3()]->getxCoordinate();
                y3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V3()]->getyCoordinate();
                z3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V3()]->getzCoordinate();
                break;
            case 1:
                x1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V5()]->getxCoordinate();
                y1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V5()]->getyCoordinate();
                z1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V5()]->getzCoordinate();
                x2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V6()]->getxCoordinate();
                y2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V6()]->getyCoordinate();
                z2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V6()]->getzCoordinate();
                x3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V7()]->getxCoordinate();
                y3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V7()]->getyCoordinate();
                z3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V7()]->getzCoordinate();
                break;
            case 2:
                x1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V1()]->getxCoordinate();
                y1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V1()]->getyCoordinate();
                z1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V1()]->getzCoordinate();
                x2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V5()]->getxCoordinate();
                y2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V5()]->getyCoordinate();
                z2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V5()]->getzCoordinate();
                x3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V6()]->getxCoordinate();
                y3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V6()]->getyCoordinate();
                z3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V6()]->getzCoordinate();
                break;
            case 3:
                x1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V4()]->getxCoordinate();
                y1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V4()]->getyCoordinate();
                z1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V4()]->getzCoordinate();
                x2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V8()]->getxCoordinate();
                y2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V8()]->getyCoordinate();
                z2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V8()]->getzCoordinate();
                x3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V7()]->getxCoordinate();
                y3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V7()]->getyCoordinate();
                z3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V7()]->getzCoordinate();
                break;
            case 4:
                x1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V1()]->getxCoordinate();
                y1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V1()]->getyCoordinate();
                z1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V1()]->getzCoordinate();
                x2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V4()]->getxCoordinate();
                y2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V4()]->getyCoordinate();
                z2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V4()]->getzCoordinate();
                x3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V8()]->getxCoordinate();
                y3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V8()]->getyCoordinate();
                z3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V8()]->getzCoordinate();
                break;
            case 5:
                x1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V2()]->getxCoordinate();
                y1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V2()]->getyCoordinate();
                z1 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V2()]->getzCoordinate();
                x2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V3()]->getxCoordinate();
                y2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V3()]->getyCoordinate();
                z2 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V3()]->getzCoordinate();
                x3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V7()]->getxCoordinate();
                y3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V7()]->getyCoordinate();
                z3 = List_Of_Vertices[List_Of_Elements[left]->getVertex_V7()]->getzCoordinate();
                break;
            default:
                std::cout << "Something went wrong in the calculation of the normals" << std::endl;
        }
        double v1x = x2 - x1;
        double v1y = y2 - y1;
        double v1z = z2 - z1;
        double v2x = x3 - x2;
        double v2y = y3 - y2;
        double v2z = z3 - z2;

        double nx = v1y*v2z-v1z*v2y;
        double ny = v1z*v2x-v1x*v2z;
        double nz = v1x*v2y-v1y*v2x;

        // Normalize
        double length_n = sqrt(nx*nx+ny*ny+nz*nz);
        nx = nx/length_n;
        ny = ny/length_n;
        nz = nz/length_n;

        //std::cout << "length_n = " << length_n << std::endl;
        //std::cout << "Boundary ID = " << (*f).getID() << ". Left El = " << left << ", Right El = " << right << ", left face = " << left_face << std::endl;
        //std::cout << nx << " " << ny << " " << nz << std::endl;

        double Area = 2.0*2.0;
        double detJacobian = length_n/Area;
        //std::cout << "detJacobian = " << detJacobian << std::endl;
        (*f)->setJacobian(detJacobian);
        (*f)->set_nx(nx);
        (*f)->set_ny(ny);
        (*f)->set_nz(nz);
    }
}
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Boundaries_Square(const std::vector<Squares2D> &List_Of_Elements, std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices)
{
    //std::cout << "Boundary Jacobian " << std::endl;
    for(auto f = List_Of_Boundaries.begin(); f < List_Of_Boundaries.end(); f++)
    {
        //std::cout << std::endl;
        int left = (*f)->getLeftElementID();
        int right = (*f)->getRightElementID();

        double left_x1 = List_Of_Vertices[List_Of_Elements[left].getVertex_V1()]->getxCoordinate();
        double left_y1 = List_Of_Vertices[List_Of_Elements[left].getVertex_V1()]->getyCoordinate();
        double left_x2 = List_Of_Vertices[List_Of_Elements[left].getVertex_V2()]->getxCoordinate();
        double left_y2 = List_Of_Vertices[List_Of_Elements[left].getVertex_V2()]->getyCoordinate();
        double left_x3 = List_Of_Vertices[List_Of_Elements[left].getVertex_V3()]->getxCoordinate();
        double left_y3 = List_Of_Vertices[List_Of_Elements[left].getVertex_V3()]->getyCoordinate();
        double left_x4 = List_Of_Vertices[List_Of_Elements[left].getVertex_V4()]->getxCoordinate();
        double left_y4 = List_Of_Vertices[List_Of_Elements[left].getVertex_V4()]->getyCoordinate();

        double right_x1 = List_Of_Vertices[List_Of_Elements[right].getVertex_V1()]->getxCoordinate();
        double right_y1 = List_Of_Vertices[List_Of_Elements[right].getVertex_V1()]->getyCoordinate();
        double right_x2 = List_Of_Vertices[List_Of_Elements[right].getVertex_V2()]->getxCoordinate();
        double right_y2 = List_Of_Vertices[List_Of_Elements[right].getVertex_V2()]->getyCoordinate();
        double right_x3 = List_Of_Vertices[List_Of_Elements[right].getVertex_V3()]->getxCoordinate();
        double right_y3 = List_Of_Vertices[List_Of_Elements[right].getVertex_V3()]->getyCoordinate();
        double right_x4 = List_Of_Vertices[List_Of_Elements[right].getVertex_V4()]->getxCoordinate();
        double right_y4 = List_Of_Vertices[List_Of_Elements[right].getVertex_V4()]->getyCoordinate();
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
                //std::cout << " yes" << std::endl;
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

        //(*f).set_Vertex_V1(ID_Vertex[0]);
        //(*f).set_Vertex_V2(ID_Vertex[1]);

        double dx = x[1]-x[0];
        double dy = y[1]-y[0];
        double Length = sqrt(dx*dx+dy*dy);
        double ReferenceLength = 2.0;
        double Jacobian = Length/ReferenceLength; // Sign Boundary Jacobian. Can this be negative? => shouldn't be for convex elements
        (*f)->setJacobian(Jacobian);
        //(*f).setJacobian(Length);
        (*f)->set_nx(dy/Length);
        (*f)->set_ny(-dx/Length);

        /*
        std::cout << "Boundary ID = " << (*f).getID() << ". Left El = " << left << ", Right El = " << right << std::endl;
        std::cout << "x[1] = " << x[1] << ", x[0] = " << x[0] << "y[1] = " << y[1] << ", y[0] = " << y[0] << std::endl;
        std::cout << "J = " << Jacobian << ". nx = " << dy/Length << ", ny = " << -dx/Length << std::endl;
        std::cout << "Length = " << Length << std::endl;
        */

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
void set_Order_Polynomials_Uniform(std::vector<std::unique_ptr<Element>> &List_Of_Elements, const unsigned int &Nx, const unsigned int &Ny, const unsigned int &Nz)
{
    for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
    {
        (*i)->set_Order_Of_Polynomials_x(Nx); //(rand() % 3) + 1
        (*i)->set_Order_Of_Polynomials_y(Ny); //(rand() % 3) + 1
        (*i)->set_Order_Of_Polynomials_z(Nz); //(rand() % 3) + 1
        unsigned int Nnodes = (Nx+1)*(Ny+1)*(Nz+1);
        (*i)->set_Number_Of_Nodes(Nnodes);
    }
}
/*--------------------------------------------------------------------------*/
extern unsigned int get_Number_Of_Nodes(std::vector<std::unique_ptr<Element>> &List_Of_Elements)
{
    unsigned int Number_Of_Nodes = 0;
    for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
    {
        (*i)->set_pos(Number_Of_Nodes);
        Number_Of_Nodes += (*i)->get_Number_Of_Nodes();
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
void set_theta_Uniform(std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries, const double &theta)
{
    for(auto f = List_Of_Boundaries.begin(); f < List_Of_Boundaries.end(); f++)
    {
        (*f)->set_theta(theta);
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
    case 60: mesh_name = "Mesh/trapezoid_60x160.msh"; break;
    case 105: mesh_name = "Mesh/trapezoid_105x280.msh"; break;
    case 142: mesh_name = "Mesh/trapezoid_142x380.msh"; break;
    default: std::cout << "***********************\n* Mesh does not exist *\n* Defaulting to 2x6   *\n***********************"; mesh_name = "Mesh/trapezoid_2x6.msh"; break;
    }

    return mesh_name;
}
/*--------------------------------------------------------------------------*/
