#include "Elements.hpp"
/*--------------------------------------------------------------------------*/
Elements::Elements(unsigned int IDg, double Jac, unsigned int IDl, unsigned int IDr, double xleft, double xright, unsigned int Order):
    ID(IDg), Jacobian(Jac), ID_Left_Boundary(IDl), ID_Right_Boundary(IDr), xCoordinateLeft(xleft), xCoordinateRight(xright), Order_Of_Polynomials(Order) {}

//Elements::Elements(const Elements &orig): ID(orig.ID), ID_Left_Boundary(orig.ID_Left_Boundary), ID_Right_Boundary(orig.ID_Right_Boundary), xCoordinateLeft(orig.xCoordinateLeft), xCoordinateRight(orig.xCoordinateRight) {}

//Elements::Elements& operator=(const Elements& that)
//{
//    ID=that.ID;
//    ID_Left_Boundary=that.ID_Left_Boundary;
//    ID_Right_Boundary=that.ID_Right_Boundary;
//    xCoordinateLeft=that.xCoordinateLeft;
//    xCoordinateRight=that.xCoordinateRight;
//}
unsigned int Elements::getID()
{
    return ID;
}
unsigned int Elements::getLeftBoundaryID()
{
    return ID_Left_Boundary;
}
unsigned int Elements::getRightBoundaryID()
{
    return ID_Right_Boundary;
}
double Elements::get_xCoordinateLeft()
{
    return xCoordinateLeft;
}
double Elements::get_xCoordinateRight()
{
    return xCoordinateRight;
}
unsigned int Elements::getOrderOfPolynomials()
{
    return Order_Of_Polynomials;
}
double Elements::getJacobian()
{
    return Jacobian;
}
void Elements::set_VertexCoordinates(double x)
{
    vertex_coordinates.push_back(x);
}
std::vector<double> Elements::get_VertexCoordinates()
{
    return vertex_coordinates;
}
void Elements::set_MassMatrix(Mat M)
{
    MatDuplicate(M, MAT_COPY_VALUES, &MassMatrix);
}
void Elements::set_MassMatrixOverRho0(Mat M1)
{
    MatDuplicate(M1, MAT_COPY_VALUES, &MassMatrixOverRho0);
}
void Elements::forget_MassMatrix()
{
    MatDestroy(&MassMatrix);
}
Mat Elements::get_MassMatrix()
{
    return MassMatrix;
}
Mat Elements::get_MassMatrixOverRho0()
{
    return MassMatrixOverRho0;
}
Elements::~Elements()
{
}
/*--------------------------------------------------------------------------*/
Boundaries::Boundaries(unsigned int IDg, int Left_Element, int Right_Element, bool isInternal, double x):
    ID(IDg), ID_Left_Element(Left_Element), ID_Right_Element(Right_Element), InternalBoundary(isInternal), xCoordinate(x) {}

unsigned int Boundaries::getID()
{
    return ID;
}
int Boundaries::getLeftElementID()
{
    return ID_Left_Element;
}
int Boundaries::getRightElementID()
{
    return ID_Right_Element;
}
bool Boundaries::isInternal()
{
    return InternalBoundary;
}
double Boundaries::getxCoordinate()
{
    return xCoordinate;
}
Boundaries::~Boundaries()
{
}
/*--------------------------------------------------------------------------*/
VertexCoordinates3D::VertexCoordinates3D(unsigned int IDg, double xC, double yC, double zC):
    ID(IDg), xCoordinate(xC), yCoordinate(yC), zCoordinate(zC) {}

unsigned int VertexCoordinates3D::getID() const
{
    return ID;
}
double VertexCoordinates3D::getxCoordinate() const
{
    return xCoordinate;
}
double VertexCoordinates3D::getyCoordinate() const
{
    return yCoordinate;
}
double VertexCoordinates3D::getzCoordinate() const
{
    return zCoordinate;
}
VertexCoordinates3D::~VertexCoordinates3D()
{
}
/*--------------------------------------------------------------------------*/
VertexCoordinates2D::VertexCoordinates2D(unsigned int IDg, double xC, double yC, bool Internal):
    ID(IDg), xCoordinate(xC), yCoordinate(yC), InternalVertex(Internal) {}

unsigned int VertexCoordinates2D::getID() const
{
    return ID;
}
double VertexCoordinates2D::getxCoordinate() const
{
    return xCoordinate;
}
double VertexCoordinates2D::getyCoordinate() const
{
    return yCoordinate;
}
bool VertexCoordinates2D::isInternal() const
{
    return InternalVertex;
}
VertexCoordinates2D::~VertexCoordinates2D()
{
}
/*--------------------------------------------------------------------------*/
Elements2D::Elements2D(unsigned int IDg, unsigned int ID_B1, unsigned int ID_B2, unsigned int ID_B3, unsigned int ID_V1, unsigned int ID_V2, unsigned int ID_V3, unsigned int N_Faces):
    ID(IDg), ID_Boundary_V1V2(ID_B1), ID_Boundary_V2V3(ID_B2), ID_Boundary_V3V1(ID_B3), ID_Vertex_V1(ID_V1), ID_Vertex_V2(ID_V2),ID_Vertex_V3(ID_V3), Number_Of_Faces(N_Faces) {}
unsigned int Elements2D::getID() const
{
    return ID;
}
unsigned int Elements2D::getVertex_V1() const
{
    return ID_Vertex_V1;
}
unsigned int Elements2D::getVertex_V2() const
{
    return ID_Vertex_V2;
}
unsigned int Elements2D::getVertex_V3() const
{
    return ID_Vertex_V3;
}
unsigned int Elements2D::getBoundary_B1() const
{
    return ID_Boundary_V1V2;
}
unsigned int Elements2D::getBoundary_B2() const
{
    return ID_Boundary_V2V3;
}
unsigned int Elements2D::getBoundary_B3() const
{
    return ID_Boundary_V3V1;
}
unsigned int Elements2D::getNumber_Of_Faces() const
{
    return Number_Of_Faces;
}
void Elements2D::setJacobian(double J)
{
    Jacobian = J;
}
void Elements2D::setArea(double A)
{
    Area = A;
}
void Elements2D::set_rx(double rx_v)
{
    rx = rx_v;
}
void Elements2D::set_ry(double ry_v)
{
    ry = ry_v;
}
void Elements2D::set_sx(double sx_v)
{
    sx = sx_v;
}
void Elements2D::set_sy(double sy_v)
{
    sy = sy_v;
}
void Elements2D::set_pos(unsigned int POS)
{
    pos = POS;
}
/*void Elements2D::computeJacobian(const std::vector<VertexCoordinates2D> &List_Of_Vertices)
{
    // ID_Vertex_V1
    // ID_Vertex_V2
    // ID_Vertex_V3
    double x1 = List_Of_Vertices[ID_Vertex_V1].getxCoordinate();
    double y1 = List_Of_Vertices[ID_Vertex_V1].getyCoordinate();
    double x2 = List_Of_Vertices[ID_Vertex_V2].getxCoordinate();
    double y2 = List_Of_Vertices[ID_Vertex_V2].getyCoordinate();
    double x3 = List_Of_Vertices[ID_Vertex_V3].getxCoordinate();
    double y3 = List_Of_Vertices[ID_Vertex_V3].getyCoordinate();

    double dxdr = (x2-x1)/2.0;
    double dydr = (y2-y1)/2.0;
    double dxds = (x3-x1)/2.0;
    double dyds = (y3-y1)/2.0;

    Jacobian = dxdr*dyds-dxds*dydr;

}*/
double Elements2D::getJacobian() const
{
    return Jacobian;
}
double Elements2D::getArea() const
{
    return Area;
}
double Elements2D::get_rx() const
{
    return rx;
}
double Elements2D::get_ry() const
{
    return ry;
}
double Elements2D::get_sx() const
{
    return sx;
}
double Elements2D::get_sy() const
{
    return sy;
}
unsigned int Elements2D::getPosition() const
{
    return pos;
}
void Elements2D::set_Order_Of_Polynomials(unsigned int N)
{
    Order_Of_Polynomials = N;
    Number_Of_Nodes = (N+1)*(N+2)/2;
}
unsigned int Elements2D::get_Order_Of_Polynomials() const
{
    return Order_Of_Polynomials;
}
unsigned int Elements2D::get_Number_Of_Nodes() const
{
    return Number_Of_Nodes;
}
void Elements2D::set_node_coordinates_x(double x)
{
    node_coordinates_x.push_back(x);
}
void Elements2D::set_node_coordinates_y(double y)
{
    node_coordinates_y.push_back(y);
}
void Elements2D::set_node_on_boundary(unsigned int type)
{
    node_on_boundary.push_back(type);
}
void Elements2D::set_node_on_boundary_1(unsigned type)
{
    node_on_boundary_1.push_back(type);
}
void Elements2D::set_node_on_boundary_2(unsigned int type)
{
    node_on_boundary_2.push_back(type);
}
void Elements2D::set_node_on_boundary_3(unsigned int type)
{
    node_on_boundary_3.push_back(type);
}
std::vector<double> Elements2D::get_node_coordinates_x() const
{
    return node_coordinates_x;
}
std::vector<double> Elements2D::get_node_coordinates_y() const
{
    return node_coordinates_y;
}
std::vector<unsigned int> Elements2D::get_node_on_boundary() const
{
    return node_on_boundary;
}
std::vector<unsigned int> Elements2D::get_node_on_boundary_1() const
{
    return node_on_boundary_1;
}
std::vector<unsigned int> Elements2D::get_node_on_boundary_2() const
{
    return node_on_boundary_2;
}
std::vector<unsigned int> Elements2D::get_node_on_boundary_3() const
{
    return node_on_boundary_3;
}
void Elements2D::setCentroidX(double setCentroidX_v)
{
    centroid_x = setCentroidX_v;
}
void Elements2D::setCentroidY(double setCentroidY_v)
{
    centroid_y = setCentroidY_v;
}
double Elements2D::get_centroid_x() const
{
    return centroid_x;
}
double Elements2D::get_centroid_y() const
{
    return centroid_y;
}
std::vector<unsigned int> Elements2D::get_nodes_on_boundary(unsigned int Type) const
{
    if (Type == 1)
    {
        return node_on_boundary_1;
    }
    else if (Type == 2)
    {
        return node_on_boundary_2;
    }
    else if (Type == 3)
    {
        return node_on_boundary_3;
    }
    else
    {
        std::cout << "Wrong Boundary Type" << std::endl;
        return {0};
    }
}
Elements2D::~Elements2D()
{
}
/*--------------------------------------------------------------------------*/
Cuboid::Cuboid(unsigned int IDg, unsigned int ID_V1, unsigned int ID_V2, unsigned int ID_V3, unsigned int ID_V4, unsigned int ID_V5, unsigned int ID_V6, unsigned int ID_V7, unsigned int ID_V8):
    ID(IDg), ID_Vertex_V1(ID_V1), ID_Vertex_V2(ID_V2), ID_Vertex_V3(ID_V3), ID_Vertex_V4(ID_V4), ID_Vertex_V5(ID_V5), ID_Vertex_V6(ID_V6), ID_Vertex_V7(ID_V7), ID_Vertex_V8(ID_V8) {}
unsigned int Cuboid::getID() const
{
    return ID;
}
unsigned int Cuboid::getVertex_V1() const
{
    return ID_Vertex_V1;
}
unsigned int Cuboid::getVertex_V2() const
{
    return ID_Vertex_V2;
}
unsigned int Cuboid::getVertex_V3() const
{
    return ID_Vertex_V3;
}
unsigned int Cuboid::getVertex_V4() const
{
    return ID_Vertex_V4;
}
unsigned int Cuboid::getVertex_V5() const
{
    return ID_Vertex_V5;
}
unsigned int Cuboid::getVertex_V6() const
{
    return ID_Vertex_V6;
}
unsigned int Cuboid::getVertex_V7() const
{
    return ID_Vertex_V7;
}
unsigned int Cuboid::getVertex_V8() const
{
    return ID_Vertex_V8;
}

Cuboid::~Cuboid()
{

}
/*--------------------------------------------------------------------------*/
InternalBoundariesCuboid::InternalBoundariesCuboid(unsigned int IDg, int ID_El_L, int ID_El_R, int Type_L, int Type_R):
    ID(IDg), ID_Left_Element(ID_El_L), ID_Right_Element(ID_El_R), Type_Left(Type_L), Type_Right(Type_R){}
unsigned int InternalBoundariesCuboid::getID() const
{
    return ID;
}
int InternalBoundariesCuboid::getLeftElementID() const
{
    return ID_Left_Element;
}
int InternalBoundariesCuboid::getRightElementID() const
{
    return ID_Right_Element;
}
//void InternalBoundariesCuboid::setJacobian(double J)
//{
//   Jacobian = J;
//}
//double InternalBoundariesCuboid::getJacobian() const
//{
//   return Jacobian;
//}
void InternalBoundariesCuboid::set_nx(double normalx)
{
    nx = normalx;
}
void InternalBoundariesCuboid::set_ny(double normaly)
{
    ny = normaly;
}
void InternalBoundariesCuboid::set_nz(double normalz)
{
    nz = normalz;
}
double InternalBoundariesCuboid::get_nx() const
{
    return nx;
}
double InternalBoundariesCuboid::get_ny() const
{
    return ny;
}
double InternalBoundariesCuboid::get_nz() const
{
    return nz;
}
unsigned int InternalBoundariesCuboid::get_Type_Left() const
{
    return Type_Left;
}
unsigned int InternalBoundariesCuboid::get_Type_Right() const
{
    return Type_Right;
}
void InternalBoundariesCuboid::set_theta(double T)
{
    theta = T;
}
double InternalBoundariesCuboid::get_theta() const
{
    return theta;
}
InternalBoundariesCuboid::~InternalBoundariesCuboid()
{}
/*--------------------------------------------------------------------------*/
Squares2D::Squares2D(unsigned int IDg, unsigned int ID_V1, unsigned int ID_V2, unsigned int ID_V3, unsigned int ID_V4):
    ID(IDg), ID_Vertex_V1(ID_V1), ID_Vertex_V2(ID_V2), ID_Vertex_V3(ID_V3), ID_Vertex_V4(ID_V4) {}
unsigned int Squares2D::getID() const
{
    return ID;
}
unsigned int Squares2D::getVertex_V1() const
{
    return ID_Vertex_V1;
}
unsigned int Squares2D::getVertex_V2() const
{
    return ID_Vertex_V2;
}
unsigned int Squares2D::getVertex_V3() const
{
    return ID_Vertex_V3;
}
unsigned int Squares2D::getVertex_V4() const
{
    return ID_Vertex_V4;
}
unsigned int Squares2D::getPosition() const
{
    return pos;
}
void Squares2D::set_pos(unsigned int POS)
{
    pos = POS;
}
void Squares2D::setJacobian(double J)
{
    Jacobian = J;
}
void Squares2D::setArea(double A)
{
    Area = A;
}
void Squares2D::set_rx(double rx_v)
{
    rx = rx_v;
}
void Squares2D::set_ry(double ry_v)
{
    ry = ry_v;
}
void Squares2D::set_sx(double sx_v)
{
    sx = sx_v;
}
void Squares2D::set_sy(double sy_v)
{
    sy = sy_v;
}
double Squares2D::getJacobian() const
{
    return Jacobian;
}
double Squares2D::getArea() const
{
    return Area;
}
double Squares2D::get_rx() const
{
    return rx;
}
double Squares2D::get_ry() const
{
    return ry;
}
double Squares2D::get_sx() const
{
    return sx;
}
double Squares2D::get_sy() const
{
    return sy;
}
void Squares2D::set_Order_Of_Polynomials(unsigned int N)
{
    Order_Of_Polynomials = N;
    Number_Of_Nodes = (N+1)*(N+1);
}
unsigned int Squares2D::get_Order_Of_Polynomials() const
{
    return Order_Of_Polynomials;
}
unsigned int Squares2D::get_Number_Of_Nodes() const
{
    return Number_Of_Nodes;
}
void Squares2D::set_node_coordinates_x(double x)
{
    node_coordinates_x.push_back(x);
}
void Squares2D::set_node_coordinates_y(double y)
{
    node_coordinates_y.push_back(y);
}
void Squares2D::set_node_on_boundary(unsigned int type)
{
    node_on_boundary.push_back(type);
}
void Squares2D::set_node_on_boundary_1(unsigned type)
{
    node_on_boundary_1.push_back(type);
}
void Squares2D::set_node_on_boundary_2(unsigned int type)
{
    node_on_boundary_2.push_back(type);
}
void Squares2D::set_node_on_boundary_3(unsigned int type)
{
    node_on_boundary_3.push_back(type);
}
void Squares2D::set_node_on_boundary_4(unsigned int type)
{
    node_on_boundary_4.push_back(type);
}
std::vector<double> Squares2D::get_node_coordinates_x() const
{
    return node_coordinates_x;
}
std::vector<double> Squares2D::get_node_coordinates_y() const
{
    return node_coordinates_y;
}
std::vector<unsigned int> Squares2D::get_node_on_boundary() const
{
    return node_on_boundary;
}
std::vector<unsigned int> Squares2D::get_node_on_boundary_1() const
{
    return node_on_boundary_1;
}
std::vector<unsigned int> Squares2D::get_node_on_boundary_2() const
{
    return node_on_boundary_2;
}
std::vector<unsigned int> Squares2D::get_node_on_boundary_3() const
{
    return node_on_boundary_3;
}
std::vector<unsigned int> Squares2D::get_node_on_boundary_4() const
{
    return node_on_boundary_4;
}

std::vector<unsigned int> Squares2D::get_nodes_on_boundary(unsigned int Type) const
{
    if (Type == 1)
    {
        return node_on_boundary_1;
    }
    else if (Type == 2)
    {
        return node_on_boundary_2;
    }
    else if (Type == 3)
    {
        return node_on_boundary_3;
    }
    else if (Type == 4)
    {
        return node_on_boundary_4;
    }
    else
    {
        std::cout << "Wrong Boundary Type" << std::endl;
        return {0};
    }
}

Squares2D::~Squares2D()
{

}
/*--------------------------------------------------------------------------*/
InternalBoundariesSquares2D::InternalBoundariesSquares2D(unsigned int IDg, int ID_El_L, int ID_El_R, int Type_L, int Type_R):
    ID(IDg), ID_Left_Element(ID_El_L), ID_Right_Element(ID_El_R), Type_Left(Type_L), Type_Right(Type_R){}
unsigned int InternalBoundariesSquares2D::getID() const
{
    return ID;
}
int InternalBoundariesSquares2D::getLeftElementID() const
{
    return ID_Left_Element;
}
int InternalBoundariesSquares2D::getRightElementID() const
{
    return ID_Right_Element;
}
void InternalBoundariesSquares2D::set_Vertex_V1(unsigned int V1)
{
    ID_Vertex_V1 = V1;
}
void InternalBoundariesSquares2D::set_Vertex_V2(unsigned int V2)
{
    ID_Vertex_V2 = V2;
}
unsigned int InternalBoundariesSquares2D::get_Vertex_V1() const
{
    return ID_Vertex_V1;
}
unsigned int InternalBoundariesSquares2D::get_Vertex_V2() const
{
    return ID_Vertex_V2;
}
void InternalBoundariesSquares2D::setJacobian(double J)
{
    Jacobian = J;
}
double InternalBoundariesSquares2D::getJacobian() const
{
    return Jacobian;
}
void InternalBoundariesSquares2D::set_nx(double normalx)
{
    nx = normalx;
}
void InternalBoundariesSquares2D::set_ny(double normaly)
{
    ny = normaly;
}
double InternalBoundariesSquares2D::get_nx() const
{
    return nx;
}
double InternalBoundariesSquares2D::get_ny() const
{
    return ny;
}
unsigned int InternalBoundariesSquares2D::get_Type_Left() const
{
    return Type_Left;
}
unsigned int InternalBoundariesSquares2D::get_Type_Right() const
{
    return Type_Right;
}
void InternalBoundariesSquares2D::set_theta(double T)
{
    theta = T;
}
double InternalBoundariesSquares2D::get_theta() const
{
    return theta;
}
InternalBoundariesSquares2D::~InternalBoundariesSquares2D()
{}


/*--------------------------------------------------------------------------*/
Boundaries2D::Boundaries2D(unsigned int IDg, bool Internal, unsigned int ID_V1, unsigned int ID_V2, int ID_El_L, int ID_El_R, int Typeg):
    ID(IDg), InternalBoundary(Internal), ID_Vertex_V1(ID_V1), ID_Vertex_V2(ID_V2), ID_Left_Element(ID_El_L), ID_Right_Element(ID_El_R), Type(Typeg) {}
unsigned int Boundaries2D::getID() const
{
    return ID;
}
int Boundaries2D::getLeftElementID()
{
    return ID_Left_Element;
}
int Boundaries2D::getRightElementID()
{
    return ID_Right_Element;
}
//void Boundaries2D::set_LeftElementID(int E_L)
//{
//    this->ID_Left_Element = E_L;
//}
//void Boundaries2D::set_RightElementID(int E_R)
//{
//    this->ID_Right_Element = E_R;
//}
void Boundaries2D::setJacobian(double J)
{
    Jacobian = J;
}
double Boundaries2D::getJacobian() const
{
    return Jacobian;
}
void Boundaries2D::setNormalX(double normalx)
{
    nx = normalx;
}
void Boundaries2D::setNormalY(double normaly)
{
    ny = normaly;
}
double Boundaries2D::get_nx() const
{
    return nx;
}
double Boundaries2D::get_ny() const
{
    return ny;
}
//void Boundaries2D::setType(double T)
//{
//    Type = T;
//}
unsigned int Boundaries2D::getType() const
{
    return Type;
}
bool Boundaries2D::isInternal()
{
    return InternalBoundary;
}
void Boundaries2D::set_node_coordinates_x(double x)
{
    node_coordinates_x.push_back(x);
}
void Boundaries2D::set_node_coordinates_y(double y)
{
    node_coordinates_y.push_back(y);
}
/*std::vector<double> Boundaries2D::get_node_coordinates_x() const
{
    return node_coordinates_x;
}
std::vector<double> Boundaries2D::get_node_coordinates_y() const
{
    return node_coordinates_y;
}
*/
unsigned int Boundaries2D::getVertex_V1() const
{
    return ID_Vertex_V1;
}
unsigned int Boundaries2D::getVertex_V2() const
{
    return ID_Vertex_V2;
}
void Boundaries2D::set_theta(double T)
{
    theta = T;
}
double Boundaries2D::get_theta() const
{
    return theta;
}
Boundaries2D::~Boundaries2D()
{
}
/*--------------------------------------------------------------------------*/
