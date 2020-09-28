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
std::vector<double> Elements2D::get_node_coordinates_x() const
{
    return node_coordinates_x;
}
std::vector<double> Elements2D::get_node_coordinates_y() const
{
    return node_coordinates_y;
}

Elements2D::~Elements2D()
{
}
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
std::vector<double> Boundaries2D::get_node_coordinates_x() const
{
    return node_coordinates_x;
}
std::vector<double> Boundaries2D::get_node_coordinates_y() const
{
    return node_coordinates_y;
}
unsigned int Boundaries2D::getVertex_V1() const
{
    return ID_Vertex_V1;
}
unsigned int Boundaries2D::getVertex_V2() const
{
    return ID_Vertex_V2;
}
Boundaries2D::~Boundaries2D()
{
}
/*--------------------------------------------------------------------------*/
