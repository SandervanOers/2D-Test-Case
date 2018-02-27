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
