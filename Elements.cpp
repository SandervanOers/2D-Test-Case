#include "Elements.hpp"
/*--------------------------------------------------------------------------*/
Vertex::Vertex(unsigned int IDg, unsigned int DIMg, double xC, double yC, double zC):
    ID(IDg), DIM(DIMg), xCoordinate(xC), yCoordinate(yC), zCoordinate(zC) {}
unsigned int Vertex::getID() const
{
    return ID;
}
unsigned int Vertex::getDIM() const
{
    return DIM;
}
double Vertex::getxCoordinate() const
{
    return xCoordinate;
}
double Vertex::getyCoordinate() const
{
    return yCoordinate;
}
double Vertex::getzCoordinate() const
{
    return zCoordinate;
}
Vertex::~Vertex()
{
}
extern void print(const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices)
{
    std::cout << "List of Vertices "  << std::endl;
    unsigned int d = List_Of_Vertices.front()->getDIM();
    switch(d)
    {
        case 1:
            std::cout << "ID : x"  << std::endl;
            for(auto i = List_Of_Vertices.begin(); i < List_Of_Vertices.end(); i++)
                std::cout << (*i)->getID() << ": " << (*i)->getxCoordinate() << std::endl;
        break;
        case 2:
            std::cout << "ID : x y"  << std::endl;
            for(auto i = List_Of_Vertices.begin(); i < List_Of_Vertices.end(); i++)
                std::cout << (*i)->getID() << ": " << (*i)->getxCoordinate() << " " << (*i)->getyCoordinate() << std::endl;
        break;
        case 3:
            std::cout << "ID : x y z"  << std::endl;
            for(auto i = List_Of_Vertices.begin(); i < List_Of_Vertices.end(); i++)
                std::cout << (*i)->getID() << ": " << (*i)->getxCoordinate() << " " << (*i)->getyCoordinate() << " " << (*i)->getzCoordinate() << std::endl;
        break;
        default:
            std::cout << "Wrong dimensions of vertices." << std::endl;
    }

}
/*--------------------------------------------------------------------------*/
Boundary::Boundary(unsigned int IDg, unsigned int DIMg, int ID_El_L, int ID_El_R, int Type_L, int Type_R):
    ID(IDg), DIM(DIMg), ID_Left_Element(ID_El_L), ID_Right_Element(ID_El_R), Type_Left(Type_L), Type_Right(Type_R) {}
unsigned int Boundary::getID() const
{
    return ID;
}
unsigned int Boundary::getDIM() const
{
    return DIM;
}
int Boundary::getLeftElementID() const
{
    return ID_Left_Element;
}
int Boundary::getRightElementID() const
{
    return ID_Right_Element;
}
unsigned int Boundary::getTypeLeft() const
{
    return Type_Left;
}
unsigned int Boundary::getTypeRight() const
{
    return Type_Right;
}
void Boundary::setJacobian(double J)
{
    Jacobian = J;
}
double Boundary::getJacobian() const
{
    return Jacobian;
}
void Boundary::set_nx(double normalx)
{
    nx = normalx;
}
void Boundary::set_ny(double normaly)
{
    ny = normaly;
}
void Boundary::set_nz(double normalz)
{
    nz = normalz;
}
double Boundary::get_nx() const
{
    return nx;
}
double Boundary::get_ny() const
{
    return ny;
}
double Boundary::get_nz() const
{
    return nz;
}
void Boundary::set_theta(double T)
{
    theta = T;
}
double Boundary::get_theta() const
{
    return theta;
}
Boundary::~Boundary()
{
}
extern void print(const std::vector<std::unique_ptr<Boundary>> &List_Of_Boundaries)
{
    std::cout << "List of Boundaries "  << std::endl;

            std::cout << "ID : Left Element ID Right Element ID Type Left Face Type Right Face"  << std::endl;
            for(auto i = List_Of_Boundaries.begin(); i < List_Of_Boundaries.end(); i++)
                std::cout << (*i)->getID() << ": " << (*i)->getLeftElementID() << " " << (*i)->getRightElementID() << " " << (*i)-> getTypeLeft() << " " << (*i)->getTypeRight() << std::endl;
}
/*--------------------------------------------------------------------------*/
Element::Element(unsigned int IDg, unsigned int DIMg, int ID_V1, int ID_V2, int ID_V3, int ID_V4, int ID_V5, int ID_V6, int ID_V7, int ID_V8):
    ID(IDg), DIM(DIMg), ID_Vertex_V1(ID_V1), ID_Vertex_V2(ID_V2), ID_Vertex_V3(ID_V3), ID_Vertex_V4(ID_V4), ID_Vertex_V5(ID_V5), ID_Vertex_V6(ID_V6), ID_Vertex_V7(ID_V7), ID_Vertex_V8(ID_V8) {}
unsigned int Element::getID() const
{
    return ID;
}
unsigned int Element::getDIM() const
{
    return DIM;
}
void Element::set_pos(unsigned int POS)
{
    pos = POS;
}
unsigned int Element::get_pos() const
{
    return pos;
}
unsigned int Element::getVertex_V1() const
{
    return ID_Vertex_V1;
}
unsigned int Element::getVertex_V2() const
{
    return ID_Vertex_V2;
}
unsigned int Element::getVertex_V3() const
{
    return ID_Vertex_V3;
}
unsigned int Element::getVertex_V4() const
{
    return ID_Vertex_V4;
}
unsigned int Element::getVertex_V5() const
{
    return ID_Vertex_V5;
}
unsigned int Element::getVertex_V6() const
{
    return ID_Vertex_V6;
}
unsigned int Element::getVertex_V7() const
{
    return ID_Vertex_V7;
}
unsigned int Element::getVertex_V8() const
{
    return ID_Vertex_V8;
}
void Element::set_Order_Of_Polynomials_x(unsigned int N)
{
    Order_Of_Polynomials_x = N;
}
void Element::set_Order_Of_Polynomials_y(unsigned int N)
{
    Order_Of_Polynomials_y = N;
}
void Element::set_Order_Of_Polynomials_z(unsigned int N)
{
    Order_Of_Polynomials_z = N;
}
unsigned int Element::get_Order_Of_Polynomials_x() const
{
    return Order_Of_Polynomials_x;
}
unsigned int Element::get_Order_Of_Polynomials_y() const
{
    return Order_Of_Polynomials_y;
}
unsigned int Element::get_Order_Of_Polynomials_z() const
{
    return Order_Of_Polynomials_z;
}
unsigned int Element::get_Number_Of_Nodes() const
{
    return Number_Of_Nodes;
}
void Element::set_Number_Of_Nodes(unsigned int Nnodes)
{
    Number_Of_Nodes = Nnodes;
}
void Element::set_node_coordinates_x(double x)
{
    node_coordinates_x.push_back(x);
}
void Element::set_node_coordinates_y(double y)
{
    node_coordinates_y.push_back(y);
}
void Element::set_node_coordinates_z(double z)
{
    node_coordinates_z.push_back(z);
}
std::vector<double> Element::get_node_coordinates_x() const
{
    return node_coordinates_x;
}
std::vector<double> Element::get_node_coordinates_y() const
{
    return node_coordinates_y;
}
std::vector<double> Element::get_node_coordinates_z() const
{
    return node_coordinates_z;
}
void Element::set_node_on_face0(unsigned int k)
{
    node_on_face0.push_back(k);
}
void Element::set_node_on_face1(unsigned int k)
{
    node_on_face1.push_back(k);
}
void Element::set_node_on_face2(unsigned int k)
{
    node_on_face2.push_back(k);
}
void Element::set_node_on_face3(unsigned int k)
{
    node_on_face3.push_back(k);
}
void Element::set_node_on_face4(unsigned int k)
{
    node_on_face4.push_back(k);
}
void Element::set_node_on_face5(unsigned int k)
{
    node_on_face5.push_back(k);
}
std::vector<unsigned int> Element::get_node_on_face0() const
{
    return node_on_face0;
}
std::vector<unsigned int> Element::get_node_on_face1() const
{
    return node_on_face1;
}
std::vector<unsigned int> Element::get_node_on_face2() const
{
    return node_on_face2;
}
std::vector<unsigned int> Element::get_node_on_face3() const
{
    return node_on_face3;
}
std::vector<unsigned int> Element::get_node_on_face4() const
{
    return node_on_face4;
}
std::vector<unsigned int> Element::get_node_on_face5() const
{
    return node_on_face5;
}
std::vector<unsigned int> Element::get_nodes_on_boundary(unsigned int Type) const
{
    switch(Type)
    {
    case 0:
        return node_on_face0;
        break;
    case 1:
        return node_on_face1;
        break;
    case 2:
        return node_on_face2;
        break;
    case 3:
        return node_on_face3;
        break;
    case 4:
        return node_on_face4;
        break;
    case 5:
        return node_on_face5;
        break;
    default:
        std::cout << "Wrong Boundary Type" << std::endl;
        return {0};
    }
}
Element::~Element()
{
}
extern void print(const std::vector<std::unique_ptr<Element>> &List_Of_Elements)
{
    std::cout << "List of Element "  << std::endl;
    unsigned int d = List_Of_Elements.front()->getDIM();
    switch(d)
    {
        case 1:
            std::cout << "ID : V1 V2"  << std::endl;
            for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
                std::cout << (*i)->getID() << ": " << (*i)->getVertex_V1() << " " << (*i)->getVertex_V2() << std::endl;
        break;
        case 2:
            std::cout << "ID : V1 V2 V3 V4"  << std::endl;
            for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
                std::cout << (*i)->getID() << ": " << (*i)->getVertex_V1() << " " << (*i)->getVertex_V2() << " " << (*i)->getVertex_V3() << " " << (*i)->getVertex_V4() << std::endl;
        break;
        case 3:
            std::cout << "ID : V1 V2 V3 V4 V5 V6 V7 V8"  << std::endl;
            for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
                std::cout << (*i)->getID() << ": " << (*i)->getVertex_V1() << " " << (*i)->getVertex_V2() << " " << (*i)->getVertex_V3() << " " << (*i)->getVertex_V4() << " " << (*i)->getVertex_V5() << " " << (*i)->getVertex_V6() << " " << (*i)->getVertex_V7() << " " << (*i)->getVertex_V8() << std::endl;
        break;
        default:
            std::cout << "Wrong dimensions of vertices." << std::endl;
    }

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
void Cuboid::set_Order_Of_Polynomials_x(unsigned int N)
{
    Order_Of_Polynomials_x = N;
}
void Cuboid::set_Order_Of_Polynomials_y(unsigned int N)
{
    Order_Of_Polynomials_y = N;
}
void Cuboid::set_Order_Of_Polynomials_z(unsigned int N)
{
    Order_Of_Polynomials_z = N;
}
unsigned int Cuboid::get_Order_Of_Polynomials_x() const
{
    return Order_Of_Polynomials_x;
}
unsigned int Cuboid::get_Order_Of_Polynomials_y() const
{
    return Order_Of_Polynomials_y;
}
unsigned int Cuboid::get_Order_Of_Polynomials_z() const
{
    return Order_Of_Polynomials_z;
}
unsigned int Cuboid::get_Number_Of_Nodes() const
{
    return Number_Of_Nodes;
}
void Cuboid::set_Number_Of_Nodes(unsigned int Nnodes)
{
    Number_Of_Nodes = Nnodes;
}
unsigned int Cuboid::get_pos() const
{
    return pos;
}
void Cuboid::set_pos(unsigned int POS)
{
    pos = POS;
}
void Cuboid::set_node_on_face0(unsigned int k)
{
    node_on_face0.push_back(k);
}
void Cuboid::set_node_on_face1(unsigned int k)
{
    node_on_face1.push_back(k);
}
void Cuboid::set_node_on_face2(unsigned int k)
{
    node_on_face2.push_back(k);
}
void Cuboid::set_node_on_face3(unsigned int k)
{
    node_on_face3.push_back(k);
}
void Cuboid::set_node_on_face4(unsigned int k)
{
    node_on_face4.push_back(k);
}
void Cuboid::set_node_on_face5(unsigned int k)
{
    node_on_face5.push_back(k);
}
std::vector<unsigned int> Cuboid::get_node_on_face0() const
{
    return node_on_face0;
}
std::vector<unsigned int> Cuboid::get_node_on_face1() const
{
    return node_on_face1;
}
std::vector<unsigned int> Cuboid::get_node_on_face2() const
{
    return node_on_face2;
}
std::vector<unsigned int> Cuboid::get_node_on_face3() const
{
    return node_on_face3;
}
std::vector<unsigned int> Cuboid::get_node_on_face4() const
{
    return node_on_face4;
}
std::vector<unsigned int> Cuboid::get_node_on_face5() const
{
    return node_on_face5;
}
std::vector<unsigned int> Cuboid::get_nodes_on_boundary(unsigned int Type) const
{
    switch(Type)
    {
    case 0:
        return node_on_face0;
        break;
    case 1:
        return node_on_face1;
        break;
    case 2:
        return node_on_face2;
        break;
    case 3:
        return node_on_face3;
        break;
    case 4:
        return node_on_face4;
        break;
    case 5:
        return node_on_face5;
        break;
    default:
        std::cout << "Wrong Boundary Type" << std::endl;
        return {0};
    }
}

Cuboid::~Cuboid()
{

}
/*--------------------------------------------------------------------------*/
/*InternalBoundariesCuboid::InternalBoundariesCuboid(unsigned int IDg, int ID_El_L, int ID_El_R, int Type_L, int Type_R):
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
void InternalBoundariesCuboid::setJacobian(double J)
{
   Jacobian = J;
}
double InternalBoundariesCuboid::getJacobian() const
{
   return Jacobian;
}
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
{}*/
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
