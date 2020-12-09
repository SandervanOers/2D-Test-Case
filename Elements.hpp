#ifndef ELEMENTS
#define ELEMENTS

#include <petscksp.h>
#include <vector>
#include <iostream>
/*--------------------------------------------------------------------------*/
class Elements
{
    unsigned int ID;
    double Jacobian;
    unsigned int ID_Left_Boundary;
    unsigned int ID_Right_Boundary;
    double xCoordinateLeft;
    double xCoordinateRight;
    unsigned int Order_Of_Polynomials;
    Mat MassMatrix;
    Mat MassMatrixOverRho0;
    std::vector<double> vertex_coordinates;
    public:
        Elements(unsigned int IDg, double Jacobian, unsigned int ID_Left_Boundary, unsigned int ID_Right_Boundary, double xCoordinateLeft, double xCoordinateRight, unsigned int Order_Of_Polynomials);
        //Elements(const Elements&);
        //Elements& operator=(const Elements& that);
        unsigned int getID();
        unsigned int getLeftBoundaryID();
        unsigned int getRightBoundaryID();
        double get_xCoordinateLeft();
        double get_xCoordinateRight();
        unsigned int getOrderOfPolynomials();
        double getJacobian();

        void set_VertexCoordinates(double x);
        std::vector<double> get_VertexCoordinates();

        void set_MassMatrix(Mat M);
        void set_MassMatrixOverRho0(Mat M1);
        void forget_MassMatrix();
        Mat get_MassMatrix();
        Mat get_MassMatrixOverRho0();

        ~Elements();
};
/*--------------------------------------------------------------------------*/
class Boundaries
{
    unsigned int ID;
    int ID_Left_Element;
    int ID_Right_Element;
    bool InternalBoundary;
    double xCoordinate;
    //Mat MassMatrix;     // G_ij = \int l_i l_j d\Gamma
    //Mat MassMatrixTimesRho0;
    //double Jacibian;
    public:
        Boundaries(unsigned int IDg, int Left_Element, int Right_Element, bool isInternal, double x);
        unsigned int getID();
        int getLeftElementID();
        int getRightElementID();
        bool isInternal();
        double getxCoordinate();

        ~Boundaries();

};
/*--------------------------------------------------------------------------*/
class VertexCoordinates2D
{
    unsigned int ID;
    double xCoordinate;
    double yCoordinate;
    bool InternalVertex;

    public:
        VertexCoordinates2D(unsigned int IDg, double xC, double yC, bool Internal);
        unsigned int getID() const;
        double getxCoordinate() const;
        double getyCoordinate() const;
        bool isInternal() const;


        ~VertexCoordinates2D();
};
/*--------------------------------------------------------------------------*/
class Elements2D
{
    unsigned int ID;
    double Jacobian;
    double Area;
    double rx;
    double ry;
    double sx;
    double sy;
    double centroid_x;
    double centroid_y;
    unsigned int ID_Boundary_V1V2;
    unsigned int ID_Boundary_V2V3;
    unsigned int ID_Boundary_V3V1;

    unsigned int ID_Vertex_V1;
    unsigned int ID_Vertex_V2;
    unsigned int ID_Vertex_V3;

    unsigned int Number_Of_Faces;
    unsigned int Order_Of_Polynomials;

    std::vector<double> node_coordinates_x;
    std::vector<double> node_coordinates_y;
    std::vector<unsigned int> node_on_boundary;     // for all nodes: on which boundary
    std::vector<unsigned int> node_on_boundary_1;   // node numbers on boundary 1
    std::vector<unsigned int> node_on_boundary_2;
    std::vector<unsigned int> node_on_boundary_3;

    unsigned int Number_Of_Nodes;
    unsigned int pos;

    public:
        Elements2D(unsigned int IDg, unsigned int ID_B1, unsigned int ID_B2, unsigned int ID_B3, unsigned int ID_V1, unsigned int ID_V2, unsigned int ID_V3, unsigned int N_Faces);
        //Elements(const Elements&);
        //Elements& operator=(const Elements& that);
        unsigned int getID() const;
        unsigned int getVertex_V1() const;
        unsigned int getVertex_V2() const;
        unsigned int getVertex_V3() const;
        unsigned int getBoundary_B1() const;
        unsigned int getBoundary_B2() const;
        unsigned int getBoundary_B3() const;
        unsigned int getNumber_Of_Faces() const;
        unsigned int getPosition() const;
        void setJacobian(double J);
        void setArea(double A);
        void set_rx(double rx_v);
        void set_ry(double ry_v);
        void set_sx(double sx_v);
        void set_sy(double sy_v);
        void setCentroidX(double setCentroidX_v);
        void setCentroidY(double setCentroidY_v);
        //void computeJacobian(const std::vector<VertexCoordinates2D> &List_Of_Vertices);
        double getJacobian() const;
        double getArea() const;
        double get_rx() const;
        double get_ry() const;
        double get_sx() const;
        double get_sy() const;
        double get_centroid_x() const;
        double get_centroid_y() const;
        void set_Order_Of_Polynomials(unsigned int N);
        unsigned int get_Order_Of_Polynomials() const;
        //void set_Number_Of_Nodes(unsigned int N);
        unsigned int get_Number_Of_Nodes() const;
        void set_pos(unsigned int POS);
        //void set_node_coordinates_x(double x);
        //void set_node_coordinates_y(double y);
        //std::vector<double> get_node_coordinates_x();
        //std::vector<double> get_node_coordinates_y();

        void set_node_coordinates_x(double x);
        void set_node_coordinates_y(double y);
        void set_node_on_boundary(unsigned int type);
        void set_node_on_boundary_1(unsigned int type);
        void set_node_on_boundary_2(unsigned int type);
        void set_node_on_boundary_3(unsigned int type);
        std::vector<double> get_node_coordinates_x() const;
        std::vector<double> get_node_coordinates_y() const;
        std::vector<unsigned int> get_node_on_boundary() const;
        std::vector<unsigned int> get_node_on_boundary_1() const;
        std::vector<unsigned int> get_node_on_boundary_2() const;
        std::vector<unsigned int> get_node_on_boundary_3() const;
        std::vector<unsigned int> get_nodes_on_boundary(unsigned int Type) const;

        //void set_MassMatrix(Mat M);
        //void set_MassMatrixOverRho0(Mat M1);
        //void forget_MassMatrix();
        //Mat get_MassMatrix();
        //Mat get_MassMatrixOverRho0();

        ~Elements2D();
};
/*--------------------------------------------------------------------------*/
class Squares2D
{
    unsigned int ID;
    double Jacobian;
    double rx;
    double ry;
    double sx;
    double sy;
    double Area;
    unsigned int ID_Boundary_V1V2;
    unsigned int ID_Boundary_V2V3;
    unsigned int ID_Boundary_V3V4;
    unsigned int ID_Boundary_V4V1;

    unsigned int ID_Vertex_V1;
    unsigned int ID_Vertex_V2;
    unsigned int ID_Vertex_V3;
    unsigned int ID_Vertex_V4;

    unsigned int Number_Of_Faces;
    unsigned int Order_Of_Polynomials;

    std::vector<double> node_coordinates_x;
    std::vector<double> node_coordinates_y;
    std::vector<unsigned int> node_on_boundary;     // for all nodes: on which boundary
    std::vector<unsigned int> node_on_boundary_1;   // node numbers on boundary 1
    std::vector<unsigned int> node_on_boundary_2;
    std::vector<unsigned int> node_on_boundary_3;
    std::vector<unsigned int> node_on_boundary_4;

    unsigned int Number_Of_Nodes;
    unsigned int pos;
  public:
    Squares2D(unsigned int IDg, unsigned int ID_V1, unsigned int ID_V2, unsigned int ID_V3, unsigned int ID_V4);
    unsigned int getID() const ;
    unsigned int getVertex_V1() const;
    unsigned int getVertex_V2() const;
    unsigned int getVertex_V3() const;
    unsigned int getVertex_V4() const;
    unsigned int getPosition() const;
    void setJacobian(double J);
    void setArea(double A);
    void set_rx(double rx_v);
    void set_ry(double ry_v);
    void set_sx(double sx_v);
    void set_sy(double sy_v);
    double getJacobian() const;
    double getArea() const;
    double get_rx() const;
    double get_ry() const;
    double get_sx() const;
    double get_sy() const;
    double get_centroid_x() const;
    double get_centroid_y() const;
    void set_Order_Of_Polynomials(unsigned int N);
    unsigned int get_Order_Of_Polynomials() const;
    unsigned int get_Number_Of_Nodes() const;
    void set_pos(unsigned int POS);

    void set_node_coordinates_x(double x);
    void set_node_coordinates_y(double y);
    void set_node_on_boundary(unsigned int type);
    void set_node_on_boundary_1(unsigned int type);
    void set_node_on_boundary_2(unsigned int type);
    void set_node_on_boundary_3(unsigned int type);
    void set_node_on_boundary_4(unsigned int type);
    std::vector<double> get_node_coordinates_x() const;
    std::vector<double> get_node_coordinates_y() const;
    std::vector<unsigned int> get_node_on_boundary() const;
    std::vector<unsigned int> get_node_on_boundary_1() const;
    std::vector<unsigned int> get_node_on_boundary_2() const;
    std::vector<unsigned int> get_node_on_boundary_3() const;
    std::vector<unsigned int> get_node_on_boundary_4() const;
    std::vector<unsigned int> get_nodes_on_boundary(unsigned int Type) const;

    ~Squares2D();
};
/*--------------------------------------------------------------------------*/
class InternalBoundariesSquares2D
{
    unsigned int ID;
    int ID_Left_Element;
    int ID_Right_Element;
    double Jacobian;
    unsigned int Type_Left;  // 1 2 3 4
    unsigned int Type_Right; // 1 2 3 4
    double theta;
    double nx;
    double ny;
    unsigned int ID_Vertex_V1;
    unsigned int ID_Vertex_V2;

    //std::vector<double> node_coordinates_x;
    //std::vector<double> node_coordinates_y;

    public:
        InternalBoundariesSquares2D(unsigned int IDg, int ID_El_L, int ID_El_R, int Type_L, int Type_R);
        unsigned int getID() const ;
        int getLeftElementID() const;
        int getRightElementID() const;
        void setJacobian(double J);
        void set_Vertex_V1(unsigned int V1);
        void set_Vertex_V2(unsigned int V2);
        unsigned int get_Vertex_V1() const;
        unsigned int get_Vertex_V2() const;
        double getJacobian() const;
        void set_nx(double normalx);
        void set_ny(double normaly);
        double get_nx() const;
        double get_ny() const;
        unsigned int get_Type_Left() const;
        unsigned int get_Type_Right() const;

        bool operator < (const InternalBoundariesSquares2D& other) const {return ID < other.ID; }

        //void set_node_coordinates_x(double x);
        //void set_node_coordinates_y(double y);
        void set_theta(double T);
        double get_theta() const;

        ~InternalBoundariesSquares2D();

};
/*--------------------------------------------------------------------------*/
class Boundaries2D
{
    unsigned int ID;
    int ID_Left_Element;
    int ID_Right_Element;
    bool InternalBoundary;
    unsigned int ID_Vertex_V1;
    unsigned int ID_Vertex_V2;
    double Jacobian;
    unsigned int Type; // 1 2 3
    double theta;
    double nx;
    double ny;

    std::vector<double> node_coordinates_x;
    std::vector<double> node_coordinates_y;

    //double xCoordinate;
    //Mat MassMatrix;     // G_ij = \int l_i l_j d\Gamma
    //Mat MassMatrixTimesRho0;
    //double Jacibian;
    public:
        //Boundaries2D(unsigned int IDg, bool Internal, unsigned int ID_V1, unsigned int ID_V2);
        Boundaries2D(unsigned int IDg, bool Internal, unsigned int ID_V1, unsigned int ID_V2, int ID_El_L, int ID_El_R, int Typeg);
        unsigned int getID() const ;
        int getLeftElementID();
        int getRightElementID();
        void setJacobian(double J);
        //void setType(double T);
        unsigned int getType() const;
        double getJacobian() const;
        void setNormalX(double normalx);
        void setNormalY(double normaly);
        double get_nx() const;
        double get_ny() const;

        //void set_LeftElementID(int E_L);
        //void set_RightElementID(int E_R);
        bool isInternal();
        //bool operator<(const Boundaries2D& lhs, const Boundaries2D& rhs) {return lhs.ID<rhs.ID;}
        bool operator < (const Boundaries2D& other) const {return ID < other.ID; }

        void set_node_coordinates_x(double x);
        void set_node_coordinates_y(double y);
        void set_theta(double T);
        //std::vector<double> get_node_coordinates_x() const;  Legendre_Gauss_Lobatto is wrong: x & y coordinates from both elements
        //std::vector<double> get_node_coordinates_y() const;
        unsigned int getVertex_V1() const;
        unsigned int getVertex_V2() const;
        double get_theta() const;

        ~Boundaries2D();

};
#endif
