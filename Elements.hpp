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
        unsigned int getID();
        double getxCoordinate();
        double getyCoordinate();
        bool isInternal();


        ~VertexCoordinates2D();
};
/*--------------------------------------------------------------------------*/
class Elements2D
{
    unsigned int ID;
    double Jacobian;
    unsigned int ID_Boundary_V1V2;
    unsigned int ID_Boundary_V2V3;
    unsigned int ID_Boundary_V3V1;

    unsigned int ID_Vertex_V1;
    unsigned int ID_Vertex_V2;
    unsigned int ID_Vertex_V3;

    unsigned int Number_Of_Faces;
    unsigned int Order_Of_Polynomials;
    Mat MassMatrix;
    Mat MassMatrixOverRho0;
    std::vector<double> node_coordinates_x; // Should we turn this into a list, can we overwrite/change this?
    std::vector<double> node_coordinates_y;
    unsigned int Number_Of_Nodes;

    public:
        Elements2D(unsigned int IDg, unsigned int ID_B1, unsigned int ID_B2, unsigned int ID_B3, unsigned int ID_V1, unsigned int ID_V2, unsigned int ID_V3, unsigned int N_Faces);
        //Elements(const Elements&);
        //Elements& operator=(const Elements& that);
        unsigned int getID();
        unsigned int getVertex_V1();
        unsigned int getVertex_V2();
        unsigned int getVertex_V3();
        unsigned int getBoundary_B1();
        unsigned int getBoundary_B2();
        unsigned int getBoundary_B3();
        unsigned int getNumber_Of_Faces();
        void setJacobian(double J);
        double getJacobian();
        void set_Order_Of_Polynomials(unsigned int N);
        unsigned int get_Order_Of_Polynomials();
        //void set_Number_Of_Nodes(unsigned int N);
        unsigned int get_Number_Of_Nodes();

        void set_node_coordinates_x(double x);
        void set_node_coordinates_y(double y);
        std::vector<double> get_node_coordinates_x();
        std::vector<double> get_node_coordinates_y();

        //void set_MassMatrix(Mat M);
        //void set_MassMatrixOverRho0(Mat M1);
        //void forget_MassMatrix();
        //Mat get_MassMatrix();
        //Mat get_MassMatrixOverRho0();

        ~Elements2D();
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
    //double xCoordinate;
    //Mat MassMatrix;     // G_ij = \int l_i l_j d\Gamma
    //Mat MassMatrixTimesRho0;
    //double Jacibian;
    public:
        //Boundaries2D(unsigned int IDg, bool Internal, unsigned int ID_V1, unsigned int ID_V2);
        Boundaries2D(unsigned int IDg, bool Internal, unsigned int ID_V1, unsigned int ID_V2, int ID_El_L, int ID_El_R);
        unsigned int getID() const ;
        int getLeftElementID();
        int getRightElementID();
        //void set_LeftElementID(int E_L);
        //void set_RightElementID(int E_R);
        bool isInternal();
        //double getxCoordinate();
        //bool operator<(const Boundaries2D& lhs, const Boundaries2D& rhs) {return lhs.ID<rhs.ID;}
        bool operator < (const Boundaries2D& other) const {return ID < other.ID; }

        ~Boundaries2D();

};
#endif
