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
#endif
