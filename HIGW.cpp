#include "HIGW.hpp"
/*--------------------------------------------------------------------------*/
void create_Matrices(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &NMat, Mat &NDerivMat)
{
    double Np = (N+1)*(N+1);
    //Mat Ex, ExT, Ey, EyT;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, 5*Np, NULL, &E);  // 2*N_Nodes x N_Nodes //number of possible nonzero blocks are 5: element and his 4 neighbours (2D)
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, 2*N_Nodes, 2*5*Np, NULL, &ET); // N_Nodes x 2*N_Nodes
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 5*Np, NULL, &Ex);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 5*Np, NULL, &ExT);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 5*Np, NULL, &Ey);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 5*Np, NULL, &EyT);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &invM);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &invM_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &M1);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M1_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NMat);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M2);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NDerivMat);

    std::cout << "Start Elemental Calculations " << std::endl;
    for (auto e = List_Of_Elements.begin(); e < List_Of_Elements.end(); e++)
    {
        unsigned int Np = (*e).get_Number_Of_Nodes();
        double J = (*e).getJacobian();
        double drdx = (*e).get_rx();
        double drdy = (*e).get_ry();
        double dsdx = (*e).get_sx();
        double dsdy = (*e).get_sy();
        unsigned int pos = (*e).getPosition();

        unsigned int Order_Polynomials = (*e).get_Order_Of_Polynomials();

        unsigned int Order_Gaussian_Quadrature = 2*Order_Polynomials+3+N_Q;//ceil(Order_Polynomials+3+N_Q); // + higher order for rho_0 term
        Order_Gaussian_Quadrature = 20  ;// 28;//std::max((uint)10, Order_Gaussian_Quadrature);

        PetscInt in[Np];
        for (unsigned int n=0;n<Np; n++)
        {
            in[n] = n+pos;
        }

        Mat M_Elemental;
        Mat invM_Elemental;
        MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &M_Elemental);

           // Quadrature: Vec JacobiGL_withWeights(const double &alpha, const double &beta, const unsigned int &N, Vec &Weights)
        //Vec R, S;
        //NodesSquares2D(unsigned int N, Vec &XX, Vec &YY)

        Vec ri;
        ri = JacobiGL(0, 0, Order_Polynomials);
        Vec Weights;
        Vec QuadraturePoints;
        QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
          //VecView(QuadraturePoints, PETSC_VIEWER_STDOUT_SELF);
        PetscScalar *w, *qp;
        VecGetArray(Weights, &w);
        VecGetArray(QuadraturePoints, &qp);

        //double diff = 0;
        //double diff2 = 0;
        //std::cout << "J = " << J << ", drdx = " << drdx << ", dsdy = " << dsdy << std::endl;
        for (unsigned int k = 1; k <= Np; k++)
        {
            unsigned int alpha = (k-1)%(N+1);
            unsigned int beta = (k-1)/(N+1);
            //std::cout << "k = " << k << "=> i = " << alpha << ", j = " << beta << std::endl;
            for (unsigned int l = 1; l <= Np; l++)
            {
                unsigned int gamma = (l-1)%(N+1);
                unsigned int delta = (l-1)/(N+1);
                double value_ex = 0.0;
                double value_ey = 0.0;
                double value_m = 0.0;
                double value_n = 0.0;
                double value_m2 = 0.0;
                double value_n_deriv = 0.0;
                for (unsigned int p = 0; p <= Order_Gaussian_Quadrature; p++)
                {
                    double L_alpha = LagrangePolynomial(ri, qp[p], alpha);
                    double L_gamma = LagrangePolynomial(ri, qp[p], gamma);
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double L_beta = LagrangePolynomial(ri, qp[q], beta);
                        double L_delta = LagrangePolynomial(ri, qp[q], delta);
                        //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dr
                        value_m += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J;
                        if (Np > 1)
                        {
                            double dL_gamma = LagrangePolynomialDeriv(ri, qp[p], gamma);
                            double dL_delta = LagrangePolynomialDeriv(ri, qp[q], delta);
                            // w_q drho_0(r_q)/dr l_i(r_q) l_j(r_q)
                            value_ex += w[p]*L_alpha*L_beta * w[q] * (drdx * dL_gamma * L_delta + dsdx * L_gamma* dL_delta) * J;
                            value_ey += w[p]*L_alpha*L_beta * w[q] * (drdy * dL_gamma * L_delta + dsdy * L_gamma* dL_delta) * J;
                        }
                    }
                }
                MatSetValue(M_Elemental, (k-1), (l-1), value_m, ADD_VALUES);
                MatSetValue(M1_small, pos+(k-1), pos+(l-1), value_m, ADD_VALUES);
                MatSetValue(NMat, pos+(k-1), pos+(l-1), value_n, ADD_VALUES);
                MatSetValue(M2, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(NDerivMat, pos+(k-1), pos+(l-1), value_n_deriv, ADD_VALUES);
                MatSetValue(M1, pos+(k-1), pos+(l-1), value_m, ADD_VALUES);
                MatSetValue(M1, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m, ADD_VALUES);
                double factor = -1.0;
                MatSetValue(E, pos+(k-1), pos+(l-1), factor*value_ex, ADD_VALUES);
                MatSetValue(E, N_Nodes+pos+(k-1), pos+(l-1), factor*value_ey, ADD_VALUES);
                value_ex = - value_ex;
                value_ey = - value_ey;
                MatSetValue(ET, pos+(l-1), pos+(k-1), factor*value_ex, ADD_VALUES);
                MatSetValue(ET, pos+(l-1), N_Nodes+pos+(k-1), factor*value_ey, ADD_VALUES);
            }
        }
        VecRestoreArray(Weights, &w);
        VecRestoreArray(QuadraturePoints, &qp);
        VecDestroy(&Weights);
        VecDestroy(&QuadraturePoints);
        VecDestroy(&ri);

        // Construction of Inverse Mass Matrix
        MatAssemblyBegin(M_Elemental, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(M_Elemental, MAT_FINAL_ASSEMBLY);
        invM_Elemental = Inverse_Matrix(M_Elemental);
        MatDestroy(&M_Elemental);
        for (PetscInt i=0; i<Np; i++)
        {
            PetscInt ncols;
            const PetscInt *cols;
            const PetscScalar *vals_invM;
            MatGetRow(invM_Elemental, i, &ncols, &cols, &vals_invM);
            const PetscInt IndexI = pos+i;
            PetscInt GlobalIndexCol[ncols];
            PetscInt GlobalIndexCol2[ncols];
            for (unsigned int j=0; j < ncols; j++)
            {
                GlobalIndexCol[j] = cols[j]+pos;
                GlobalIndexCol2[j] = N_Nodes+cols[j]+pos;
            }
            PetscInt GlobalIndex[1] = {i + pos};
            PetscInt GlobalIndex2[1] = {N_Nodes+i + pos};

            MatSetValues(invM_small,1, GlobalIndex, ncols, GlobalIndexCol,vals_invM,ADD_VALUES);
            MatSetValues(invM,1, GlobalIndex, ncols, GlobalIndexCol,vals_invM,ADD_VALUES);
            MatSetValues(invM,1, GlobalIndex2, ncols, GlobalIndexCol2,vals_invM,ADD_VALUES);

            MatRestoreRow(invM_Elemental, i, &ncols, &cols, &vals_invM);
        }
        MatDestroy(&invM_Elemental);

                    //auto t0 = std::chrono::high_resolution_clock::now();
                    //double Lj = LagrangePolynomial(ri, qp[q], j);
                    //double Li2 = LagrangePolynomial_Test(ri, qp[q], i);
                    //std::cout << "i = " << i << ", q = " << q << ", qp[q] = " << qp[q] << ", Li = " << Li << ", Li2 = " << Li2 << ", diff = " << Li2-Li << std::endl;
                    //diff += abs(Li2-Li);

                    /// Order 5, 9, 13, 17 give a nonzero difference
                    /// Independent of Order_Gaussian_Quadrature
                    //if (Np > 1)
                    //{
                        //value_ex += cubW_a[i]*Lk[i]*dLdxl[i]*Rho0[i]*J;
                        //value_ey += cubW_a[i]*Lk[i]*dLdyl[i]*Rho0[i]*J;
                    //}


                    //double dLi = LagrangePolynomialDeriv(ri, qp[q], i);
                    //double dLi2 = LagrangePolynomialDeriv_Test(ri, qp[q], i);
                    //std::cout << "i = " << i << ", q = " << q << ", qp[q] = " << qp[q] << ", dLi = " << dLi << ", dLi2 = " << dLi2 << ", diff = " << dLi2-dLi << std::endl;
                    //diff2 += abs(dLi2-dLi);

        //std::cout << "ri = " << std::endl;
        //VecView(ri, PETSC_VIEWER_STDOUT_SELF);

    }
    MatAssemblyBegin(M1, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(M2, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M2, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(NMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(NMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(NDerivMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(NDerivMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(M1_small, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1_small, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(invM, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(invM, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(invM_small, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(invM_small, MAT_FINAL_ASSEMBLY);

    for (auto f = List_Of_Boundaries.begin(); f < List_Of_Boundaries.end(); f++)
    {
            //double Jacobian = (*f).getJacobian();
            unsigned int Type_Boundary_Left = (*f).get_Type_Left();
            unsigned int Type_Boundary_Right = (*f).get_Type_Right();
            double theta = (*f).get_theta();
            double nx = (*f).get_nx();
            double ny = 0;//(*f).get_ny();

            int left = (*f).getLeftElementID();
            int right = (*f).getRightElementID();

            unsigned int Np_left = List_Of_Elements[left].get_Number_Of_Nodes();
            unsigned int Np_right = List_Of_Elements[right].get_Number_Of_Nodes();

            unsigned int posL = List_Of_Elements[left].getPosition();
            unsigned int posR = List_Of_Elements[right].getPosition();

            unsigned int Order_Polynomials_left = (List_Of_Elements[left]).get_Order_Of_Polynomials();
            unsigned int Order_Polynomials_right = (List_Of_Elements[right]).get_Order_Of_Polynomials();

            std::vector<unsigned int> Node_Numbers_On_Boundary_Left = List_Of_Elements[left].get_nodes_on_boundary(Type_Boundary_Left);
            std::vector<unsigned int> Node_Numbers_On_Boundary_Right = List_Of_Elements[right].get_nodes_on_boundary(Type_Boundary_Right);


            double Jacobian = List_Of_Elements[left].getJacobian();


            //std::cout << " Type_Boundary_Left = " << Type_Boundary_Left << std::endl;
            //std::cout << " Type_Boundary_Right = " << Type_Boundary_Right << std::endl;

            //if (Type_Boundary_Left == 4)
            //{
            //std::cout << "Node Numbers left: " << std::endl;
            //for (auto i = Node_Numbers_On_Boundary_Left.begin(); i < Node_Numbers_On_Boundary_Left.end(); i++)
            //{
            //    std::cout << (*i) << std::endl;
            //}
            //}
            //if (Type_Boundary_Left == 4)
            //{
            //std::cout << "Node Numbers right: " << std::endl;
            //for (auto i = Node_Numbers_On_Boundary_Right.begin(); i < Node_Numbers_On_Boundary_Right.end(); i++)
            //{
            //    std::cout << (*i) << std::endl;
            //}
            //}
            /// Reverse Nodes on Faces 3 & 4?

            /// Or use two different gaussian quadratures
            // unsigned int Order_Gaussian_Quadrature_L
            // unsigned int Order_Gaussian_Quadrature_R
            unsigned int Order_Gaussian_Quadrature  = ceil(std::max(2*Order_Polynomials_left, 2*Order_Polynomials_right)+3+N_Q);
            Order_Gaussian_Quadrature = 10;//std::max((uint)10, Order_Gaussian_Quadrature);
            // Order_Gaussian_Quadrature+1 = Number of Points

            Vec Weights;
            Vec QuadraturePoints;
            QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
            PetscScalar *w_a, *r_a;
            VecGetArray(QuadraturePoints, &r_a);
            VecGetArray(Weights, &w_a);

            Vec ri_left, ri_right;
            /// Does this not determine the ordering?
            ri_left = JacobiGL(0, 0, Order_Polynomials_left);
            ri_right = JacobiGL(0, 0, Order_Polynomials_right);

            double factor = 1.0;

            //std::cout << "Jacobian = " << Jacobian << std::endl;
            // GLL
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++) // Should be N_Left + 1, etc.
            {
                for (unsigned int j = 0; j <= Order_Polynomials_left; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_L
                    {
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += (1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor;

                    }
                    //std::cout << "value_e = " << value_e << std::endl;
                    //if (abs(ny) > 1e-5)
                    {
                    //std::cout << "Boundary ID = " << (*f).getID() << " Boundary Types: " << Type_Boundary_Left << ", " << Type_Boundary_Right << ". Left El = " << left << ", Right El = " << right << std::endl;
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Left[j] << std::endl;
                    }
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    //}
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);
                }
            }
            // GLR
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -(1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor;
                    }
                    //if (abs(nx) > 1e-5)
                    //{
                    //std::cout << "Boundary ID = " << (*f).getID() << " Boundary Types: " << Type_Boundary_Left << ", " << Type_Boundary_Right << ". Left El = " << left << ", Right El = " << right << std::endl;
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ". posL = " << posL << ", posR = " << posL <<  ". Node Numbers =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    //}
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    //std::cout << posL+Node_Numbers_On_Boundary_Left[i] << " " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                }
            }
            // GRL
            for (unsigned int i = 0; i <= Order_Polynomials_right; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_left; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += theta*w_a[q]*Li*Lj*Jacobian*factor;
                    }
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);

                    //std::cout << posR+Node_Numbers_On_Boundary_Right[i] << " " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                }
            }
            // GRR
            for (unsigned int i = 0; i <= Order_Polynomials_right; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_R
                    {
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -theta*w_a[q]*Li*Lj*Jacobian*factor;
                    }
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);

                    //std::cout << posR+Node_Numbers_On_Boundary_Right[i] << " " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                }
            }


            VecDestroy(&ri_left);
            VecDestroy(&ri_right);
            VecRestoreArray(QuadraturePoints, &r_a);
            VecRestoreArray(Weights, &w_a);
            VecDestroy(&Weights);
            VecDestroy(&QuadraturePoints);
    }
    MatAssemblyBegin(E, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(E, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY);

}
/*--------------------------------------------------------------------------*/
extern void create_Compressible_System_MidPoint(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, Mat &A, Mat &B)
{
    double Np = (N+1)*(N+1);
    std::cout << "Start Global Matrices Construction" << std::endl;

    double fillBF = 1;

    Mat BF1, BF1_TEMP1, BF1_TEMP2;
    MatMatMult(E, invM_small, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP1);
    MatMatMult(BF1_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP2);
    MatMatMult(invM, BF1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF1);
    MatDestroy(&BF1_TEMP1);
    MatDestroy(&BF1_TEMP2);
/*
    Mat BF2, BF2_TEMP1, BF2_TEMP2;
    MatMatMult(E, invM, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP1);
    MatMatMult(BF2_TEMP1, M2, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP2);
    MatMatMult(invM, BF2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF2);
    MatDestroy(&BF2_TEMP1);
    MatDestroy(&BF2_TEMP2);
*/
    Mat DIV, DIV_TEMP1, DIV_TEMP2;
    MatMatMult(ET, invM, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP1);
    MatMatMult(DIV_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP2);
    MatMatMult(invM_small, DIV_TEMP2, MAT_INITIAL_MATRIX, fillBF, &DIV);
    MatDestroy(&DIV_TEMP1);
    MatDestroy(&DIV_TEMP2);
/*
    Mat C, C_TEMP1, C_TEMP2;
    MatMatMult(NMat, invM, MAT_INITIAL_MATRIX, fillBF, &C_TEMP1);
    MatMatMult(C_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &C_TEMP2);
    MatMatMult(invM, C_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C);
    MatDestroy(&C_TEMP1);
    MatDestroy(&C_TEMP2);

    Mat C2, C2_TEMP1, C2_TEMP2;
    MatMatMult(NMat, invM, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP1);
    MatMatMult(C2_TEMP1, M2, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP2);
    MatMatMult(invM, C2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C2);
    MatDestroy(&C2_TEMP1);
    MatDestroy(&C2_TEMP2);

    Mat D1, D1_TEMP1, D1_TEMP2;
    MatMatMult(NDerivMat, invM, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP1);
    MatMatMult(D1_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP2);
    MatMatMult(invM, D1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D1);
    MatDestroy(&D1_TEMP1);
    MatDestroy(&D1_TEMP2);

    Mat D2, D2_TEMP1, D2_TEMP2;
    MatMatMult(NDerivMat, invM, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP1);
    MatMatMult(D2_TEMP1, M2, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP2);
    MatMatMult(invM, D2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D2);
    MatDestroy(&D2_TEMP1);
    MatDestroy(&D2_TEMP2);

    MatDestroy(&E);
    MatDestroy(&ET);
    MatDestroy(&invM);
    MatDestroy(&NMat);
    MatDestroy(&NDerivMat);

    Mat Laplacian;
	MatMatMult(BF1, DIV, MAT_INITIAL_MATRIX, 1, &Laplacian);
	MatScale(Laplacian, nu);

	//MatView(Laplacian, viewer_info);
    */

    //std::cout << "DIV = " << std::endl;
    //MatView(DIV, PETSC_VIEWER_STDOUT_SELF);

       // factor 4 = Number of Variables
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+3*Np*Np,  NULL, &A); // Check Factor 10 // Change from Triangles to Quadrilaterals
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+3*Np*Np,  NULL, &B);
    std::cout << " Global Matrices Preallocated" << std::endl;
    for (unsigned int i = 0; i < 2*N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        MatGetRow(BF1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            //if (dummy!=0.0)
            if (abs(dummy)>1e-10)
            {
                MatSetValue(A, 	i,     3*N_Nodes+cols[j],    -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     3*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	    ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);
    }
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;

        // Fill Diagonals
        MatSetValue(A, i, i, 1.0, ADD_VALUES);
        MatSetValue(A, N_Nodes+i, N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(A, 2*N_Nodes+i, 2*N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(A, 3*N_Nodes+i, 3*N_Nodes+i, 1.0, ADD_VALUES);

        MatSetValue(B, i, i, 1.0, ADD_VALUES);
        MatSetValue(B, N_Nodes+i, N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(B, 2*N_Nodes+i, 2*N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(B, 3*N_Nodes+i, 3*N_Nodes+i, 1.0, ADD_VALUES);

        /*
        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        */

        /*
        MatGetRow(BF2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     2*Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     2*Np*Number_Of_Elements+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	i,     Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     Np*Number_Of_Elements+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF2, i, &numberOfNonZeros, &cols, &values);
        * /
        MatGetRow(C, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     2*Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     2*Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	2*Np*Number_Of_Elements+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	2*Np*Number_Of_Elements+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(C2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	i,     2*Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     2*Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C2, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(D1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	Np*Number_Of_Elements+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	Np*Number_Of_Elements+i,     cols[j],    -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D1, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	i,     2*Np*Number_Of_Elements+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     2*Np*Number_Of_Elements+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);
        */
        MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            //if (dummy!=0.0)
            if (abs(dummy)>1e-10)
            {
                MatSetValue(A, 	3*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                //MatSetValue(A, 	Np*Number_Of_Elements+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	Np*Number_Of_Elements+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);

    }
    MatDestroy(&BF1);
    MatDestroy(&DIV);

    std::cout << " Global Matrices Constructed" << std::endl;
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
}
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecU, Vec &VecW, Vec &VecR, Vec &VecP)
{
    PetscScalar   sigma;
    sigma = calculate_sigma_2D_system1(rho_0_Deriv, kxmode, kzmode);
    PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);
    std::cout << "Computing Initial Condition " << std::endl;
    // Initial Condition
    // Size = sum_i Np_i, i = 1 .. Nel
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecU);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecW);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecR);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecP);

    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int pos = (*k).getPosition();
        std::vector<double> xCoor, zCoor;
        xCoor = (*k).get_node_coordinates_x();
        zCoor = (*k).get_node_coordinates_y();
        double t = 0;
        for (unsigned int n = 0; n < Np; n++)
        {
            double value = Exact_Solution_mx_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecU, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, pos + n, value, INSERT_VALUES);

            value = Exact_Solution_mz_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecW, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, N_Nodes+ pos + n, value, INSERT_VALUES);

            value = Exact_Solution_r_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecR, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, 2*N_Nodes+ pos + n, value, INSERT_VALUES);

            value = Exact_Solution_p_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecP, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, 3*N_Nodes+ pos + n, value, INSERT_VALUES);
            //std::cout << "n = " << n << ", pos = " << pos << ", xCoor[n] = " << xCoor[n] << ", zCoor[n] = " << zCoor[n] << ", p = " << value << std::endl;
        }

    }
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    VecAssemblyBegin(VecU);
    VecAssemblyEnd(VecU);
    VecAssemblyBegin(VecW);
    VecAssemblyEnd(VecW);
    VecAssemblyBegin(VecR);
    VecAssemblyEnd(VecR);
    VecAssemblyBegin(VecP);
    VecAssemblyEnd(VecP);
}
/*--------------------------------------------------------------------------*/
extern void Simulate(const Mat &A, const Mat &B, const Mat &M1_small, const Vec &Initial_Condition, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol)
{

    double H0 = calculate_Hamiltonian2D(M1_small, Initial_Condition, List_Of_Elements, N_Nodes);
    std::cout << "Initial Energy = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Start Simulations " << std::endl;

    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    KSPSetUp(ksp);
    KSPSetTolerances(ksp, 1e-12, 1e-12, 1e30, PETSC_DEFAULT);

    KSPSetType(ksp,KSPPREONLY);
    //KSPSetType(ksp,KSPCG);
    //KSPSetType(ksp,KSPGMRES);
    //KSPSetType(ksp,KSPBCGS);
    //KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    //PCSetType(pc,PCILU);
    //PCSetType(pc,PCNONE);
    //PCSetType(pc,PCSOR);

    KSPSetFromOptions(ksp);

    Vec QX;
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes, &Sol);
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes, &QX);
    VecCopy(Initial_Condition, Sol);
    double H1 = 0.0;

    //PetscPrintf(PETSC_COMM_SELF,"Size Global Matrices %6.4e\n",(double)sigma);
    // MatView(A, viewer_info);

    char szFileName[255] = {0};
    FILE *f = fopen("Energy.txt", "w");
    // Solve Linear System
    std::cout << "Start Time Stepping" << std::endl;
    double time = 0.0;

    double Hold = H0;
    for (unsigned int t = 0; t < Number_Of_TimeSteps; t++) //
    {

        fprintf(f, "%1.16e \t %1.16e\n", time, H1);
        time = (t+1)*DeltaT;
            MatMult(B, Sol, QX);
            KSPSolve(ksp, QX, Sol);

            H1 = calculate_Hamiltonian2D(M1_small, Sol, List_Of_Elements, N_Nodes);

            //std::cout << "Energy Diff= " << std::setprecision(16) << H1-Hold <<std::endl;
            Hold = H1;

           // std::cout << "QX = " << std::endl;
            //VecView(QX, PETSC_VIEWER_STDOUT_SELF);
            //std::cout << "Solution = " << std::endl;
            //VecView(Sol, PETSC_VIEWER_STDOUT_SELF);
        /*
        PetscViewer viewer2;
        sprintf(szFileName, "solution%d.txt", t);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
        VecView(Sol, viewer2);
        PetscViewerDestroy(&viewer2);
        */

    }
    fprintf(f, "%1.16e \t %1.16e\n", time, H1);
    fclose(f);
    std::cout << "Final Time = " << time << std::endl;
    std::cout << "End Time Stepping" << std::endl;
    KSPDestroy(&ksp);
    VecDestroy(&QX);
    std::cout << "Initial Energy    = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Final Energy      = " << std::setprecision(16) << H1 << std::endl;
    std::cout << "Difference Energy = " << std::setprecision(16) << H1-H0 << std::endl;
}
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix_local(const Mat &V)
{
    Mat Product;
    // Calculate inverse Vandermonde Matrix
    Mat A, B, X;
    MatDuplicate(V,MAT_COPY_VALUES,&A);
    MatDuplicate(V,MAT_DO_NOT_COPY_VALUES,&X);
    MatDuplicate(V,MAT_DO_NOT_COPY_VALUES,&B);
    MatConvert(B, MATSEQDENSE, MAT_INPLACE_MATRIX, &B);
    MatConvert(X, MATSEQDENSE, MAT_INPLACE_MATRIX, &X);

    MatShift(B, 1.0);
    MatOrderingType rtype = MATORDERINGNATURAL;
    IS row, col;
    MatGetOrdering(A, rtype, &row, &col);
    MatFactorInfo info;
    MatFactorInfoInitialize(&info);
    info.fill=1.0;
    MatLUFactor(A, row, col, &info);
    MatMatSolve(A, B, X);
    //MatView(X, PETSC_VIEWER_STDOUT_SELF);
    // X is the inverse

    //M = inv(V^T*V)
    MatTransposeMatMult(X, X, MAT_INITIAL_MATRIX, 1.0, &Product);
    MatDestroy(&A);
    MatDestroy(&B);
    MatDestroy(&X);
    ISDestroy(&row);
    ISDestroy(&col);

    return Product;
}
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix_inverse_local(const Mat &V)
{
    Mat M;
    MatMatTransposeMult(V, V, MAT_INITIAL_MATRIX, 1, &M);
    //MatScale(M, 2.0/deltaX);
    MatConvert(M, MATSEQDENSE, MAT_INPLACE_MATRIX, &M);
    return M;
}
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian2D(const Mat &M1, const Vec &Solution, const std::vector<Elements2D> &List_Of_Elements2D, const unsigned int &N_Nodes)
{
    double H = 0.0;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (auto k = List_Of_Elements2D.begin(); k < List_Of_Elements2D.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int pos = (*k).getPosition();

        // Get Submatrix from M1
        Mat M1_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, pos, 1, &isrow);
        MatGetSubMatrix(M1,isrow, isrow,MAT_INITIAL_MATRIX, &M1_Elemental);

        for (unsigned int i = 0; i < Np; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                double massij;
                const PetscInt idxm[1] = {i};
                const PetscInt idxn[1] = {j};

                MatGetValues(M1_Elemental, 1, idxm, 1,  idxn, &massij);
                H += 0.5*massij*(XTemp[pos+i]*XTemp[pos+j]+XTemp[N_Nodes+pos+i]*XTemp[N_Nodes+pos+j]+XTemp[3*N_Nodes+pos+i]*XTemp[3*N_Nodes+pos+j]);

            }
        }

        MatDestroy(&M1_Elemental);
        ISDestroy(&isrow);
    }
    VecRestoreArray(Solution, &XTemp);
    return H;
}
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian2D(const Mat &M1, const Vec &Solution, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes)
{
    double H = 0.0;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int pos = (*k).getPosition();

        // Get Submatrix from M1
        Mat M1_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, pos, 1, &isrow);
        MatGetSubMatrix(M1,isrow, isrow,MAT_INITIAL_MATRIX, &M1_Elemental);

        for (unsigned int i = 0; i < Np; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                double massij;
                const PetscInt idxm[1] = {i};
                const PetscInt idxn[1] = {j};

                MatGetValues(M1_Elemental, 1, idxm, 1,  idxn, &massij);
                H += 0.5*massij*(XTemp[pos+i]*XTemp[pos+j]+XTemp[N_Nodes+pos+i]*XTemp[N_Nodes+pos+j]+XTemp[3*N_Nodes+pos+i]*XTemp[3*N_Nodes+pos+j]);

            }
        }

        MatDestroy(&M1_Elemental);
        ISDestroy(&isrow);
    }
    VecRestoreArray(Solution, &XTemp);
    return H;
}
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian(const Mat &M1, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np)
{

    double H = 0.0;
    //PetscInt size_s = Number_Of_Elements*Np;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (unsigned int k = 0; k < Number_Of_Elements; k++)
    {
        // Get Submatrix from M1
        Mat M1_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, k*Np, 1, &isrow);
        MatGetSubMatrix(M1,isrow, isrow,MAT_INITIAL_MATRIX, &M1_Elemental);
        for (unsigned int i = 0; i < Np; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                double massij;
                const PetscInt idxm[1] = {i};
                const PetscInt idxn[1] = {j};

                MatGetValues(M1_Elemental, 1, idxm, 1,  idxn, &massij);
                //massij = 1;
                H += 0.5*massij*(XTemp[k*Np+i]*XTemp[k*Np+j]+XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[Np*Number_Of_Elements+k*Np+j]);
            }
        }
        MatDestroy(&M1_Elemental);
        ISDestroy(&isrow);
    }

    VecRestoreArray(Solution, &XTemp);
    return H;

}
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian_pres(const Mat &M1, const Mat &M2, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np)
{

    double H = 0.0;
    //PetscInt size_s = Number_Of_Elements*Np;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (unsigned int k = 0; k < Number_Of_Elements; k++)
    {
        // Get Submatrix from M1 and M2
        Mat M1_Elemental, M2_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, k*Np, 1, &isrow);
        MatGetSubMatrix(M1,isrow, isrow,MAT_INITIAL_MATRIX, &M1_Elemental);
        MatGetSubMatrix(M2,isrow, isrow,MAT_INITIAL_MATRIX, &M2_Elemental);

        for (unsigned int i = 0; i < Np; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                double massij1, massij2;
                const PetscInt idxm[1] = {i};
                const PetscInt idxn[1] = {j};

                MatGetValues(M1_Elemental, 1, idxm, 1,  idxn, &massij1);
                MatGetValues(M2_Elemental, 1, idxm, 1,  idxn, &massij2);
                H += 0.5*massij1*(XTemp[k*Np+i]*XTemp[k*Np+j]+XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[Np*Number_Of_Elements+k*Np+j])+0.5*massij2*(XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[Np*Number_Of_Elements+k*Np+j]);
            }
        }
        MatDestroy(&M1_Elemental);
        MatDestroy(&M2_Elemental);
        ISDestroy(&isrow);
    }

    VecRestoreArray(Solution, &XTemp);
    return H;

}
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian_comp(const Mat &M1, const Mat &M2, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np)
{

    double H = 0.0;
    //PetscInt size_s = Number_Of_Elements*Np;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (unsigned int k = 0; k < Number_Of_Elements; k++)
    {
        // Get Submatrix from M1 and M2
        Mat M1_Elemental, M2_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, k*Np, 1, &isrow);
        MatGetSubMatrix(M1,isrow, isrow,MAT_INITIAL_MATRIX, &M1_Elemental);
        MatGetSubMatrix(M2,isrow, isrow,MAT_INITIAL_MATRIX, &M2_Elemental);

        for (unsigned int i = 0; i < Np; i++)
        {
            for (unsigned int j = 0; j < Np; j++)
            {
                double massij1, massij2;
                const PetscInt idxm[1] = {i};
                const PetscInt idxn[1] = {j};

                MatGetValues(M1_Elemental, 1, idxm, 1,  idxn, &massij1);
                MatGetValues(M2_Elemental, 1, idxm, 1,  idxn, &massij2);
                H += 0.5*massij1*(XTemp[k*Np+i]*XTemp[k*Np+j]+XTemp[2*Np*Number_Of_Elements+k*Np+i]*XTemp[2*Np*Number_Of_Elements+k*Np+j])
                    +0.5*massij2*(XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[Np*Number_Of_Elements+k*Np+j]-XTemp[Np*Number_Of_Elements+k*Np+i]*XTemp[2*Np*Number_Of_Elements+k*Np+j]
                        -XTemp[Np*Number_Of_Elements+k*Np+j]*XTemp[2*Np*Number_Of_Elements+k*Np+i]+XTemp[2*Np*Number_Of_Elements+k*Np+i]*XTemp[2*Np*Number_Of_Elements+k*Np+j]);
            }
        }
        MatDestroy(&M1_Elemental);
        MatDestroy(&M2_Elemental);
        ISDestroy(&isrow);
    }

    VecRestoreArray(Solution, &XTemp);
    return H;

}
/*--------------------------------------------------------------------------*/
extern double calculate_Error(const Vec &Exact, const Vec &Solution, const unsigned int &Number_Of_Elements, const unsigned int &Np, const double &DeltaX)
{
    PetscInt size_r;
    VecGetSize(Solution, &size_r);
    PetscScalar *ExactTemp, *SolutionTemp;
    VecGetArray(Exact, &ExactTemp);
    VecGetArray(Solution, &SolutionTemp);

    double error = 0.0;
    for (unsigned int i=0; i < size_r; i++)
    {
        if (Np > 1 && i > 0 && i < size_r-1 && (i+1)%(Np*Number_Of_Elements)!=0   && (i+1)%Np == 0 && i%Np==Np-1)
        {
            // Average over Interface
            SolutionTemp[i+1] = 0.5*(SolutionTemp[i]+SolutionTemp[i+1]);
        }
        else
        {
            error += (SolutionTemp[i]-ExactTemp[i])*(SolutionTemp[i]-ExactTemp[i]);
        }
    }
    error = sqrt(error);
    error *= sqrt(DeltaX);


    VecRestoreArray(Exact, &ExactTemp);
    VecRestoreArray(Solution, &SolutionTemp);

    return error;
}
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const double &DeltaX, const double &DeltaY, const unsigned int &Np)
{
    Vec Difference;
    VecDuplicate(Exact, &Difference);

    VecWAXPY(Difference, -1.0, Exact, Solution);

    PetscReal error;
    if (Norm_Type == 1)
    {
        VecNorm(Difference, NORM_1, &error);
    }
    else if (Norm_Type == 2)
    {
        VecNorm(Difference, NORM_2, &error);
    }
    else
    {
        VecNorm(Difference, NORM_INFINITY, &error);
    }

    VecDestroy(&Difference);
    error /= sqrt(Np);
    error *= sqrt(0.5*DeltaX*DeltaY);
    return error;

}
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes)
{
    std::cout << "Start Elemental Calculations " << std::endl;
    PetscScalar *exact_a, *sol_a;
    VecGetArray(Exact, &exact_a);
    VecGetArray(Solution, &sol_a);
    double error_u = 0, error_w = 0, error_r = 0, error_p = 0;

    //std::cout << "Exact = " << std::endl;
    //VecView(Exact , PETSC_VIEWER_STDOUT_SELF);
    //std::cout << "Sol = " << std::endl;
    //VecView(Solution , PETSC_VIEWER_STDOUT_SELF);
    for (auto e = List_Of_Elements.begin(); e < List_Of_Elements.end(); e++)
    {
        unsigned int Np = (*e).get_Number_Of_Nodes();
        double J = (*e).getJacobian();
        double drdx = (*e).get_rx();
        double drdy = (*e).get_ry();
        double dsdx = (*e).get_sx();
        double dsdy = (*e).get_sy();
        unsigned int pos = (*e).getPosition();

        unsigned int N = (*e).get_Order_Of_Polynomials();
        std::vector<double> xCoor, zCoor;
        xCoor = (*e).get_node_coordinates_x();
        zCoor = (*e).get_node_coordinates_y();

        Vec Weights;
        Vec QuadraturePoints;
        QuadraturePoints = JacobiGL_withWeights(0, 0, N, Weights);

        PetscScalar *w, *qp;
        VecGetArray(Weights, &w);
        VecGetArray(QuadraturePoints, &qp);

        double diff_u = 0;
        double exact_u = 0;
        double diff_w = 0;
        double exact_w = 0;
        double diff_r = 0;
        double exact_r = 0;
        double diff_p = 0;
        double exact_p = 0;

        for (unsigned int i = 0; i < (N+1); i++)
        {
            for (unsigned int j = 0; j < (N+1); j++)
            {
                double ex_u = exact_a[pos+i*(N+1)+j];
                double sol_u = sol_a[pos+i*(N+1)+j];
                double ex_w = exact_a[N_Nodes+pos+i*(N+1)+j];
                double sol_w = sol_a[N_Nodes+pos+i*(N+1)+j];
                double ex_r = exact_a[2*N_Nodes+pos+i*(N+1)+j];
                double sol_r = sol_a[2*N_Nodes+pos+i*(N+1)+j];
                double ex_p = exact_a[3*N_Nodes+pos+i*(N+1)+j];
                double sol_p = sol_a[3*N_Nodes+pos+i*(N+1)+j];

                //std::cout << "n = " << (*e).getID() << ", J =  " << J << ", pos = " << pos << ", xCoor[n] = " << xCoor[i*(N+1)+j] << ", zCoor[n] = " << zCoor[i*(N+1)+j] << ", ex p = " << ex_p << ", sol p = " << sol_p << std::endl;

                exact_u += w[i]*w[j]*ex_u*ex_u*J;
                diff_u += w[i]*w[j]*(ex_u-sol_u)*(ex_u-sol_u)*J;
                exact_w += w[i]*w[j]*ex_w*ex_w*J;
                diff_w += w[i]*w[j]*(ex_w-sol_w)*(ex_w-sol_w)*J;
                exact_r += w[i]*w[j]*ex_r*ex_r*J;
                diff_r += w[i]*w[j]*(ex_r-sol_r)*(ex_r-sol_r)*J;
                exact_p += w[i]*w[j]*ex_p*ex_p*J;
                diff_p += w[i]*w[j]*(ex_p-sol_p)*(ex_p-sol_p)*J;
            }
        }
        /*
        std::cout << "diff_u = " << diff_u << std::endl;
        std::cout << "exact_u = " << exact_u << std::endl;
        std::cout << "diff_w = " << diff_w << std::endl;
        std::cout << "exact_w = " << exact_w << std::endl;
        std::cout << "diff_r = " << diff_r << std::endl;
        std::cout << "exact_r = " << exact_r << std::endl;
        std::cout << "diff_p = " << diff_p << std::endl;
        std::cout << "exact_p = " << exact_p << std::endl;
        */
        //error_u += sqrt(diff_u/exact_u);
        //error_w += sqrt(diff_w/exact_w);
        //error_r += sqrt(diff_r/exact_r);
        //error_p += sqrt(diff_p/exact_p);

        error_u += sqrt(diff_u);
        error_w += sqrt(diff_w);
        error_r += sqrt(diff_r);
        error_p += sqrt(diff_p);


        VecRestoreArray(Weights, &w);
        VecRestoreArray(QuadraturePoints, &qp);
        VecDestroy(&Weights);
        VecDestroy(&QuadraturePoints);
    }
    VecRestoreArray(Exact, &exact_a);
    VecRestoreArray(Solution, &sol_a);
    //std::cout << "Error u = " << error_u << std::endl;
    //std::cout << "Error w = " << error_w << std::endl;
    //std::cout << "Error r = " << error_r << std::endl;
    //std::cout << "Error p = " << error_p << std::endl;

    //return sqrt(error_u+error_w+error_r+error_p);
    return error_p;

}
/*--------------------------------------------------------------------------*/
