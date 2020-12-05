#include "HIGW.hpp"
/*--------------------------------------------------------------------------*/
void printFullMatrixInfo(Mat& matrix, const std::string& name)
{
    PetscInt m=0,n=0;
    MatGetSize(matrix,&m,&n);

    MatInfo info;
    MatGetInfo(matrix,MAT_LOCAL, &info);
    std::cout<<name<<std::endl;
    //printf("N = %d, N = %d, block_size = %d, memory = %d, assemblies = %d, mallocs = %d, matrix nonzeros (SeqBAIJ format) = %d, allocated nonzeros= %d\n", m, n, (int)info.block_size, (int)info.memory, (int)info.assemblies, (int)info.mallocs,(int)info.nz_used,(int)info.nz_allocated);
	printf("N = %d, N = %d, block_size = %d, memory = %u, assemblies = %d, mallocs = %d, matrix nonzeros (SeqBAIJ format) = %d, allocated nonzeros= %d\n", m, n, (int)info.block_size, (unsigned int)info.memory, (int)info.assemblies, (int)info.mallocs,(int)info.nz_used,(int)info.nz_allocated);
}
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
        //std::cout << "J = " << J << ", drdx = " << drdx << ", drdy = " << drdy << ", dsdy = " << dsdy << ", dsdx = " << dsdx  << std::endl;
        for (unsigned int k = 1; k <= Np; k++)
        {
            unsigned int alpha = (k-1)%(N+1);
            unsigned int beta = (k-1)/(N+1);
            for (unsigned int l = 1; l <= Np; l++)
            {
                unsigned int gamma = (l-1)%(N+1);
                unsigned int delta = (l-1)/(N+1);
                double value_ex = 0.0;
                double value_ey = 0.0;
                double value_m = 0.0;
                double value_m1 = 0.0;
                double value_n = 0.0;
                double value_m2 = 0.0;
                double value_n_deriv = 0.0;

                //double sum = 0;
                for (unsigned int p = 0; p <= Order_Gaussian_Quadrature; p++)
                {
                    double L_alpha = LagrangePolynomial(ri, qp[p], alpha);
                    double L_gamma = LagrangePolynomial(ri, qp[p], gamma);
                    //std::cout << "w[p] = " << w[p] << std::endl;
                    //sum += w[p];
                    //std::cout << "sum w[p] = " << sum << std::endl;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double L_beta = LagrangePolynomial(ri, qp[q], beta);
                        double L_delta = LagrangePolynomial(ri, qp[q], delta);
                        //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dr
                        value_m += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J;

                        // value_ey += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0deriv;
                        // value_m1 += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J / rho0 ;
                        if (Np > 1)
                        {
                            double dL_gamma = LagrangePolynomialDeriv(ri, qp[p], gamma);
                            double dL_delta = LagrangePolynomialDeriv(ri, qp[q], delta);
                            //std::cout << dL_gamma << " ";
                            // w_q drho_0(r_q)/dr l_i(r_q) l_j(r_q)
                            value_ex += w[p]*L_alpha*L_beta * w[q] * (drdx * dL_gamma * L_delta + dsdx * L_gamma* dL_delta) * J;
                            //value_ex += w[p]*L_alpha*L_beta * w[q]*dL_gamma*L_delta;
                            value_ey += w[p]*L_alpha*L_beta * w[q] * (drdy * dL_gamma * L_delta + dsdy * L_gamma* dL_delta) * J;
                        }
                    }
                }
                //std::cout << std::endl;
                //std::cout << "k = " << k << "=> i = " << alpha << ", j = " << beta << ", value_ex = " << value_ex << std::endl;
                MatSetValue(M_Elemental, (k-1), (l-1), value_m, ADD_VALUES);
                MatSetValue(M1_small, pos+(k-1), pos+(l-1), value_m, ADD_VALUES);
                MatSetValue(NMat, pos+(k-1), pos+(l-1), value_n, ADD_VALUES);
                MatSetValue(M2, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(NDerivMat, pos+(k-1), pos+(l-1), value_n_deriv, ADD_VALUES);
                MatSetValue(M1, pos+(k-1), pos+(l-1), value_m, ADD_VALUES);
                MatSetValue(M1, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m, ADD_VALUES);
                /****************************************************/
                double factor = -1.0;
                /****************************************************/
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
            double Jacobian = (*f).getJacobian();
            unsigned int Type_Boundary_Left = (*f).get_Type_Left();
            unsigned int Type_Boundary_Right = (*f).get_Type_Right();
            double theta = (*f).get_theta();
            double nx = (*f).get_nx();
            double ny = (*f).get_ny();

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

            //std::cout << "left = " << left;
            //std::cout << ", right = " << right << std::endl;
            //double Jacobian = List_Of_Elements[left].getJacobian();

            //std::cout << "Element Jacobian = " << List_Of_Elements[left].getJacobian() << std::endl;
            //std::cout << "Face Jacobian = " << (*f).getJacobian() << std::endl;
            //std::cout << "Element rx = " << List_Of_Elements[left].get_rx() << std::endl;
            //std::cout << "Element sy = " << List_Of_Elements[left].get_sy() << std::endl;

            //nx = nx / List_Of_Elements[left].get_rx();
            //ny = ny / List_Of_Elements[left].get_sy();


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
            Order_Gaussian_Quadrature = std::max((uint)10, Order_Gaussian_Quadrature);
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
            //std::cout << " *---* " << std::endl;
            //std::cout << "Boundary ID = " << (*f).getID() << " Boundary Types: " << Type_Boundary_Left << ", " << Type_Boundary_Right << ". Left El = " << left << ", Right El = " << right << std::endl;
            //std::cout << " *---* " << std::endl;
            //std::cout << "GLL: " << std::endl;


            /*
        for (int j = 0; j < numberOfDegreesOfFreedom; ++j) {
            bFL = left->basisFunction(j, pL);
            bFR = right->basisFunction(j, pR);

            LinearAlgebra::MiddleSizeMatrix& leftReturnData = ret.left_[j];

            LinearAlgebra::MiddleSizeMatrix& rightReturnData = ret.right_[j];

            for (unsigned int i = 0; i < numberOfDegreesOfFreedom; ++i) {

                BFevalL = left->basisFunction(i, pL);

                BFevalR = right->basisFunction(i, pR);



            }
        }

            for (unsigned int j = 0; j < nb; ++j) {
                for (unsigned int i = 0; i < nb; ++i) {




                leftReturnData(0, i) = BFevalL * nx * (1 - theta) * bFL;
                leftReturnData(3, i) = BFevalR * nx * (theta)*bFL;
                leftReturnData(6, i) = bFL * nx * (1 - theta) * BFevalL;
                leftReturnData(9, i) = bFL * nx * (theta)*BFevalR;
                rightReturnData(0, i) = BFevalL * nx * (1 - theta) * bFR;
                rightReturnData(3, i) = BFevalR * nx * (theta)*bFR;
                rightReturnData(6, i) = bFR * nx * (1 - theta) * BFevalL;
                rightReturnData(9, i) = bFR * nx * (theta)*BFevalR;
                    MatSetValue(BF, posR + j, posL + i, fData.right_[j](0, i), ADD_VALUES);  // U coefficient
                    MatSetValue(BF, posR + j, posR + i, -fData.right_[j](3, i), ADD_VALUES);  // U coefficient
                    MatSetValue(BF, posL + j, posL + i, fData.left_[j](0, i), ADD_VALUES);  // U coefficient
                    MatSetValue(BF, posL + j, posR + i, -fData.left_[j](3, i),  ADD_VALUES);  // U coefficient

                    MatSetValue(globalDIV, posR + j, posL + i, fData.right_[j](6, i),  ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posL + j, posL + i,  -fData.left_[j](6, i), ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posL + j, posR + i,  -fData.left_[j](9, i), ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posR + j, posR + i, fData.right_[j](9, i), ADD_VALUES);  // U coefficient


                            GRL
                    rightReturnData(1, i) = BFevalL * ny * (1 - theta) * bFR;
                    rightReturnData(7, i) = bFR * ny * (1 - theta) * BFevalL;
                    MatSetValue(BF, Nu + posR + j, posL + i, fData.right_[j](1, i), ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posR + j, Nu + posL + i, fData.right_[j](7, i),  ADD_VALUES);  // V coefficient

                            GLR
                leftReturnData(4, i) = BFevalR * ny * (theta)*bFL;
                leftReturnData(10, i) = bFL * ny * (theta)*BFevalR;
                    MatSetValue(BF, Nu + posL + j, posR + i, -fData.left_[j](4, i),  ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posL + j, Nu + posR + i, -fData.left_[j](10, i),  ADD_VALUES);  // V coefficient

                                    GLL
                leftReturnData(1, i) = BFevalL * ny * (1 - theta) * bFL;
                leftReturnData(7, i) = bFL * ny * (1 - theta) * BFevalL;
                    MatSetValue(BF, Nu + posL + j, posL + i, fData.left_[j](1, i), ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posL + j, Nu + posL + i, -fData.left_[j](7, i),  ADD_VALUES);  // V coefficient

                                    GRR
                    rightReturnData(4, i) = BFevalR * ny * (theta)*bFR;
                    rightReturnData(10, i) = bFR * ny * (theta)*BFevalR;
                    MatSetValue(BF, Nu + posR + j, posR + i, -fData.right_[j](4, i), ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posR + j, Nu + posR + i, fData.right_[j](10, i),  ADD_VALUES);  // V coefficient





                leftReturnData(2, i) = BFevalL * nz * (1 - theta) * bFL;
                leftReturnData(5, i) = BFevalR * nz * (theta)*bFL;
                leftReturnData(8, i) = bFL * nz * (1 - theta) * BFevalL;
                leftReturnData(11, i) = bFL * nz * (theta)*BFevalR;
                rightReturnData(2, i) = BFevalL * nz * (1 - theta) * bFR;
                rightReturnData(5, i) = BFevalR * nz * (theta)*bFR;
                rightReturnData(8, i) = bFR * nz * (1 - theta) * BFevalL;
                rightReturnData(11, i) = bFR * nz * (theta)*BFevalR;
                    MatSetValue(BF, Nu + Nv + posL + j, posR + i,  -fData.left_[j](5, i),  ADD_VALUES);  // W coefficient
                    MatSetValue(BF, Nu + Nv + posL + j, posL + i, fData.left_[j](2, i), ADD_VALUES);  // W coefficient
                    MatSetValue(BF, Nu + Nv + posR + j, posR + i,  -fData.right_[j](5, i),  ADD_VALUES);  // W coefficient
                    MatSetValue(BF, Nu + Nv + posR + j, posL + i, fData.right_[j](2, i), ADD_VALUES);  // W coefficient

                    MatSetValue(globalDIV, posL + j, Nu + Nv + posL + i,  -fData.left_[j](8, i),  ADD_VALUES);  // W coefficient
                    MatSetValue(globalDIV, posR + j, Nu + Nv + posL + i, fData.right_[j](8, i),  ADD_VALUES);  // W coefficient
                    MatSetValue(globalDIV, posL + j, Nu + Nv + posR + i, -fData.left_[j](11, i), ADD_VALUES);  // W coefficient
                    MatSetValue(globalDIV, posR + j, Nu + Nv + posR + i, fData.right_[j](11, i), ADD_VALUES);  // W coefficient



                    leftReturnData(0, i) = BFevalL * nx * (1 - theta) * bFL;
                    leftReturnData(6, i) = bFL * nx * (1 - theta) * BFevalL;
                    MatSetValue(BF, posL + j, posL + i, fData.left_[j](0, i), ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posL + j, posL + i, -fData.left_[j](6, i), ADD_VALUES);  // U coefficient

                    rightReturnData(0, i) = BFevalL * nx * (1 - theta) * bFR;
                    rightReturnData(6, i) = bFR * nx * (1 - theta) * BFevalL;
                    MatSetValue(BF, posR + j, posL + i, fData.right_[j](0, i), ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posR + j, posL + i, fData.right_[j](6, i), ADD_VALUES);  // U coefficient

                    leftReturnData(3, i) = BFevalR * nx * (theta)*bFL;
                    leftReturnData(9, i) = bFL * nx * (theta)*BFevalR;
                    MatSetValue(BF, posL + j, posR + i, -fData.left_[j](3, i), ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posL + j, posR + i, -fData.left_[j](9, i),ADD_VALUES);  // U coefficient

                    rightReturnData(3, i) = BFevalR * nx * (theta)*bFR;
                    rightReturnData(9, i) = bFR * nx * (theta)*BFevalR;
                    MatSetValue(BF, posR + j, posR + i, -fData.right_[j](3, i),  ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posR + j, posR + i, fData.right_[j](9, i),ADD_VALUES);  // U coefficient

                }
            }
*/
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++) // Should be N_Left + 1, etc.
            {
                for (unsigned int j = 0; j <= Order_Polynomials_left; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    // GLL
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_L
                    {
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += (1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor;

                    }
                    //std::cout << "value_e = " << value_e << std::endl;
                    //if (abs(ny) > 1e-5)
                    {
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Left[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    }
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], -nx*value_e, ADD_VALUES);
                    //}
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                }
            }
            //std::cout << "GLR: " << std::endl;
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
                        //value_e += (1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor;
                    }
                    //if (abs(nx) > 1e-5)
                    //{
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    //MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[j],  posL+Node_Numbers_On_Boundary_Left[i] , nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i],   posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i],  ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j],  -ny*value_e, ADD_VALUES);

                    //std::cout << posL+Node_Numbers_On_Boundary_Left[i] << " " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                }
            }
/*

                            GRL
                    rightReturnData(1, i) = BFevalL * ny * (1 - theta) * bFR;
                    rightReturnData(7, i) = bFR * ny * (1 - theta) * BFevalL;
                    MatSetValue(BF, Nu + posR + j, posL + i, fData.right_[j](1, i), ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posR + j, Nu + posL + i, fData.right_[j](7, i),  ADD_VALUES);  // V coefficient

                            GLR
                leftReturnData(4, i) = BFevalR * ny * (theta)*bFL;
                leftReturnData(10, i) = bFL * ny * (theta)*BFevalR;
                    MatSetValue(BF, Nu + posL + j, posR + i, -fData.left_[j](4, i),  ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posL + j, Nu + posR + i, -fData.left_[j](10, i),  ADD_VALUES);  // V coefficient
                                */
            //std::cout << "GRL: " << std::endl;
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
                        //value_e += -theta*w_a[q]*Li*Lj*Jacobian*factor;
                    }
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    //MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i],  nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j],   -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i], ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                    //std::cout << posR+Node_Numbers_On_Boundary_Right[i] << " " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                }
            }
            //std::cout << "GRR: " << std::endl;
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
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], -ny*value_e, ADD_VALUES);

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
/*--------------------------------------------------------------------------*//*--------------------------------------------------------------------------*/
void create_Matrices_Quads(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &NMat, Mat &NDerivMat)
{
    double Np = (N+1)*(N+1);
    //Mat Ex, ExT, Ey, EyT;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, 5*Np, NULL, &E);  // 2*N_Nodes x N_Nodes //number of possible nonzero blocks are 5: element and his 4 neighbours (2D)
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, 2*N_Nodes, 2*5*Np, NULL, &ET); // N_Nodes x 2*N_Nodes
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

        Vec ri;
        ri = JacobiGL(0, 0, Order_Polynomials);
        Vec Weights;
        Vec QuadraturePoints;
        QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
        PetscScalar *w, *qp;
        VecGetArray(Weights, &w);
        VecGetArray(QuadraturePoints, &qp);

        for (unsigned int k = 1; k <= Np; k++)
        {
            unsigned int alpha = (k-1)%(N+1);
            unsigned int beta = (k-1)/(N+1);
            for (unsigned int l = 1; l <= Np; l++)
            {
                unsigned int gamma = (l-1)%(N+1);
                unsigned int delta = (l-1)/(N+1);
                double value_ex = 0.0;
                double value_ey = 0.0;
                double value_m = 0.0;
                double value_m1 = 0.0;
                double value_n = 0.0;
                double value_m2 = 0.0;
                double value_n_deriv = 0.0;

                double value_area = 0.0;

                //double sum = 0;
                for (unsigned int p = 0; p <= Order_Gaussian_Quadrature; p++)
                {
                    double L_alpha = LagrangePolynomial(ri, qp[p], alpha);
                    double L_gamma = LagrangePolynomial(ri, qp[p], gamma);
                    //std::cout << "w[p] = " << w[p] << std::endl;
                    //sum += w[p];
                    //std::cout << "sum w[p] = " << sum << std::endl;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {

                        double J, drdx, drdy, dsdx, dsdy, x, y;
                        Calculate_Jacobian_Quadrilateral((*e), List_Of_Vertices, qp[p], qp[q], J, drdx, drdy, dsdx, dsdy, x, y);
                        //std::cout << "J = " << J << ", drdx = " << drdx << ", drdy = " << drdy << ", dsdy = " << dsdy << ", dsdx = " << dsdx  << std::endl;

                        double L_beta = LagrangePolynomial(ri, qp[q], beta);
                        double L_delta = LagrangePolynomial(ri, qp[q], delta);
                        //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dr
                        value_m += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J;
                        value_area += w[p]*w[q]*J;

                        // value_ey += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0deriv;
                        // value_m1 += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J / rho0 ;
                        if (Np > 1)
                        {
                            double dL_gamma = LagrangePolynomialDeriv(ri, qp[p], gamma);
                            double dL_delta = LagrangePolynomialDeriv(ri, qp[q], delta);
                            //std::cout << dL_gamma << " ";
                            // w_q drho_0(r_q)/dr l_i(r_q) l_j(r_q)
                            value_ex += w[p]*L_alpha*L_beta * w[q] * (drdx * dL_gamma * L_delta + dsdx * L_gamma* dL_delta) * J;
                            //value_ex += w[p]*L_alpha*L_beta * w[q]*dL_gamma*L_delta;
                            value_ey += w[p]*L_alpha*L_beta * w[q] * (drdy * dL_gamma * L_delta + dsdy * L_gamma* dL_delta) * J;
                        }
                    }
                }
                //std::cout << std::endl;
                //std::cout << "ID = " << (*e).getID() << ", Area = " << value_area << std::endl;
                //std::cout << std::endl;
                //std::cout << "k = " << k << "=> i = " << alpha << ", j = " << beta << ", value_ex = " << value_ex << std::endl;
                MatSetValue(M_Elemental, (k-1), (l-1), value_m, ADD_VALUES);
                MatSetValue(M1_small, pos+(k-1), pos+(l-1), value_m, ADD_VALUES);
                MatSetValue(NMat, pos+(k-1), pos+(l-1), value_n, ADD_VALUES);
                MatSetValue(M2, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(NDerivMat, pos+(k-1), pos+(l-1), value_n_deriv, ADD_VALUES);
                MatSetValue(M1, pos+(k-1), pos+(l-1), value_m, ADD_VALUES);
                MatSetValue(M1, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m, ADD_VALUES);
                /****************************************************/
                double factor = -1.0;
                /****************************************************/
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
            double Jacobian = (*f).getJacobian();
            unsigned int Type_Boundary_Left = (*f).get_Type_Left();
            unsigned int Type_Boundary_Right = (*f).get_Type_Right();
            double theta = (*f).get_theta();
            double nx = (*f).get_nx();
            double ny = (*f).get_ny();


            //std::cout << "J = " << Jacobian << std::endl;
            //std::cout << "nx = " << nx << std::endl;
            //std::cout << "ny = " << ny << std::endl;
            //std::cout << "nx^2 + ny^2 = " << nx*nx+ny*ny << std::endl;
            //std::cout << std::endl;


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

            //std::cout << "left = " << left;
            //std::cout << ", right = " << right << std::endl;
            //double Jacobian = List_Of_Elements[left].getJacobian();

            //std::cout << "Element Jacobian = " << List_Of_Elements[left].getJacobian() << std::endl;
            //std::cout << "Face Jacobian = " << (*f).getJacobian() << std::endl;
            //std::cout << "Element rx = " << List_Of_Elements[left].get_rx() << std::endl;
            //std::cout << "Element sy = " << List_Of_Elements[left].get_sy() << std::endl;

            //nx = nx / List_Of_Elements[left].get_rx();
            //ny = ny / List_Of_Elements[left].get_sy();


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
            Order_Gaussian_Quadrature = std::max((uint)10, Order_Gaussian_Quadrature);
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
            //std::cout << " *---* " << std::endl;
            //std::cout << "Boundary ID = " << (*f).getID() << " Boundary Types: " << Type_Boundary_Left << ", " << Type_Boundary_Right << ". Left El = " << left << ", Right El = " << right << std::endl;
            //std::cout << " *---* " << std::endl;
            //std::cout << "GLL: " << std::endl;


            /*
        for (int j = 0; j < numberOfDegreesOfFreedom; ++j) {
            bFL = left->basisFunction(j, pL);
            bFR = right->basisFunction(j, pR);

            LinearAlgebra::MiddleSizeMatrix& leftReturnData = ret.left_[j];

            LinearAlgebra::MiddleSizeMatrix& rightReturnData = ret.right_[j];

            for (unsigned int i = 0; i < numberOfDegreesOfFreedom; ++i) {

                BFevalL = left->basisFunction(i, pL);

                BFevalR = right->basisFunction(i, pR);



            }
        }

            for (unsigned int j = 0; j < nb; ++j) {
                for (unsigned int i = 0; i < nb; ++i) {




                leftReturnData(0, i) = BFevalL * nx * (1 - theta) * bFL;
                leftReturnData(3, i) = BFevalR * nx * (theta)*bFL;
                leftReturnData(6, i) = bFL * nx * (1 - theta) * BFevalL;
                leftReturnData(9, i) = bFL * nx * (theta)*BFevalR;
                rightReturnData(0, i) = BFevalL * nx * (1 - theta) * bFR;
                rightReturnData(3, i) = BFevalR * nx * (theta)*bFR;
                rightReturnData(6, i) = bFR * nx * (1 - theta) * BFevalL;
                rightReturnData(9, i) = bFR * nx * (theta)*BFevalR;
                    MatSetValue(BF, posR + j, posL + i, fData.right_[j](0, i), ADD_VALUES);  // U coefficient
                    MatSetValue(BF, posR + j, posR + i, -fData.right_[j](3, i), ADD_VALUES);  // U coefficient
                    MatSetValue(BF, posL + j, posL + i, fData.left_[j](0, i), ADD_VALUES);  // U coefficient
                    MatSetValue(BF, posL + j, posR + i, -fData.left_[j](3, i),  ADD_VALUES);  // U coefficient

                    MatSetValue(globalDIV, posR + j, posL + i, fData.right_[j](6, i),  ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posL + j, posL + i,  -fData.left_[j](6, i), ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posL + j, posR + i,  -fData.left_[j](9, i), ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posR + j, posR + i, fData.right_[j](9, i), ADD_VALUES);  // U coefficient


                            GRL
                    rightReturnData(1, i) = BFevalL * ny * (1 - theta) * bFR;
                    rightReturnData(7, i) = bFR * ny * (1 - theta) * BFevalL;
                    MatSetValue(BF, Nu + posR + j, posL + i, fData.right_[j](1, i), ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posR + j, Nu + posL + i, fData.right_[j](7, i),  ADD_VALUES);  // V coefficient

                            GLR
                leftReturnData(4, i) = BFevalR * ny * (theta)*bFL;
                leftReturnData(10, i) = bFL * ny * (theta)*BFevalR;
                    MatSetValue(BF, Nu + posL + j, posR + i, -fData.left_[j](4, i),  ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posL + j, Nu + posR + i, -fData.left_[j](10, i),  ADD_VALUES);  // V coefficient

                                    GLL
                leftReturnData(1, i) = BFevalL * ny * (1 - theta) * bFL;
                leftReturnData(7, i) = bFL * ny * (1 - theta) * BFevalL;
                    MatSetValue(BF, Nu + posL + j, posL + i, fData.left_[j](1, i), ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posL + j, Nu + posL + i, -fData.left_[j](7, i),  ADD_VALUES);  // V coefficient

                                    GRR
                    rightReturnData(4, i) = BFevalR * ny * (theta)*bFR;
                    rightReturnData(10, i) = bFR * ny * (theta)*BFevalR;
                    MatSetValue(BF, Nu + posR + j, posR + i, -fData.right_[j](4, i), ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posR + j, Nu + posR + i, fData.right_[j](10, i),  ADD_VALUES);  // V coefficient





                leftReturnData(2, i) = BFevalL * nz * (1 - theta) * bFL;
                leftReturnData(5, i) = BFevalR * nz * (theta)*bFL;
                leftReturnData(8, i) = bFL * nz * (1 - theta) * BFevalL;
                leftReturnData(11, i) = bFL * nz * (theta)*BFevalR;
                rightReturnData(2, i) = BFevalL * nz * (1 - theta) * bFR;
                rightReturnData(5, i) = BFevalR * nz * (theta)*bFR;
                rightReturnData(8, i) = bFR * nz * (1 - theta) * BFevalL;
                rightReturnData(11, i) = bFR * nz * (theta)*BFevalR;
                    MatSetValue(BF, Nu + Nv + posL + j, posR + i,  -fData.left_[j](5, i),  ADD_VALUES);  // W coefficient
                    MatSetValue(BF, Nu + Nv + posL + j, posL + i, fData.left_[j](2, i), ADD_VALUES);  // W coefficient
                    MatSetValue(BF, Nu + Nv + posR + j, posR + i,  -fData.right_[j](5, i),  ADD_VALUES);  // W coefficient
                    MatSetValue(BF, Nu + Nv + posR + j, posL + i, fData.right_[j](2, i), ADD_VALUES);  // W coefficient

                    MatSetValue(globalDIV, posL + j, Nu + Nv + posL + i,  -fData.left_[j](8, i),  ADD_VALUES);  // W coefficient
                    MatSetValue(globalDIV, posR + j, Nu + Nv + posL + i, fData.right_[j](8, i),  ADD_VALUES);  // W coefficient
                    MatSetValue(globalDIV, posL + j, Nu + Nv + posR + i, -fData.left_[j](11, i), ADD_VALUES);  // W coefficient
                    MatSetValue(globalDIV, posR + j, Nu + Nv + posR + i, fData.right_[j](11, i), ADD_VALUES);  // W coefficient



                    leftReturnData(0, i) = BFevalL * nx * (1 - theta) * bFL;
                    leftReturnData(6, i) = bFL * nx * (1 - theta) * BFevalL;
                    MatSetValue(BF, posL + j, posL + i, fData.left_[j](0, i), ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posL + j, posL + i, -fData.left_[j](6, i), ADD_VALUES);  // U coefficient

                    rightReturnData(0, i) = BFevalL * nx * (1 - theta) * bFR;
                    rightReturnData(6, i) = bFR * nx * (1 - theta) * BFevalL;
                    MatSetValue(BF, posR + j, posL + i, fData.right_[j](0, i), ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posR + j, posL + i, fData.right_[j](6, i), ADD_VALUES);  // U coefficient

                    leftReturnData(3, i) = BFevalR * nx * (theta)*bFL;
                    leftReturnData(9, i) = bFL * nx * (theta)*BFevalR;
                    MatSetValue(BF, posL + j, posR + i, -fData.left_[j](3, i), ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posL + j, posR + i, -fData.left_[j](9, i),ADD_VALUES);  // U coefficient

                    rightReturnData(3, i) = BFevalR * nx * (theta)*bFR;
                    rightReturnData(9, i) = bFR * nx * (theta)*BFevalR;
                    MatSetValue(BF, posR + j, posR + i, -fData.right_[j](3, i),  ADD_VALUES);  // U coefficient
                    MatSetValue(globalDIV, posR + j, posR + i, fData.right_[j](9, i),ADD_VALUES);  // U coefficient

                }
            }
*/
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++) // Should be N_Left + 1, etc.
            {
                for (unsigned int j = 0; j <= Order_Polynomials_left; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    // GLL
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_L
                    {
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += (1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor;

                    }
                    //std::cout << "value_e = " << value_e << std::endl;
                    //if (abs(ny) > 1e-5)
                    {
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Left[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    }
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], -nx*value_e, ADD_VALUES);
                    //}
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                }
            }
            //std::cout << "GLR: " << std::endl;
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
                        //value_e += (1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor;
                    }
                    //if (abs(nx) > 1e-5)
                    //{
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    //MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[j],  posL+Node_Numbers_On_Boundary_Left[i] , nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i],   posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i],  ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j],  -ny*value_e, ADD_VALUES);

                    //std::cout << posL+Node_Numbers_On_Boundary_Left[i] << " " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                }
            }
/*

                            GRL
                    rightReturnData(1, i) = BFevalL * ny * (1 - theta) * bFR;
                    rightReturnData(7, i) = bFR * ny * (1 - theta) * BFevalL;
                    MatSetValue(BF, Nu + posR + j, posL + i, fData.right_[j](1, i), ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posR + j, Nu + posL + i, fData.right_[j](7, i),  ADD_VALUES);  // V coefficient

                            GLR
                leftReturnData(4, i) = BFevalR * ny * (theta)*bFL;
                leftReturnData(10, i) = bFL * ny * (theta)*BFevalR;
                    MatSetValue(BF, Nu + posL + j, posR + i, -fData.left_[j](4, i),  ADD_VALUES);  // V coefficient
                    MatSetValue(globalDIV, posL + j, Nu + posR + i, -fData.left_[j](10, i),  ADD_VALUES);  // V coefficient
                                */
            //std::cout << "GRL: " << std::endl;
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
                        //value_e += -theta*w_a[q]*Li*Lj*Jacobian*factor;
                    }
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    //MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i],  nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j],   -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i], ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                    //std::cout << posR+Node_Numbers_On_Boundary_Right[i] << " " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                }
            }
            //std::cout << "GRR: " << std::endl;
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
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], -ny*value_e, ADD_VALUES);

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
/*--------------------------------------------------------------------------*//*--------------------------------------------------------------------------*/
void create_Matrices_Quads_Full(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, const double &rho_0_Deriv, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &M2_small, Mat &NMat, Mat &NDerivMat)
{
    double Np = (N+1)*(N+1);
    //Mat Ex, ExT, Ey, EyT;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, 5*Np, NULL, &E);   //number of possible nonzero blocks are 5: element and his 4 neighbours (2D)
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, 2*N_Nodes, 2*5*Np, NULL, &ET);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &invM);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &invM_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &M1);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M1_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M2_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NMat);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &M2);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NDerivMat);

    std::cout << "Start Elemental Calculations " << std::endl;
    for (auto e = List_Of_Elements.begin(); e < List_Of_Elements.end(); e++)
    {
        unsigned int Np = (*e).get_Number_Of_Nodes();
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

        Vec ri;
        ri = JacobiGL(0, 0, Order_Polynomials);
        Vec Weights;
        Vec QuadraturePoints;
        //QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
        QuadraturePoints = JacobiGL_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
        PetscScalar *w, *qp;
        VecGetArray(Weights, &w);
        VecGetArray(QuadraturePoints, &qp);

        for (unsigned int k = 1; k <= Np; k++)
        {
            unsigned int alpha = (k-1)%(N+1);
            unsigned int beta = (k-1)/(N+1);
            for (unsigned int l = 1; l <= Np; l++)
            {
                unsigned int gamma = (l-1)%(N+1);
                unsigned int delta = (l-1)/(N+1);
                double value_ex = 0.0;
                double value_ey = 0.0;
                double value_m = 0.0;
                double value_m1 = 0.0;
                double value_n = 0.0;
                double value_m2 = 0.0;
                double value_n_deriv = 0.0;

                double value_area = 0.0;
                for (unsigned int p = 0; p <= Order_Gaussian_Quadrature; p++)
                {
                    double L_alpha = LagrangePolynomial(ri, qp[p], alpha);
                    double L_gamma = LagrangePolynomial(ri, qp[p], gamma);
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double J, drdx, drdy, dsdx, dsdy, x, y;
                        Calculate_Jacobian_Quadrilateral((*e), List_Of_Vertices, qp[p], qp[q], J, drdx, drdy, dsdx, dsdy, x, y);
                        //std::cout << "J = " << J << ", drdx = " << drdx << ", drdy = " << drdy << ", dsdy = " << dsdy << ", dsdx = " << dsdx  << std::endl;
                        //std::cout << "qp[p] = " << qp[p] << ", qp[q] = " << qp[q] << std::endl;
                        //std::cout << "ID = " << (*e).getID() << ", x = " << x << ", y = " << y << ", rho_0 = " << rho_0_2D_system2(y, rho_0_Deriv) << std::endl;
                        double rho0 = rho_0_2D_system2(y, rho_0_Deriv);
                        double rho0deriv = rho_0_deriv_2D_system2(y, rho_0_Deriv);
                        double N2 = N_2_2D_system2(y, rho_0_Deriv);

                        double L_beta = LagrangePolynomial(ri, qp[q], beta);
                        double L_delta = LagrangePolynomial(ri, qp[q], delta);

                        value_m += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J;
                        value_m1 += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J / rho0;
                        //if (N2 == 0)
                        //{
                        //    std::cout << "N2 is zero " << std::endl;
                        //}
                        value_m2 += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J / rho0 / N2;

                        value_n += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0;
                        value_n_deriv += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0deriv;
                        value_area += w[p]*w[q]*J;

                        // w_q drho_0(r_q)/dr l_i(r_q) l_j(r_q)
                        value_ey += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0deriv;
                        if (Np > 1)
                        {
                            double dL_gamma = LagrangePolynomialDeriv(ri, qp[p], gamma);
                            double dL_delta = LagrangePolynomialDeriv(ri, qp[q], delta);
                            //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dx
                            value_ex += w[p]*L_alpha*L_beta * w[q] * (drdx * dL_gamma * L_delta + dsdx * L_gamma* dL_delta) * J * rho0;
                            //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dy
                            value_ey += w[p]*L_alpha*L_beta * w[q] * (drdy * dL_gamma * L_delta + dsdy * L_gamma* dL_delta) * J * rho0;
                        }
                    }
                }
                //std::cout << std::endl;
                //std::cout << "ID = " << (*e).getID() << ", Area = " << value_area << std::endl;
                //std::cout << std::endl;
                //std::cout << "k = " << k << "=> i = " << alpha << ", j = " << beta << ", value_ex = " << value_ex << std::endl;
                MatSetValue(M_Elemental, (k-1), (l-1), value_m, ADD_VALUES);
                MatSetValue(M1_small, pos+(k-1), pos+(l-1), value_m1, ADD_VALUES);
                MatSetValue(NMat, pos+(k-1), pos+(l-1), value_n, ADD_VALUES);
                MatSetValue(M2_small, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(M2, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(M2, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(NDerivMat, pos+(k-1), pos+(l-1), value_n_deriv, ADD_VALUES);
                MatSetValue(M1, pos+(k-1), pos+(l-1), value_m1, ADD_VALUES);
                MatSetValue(M1, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m1, ADD_VALUES);
                /****************************************************/
                double factor = -1.0;
                /****************************************************/
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
    MatAssemblyBegin(M2_small, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M2_small, MAT_FINAL_ASSEMBLY);
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
            double Jacobian = (*f).getJacobian();
            unsigned int Type_Boundary_Left = (*f).get_Type_Left();
            unsigned int Type_Boundary_Right = (*f).get_Type_Right();
            double theta = (*f).get_theta();
            double nx = (*f).get_nx();
            double ny = (*f).get_ny();

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

            /// Or use two different gaussian quadratures
            // unsigned int Order_Gaussian_Quadrature_L
            // unsigned int Order_Gaussian_Quadrature_R
            unsigned int Order_Gaussian_Quadrature  = ceil(std::max(2*Order_Polynomials_left, 2*Order_Polynomials_right)+3+N_Q);
            Order_Gaussian_Quadrature = std::max((uint)10, Order_Gaussian_Quadrature);
            // Order_Gaussian_Quadrature+1 = Number of Points

            Vec Weights;
            Vec QuadraturePoints;
            //QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
            QuadraturePoints = JacobiGL_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
            PetscScalar *w_a, *r_a;
            VecGetArray(QuadraturePoints, &r_a);
            VecGetArray(Weights, &w_a);

            Vec ri_left, ri_right;
            ri_left = JacobiGL(0, 0, Order_Polynomials_left);
            ri_right = JacobiGL(0, 0, Order_Polynomials_right);

            double x1 = List_Of_Vertices[(*f).get_Vertex_V1()].getxCoordinate();
            double y1 = List_Of_Vertices[(*f).get_Vertex_V1()].getyCoordinate();
            double x2 = List_Of_Vertices[(*f).get_Vertex_V2()].getxCoordinate();
            double y2 = List_Of_Vertices[(*f).get_Vertex_V2()].getyCoordinate();

            double dx = (x2-x1);
            double dy = (y2-y1);

            double factor = 1.0;
            // GLL
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++) // Should be N_Left + 1, etc.
            {
                for (unsigned int j = 0; j <= Order_Polynomials_left; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    // GLL
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_L
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2D_system2(y, rho_0_Deriv);
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += (1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor*rho0;

                    }
                    //if (abs(ny) > 1e-5)
                    {
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Left[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    }
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], -nx*value_e, ADD_VALUES);
                    //}
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                }
            }
            //std::cout << "GLR: " << std::endl;
            // GLR
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2D_system2(y, rho_0_Deriv);
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -(1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //if (abs(nx) > 1e-5)
                    //{
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    //MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[j],  posL+Node_Numbers_On_Boundary_Left[i] , nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i],   posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i],  ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j],  -ny*value_e, ADD_VALUES);
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
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2D_system2(y, rho_0_Deriv);
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += theta*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    //MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i],  nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j],   -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i], ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                }
            }
            //std::cout << "GRR: " << std::endl;
            // GRR
            for (unsigned int i = 0; i <= Order_Polynomials_right; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_R
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2D_system2(y, rho_0_Deriv);
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -theta*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], -ny*value_e, ADD_VALUES);
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
/*--------------------------------------------------------------------------*//*--------------------------------------------------------------------------*/
void create_Matrices_Quads_Full_IC(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, const double &rho_0_Deriv, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &M2_small, Mat &NMat, Mat &NDerivMat)
{
    double Np = (N+1)*(N+1);
    //Mat Ex, ExT, Ey, EyT;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, 5*Np, NULL, &E);   //number of possible nonzero blocks are 5: element and his 4 neighbours (2D)
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, 2*N_Nodes, 2*5*Np, NULL, &ET);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &invM);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &invM_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &M1);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M1_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M2_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NMat);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &M2);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NDerivMat);

    std::cout << "Start Elemental Calculations " << std::endl;
    for (auto e = List_Of_Elements.begin(); e < List_Of_Elements.end(); e++)
    {
        unsigned int Np = (*e).get_Number_Of_Nodes();
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

        Vec ri;
        ri = JacobiGL(0, 0, Order_Polynomials);
        Vec Weights;
        Vec QuadraturePoints;
        //QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
        QuadraturePoints = JacobiGL_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
        PetscScalar *w, *qp;
        VecGetArray(Weights, &w);
        VecGetArray(QuadraturePoints, &qp);

        for (unsigned int k = 1; k <= Np; k++)
        {
            unsigned int alpha = (k-1)%(N+1);
            unsigned int beta = (k-1)/(N+1);
            for (unsigned int l = 1; l <= Np; l++)
            {
                unsigned int gamma = (l-1)%(N+1);
                unsigned int delta = (l-1)/(N+1);
                double value_ex = 0.0;
                double value_ey = 0.0;
                double value_m = 0.0;
                double value_m1 = 0.0;
                double value_n = 0.0;
                double value_m2 = 0.0;
                double value_n_deriv = 0.0;

                double value_area = 0.0;
                for (unsigned int p = 0; p <= Order_Gaussian_Quadrature; p++)
                {
                    double L_alpha = LagrangePolynomial(ri, qp[p], alpha);
                    double L_gamma = LagrangePolynomial(ri, qp[p], gamma);
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double J, drdx, drdy, dsdx, dsdy, x, y;
                        Calculate_Jacobian_Quadrilateral((*e), List_Of_Vertices, qp[p], qp[q], J, drdx, drdy, dsdx, dsdy, x, y);
                        //std::cout << "J = " << J << ", drdx = " << drdx << ", drdy = " << drdy << ", dsdy = " << dsdy << ", dsdx = " << dsdx  << std::endl;
                        //std::cout << "qp[p] = " << qp[p] << ", qp[q] = " << qp[q] << std::endl;
                        //std::cout << "ID = " << (*e).getID() << ", x = " << x << ", y = " << y << ", rho_0 = " << rho_0_2D_system2(y, rho_0_Deriv) << std::endl;
                        double rho0 = rho_0_2DIC(y, rho_0_Deriv);
                        double rho0deriv = rho_0_deriv_2DIC(y, rho_0_Deriv);
                        double N2 = N_2_2DIC(y, rho_0_Deriv);

                        double L_beta = LagrangePolynomial(ri, qp[q], beta);
                        double L_delta = LagrangePolynomial(ri, qp[q], delta);

                        value_m += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J;
                        value_m1 += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J / rho0;
                        //if (N2 == 0)
                        //{
                        //    std::cout << "N2 is zero " << std::endl;
                        //}
                        value_m2 += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J / rho0 / N2;

                        value_n += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0;
                        value_n_deriv += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0deriv;
                        value_area += w[p]*w[q]*J;

                        // w_q drho_0(r_q)/dr l_i(r_q) l_j(r_q)
                        value_ey += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0deriv;
                        if (Np > 1)
                        {
                            double dL_gamma = LagrangePolynomialDeriv(ri, qp[p], gamma);
                            double dL_delta = LagrangePolynomialDeriv(ri, qp[q], delta);
                            //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dx
                            value_ex += w[p]*L_alpha*L_beta * w[q] * (drdx * dL_gamma * L_delta + dsdx * L_gamma* dL_delta) * J * rho0;
                            //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dy
                            value_ey += w[p]*L_alpha*L_beta * w[q] * (drdy * dL_gamma * L_delta + dsdy * L_gamma* dL_delta) * J * rho0;
                        }
                    }
                }
                //std::cout << std::endl;
                //std::cout << "ID = " << (*e).getID() << ", Area = " << value_area << std::endl;
                //std::cout << std::endl;
                //std::cout << "k = " << k << "=> i = " << alpha << ", j = " << beta << ", value_ex = " << value_ex << std::endl;
                MatSetValue(M_Elemental, (k-1), (l-1), value_m, ADD_VALUES);
                MatSetValue(M1_small, pos+(k-1), pos+(l-1), value_m1, ADD_VALUES);
                MatSetValue(NMat, pos+(k-1), pos+(l-1), value_n, ADD_VALUES);
                MatSetValue(M2_small, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(M2, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(M2, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(NDerivMat, pos+(k-1), pos+(l-1), value_n_deriv, ADD_VALUES);
                MatSetValue(M1, pos+(k-1), pos+(l-1), value_m1, ADD_VALUES);
                MatSetValue(M1, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m1, ADD_VALUES);
                /****************************************************/
                double factor = -1.0;
                /****************************************************/
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
    MatAssemblyBegin(M2_small, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M2_small, MAT_FINAL_ASSEMBLY);
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
            double Jacobian = (*f).getJacobian();
            unsigned int Type_Boundary_Left = (*f).get_Type_Left();
            unsigned int Type_Boundary_Right = (*f).get_Type_Right();
            double theta = (*f).get_theta();
            double nx = (*f).get_nx();
            double ny = (*f).get_ny();

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

            /// Or use two different gaussian quadratures
            // unsigned int Order_Gaussian_Quadrature_L
            // unsigned int Order_Gaussian_Quadrature_R
            unsigned int Order_Gaussian_Quadrature  = ceil(std::max(2*Order_Polynomials_left, 2*Order_Polynomials_right)+3+N_Q);
            Order_Gaussian_Quadrature = std::max((uint)10, Order_Gaussian_Quadrature);
            // Order_Gaussian_Quadrature+1 = Number of Points

            Vec Weights;
            Vec QuadraturePoints;
            //QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
            QuadraturePoints = JacobiGL_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
            PetscScalar *w_a, *r_a;
            VecGetArray(QuadraturePoints, &r_a);
            VecGetArray(Weights, &w_a);

            Vec ri_left, ri_right;
            ri_left = JacobiGL(0, 0, Order_Polynomials_left);
            ri_right = JacobiGL(0, 0, Order_Polynomials_right);

            double x1 = List_Of_Vertices[(*f).get_Vertex_V1()].getxCoordinate();
            double y1 = List_Of_Vertices[(*f).get_Vertex_V1()].getyCoordinate();
            double x2 = List_Of_Vertices[(*f).get_Vertex_V2()].getxCoordinate();
            double y2 = List_Of_Vertices[(*f).get_Vertex_V2()].getyCoordinate();

            double dx = (x2-x1);
            double dy = (y2-y1);

            double factor = 1.0;
            // GLL
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++) // Should be N_Left + 1, etc.
            {
                for (unsigned int j = 0; j <= Order_Polynomials_left; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    // GLL
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_L
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2DIC(y, rho_0_Deriv);
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += (1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor*rho0;

                    }
                    //if (abs(ny) > 1e-5)
                    {
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Left[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    }
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], -nx*value_e, ADD_VALUES);
                    //}
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                }
            }
            //std::cout << "GLR: " << std::endl;
            // GLR
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2DIC(y, rho_0_Deriv);
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -(1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //if (abs(nx) > 1e-5)
                    //{
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    //MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[j],  posL+Node_Numbers_On_Boundary_Left[i] , nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i],   posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i],  ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j],  -ny*value_e, ADD_VALUES);
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
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2DIC(y, rho_0_Deriv);
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += theta*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    //MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i],  nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j],   -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i], ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                }
            }
            //std::cout << "GRR: " << std::endl;
            // GRR
            for (unsigned int i = 0; i <= Order_Polynomials_right; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_R
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2DIC(y, rho_0_Deriv);
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -theta*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], -ny*value_e, ADD_VALUES);
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
void create_Matrices_Quads_Full_WA(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, const double &rho_0_Deriv, const double &Fr, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &M2_small, Mat &NMat, Mat &NDerivMat)
{
    double Np = (N+1)*(N+1);
    //Mat Ex, ExT, Ey, EyT;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, 5*Np, NULL, &E);   //number of possible nonzero blocks are 5: element and his 4 neighbours (2D)
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, 2*N_Nodes, 2*5*Np, NULL, &ET);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &invM);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &invM_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &M1);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M1_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M2_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NMat);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &M2);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NDerivMat);

    std::cout << "Start Elemental Calculations " << std::endl;
    for (auto e = List_Of_Elements.begin(); e < List_Of_Elements.end(); e++)
    {
        unsigned int Np = (*e).get_Number_Of_Nodes();
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

        Vec ri;
        ri = JacobiGL(0, 0, Order_Polynomials);
        Vec Weights;
        Vec QuadraturePoints;
        //QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
        QuadraturePoints = JacobiGL_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
        PetscScalar *w, *qp;
        VecGetArray(Weights, &w);
        VecGetArray(QuadraturePoints, &qp);

        for (unsigned int k = 1; k <= Np; k++)
        {
            unsigned int alpha = (k-1)%(N+1);
            unsigned int beta = (k-1)/(N+1);
            for (unsigned int l = 1; l <= Np; l++)
            {
                unsigned int gamma = (l-1)%(N+1);
                unsigned int delta = (l-1)/(N+1);
                double value_ex = 0.0;
                double value_ey = 0.0;
                double value_m = 0.0;
                double value_m1 = 0.0;
                double value_n = 0.0;
                double value_m2 = 0.0;
                double value_n_deriv = 0.0;

                double value_area = 0.0;
                for (unsigned int p = 0; p <= Order_Gaussian_Quadrature; p++)
                {
                    double L_alpha = LagrangePolynomial(ri, qp[p], alpha);
                    double L_gamma = LagrangePolynomial(ri, qp[p], gamma);
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double J, drdx, drdy, dsdx, dsdy, x, y;
                        Calculate_Jacobian_Quadrilateral((*e), List_Of_Vertices, qp[p], qp[q], J, drdx, drdy, dsdx, dsdy, x, y);
                        //std::cout << "J = " << J << ", drdx = " << drdx << ", drdy = " << drdy << ", dsdy = " << dsdy << ", dsdx = " << dsdx  << std::endl;
                        //std::cout << "qp[p] = " << qp[p] << ", qp[q] = " << qp[q] << std::endl;
                        //std::cout << "ID = " << (*e).getID() << ", x = " << x << ", y = " << y << ", rho_0 = " << rho_0_2D_system2(y, rho_0_Deriv) << std::endl;
                        double rho0 = rho_0_2DWA(y, rho_0_Deriv, Fr);
                        double rho0deriv = rho_0_deriv_2DWA(y, rho_0_Deriv, Fr);
                        double N2 = N_2_2DWA(y, rho_0_Deriv, Fr);

                        double L_beta = LagrangePolynomial(ri, qp[q], beta);
                        double L_delta = LagrangePolynomial(ri, qp[q], delta);

                        value_m += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J;
                        value_m1 += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J / rho0;
                        //if (N2 == 0)
                        //{
                        //    std::cout << "N2 is zero " << std::endl;
                        //}
                        value_m2 += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J / rho0 / N2;

                        value_n += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0;
                        value_n_deriv += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0deriv;
                        value_area += w[p]*w[q]*J;

                        // w_q drho_0(r_q)/dr l_i(r_q) l_j(r_q)
                        value_ey += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0deriv;
                        if (Np > 1)
                        {
                            double dL_gamma = LagrangePolynomialDeriv(ri, qp[p], gamma);
                            double dL_delta = LagrangePolynomialDeriv(ri, qp[q], delta);
                            //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dx
                            value_ex += w[p]*L_alpha*L_beta * w[q] * (drdx * dL_gamma * L_delta + dsdx * L_gamma* dL_delta) * J * rho0;
                            //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dy
                            value_ey += w[p]*L_alpha*L_beta * w[q] * (drdy * dL_gamma * L_delta + dsdy * L_gamma* dL_delta) * J * rho0;
                        }
                    }
                }
                //std::cout << std::endl;
                //std::cout << "ID = " << (*e).getID() << ", Area = " << value_area << std::endl;
                //std::cout << std::endl;
                //std::cout << "k = " << k << "=> i = " << alpha << ", j = " << beta << ", value_ex = " << value_ex << std::endl;
                MatSetValue(M_Elemental, (k-1), (l-1), value_m, ADD_VALUES);
                MatSetValue(M1_small, pos+(k-1), pos+(l-1), value_m1, ADD_VALUES);
                MatSetValue(NMat, pos+(k-1), pos+(l-1), value_n, ADD_VALUES);
                MatSetValue(M2_small, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(M2, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(M2, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(NDerivMat, pos+(k-1), pos+(l-1), value_n_deriv, ADD_VALUES);
                MatSetValue(M1, pos+(k-1), pos+(l-1), value_m1, ADD_VALUES);
                MatSetValue(M1, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m1, ADD_VALUES);
                /****************************************************/
                double factor = -1.0;
                /****************************************************/
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
    MatAssemblyBegin(M2_small, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M2_small, MAT_FINAL_ASSEMBLY);
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
            double Jacobian = (*f).getJacobian();
            unsigned int Type_Boundary_Left = (*f).get_Type_Left();
            unsigned int Type_Boundary_Right = (*f).get_Type_Right();
            double theta = (*f).get_theta();
            double nx = (*f).get_nx();
            double ny = (*f).get_ny();

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

            /// Or use two different gaussian quadratures
            // unsigned int Order_Gaussian_Quadrature_L
            // unsigned int Order_Gaussian_Quadrature_R
            unsigned int Order_Gaussian_Quadrature  = ceil(std::max(2*Order_Polynomials_left, 2*Order_Polynomials_right)+3+N_Q);
            Order_Gaussian_Quadrature = std::max((uint)10, Order_Gaussian_Quadrature);
            // Order_Gaussian_Quadrature+1 = Number of Points

            Vec Weights;
            Vec QuadraturePoints;
            //QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
            QuadraturePoints = JacobiGL_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
            PetscScalar *w_a, *r_a;
            VecGetArray(QuadraturePoints, &r_a);
            VecGetArray(Weights, &w_a);

            Vec ri_left, ri_right;
            ri_left = JacobiGL(0, 0, Order_Polynomials_left);
            ri_right = JacobiGL(0, 0, Order_Polynomials_right);

            double x1 = List_Of_Vertices[(*f).get_Vertex_V1()].getxCoordinate();
            double y1 = List_Of_Vertices[(*f).get_Vertex_V1()].getyCoordinate();
            double x2 = List_Of_Vertices[(*f).get_Vertex_V2()].getxCoordinate();
            double y2 = List_Of_Vertices[(*f).get_Vertex_V2()].getyCoordinate();

            double dx = (x2-x1);
            double dy = (y2-y1);

            double factor = 1.0;
            // GLL
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++) // Should be N_Left + 1, etc.
            {
                for (unsigned int j = 0; j <= Order_Polynomials_left; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    // GLL
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_L
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2DWA(y, rho_0_Deriv, Fr);
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += (1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor*rho0;

                    }
                    //if (abs(ny) > 1e-5)
                    {
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Left[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    }
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], -nx*value_e, ADD_VALUES);
                    //}
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                }
            }
            //std::cout << "GLR: " << std::endl;
            // GLR
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2DWA(y, rho_0_Deriv, Fr);
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -(1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //if (abs(nx) > 1e-5)
                    //{
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    //MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[j],  posL+Node_Numbers_On_Boundary_Left[i] , nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i],   posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i],  ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j],  -ny*value_e, ADD_VALUES);
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
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2DWA(y, rho_0_Deriv, Fr);
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += theta*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    //MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i],  nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j],   -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i], ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                }
            }
            //std::cout << "GRR: " << std::endl;
            // GRR
            for (unsigned int i = 0; i <= Order_Polynomials_right; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_R
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double rho0 = rho_0_2DWA(y, rho_0_Deriv, Fr);
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -theta*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], -ny*value_e, ADD_VALUES);
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
void create_Matrices_Quads_EB(const std::vector<VertexCoordinates2D> &List_Of_Vertices, const std::vector<InternalBoundariesSquares2D> &List_Of_Boundaries, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &N, const unsigned int &N_Q, const double &rho_0_Deriv, const double &Fr, Mat &E, Mat &ET, Mat &invM, Mat &invM_small, Mat &M1, Mat &M1_small, Mat &M2, Mat &M2_small, Mat &NMat, Mat &NDerivMat)
{
    double Np = (N+1)*(N+1);
    //Mat Ex, ExT, Ey, EyT;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, 5*Np, NULL, &E);   //number of possible nonzero blocks are 5: element and his 4 neighbours (2D)
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, 2*N_Nodes, 2*5*Np, NULL, &ET);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &invM);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &invM_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &M1);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M1_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &M2_small);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NMat);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 2*Np, NULL, &M2);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, N_Nodes, N_Nodes, 2*Np, NULL, &NDerivMat);

    std::cout << "Start Elemental Calculations " << std::endl;
    for (auto e = List_Of_Elements.begin(); e < List_Of_Elements.end(); e++)
    {
        unsigned int Np = (*e).get_Number_Of_Nodes();
        unsigned int pos = (*e).getPosition();

        unsigned int Order_Polynomials = (*e).get_Order_Of_Polynomials();

        unsigned int Order_Gaussian_Quadrature = 2*Order_Polynomials+3+N_Q;//ceil(Order_Polynomials+3+N_Q); // + higher order for rho_0 term
        Order_Gaussian_Quadrature = 1;// 28;//std::max((uint)10, Order_Gaussian_Quadrature);

        PetscInt in[Np];
        for (unsigned int n=0;n<Np; n++)
        {
            in[n] = n+pos;
        }

        Mat M_Elemental;
        Mat invM_Elemental;
        MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Np, Np, NULL, &M_Elemental);

        Vec ri;
        ri = JacobiGL(0, 0, Order_Polynomials);
        Vec Weights;
        Vec QuadraturePoints;
        //QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
        QuadraturePoints = JacobiGL_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
        PetscScalar *w, *qp;
        VecGetArray(Weights, &w);
        VecGetArray(QuadraturePoints, &qp);
        double rho0 = 1.0;

        for (unsigned int k = 1; k <= Np; k++)
        {
            unsigned int alpha = (k-1)%(N+1);
            unsigned int beta = (k-1)/(N+1);
            for (unsigned int l = 1; l <= Np; l++)
            {
                unsigned int gamma = (l-1)%(N+1);
                unsigned int delta = (l-1)/(N+1);
                double value_ex = 0.0;
                double value_ey = 0.0;
                double value_m = 0.0;
                double value_m1 = 0.0;
                double value_n = 0.0;
                double value_m2 = 0.0;
                double value_n_deriv = 0.0;

                double value_area = 0.0;
                for (unsigned int p = 0; p <= Order_Gaussian_Quadrature; p++)
                {
                    double L_alpha = LagrangePolynomial(ri, qp[p], alpha);
                    double L_gamma = LagrangePolynomial(ri, qp[p], gamma);
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double J, drdx, drdy, dsdx, dsdy, x, y;
                        Calculate_Jacobian_Quadrilateral((*e), List_Of_Vertices, qp[p], qp[q], J, drdx, drdy, dsdx, dsdy, x, y);
                        //std::cout << "J = " << J << ", drdx = " << drdx << ", drdy = " << drdy << ", dsdy = " << dsdy << ", dsdx = " << dsdx  << std::endl;
                        //std::cout << "qp[p] = " << qp[p] << ", qp[q] = " << qp[q] << std::endl;
                        //std::cout << "ID = " << (*e).getID() << ", x = " << x << ", y = " << y << ", rho_0 = " << rho_0_2D_system2(y, rho_0_Deriv) << std::endl;
                        double N2 = N_2_2DEB(y, rho_0_Deriv, Fr);
                        double rho0deriv = rho_0_deriv_2DEB(y, rho_0_Deriv, Fr); // = - rho0 * N2 / g

                        double L_beta = LagrangePolynomial(ri, qp[q], beta);
                        double L_delta = LagrangePolynomial(ri, qp[q], delta);

                        value_m += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J;
                        value_m1 += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J / rho0;
                        //if (N2 == 0)
                        //{
                        //    std::cout << "N2 is zero " << std::endl;
                        //}
                        value_m2 += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J / rho0 / N2 /Fr/Fr;

                        value_n += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0;
                        value_n_deriv += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0deriv;
                        value_area += w[p]*w[q]*J;

                        // w_q drho_0(r_q)/dr l_i(r_q) l_j(r_q)
                        //value_ey += w[p]*L_alpha*L_beta * w[q]*L_gamma*L_delta * J * rho0deriv;
                        if (Np > 1)
                        {
                            double dL_gamma = LagrangePolynomialDeriv(ri, qp[p], gamma);
                            double dL_delta = LagrangePolynomialDeriv(ri, qp[q], delta);
                            //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dx
                            value_ex += w[p]*L_alpha*L_beta * w[q] * (drdx * dL_gamma * L_delta + dsdx * L_gamma* dL_delta) * J * rho0;
                            //w_q rho_0(r_q) l_i(r_q) dl_j(r_q)/dy
                            value_ey += w[p]*L_alpha*L_beta * w[q] * (drdy * dL_gamma * L_delta + dsdy * L_gamma* dL_delta) * J * rho0;
                        }
                    }
                }
                //std::cout << std::endl;
                //std::cout << "ID = " << (*e).getID() << ", Area = " << value_area << std::endl;
                //std::cout << std::endl;
                //std::cout << "k = " << k << "=> i = " << alpha << ", j = " << beta << ", value_ex = " << value_ex << std::endl;
                MatSetValue(M_Elemental, (k-1), (l-1), value_m, ADD_VALUES);
                MatSetValue(M1_small, pos+(k-1), pos+(l-1), value_m1, ADD_VALUES);
                MatSetValue(NMat, pos+(k-1), pos+(l-1), value_n, ADD_VALUES);
                MatSetValue(M2_small, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(M2, pos+(k-1), pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(M2, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m2, ADD_VALUES);
                MatSetValue(NDerivMat, pos+(k-1), pos+(l-1), value_n_deriv, ADD_VALUES);
                MatSetValue(M1, pos+(k-1), pos+(l-1), value_m1, ADD_VALUES);
                MatSetValue(M1, N_Nodes+pos+(k-1), N_Nodes+pos+(l-1), value_m1, ADD_VALUES);
                /****************************************************/
                double factor = -1.0;
                /****************************************************/
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
                    /// Order 5, 9, 13, 17 give a nonzero difference
                    /// Independent of Order_Gaussian_Quadrature

    }
    MatAssemblyBegin(M1, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M1, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(M2, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M2, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(M2_small, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M2_small, MAT_FINAL_ASSEMBLY);
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
            double Jacobian = (*f).getJacobian();
            unsigned int Type_Boundary_Left = (*f).get_Type_Left();
            unsigned int Type_Boundary_Right = (*f).get_Type_Right();
            double theta = (*f).get_theta();
            double nx = (*f).get_nx();
            double ny = (*f).get_ny();

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

            /// Or use two different gaussian quadratures
            // unsigned int Order_Gaussian_Quadrature_L
            // unsigned int Order_Gaussian_Quadrature_R
            unsigned int Order_Gaussian_Quadrature  = ceil(std::max(2*Order_Polynomials_left, 2*Order_Polynomials_right)+3+N_Q);
            Order_Gaussian_Quadrature = std::max((uint)10, Order_Gaussian_Quadrature);
            // Order_Gaussian_Quadrature+1 = Number of Points

            Vec Weights;
            Vec QuadraturePoints;
            //QuadraturePoints = JacobiGQ_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
            QuadraturePoints = JacobiGL_withWeights(0, 0, Order_Gaussian_Quadrature, Weights);
            PetscScalar *w_a, *r_a;
            VecGetArray(QuadraturePoints, &r_a);
            VecGetArray(Weights, &w_a);

            Vec ri_left, ri_right;
            ri_left = JacobiGL(0, 0, Order_Polynomials_left);
            ri_right = JacobiGL(0, 0, Order_Polynomials_right);

            double x1 = List_Of_Vertices[(*f).get_Vertex_V1()].getxCoordinate();
            double y1 = List_Of_Vertices[(*f).get_Vertex_V1()].getyCoordinate();
            double x2 = List_Of_Vertices[(*f).get_Vertex_V2()].getxCoordinate();
            double y2 = List_Of_Vertices[(*f).get_Vertex_V2()].getyCoordinate();

            double dx = (x2-x1);
            double dy = (y2-y1);

            double factor = 1.0;
            double rho0 = 1.0;
            // GLL
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++) // Should be N_Left + 1, etc.
            {
                for (unsigned int j = 0; j <= Order_Polynomials_left; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    // GLL
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_L
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += (1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor*rho0;

                    }
                    //if (abs(ny) > 1e-5)
                    {
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Left[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    }
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], -nx*value_e, ADD_VALUES);
                    //}
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j], N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                }
            }
            //std::cout << "GLR: " << std::endl;
            // GLR
            for (unsigned int i = 0; i <= Order_Polynomials_left; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++)
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double Li = LagrangePolynomial(ri_left, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -(1.0-theta)*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //if (abs(nx) > 1e-5)
                    //{
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posL+Node_Numbers_On_Boundary_Left[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i],  posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[i], -ny*value_e, ADD_VALUES);

                    //MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[j],  posL+Node_Numbers_On_Boundary_Left[i] , nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i],   posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], posL+Node_Numbers_On_Boundary_Left[i],  ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j],  -ny*value_e, ADD_VALUES);
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
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_left, r_a[q], j);
                        value_e += theta*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posL+Node_Numbers_On_Boundary_Left[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posL+Node_Numbers_On_Boundary_Left[j],  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    //MatSetValue(E,  posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i],  nx*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posL+Node_Numbers_On_Boundary_Left[j],   -nx*value_e, ADD_VALUES);
                    //MatSetValue(E,  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], posR+Node_Numbers_On_Boundary_Right[i], ny*value_e, ADD_VALUES);
                    //MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i],  N_Nodes+posL+Node_Numbers_On_Boundary_Left[j], -ny*value_e, ADD_VALUES);
                }
            }
            //std::cout << "GRR: " << std::endl;
            // GRR
            for (unsigned int i = 0; i <= Order_Polynomials_right; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials_right; j++)
                {
                    // E Matrix
                    double value_e = 0.0;
                    for (unsigned int q = 0; q <= Order_Gaussian_Quadrature; q++) //Order_Gaussian_Quadrature_R
                    {
                        double y = y1 + (r_a[q]+1.0)*dy;
                        double Li = LagrangePolynomial(ri_right, r_a[q], i);
                        double Lj = LagrangePolynomial(ri_right, r_a[q], j);
                        value_e += -theta*w_a[q]*Li*Lj*Jacobian*factor*rho0;
                    }
                    //std::cout << "value_e = " << value_e << ", nx = " << nx << ", ny = " << ny << ". posL = " << posL << ". Node Numbers Local =  " << Node_Numbers_On_Boundary_Left[i] << ", " << Node_Numbers_On_Boundary_Right[j] << std::endl;
                    //std::cout << ". Node Numbers Global =  " << posR+Node_Numbers_On_Boundary_Right[i] << ", " << posR+Node_Numbers_On_Boundary_Right[j] << std::endl;
                    MatSetValue(E,  posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], nx*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], posR+Node_Numbers_On_Boundary_Right[i], -nx*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], -nx*value_e, ADD_VALUES);
                    MatSetValue(E,  N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], posR+Node_Numbers_On_Boundary_Right[j], ny*value_e, ADD_VALUES);
                    MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[j], N_Nodes+posR+Node_Numbers_On_Boundary_Right[i], -ny*value_e, ADD_VALUES);
                    ///MatSetValue(ET, posR+Node_Numbers_On_Boundary_Right[i], N_Nodes+posR+Node_Numbers_On_Boundary_Right[j], -ny*value_e, ADD_VALUES);
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
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+3*Np*Np,  NULL, &A); // Check Factor 10 // Change from Triangles to Quadrilaterals
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 3*N_Nodes, 3*N_Nodes, 10*6*Np+1+3*Np*Np,  NULL, &A); // Check Factor 10 // Change from Triangles to Quadrilaterals
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+3*Np*Np,  NULL, &B);
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 3*N_Nodes, 3*N_Nodes, 10*6*Np+1+3*Np*Np,  NULL, &B);
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
                //MatSetValue(A, 	i,     3*N_Nodes+cols[j],    -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(A, 	i,     2*N_Nodes+cols[j],    -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	i,     3*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	    ADD_VALUES);
                MatSetValue(B, 	i,     2*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	    ADD_VALUES);
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
        //MatSetValue(A, 3*N_Nodes+i, 3*N_Nodes+i, 1.0, ADD_VALUES);

        MatSetValue(B, i, i, 1.0, ADD_VALUES);
        MatSetValue(B, N_Nodes+i, N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(B, 2*N_Nodes+i, 2*N_Nodes+i, 1.0, ADD_VALUES);
        //MatSetValue(B, 3*N_Nodes+i, 3*N_Nodes+i, 1.0, ADD_VALUES);

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
                //MatSetValue(A, 	3*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(A, 	2*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	3*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	2*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

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
}/*--------------------------------------------------------------------------*/
extern void create_Compressible_System_MidPoint_Full(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B)
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

    Mat BF2, BF2_TEMP1, BF2_TEMP2;
    MatMatMult(E, invM_small, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP1);
    MatMatMult(BF2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP2);
    MatMatMult(invM, BF2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF2);
    MatDestroy(&BF2_TEMP1);
    MatDestroy(&BF2_TEMP2);

    Mat DIV, DIV_TEMP1, DIV_TEMP2;
    MatMatMult(ET, invM, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP1);
    MatMatMult(DIV_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP2);
    MatMatMult(invM_small, DIV_TEMP2, MAT_INITIAL_MATRIX, fillBF, &DIV);
    MatDestroy(&DIV_TEMP1);
    MatDestroy(&DIV_TEMP2);

    Mat C, C_TEMP1, C_TEMP2;
    MatMatMult(NMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &C_TEMP1);
    MatMatMult(C_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &C_TEMP2);
    MatMatMult(invM_small, C_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C);
    MatDestroy(&C_TEMP1);
    MatDestroy(&C_TEMP2);

    Mat C2, C2_TEMP1, C2_TEMP2;
    MatMatMult(NMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP1);
    MatMatMult(C2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP2);
    MatMatMult(invM_small, C2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C2);
    MatDestroy(&C2_TEMP1);
    MatDestroy(&C2_TEMP2);

    Mat D1, D1_TEMP1, D1_TEMP2;
    MatMatMult(NDerivMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP1);
    MatMatMult(D1_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP2);
    MatMatMult(invM_small, D1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D1);
    MatDestroy(&D1_TEMP1);
    MatDestroy(&D1_TEMP2);

    Mat D2, D2_TEMP1, D2_TEMP2;
    MatMatMult(NDerivMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP1);
    MatMatMult(D2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP2);
    MatMatMult(invM_small, D2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D2);
    MatDestroy(&D2_TEMP1);
    MatDestroy(&D2_TEMP2);

    Mat Laplacian;
	MatMatMult(BF1, DIV, MAT_INITIAL_MATRIX, 1, &Laplacian);
	MatScale(Laplacian, nu);



    //std::cout << "Laplacian = " << std::endl;
	//MatView(Laplacian, PETSC_VIEWER_STDOUT_SELF);

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
        /*
        MatGetRow(BF2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,    3*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,    2*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF2, i, &numberOfNonZeros, &cols, &values);
        */
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

        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(A, 	N_Nodes+i,     N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(C, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	3*N_Nodes+i,     N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(C2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C2, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(D1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	2*N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	2*N_Nodes+i,     N_Nodes+cols[j],    -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D1, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            //if (dummy!=0.0)
            if (abs(dummy)>1e-10)
            {
                MatSetValue(A, 	2*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	2*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	3*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);

    }
    MatDestroy(&BF1);
    MatDestroy(&BF2);
    MatDestroy(&DIV);
    MatDestroy(&C);
    MatDestroy(&C2);
    MatDestroy(&D1);
    MatDestroy(&D2);
    MatDestroy(&Laplacian);

    std::cout << " Global Matrices Constructed" << std::endl;
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
}
/*--------------------------------------------------------------------------*/
extern void create_Incompressible_System_MidPoint_Full(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B, Mat &DIV)
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
    MatMatMult(E, invM_small, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP1);
    MatMatMult(BF2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP2);
    MatMatMult(invM, BF2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF2);
    MatDestroy(&BF2_TEMP1);
    MatDestroy(&BF2_TEMP2);
    MatDestroy(&BF2);
    */

    // Not needed since we scale both sides of the equation (= constructed incompressibility)
    Mat DIV_TEMP1, DIV_TEMP2;
    MatMatMult(ET, invM, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP1);
    MatMatMult(DIV_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP2);
    MatMatMult(invM_small, DIV_TEMP2, MAT_INITIAL_MATRIX, fillBF, &DIV);
    MatDestroy(&DIV_TEMP1);

    Mat C, C_TEMP1, C_TEMP2;
    MatMatMult(NMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &C_TEMP1);
    MatMatMult(C_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &C_TEMP2);
    MatMatMult(invM_small, C_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C);
    MatDestroy(&C_TEMP1);
    MatDestroy(&C_TEMP2);

    Mat C2, C2_TEMP1, C2_TEMP2;
    MatMatMult(NMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP1);
    MatMatMult(C2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP2);
    MatMatMult(invM_small, C2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C2);
    MatDestroy(&C2_TEMP1);
    MatDestroy(&C2_TEMP2);

    Mat D1, D1_TEMP1, D1_TEMP2;
    MatMatMult(NDerivMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP1);
    MatMatMult(D1_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP2);
    MatMatMult(invM_small, D1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D1);
    MatDestroy(&D1_TEMP1);
    MatDestroy(&D1_TEMP2);

    Mat D2, D2_TEMP1, D2_TEMP2;
    MatMatMult(NDerivMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP1);
    MatMatMult(D2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP2);
    MatMatMult(invM_small, D2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D2);
    MatDestroy(&D2_TEMP1);
    MatDestroy(&D2_TEMP2);

    Mat Identity3;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, Np, NULL, &Identity3);
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(Identity3, 	N_Nodes+i,     cols[j],     dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);
    }
    MatAssemblyBegin(Identity3,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Identity3,MAT_FINAL_ASSEMBLY);


    Mat Laplacian;
	MatMatMult(DIV, BF1, MAT_INITIAL_MATRIX, 1, &Laplacian);
	MatScale(Laplacian, nu);

    //std::cout << "Laplacian = " << std::endl;
    //MatView(Laplacian, PETSC_VIEWER_STDOUT_SELF);

    Mat VectorLaplacian;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 10*6*Np+1+3*Np*Np, NULL, &VectorLaplacian);
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(VectorLaplacian, 	i,     cols[j],     dummy, 	ADD_VALUES);
                MatSetValue(VectorLaplacian, 	N_Nodes+i,     N_Nodes+cols[j],     dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
    }
    MatAssemblyBegin(VectorLaplacian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(VectorLaplacian,MAT_FINAL_ASSEMBLY);

    //std::cout << "VectorLaplacian = " << std::endl;
    //MatView(VectorLaplacian, PETSC_VIEWER_STDOUT_SELF);

    // invM_small Not needed since we scale both sides of the equation (= constructed incompressibility)
    Mat ALaplacian, ALaplacian_h, ADivLaplacian;
	MatMatMult(DIV_TEMP2, VectorLaplacian, MAT_INITIAL_MATRIX, 1, &ADivLaplacian);
	MatMatMult(DIV_TEMP2, BF1, MAT_INITIAL_MATRIX, 1, &ALaplacian);
	MatMatMult(DIV_TEMP2, Identity3, MAT_INITIAL_MATRIX, 1, &ALaplacian_h);
    MatDestroy(&DIV_TEMP2);
    MatDestroy(&Identity3);
    MatDestroy(&VectorLaplacian);

    //std::cout << "ADivLaplacian = " << std::endl;
    //MatView(ADivLaplacian, PETSC_VIEWER_STDOUT_SELF);

    /// Pass DIV_TEMP2 instead of DIV to compute L infinity norm of divergence

    // factor 4 = Number of Variables
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+9*Np*Np,  NULL, &A); // Check Factor 10 // Change from Triangles to Quadrilaterals
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+9*Np*Np,  NULL, &B);
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
            //if (abs(dummy)>1e-10)
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     3*N_Nodes+cols[j],    -DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);
        /*
        MatGetRow(BF2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,    3*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,    2*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF2, i, &numberOfNonZeros, &cols, &values);
        */
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

        MatSetValue(B, i, i, 1.0, ADD_VALUES);
        MatSetValue(B, N_Nodes+i, N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(B, 2*N_Nodes+i, 2*N_Nodes+i, 1.0, ADD_VALUES);


        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(A, 	N_Nodes+i,     N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);

        /*
        MatGetRow(C, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	3*N_Nodes+i,     N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(C2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C2, i, &numberOfNonZeros, &cols, &values);
        */
        MatGetRow(D1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	2*N_Nodes+i,     N_Nodes+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	2*N_Nodes+i,     N_Nodes+cols[j],    -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D1, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                //MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);

        /*
        MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            //if (dummy!=0.0)
            if (abs(dummy)>1e-10)
            {
                //MatSetValue(A, 	2*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	2*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                //MatSetValue(A, 	3*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	3*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);
        */
        MatGetRow(ALaplacian_h, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     2*N_Nodes+cols[j],     0.5*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(ALaplacian_h, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(ALaplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     	3*N_Nodes+cols[j],     dummy, 		ADD_VALUES);
            }
        }
        MatRestoreRow(ALaplacian, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(ADivLaplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     	cols[j],     0.5*dummy, 		ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     	cols[j],     -0.5*dummy, 		ADD_VALUES);
            }
        }
        MatRestoreRow(ADivLaplacian, i, &numberOfNonZeros, &cols, &values);

    }
    MatDestroy(&BF1);
    MatDestroy(&C);
    MatDestroy(&C2);
    MatDestroy(&D1);
    MatDestroy(&D2);
    MatDestroy(&Laplacian);
    MatDestroy(&ADivLaplacian);
    MatDestroy(&ALaplacian);
    MatDestroy(&ALaplacian_h);

    std::cout << " Global Matrices Constructed" << std::endl;
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
}
/*--------------------------------------------------------------------------*/
extern void create_WA_System_MidPoint(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B, Mat &DIV, const double &Re, const double &Fr)
{
    /*
    [I 0 0 dx] [U^n+1]   [I 0 0  0] [U^n]
    [0 I I dz] [W^n+1] = [0 I I  0] [W^n]
    [0 I I  0] [R^n+1]   [0 I I  0] [R^n]
    [0 0 dz L] [P^n+1]   [0 0 dz 0] [P^n] <= P^n not used
    */
    double Np = (N+1)*(N+1);
    std::cout << "Start Global Matrices Construction" << std::endl;

    double fillBF = 1;
    double gamma = PETSC_PI/20.0;

    Mat BF1, BF1_TEMP1, BF1_TEMP2;
    MatMatMult(E, invM_small, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP1);
    MatMatMult(BF1_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP2);
    MatMatMult(invM, BF1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF1);
    MatDestroy(&BF1_TEMP1);
    MatDestroy(&BF1_TEMP2);

    /*
    Mat BF2, BF2_TEMP1, BF2_TEMP2;
    MatMatMult(E, invM_small, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP1);
    MatMatMult(BF2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP2);
    MatMatMult(invM, BF2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF2);
    MatDestroy(&BF2_TEMP1);
    MatDestroy(&BF2_TEMP2);
    MatDestroy(&BF2);
    */

    // Not needed since we scale both sides of the equation (= constructed incompressibility)
    Mat DIV_TEMP1, DIV_TEMP2;
    MatMatMult(ET, invM, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP1);
    MatMatMult(DIV_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP2);
    MatMatMult(invM_small, DIV_TEMP2, MAT_INITIAL_MATRIX, fillBF, &DIV);
    MatDestroy(&DIV_TEMP1);

    Mat C, C_TEMP1, C_TEMP2;
    MatMatMult(NMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &C_TEMP1);
    MatMatMult(C_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &C_TEMP2);
    MatMatMult(invM_small, C_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C);
    MatDestroy(&C_TEMP1);
    MatDestroy(&C_TEMP2);

    Mat C2, C2_TEMP1, C2_TEMP2;
    MatMatMult(NMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP1);
    MatMatMult(C2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP2);
    MatMatMult(invM_small, C2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C2);
    MatDestroy(&C2_TEMP1);
    MatDestroy(&C2_TEMP2);

    Mat D1, D1_TEMP1, D1_TEMP2;
    MatMatMult(NDerivMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP1);
    MatMatMult(D1_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP2);
    MatMatMult(invM_small, D1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D1);
    MatDestroy(&D1_TEMP1);
    MatDestroy(&D1_TEMP2);

    Mat D2, D2_TEMP1, D2_TEMP2;
    MatMatMult(NDerivMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP1);
    MatMatMult(D2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP2);
    MatMatMult(invM_small, D2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D2);
    MatDestroy(&D2_TEMP1);
    MatDestroy(&D2_TEMP2);
	//MatScale(D2, 1.0/Fr/Fr);

    Mat Identity3;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, Np, NULL, &Identity3);
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(Identity3, 	i,     cols[j],     dummy*sin(gamma), 	ADD_VALUES);
                MatSetValue(Identity3, 	N_Nodes+i,     cols[j],     dummy*cos(gamma), 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);
    }
    MatAssemblyBegin(Identity3,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Identity3,MAT_FINAL_ASSEMBLY);


    Mat Laplacian;
	MatMatMult(DIV, BF1, MAT_INITIAL_MATRIX, 1, &Laplacian);
	MatScale(Laplacian, nu*1.0/Re); // nu = 0 or 1

    //std::cout << "Laplacian = " << std::endl;
    //MatView(Laplacian, PETSC_VIEWER_STDOUT_SELF);

    Mat VectorLaplacian;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 10*6*Np+1+3*Np*Np, NULL, &VectorLaplacian);
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(VectorLaplacian, 	i,     cols[j],     dummy, 	ADD_VALUES);
                MatSetValue(VectorLaplacian, 	N_Nodes+i,     N_Nodes+cols[j],     dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
    }
    MatAssemblyBegin(VectorLaplacian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(VectorLaplacian,MAT_FINAL_ASSEMBLY);

    //std::cout << "VectorLaplacian = " << std::endl;
    //MatView(VectorLaplacian, PETSC_VIEWER_STDOUT_SELF);

    // invM_small Not needed since we scale both sides of the equation (= constructed incompressibility)
    Mat ALaplacian, ALaplacian_h, ADivLaplacian;
	MatMatMult(DIV_TEMP2, VectorLaplacian, MAT_INITIAL_MATRIX, 1, &ADivLaplacian);
	MatMatMult(DIV_TEMP2, BF1, MAT_INITIAL_MATRIX, 1, &ALaplacian);
	MatMatMult(DIV_TEMP2, Identity3, MAT_INITIAL_MATRIX, 1, &ALaplacian_h);
    MatDestroy(&DIV_TEMP2);
    MatDestroy(&Identity3);
    MatDestroy(&VectorLaplacian);

    //std::cout << "ADivLaplacian = " << std::endl;
    //MatView(ADivLaplacian, PETSC_VIEWER_STDOUT_SELF);

    /// Pass DIV_TEMP2 instead of DIV to compute L infinity norm of divergence


    // factor 4 = Number of Variables
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+9*Np*Np,  NULL, &A); // Check Factor 10 // Change from Triangles to Quadrilaterals
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+9*Np*Np,  NULL, &B);
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
            //if (abs(dummy)>1e-10)
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     3*N_Nodes+cols[j],    -DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);
        /*
        MatGetRow(BF2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,    3*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,    2*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF2, i, &numberOfNonZeros, &cols, &values);
        */
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

        MatSetValue(B, i, i, 1.0, ADD_VALUES);
        MatSetValue(B, N_Nodes+i, N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(B, 2*N_Nodes+i, 2*N_Nodes+i, 1.0, ADD_VALUES);


        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(A, 	N_Nodes+i,     N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);

        /*
        MatGetRow(C, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	3*N_Nodes+i,     N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(C2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C2, i, &numberOfNonZeros, &cols, &values);
        */
        MatGetRow(D1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	2*N_Nodes+i,    cols[j],    0.5*DeltaT*dummy*sin(gamma), 	ADD_VALUES);
                MatSetValue(B, 	2*N_Nodes+i,     cols[j],    -0.5*DeltaT*dummy*sin(gamma), 	ADD_VALUES);
                MatSetValue(A, 	2*N_Nodes+i,     N_Nodes+cols[j],    0.5*DeltaT*dummy*cos(gamma), 	ADD_VALUES);
                MatSetValue(B, 	2*N_Nodes+i,     N_Nodes+cols[j],    -0.5*DeltaT*dummy*cos(gamma), 	ADD_VALUES);
            }
        }
        MatRestoreRow(D1, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy*sin(gamma), 	ADD_VALUES);
                MatSetValue(B, 	i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy*sin(gamma), 	ADD_VALUES);
                MatSetValue(A, 	N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy*cos(gamma), 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy*cos(gamma), 	ADD_VALUES);

                //MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);

        /*
        MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            //if (dummy!=0.0)
            if (abs(dummy)>1e-10)
            {
                //MatSetValue(A, 	2*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	2*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                //MatSetValue(A, 	3*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	3*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);
        */
        MatGetRow(ALaplacian_h, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     2*N_Nodes+cols[j],     0.5*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(ALaplacian_h, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(ALaplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     	3*N_Nodes+cols[j],     dummy, 		ADD_VALUES);
            }
        }
        MatRestoreRow(ALaplacian, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(ADivLaplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     	cols[j],     0.5*dummy, 		ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     	cols[j],     -0.5*dummy, 		ADD_VALUES);
            }
        }
        MatRestoreRow(ADivLaplacian, i, &numberOfNonZeros, &cols, &values);

    }
    MatDestroy(&BF1);
    MatDestroy(&C);
    MatDestroy(&C2);
    MatDestroy(&D1);
    MatDestroy(&D2);
    MatDestroy(&Laplacian);
    MatDestroy(&ADivLaplacian);
    MatDestroy(&ALaplacian);
    MatDestroy(&ALaplacian_h);

    std::cout << " Global Matrices Constructed" << std::endl;
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
}
/*--------------------------------------------------------------------------*/
extern void create_WA_System_Forced_MidPoint(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const Vec &Forcing_a, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B, Mat &DIV, const double &Re, const double &Fr)
{
    double Np = (N+1)*(N+1);
    std::cout << "Start Global Matrices Construction" << std::endl;

    /*
    [I 0 0 dx] [U^n+1]   [I 0 0   ax  ] [U^n]
    [0 I I dz] [W^n+1] = [0 I I   az  ] [W^n]
    [0 I I  0] [R^n+1]   [0 I I    0  ] [R^n]
    [0 0 dz L] [P^n+1]   [0 0 dz dx+dz] [F^n+1/2]
    */
    double fillBF = 1;
    double gamma = 0.0;//PETSC_PI/20.0;


    Mat F;
    MatCreateSeqAIJ(PETSC_COMM_SELF,2*N_Nodes,N_Nodes,1, NULL, &F);

    PetscScalar *F_a;
    VecGetArray(Forcing_a, &F_a);
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        MatSetValue(F, i, i, F_a[i], INSERT_VALUES);
        MatSetValue(F, N_Nodes+i, i, F_a[N_Nodes+i], INSERT_VALUES);
    }
    VecRestoreArray(Forcing_a, &F_a);
    MatAssemblyBegin(F,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(F,MAT_FINAL_ASSEMBLY);

    Mat F1;//, F_TEMP1, F_TEMP2;
    //MatMatMult(F, invM_small, MAT_INITIAL_MATRIX, fillBF, &F_TEMP1);
    //MatMatMult(F_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &F_TEMP2);
    //MatMatMult(invM, F_TEMP2, MAT_INITIAL_MATRIX, fillBF, &F1);
    MatMatMult(invM, F, MAT_INITIAL_MATRIX, fillBF, &F1);
    //MatDestroy(&F_TEMP1);
    //MatDestroy(&F_TEMP2);
    MatDestroy(&F);

    Mat BF1, BF1_TEMP1, BF1_TEMP2;
    MatMatMult(E, invM_small, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP1);
    MatMatMult(BF1_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &BF1_TEMP2);
    MatMatMult(invM, BF1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF1);
    MatDestroy(&BF1_TEMP1);
    MatDestroy(&BF1_TEMP2);

    /*
    Mat BF2, BF2_TEMP1, BF2_TEMP2;
    MatMatMult(E, invM_small, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP1);
    MatMatMult(BF2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP2);
    MatMatMult(invM, BF2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF2);
    MatDestroy(&BF2_TEMP1);
    MatDestroy(&BF2_TEMP2);
    MatDestroy(&BF2);
    */

    // Not needed since we scale both sides of the equation (= constructed incompressibility)
    Mat DIV_TEMP1, DIV_TEMP2;
    MatMatMult(ET, invM, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP1);
    MatMatMult(DIV_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP2);
    MatMatMult(invM_small, DIV_TEMP2, MAT_INITIAL_MATRIX, fillBF, &DIV);
    MatDestroy(&DIV_TEMP1);

    Mat C, C_TEMP1, C_TEMP2;
    MatMatMult(NMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &C_TEMP1);
    MatMatMult(C_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &C_TEMP2);
    MatMatMult(invM_small, C_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C);
    MatDestroy(&C_TEMP1);
    MatDestroy(&C_TEMP2);

    Mat C2, C2_TEMP1, C2_TEMP2;
    MatMatMult(NMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP1);
    MatMatMult(C2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP2);
    MatMatMult(invM_small, C2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C2);
    MatDestroy(&C2_TEMP1);
    MatDestroy(&C2_TEMP2);

    Mat D1, D1_TEMP1, D1_TEMP2;
    MatMatMult(NDerivMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP1);
    MatMatMult(D1_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP2);
    MatMatMult(invM_small, D1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D1);
    MatDestroy(&D1_TEMP1);
    MatDestroy(&D1_TEMP2);

    Mat D2, D2_TEMP1, D2_TEMP2;
    MatMatMult(NDerivMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP1);
    MatMatMult(D2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP2);
    MatMatMult(invM_small, D2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D2);
    MatDestroy(&D2_TEMP1);
    MatDestroy(&D2_TEMP2);
	//MatScale(D2, 1.0/Fr/Fr);

    Mat Identity3;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, Np, NULL, &Identity3);
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(Identity3, 	i,     cols[j],     dummy*sin(gamma), 	ADD_VALUES);
                MatSetValue(Identity3, 	N_Nodes+i,     cols[j],     dummy*cos(gamma), 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);
    }
    MatAssemblyBegin(Identity3,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Identity3,MAT_FINAL_ASSEMBLY);


    Mat Laplacian;
	MatMatMult(DIV, BF1, MAT_INITIAL_MATRIX, 1, &Laplacian);
	MatScale(Laplacian, nu*1.0/Re); // nu = 0 or 1

    //std::cout << "Laplacian = " << std::endl;
    //MatView(Laplacian, PETSC_VIEWER_STDOUT_SELF);

    Mat VectorLaplacian;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 10*6*Np+1+3*Np*Np, NULL, &VectorLaplacian);
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(VectorLaplacian, 	i,     cols[j],     dummy, 	ADD_VALUES);
                MatSetValue(VectorLaplacian, 	N_Nodes+i,     N_Nodes+cols[j],     dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
    }
    MatAssemblyBegin(VectorLaplacian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(VectorLaplacian,MAT_FINAL_ASSEMBLY);

    //std::cout << "VectorLaplacian = " << std::endl;
    //MatView(VectorLaplacian, PETSC_VIEWER_STDOUT_SELF);

    // invM_small Not needed since we scale both sides of the equation (= constructed incompressibility)
    Mat ALaplacian, ALaplacian_h, ADivLaplacian, AF;
	MatMatMult(DIV_TEMP2, VectorLaplacian, MAT_INITIAL_MATRIX, 1, &ADivLaplacian);
	MatMatMult(DIV_TEMP2, BF1, MAT_INITIAL_MATRIX, 1, &ALaplacian);
	MatMatMult(DIV_TEMP2, F1, MAT_INITIAL_MATRIX, 1, &AF);
	MatMatMult(DIV_TEMP2, Identity3, MAT_INITIAL_MATRIX, 1, &ALaplacian_h);
    MatDestroy(&DIV_TEMP2);
    MatDestroy(&Identity3);
    MatDestroy(&VectorLaplacian);

    //std::cout << "ADivLaplacian = " << std::endl;
    //MatView(ADivLaplacian, PETSC_VIEWER_STDOUT_SELF);

    /// Pass DIV_TEMP2 instead of DIV to compute L infinity norm of divergence


    // factor 4 = Number of Variables
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+9*Np*Np,  NULL, &A); // Check Factor 10 // Change from Triangles to Quadrilaterals
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+9*Np*Np,  NULL, &B);
    std::cout << "Global Matrices Preallocated" << std::endl;
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
            //if (abs(dummy)>1e-10)
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     3*N_Nodes+cols[j],    -DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(F1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(B, 	i,    3*N_Nodes+cols[j],    DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(F1, i, &numberOfNonZeros, &cols, &values);
    }
    std::cout <<"BF1 done" << std::endl;
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

        MatSetValue(B, i, i, 1.0, ADD_VALUES);
        MatSetValue(B, N_Nodes+i, N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(B, 2*N_Nodes+i, 2*N_Nodes+i, 1.0, ADD_VALUES);


        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(A, 	N_Nodes+i,     N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);

        /*
        MatGetRow(C, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	3*N_Nodes+i,     N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(C2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C2, i, &numberOfNonZeros, &cols, &values);
        */
        MatGetRow(D1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	2*N_Nodes+i,    cols[j],    0.5*DeltaT*dummy*sin(gamma), 	ADD_VALUES);
                MatSetValue(B, 	2*N_Nodes+i,     cols[j],    -0.5*DeltaT*dummy*sin(gamma), 	ADD_VALUES);
                MatSetValue(A, 	2*N_Nodes+i,     N_Nodes+cols[j],    0.5*DeltaT*dummy*cos(gamma), 	ADD_VALUES);
                MatSetValue(B, 	2*N_Nodes+i,     N_Nodes+cols[j],    -0.5*DeltaT*dummy*cos(gamma), 	ADD_VALUES);
            }
        }
        MatRestoreRow(D1, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy*sin(gamma), 	ADD_VALUES);
                MatSetValue(B, 	i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy*sin(gamma), 	ADD_VALUES);
                MatSetValue(A, 	N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy*cos(gamma), 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy*cos(gamma), 	ADD_VALUES);

                //MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);

        /*
        MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            //if (dummy!=0.0)
            if (abs(dummy)>1e-10)
            {
                //MatSetValue(A, 	2*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	2*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                //MatSetValue(A, 	3*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	3*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);
        */
        MatGetRow(ALaplacian_h, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     2*N_Nodes+cols[j],     0.5*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(ALaplacian_h, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(ALaplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     	3*N_Nodes+cols[j],     dummy, 		ADD_VALUES);
            }
        }
        MatRestoreRow(ALaplacian, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(AF, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(B, 	3*N_Nodes+i,     	3*N_Nodes+cols[j],     dummy, 		ADD_VALUES);
            }
        }
        MatRestoreRow(AF, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(ADivLaplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     	cols[j],     0.5*dummy, 		ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     	cols[j],     -0.5*dummy, 		ADD_VALUES);
            }
        }
        MatRestoreRow(ADivLaplacian, i, &numberOfNonZeros, &cols, &values);

    }
    MatDestroy(&BF1);
    MatDestroy(&F1);
    MatDestroy(&C);
    MatDestroy(&C2);
    MatDestroy(&D1);
    MatDestroy(&D2);
    MatDestroy(&Laplacian);
    MatDestroy(&ADivLaplacian);
    MatDestroy(&ALaplacian);
    MatDestroy(&ALaplacian_h);
    MatDestroy(&AF);

    std::cout << " Global Matrices Constructed" << std::endl;
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
}
/*--------------------------------------------------------------------------*/
extern void create_EB_System_MidPoint(const Mat &E, const Mat &ET, const Mat &invM, const Mat &invM_small, const Mat &M1, const Mat &M1_small, const Mat &M2, const Mat &M2_small, const Mat &NMat, const Mat &NDerivMat, const unsigned int &N_Nodes, const unsigned int &N, const double &DeltaT, const double &nu, Mat &A, Mat &B, Mat &DIV, const double &Re, const double &Fr)
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
    MatMatMult(E, invM_small, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP1);
    MatMatMult(BF2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &BF2_TEMP2);
    MatMatMult(invM, BF2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &BF2);
    MatDestroy(&BF2_TEMP1);
    MatDestroy(&BF2_TEMP2);
    MatDestroy(&BF2);
    */

    // Not needed since we scale both sides of the equation (= constructed incompressibility)
    Mat DIV_TEMP1, DIV_TEMP2;
    MatMatMult(ET, invM, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP1);
    MatMatMult(DIV_TEMP1, M1, MAT_INITIAL_MATRIX, fillBF, &DIV_TEMP2);
    MatMatMult(invM_small, DIV_TEMP2, MAT_INITIAL_MATRIX, fillBF, &DIV);
    MatDestroy(&DIV_TEMP1);

    Mat C, C_TEMP1, C_TEMP2;
    MatMatMult(NMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &C_TEMP1);
    MatMatMult(C_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &C_TEMP2);
    MatMatMult(invM_small, C_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C);
    MatDestroy(&C_TEMP1);
    MatDestroy(&C_TEMP2);

    Mat C2, C2_TEMP1, C2_TEMP2;
    MatMatMult(NMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP1);
    MatMatMult(C2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &C2_TEMP2);
    MatMatMult(invM_small, C2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &C2);
    MatDestroy(&C2_TEMP1);
    MatDestroy(&C2_TEMP2);

    Mat D1, D1_TEMP1, D1_TEMP2;
    MatMatMult(NDerivMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP1);
    MatMatMult(D1_TEMP1, M1_small, MAT_INITIAL_MATRIX, fillBF, &D1_TEMP2);
    MatMatMult(invM_small, D1_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D1);
    MatDestroy(&D1_TEMP1);
    MatDestroy(&D1_TEMP2);

    Mat D2, D2_TEMP1, D2_TEMP2;
    MatMatMult(NDerivMat, invM_small, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP1);
    MatMatMult(D2_TEMP1, M2_small, MAT_INITIAL_MATRIX, fillBF, &D2_TEMP2);
    MatMatMult(invM_small, D2_TEMP2, MAT_INITIAL_MATRIX, fillBF, &D2);
    MatDestroy(&D2_TEMP1);
    MatDestroy(&D2_TEMP2);
	//MatScale(D2, 1.0/Fr/Fr);

    Mat Identity3;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, N_Nodes, Np, NULL, &Identity3);
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(Identity3, 	N_Nodes+i,     cols[j],     dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);
    }
    MatAssemblyBegin(Identity3,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Identity3,MAT_FINAL_ASSEMBLY);


    Mat Laplacian;
	MatMatMult(DIV, BF1, MAT_INITIAL_MATRIX, 1, &Laplacian);
	MatScale(Laplacian, nu*1.0/Re); // nu = 0 or 1

    //std::cout << "Laplacian = " << std::endl;
    //MatView(Laplacian, PETSC_VIEWER_STDOUT_SELF);

    Mat VectorLaplacian;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*N_Nodes, 2*N_Nodes, 10*6*Np+1+3*Np*Np, NULL, &VectorLaplacian);
    for (unsigned int i = 0; i < N_Nodes; i++)
    {
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(VectorLaplacian, 	i,     cols[j],     dummy, 	ADD_VALUES);
                MatSetValue(VectorLaplacian, 	N_Nodes+i,     N_Nodes+cols[j],     dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
    }
    MatAssemblyBegin(VectorLaplacian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(VectorLaplacian,MAT_FINAL_ASSEMBLY);

    //std::cout << "VectorLaplacian = " << std::endl;
    //MatView(VectorLaplacian, PETSC_VIEWER_STDOUT_SELF);

    // invM_small Not needed since we scale both sides of the equation (= constructed incompressibility)
    Mat ALaplacian, ALaplacian_h, ADivLaplacian;
	MatMatMult(DIV_TEMP2, VectorLaplacian, MAT_INITIAL_MATRIX, 1, &ADivLaplacian);
	MatMatMult(DIV_TEMP2, BF1, MAT_INITIAL_MATRIX, 1, &ALaplacian);
	MatMatMult(DIV_TEMP2, Identity3, MAT_INITIAL_MATRIX, 1, &ALaplacian_h);
    MatDestroy(&DIV_TEMP2);
    MatDestroy(&Identity3);
    MatDestroy(&VectorLaplacian);

    //std::cout << "ADivLaplacian = " << std::endl;
    //MatView(ADivLaplacian, PETSC_VIEWER_STDOUT_SELF);

    /// Pass DIV_TEMP2 instead of DIV to compute L infinity norm of divergence


    // factor 4 = Number of Variables
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+9*Np*Np,  NULL, &A); // Check Factor 10 // Change from Triangles to Quadrilaterals
    MatCreateSeqAIJ(PETSC_COMM_WORLD, 4*N_Nodes, 4*N_Nodes, 10*6*Np+1+9*Np*Np,  NULL, &B);
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
            //if (abs(dummy)>1e-10)
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     3*N_Nodes+cols[j],    -DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);
        /*
        MatGetRow(BF2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,    3*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,    2*N_Nodes+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF2, i, &numberOfNonZeros, &cols, &values);
        */
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

        MatSetValue(B, i, i, 1.0, ADD_VALUES);
        MatSetValue(B, N_Nodes+i, N_Nodes+i, 1.0, ADD_VALUES);
        MatSetValue(B, 2*N_Nodes+i, 2*N_Nodes+i, 1.0, ADD_VALUES);


        MatGetRow(Laplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(A, 	N_Nodes+i,     N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(Laplacian, i, &numberOfNonZeros, &cols, &values);

        /*
        MatGetRow(C, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	3*N_Nodes+i,     N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(C2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C2, i, &numberOfNonZeros, &cols, &values);
        */
        MatGetRow(D1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	2*N_Nodes+i,     N_Nodes+cols[j],    0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	2*N_Nodes+i,     N_Nodes+cols[j],    -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D1, i, &numberOfNonZeros, &cols, &values);
        MatGetRow(D2, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                MatSetValue(A, 	N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                MatSetValue(B, 	N_Nodes+i,     2*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                //MatSetValue(A, 	N_Nodes+i,     3*N_Nodes+cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	N_Nodes+i,     3*N_Nodes+cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(D2, i, &numberOfNonZeros, &cols, &values);

        /*
        MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            //if (dummy!=0.0)
            if (abs(dummy)>1e-10)
            {
                //MatSetValue(A, 	2*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	2*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);

                //MatSetValue(A, 	3*N_Nodes+i,     cols[j],     -0.5*DeltaT*dummy, 	ADD_VALUES);
                //MatSetValue(B, 	3*N_Nodes+i,     cols[j],     0.5*DeltaT*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);
        */
        MatGetRow(ALaplacian_h, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     2*N_Nodes+cols[j],     0.5*dummy, 	ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     2*N_Nodes+cols[j],     -0.5*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(ALaplacian_h, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(ALaplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     	3*N_Nodes+cols[j],     dummy, 		ADD_VALUES);
            }
        }
        MatRestoreRow(ALaplacian, i, &numberOfNonZeros, &cols, &values);

        MatGetRow(ADivLaplacian, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(A, 	3*N_Nodes+i,     	cols[j],     0.5*dummy, 		ADD_VALUES);
                MatSetValue(B, 	3*N_Nodes+i,     	cols[j],     -0.5*dummy, 		ADD_VALUES);
            }
        }
        MatRestoreRow(ADivLaplacian, i, &numberOfNonZeros, &cols, &values);

    }
    MatDestroy(&BF1);
    MatDestroy(&C);
    MatDestroy(&C2);
    MatDestroy(&D1);
    MatDestroy(&D2);
    MatDestroy(&Laplacian);
    MatDestroy(&ADivLaplacian);
    MatDestroy(&ALaplacian);
    MatDestroy(&ALaplacian_h);

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
    //PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);
    std::cout << "Computing Initial Condition " << std::endl;
    // Initial Condition
    // Size = sum_i Np_i, i = 1 .. Nel
    //VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, 3*N_Nodes,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecU);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecW);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecR);
    VecCreateSeq(PETSC_COMM_WORLD, N_Nodes, &VecP);

    char szFileName[255] = {0};
    FILE *f = fopen("Solution/Coordinates.txt", "w");
    fprintf(f, "n \t pos \t xCoor \t zCoor \t p value \n");
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int ID = (*k).getID();
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
            VecSetValue(Initial_Condition, N_Nodes + pos + n, value, INSERT_VALUES);

            //value = Exact_Solution_r_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecR, pos + n, value, INSERT_VALUES);
            //VecSetValue(Initial_Condition, 2*N_Nodes + pos + n, value, INSERT_VALUES);

            value = Exact_Solution_p_2D_system1(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecP, pos + n, value, INSERT_VALUES);
            //VecSetValue(Initial_Condition, 3*N_Nodes + pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, 2*N_Nodes + pos + n, value, INSERT_VALUES);
            //std::cout << "ID = " << ID << ", n = " << n << ", pos = " << pos << ", xCoor[n] = " << xCoor[n] << ", zCoor[n] = " << zCoor[n] << ", p = " << value << std::endl;

            //fprintf(f, "%1u \t %u \t %1.16e \t %1.16e \t %1.16e \n", n, pos, xCoor[n], zCoor[n], value);
        }

    }
    fclose(f);
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
extern void compute_InitialCondition_system2(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecU, Vec &VecW, Vec &VecR, Vec &VecP, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc)
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

    char szFileName[255] = {0};


    std::string store_solution = "Solution/Solutions/Coordinates_n"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+"N"+std::to_string(N_Petsc)+"Ts"+std::to_string(Number_Of_TimeSteps_In_One_Period)+".txt";
    const char *store_solution_char = store_solution.c_str();
    //FILE *f = fopen("Solution/Coordinates.txt", "w");
    FILE *f = fopen(store_solution_char, "w");
    fprintf(f, "n \t pos \t xCoor \t zCoor \t p value \n");
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int ID = (*k).getID();
        unsigned int pos = (*k).getPosition();
        std::vector<double> xCoor, zCoor;
        xCoor = (*k).get_node_coordinates_x();
        zCoor = (*k).get_node_coordinates_y();
        double t = 0;
        for (unsigned int n = 0; n < Np; n++)
        {
            double value = Exact_Solution_mx_2D_system2(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecU, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, pos + n, value, INSERT_VALUES);

            value = Exact_Solution_mz_2D_system2(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecW, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, N_Nodes + pos + n, value, INSERT_VALUES);

            value = Exact_Solution_r_2D_system2(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecR, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, 2*N_Nodes + pos + n, value, INSERT_VALUES);

            value = Exact_Solution_p_2D_system2(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(VecP, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, 3*N_Nodes + pos + n, value, INSERT_VALUES);
            //std::cout << "ID = " << ID << ", n = " << n << ", pos = " << pos << ", xCoor[n] = " << xCoor[n] << ", zCoor[n] = " << zCoor[n] << ", p = " << value << std::endl;

            fprintf(f, "%1u \t %u \t %1.16e \t %1.16e \t %1.16e \n", n, pos, xCoor[n], zCoor[n], value);
        }

    }
    fclose(f);
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
/*--------------------------------------------------------------------------*//*--------------------------------------------------------------------------*/
extern void compute_InitialCondition_Incompressible(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecNodes, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc, const Mat &DIV)
{
    PetscScalar   sigma;
    sigma = calculate_sigma_2DIC(rho_0_Deriv, kxmode, kzmode);
    PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);
    std::cout << "Computing Initial Condition " << std::endl;
    // Initial Condition
    // Size = sum_i Np_i, i = 1 .. Nel
    Vec VecNodesVel;
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&VecNodes);
    VecCreateSeq(PETSC_COMM_WORLD, 2*N_Nodes,&VecNodesVel);
    Vec Velocity;
    VecCreateSeq(PETSC_COMM_WORLD, 2*N_Nodes, &Velocity);

    char szFileName[255] = {0};


    std::string store_solution = "Solution/Solutions/Coordinates_n"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+"N"+std::to_string(N_Petsc)+"Ts"+std::to_string(Number_Of_TimeSteps_In_One_Period)+".txt";
    const char *store_solution_char = store_solution.c_str();
    //FILE *f = fopen("Solution/Coordinates.txt", "w");
    FILE *f = fopen(store_solution_char, "w");
    fprintf(f, "n \t pos \t xCoor \t zCoor \t p value \n");
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int ID = (*k).getID();
        unsigned int pos = (*k).getPosition();
        std::vector<double> xCoor, zCoor;
        xCoor = (*k).get_node_coordinates_x();
        zCoor = (*k).get_node_coordinates_y();
        double t = 0;
        for (unsigned int n = 0; n < Np; n++)
        {
            VecSetValue(VecNodesVel, pos + n, -zCoor[n],  INSERT_VALUES);
            VecSetValue(VecNodesVel, N_Nodes + pos + n, xCoor[n], INSERT_VALUES);
            VecSetValue(VecNodes, pos + n, -zCoor[n],  INSERT_VALUES);
            VecSetValue(VecNodes, N_Nodes + pos + n, xCoor[n], INSERT_VALUES);
            VecSetValue(VecNodes, 2*N_Nodes + pos + n, 0, INSERT_VALUES);
            VecSetValue(VecNodes, 3*N_Nodes + pos + n, 0, INSERT_VALUES);

            double value = 0;//Exact_Solution_mx_2DIC(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(Velocity, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, pos + n, value, INSERT_VALUES);

            //value = Exact_Solution_mz_2DIC(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(Velocity, N_Nodes+pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, N_Nodes + pos + n, value, INSERT_VALUES);

            //value = Exact_Solution_r_2DIC(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(Initial_Condition, 2*N_Nodes + pos + n, value, INSERT_VALUES);

            //value = Exact_Solution_p_2DIC(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode);
            VecSetValue(Initial_Condition, 3*N_Nodes + pos + n, value, INSERT_VALUES);
            //std::cout << "ID = " << ID << ", n = " << n << ", pos = " << pos << ", xCoor[n] = " << xCoor[n] << ", zCoor[n] = " << zCoor[n] << ", p = " << value << std::endl;

            fprintf(f, "%1u \t %u \t %1.16e \t %1.16e \t %1.16e \n", n, pos, xCoor[n], zCoor[n], value);
        }

    }
    fclose(f);
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    VecAssemblyBegin(Velocity);
    VecAssemblyEnd(Velocity);
    VecAssemblyBegin(VecNodes);
    VecAssemblyEnd(VecNodes);
    VecAssemblyBegin(VecNodesVel);
    VecAssemblyEnd(VecNodesVel);

    compute_Divergence_Velocity(Initial_Condition, N_Nodes, DIV);
    correctInitialProjectionOfVelocity(N_Nodes, Velocity, DIV);

	PetscScalar* XTEMP;
	VecGetArray(Velocity, &XTEMP);
    PetscInt ix[2*N_Nodes];
    for (unsigned int k=0;k<2*N_Nodes; k++)
    {
        ix[k] = k;
    }
    VecSetValues(Initial_Condition, 2*N_Nodes, ix, XTEMP , INSERT_VALUES);
	VecRestoreArray(Velocity, &XTEMP);
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    compute_Divergence_Velocity(Initial_Condition, N_Nodes, DIV);
    VecDestroy(&Velocity);

    std::cout << "Forcing" << std::endl;
    compute_Divergence_Velocity(VecNodes, N_Nodes, DIV);
    correctInitialProjectionOfVelocity(N_Nodes, VecNodesVel, DIV);
	PetscScalar* XTEMP2;
	VecGetArray(VecNodesVel, &XTEMP2);
    VecSetValues(VecNodes, 2*N_Nodes, ix, XTEMP2, INSERT_VALUES);
	VecRestoreArray(VecNodesVel, &XTEMP2);
    VecAssemblyBegin(VecNodes);
    VecAssemblyEnd(VecNodes);
    compute_Divergence_Velocity(VecNodes, N_Nodes, DIV);
    VecDestroy(&VecNodesVel);
}
/*--------------------------------------------------------------------------*/
extern void ComputeNodesForcing(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, Vec &VecNodes)
{
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&VecNodes);
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int ID = (*k).getID();
        unsigned int pos = (*k).getPosition();
        std::vector<double> xCoor, zCoor;
        xCoor = (*k).get_node_coordinates_x();
        zCoor = (*k).get_node_coordinates_y();
        for (unsigned int n = 0; n < Np; n++)
        {
            VecSetValue(VecNodes, pos + n, zCoor[n],  INSERT_VALUES);
            VecSetValue(VecNodes, N_Nodes + pos + n, -xCoor[n], INSERT_VALUES);
            VecSetValue(VecNodes, 2*N_Nodes + pos + n, 0.0, INSERT_VALUES);
            VecSetValue(VecNodes, 3*N_Nodes + pos + n, 0, INSERT_VALUES);
        }
    }

    VecAssemblyBegin(VecNodes);
    VecAssemblyEnd(VecNodes);
}/*--------------------------------------------------------------------------*/
extern void ComputeForcing(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, Vec &Forcing_a)
{
    VecCreateSeq(PETSC_COMM_WORLD, 2*N_Nodes,&Forcing_a);
    double x0 = 1;
    double z0 = 0.75/2.0;
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int ID = (*k).getID();
        unsigned int pos = (*k).getPosition();
        std::vector<double> xCoor, zCoor;
        xCoor = (*k).get_node_coordinates_x();
        zCoor = (*k).get_node_coordinates_y();
        for (unsigned int n = 0; n < Np; n++)
        {
            VecSetValue(Forcing_a, pos + n, (zCoor[n])/0.75,  INSERT_VALUES);
            VecSetValue(Forcing_a, N_Nodes + pos + n, - (xCoor[n]-x0), INSERT_VALUES);
        }
    }

    VecAssemblyBegin(Forcing_a);
    VecAssemblyEnd(Forcing_a);
}
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition_WA(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &Fr, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecNodes, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc, const Mat &DIV)
{
    PetscScalar   sigma;
    sigma = calculate_sigma_2DWA(rho_0_Deriv, kxmode, kzmode, Fr);
    //PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);
    std::cout << "Computing Initial Condition " << std::endl;
    // Initial Condition
    // Size = sum_i Np_i, i = 1 .. Nel
    Vec VecNodesVel;
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&VecNodes);
    VecCreateSeq(PETSC_COMM_WORLD, 2*N_Nodes,&VecNodesVel);
    Vec Velocity;
    VecCreateSeq(PETSC_COMM_WORLD, 2*N_Nodes, &Velocity);

    char szFileName[255] = {0};


    std::string store_solution = "Solution/Coordinates_n"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+"N"+std::to_string(N_Petsc)+"Ts"+std::to_string(Number_Of_TimeSteps_In_One_Period)+".txt";
    const char *store_solution_char = store_solution.c_str();
    //FILE *f = fopen("Solution/Coordinates.txt", "w");
    FILE *f = fopen(store_solution_char, "w");
    fprintf(f, "n \t pos \t xCoor \t zCoor \t p value \n");
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int ID = (*k).getID();
        unsigned int pos = (*k).getPosition();
        std::vector<double> xCoor, zCoor;
        xCoor = (*k).get_node_coordinates_x();
        zCoor = (*k).get_node_coordinates_y();
        double t = 0;
        /*
        double xleft = 0.0;
        if (Np > 1)
        {
            xleft = 1.0e-12;
        }
        else
        {
            xleft = 0.1;
        }
        */
        for (unsigned int n = 0; n < Np; n++)
        {
            VecSetValue(VecNodesVel, pos + n, 0.0,  INSERT_VALUES);
            VecSetValue(VecNodesVel, N_Nodes + pos + n, 0.0, INSERT_VALUES);
            //VecSetValue(VecNodes, pos + n, 0.0,  INSERT_VALUES);
            //VecSetValue(VecNodes, N_Nodes + pos + n, 0.0, INSERT_VALUES);
            VecSetValue(VecNodes, pos + n, 1.0,  INSERT_VALUES);
            VecSetValue(VecNodes, N_Nodes + pos + n, 0, INSERT_VALUES);
            //std::cout << xCoor[n] << std::endl;

            VecSetValue(VecNodes, 2*N_Nodes + pos + n, 0.0, INSERT_VALUES);
            /*
            if (xCoor[n] < xleft)
            {
                VecSetValue(VecNodes, 2*N_Nodes + pos + n, sin(2.0*PETSC_PI*zCoor[n]/0.75), INSERT_VALUES);
            }
            else
            {
                VecSetValue(VecNodes, 2*N_Nodes + pos + n, 0, INSERT_VALUES);
            }
            */
            VecSetValue(VecNodes, 3*N_Nodes + pos + n, 0, INSERT_VALUES);


            double value = 0.0;//Exact_Solution_mx_2DWA(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode, Fr);
            VecSetValue(Velocity, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, pos + n, value, INSERT_VALUES);

            //value = Exact_Solution_mz_2DWA(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode, Fr);
            VecSetValue(Velocity, N_Nodes+pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, N_Nodes + pos + n, value, INSERT_VALUES);

            //value = Exact_Solution_r_2DWA(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode, Fr);
            VecSetValue(Initial_Condition, 2*N_Nodes + pos + n, value, INSERT_VALUES);

            //value = Exact_Solution_p_2DWA(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode, Fr);
            VecSetValue(Initial_Condition, 3*N_Nodes + pos + n, value, INSERT_VALUES);
            //std::cout << "ID = " << ID << ", n = " << n << ", pos = " << pos << ", xCoor[n] = " << xCoor[n] << ", zCoor[n] = " << zCoor[n] << ", p = " << value << std::endl;

            fprintf(f, "%1u \t %u \t %1.16e \t %1.16e \t %1.16e \n", n, pos, xCoor[n], zCoor[n], value);
        }

    }
    fclose(f);
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    VecAssemblyBegin(Velocity);
    VecAssemblyEnd(Velocity);
    VecAssemblyBegin(VecNodes);
    VecAssemblyEnd(VecNodes);
    VecAssemblyBegin(VecNodesVel);
    VecAssemblyEnd(VecNodesVel);

    compute_Divergence_Velocity(Initial_Condition, N_Nodes, DIV);
    correctInitialProjectionOfVelocity(N_Nodes, Velocity, DIV);

	PetscScalar* XTEMP;
	VecGetArray(Velocity, &XTEMP);
    PetscInt ix[2*N_Nodes];
    for (unsigned int k=0;k<2*N_Nodes; k++)
    {
        ix[k] = k;
    }
    VecSetValues(Initial_Condition, 2*N_Nodes, ix, XTEMP , INSERT_VALUES);
	VecRestoreArray(Velocity, &XTEMP);
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    compute_Divergence_Velocity(Initial_Condition, N_Nodes, DIV);
    VecDestroy(&Velocity);
    /*
    std::cout << "Forcing" << std::endl;
    compute_Divergence_Velocity(VecNodes, N_Nodes, DIV);
    correctInitialProjectionOfVelocity(N_Nodes, VecNodesVel, DIV);
	PetscScalar* XTEMP2;
	VecGetArray(VecNodesVel, &XTEMP2);
    VecSetValues(VecNodes, 2*N_Nodes, ix, XTEMP2, INSERT_VALUES);
	VecRestoreArray(VecNodesVel, &XTEMP2);
    VecAssemblyBegin(VecNodes);
    VecAssemblyEnd(VecNodes);
    compute_Divergence_Velocity(VecNodes, N_Nodes, DIV);
    */
    VecDestroy(&VecNodesVel);
}
/*--------------------------------------------------------------------------*/
extern void compute_InitialCondition_EB(const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const double &rho_0_Deriv, const double &Fr, const double &kxmode, const double &kzmode, Vec &Initial_Condition, Vec &VecNodes, const unsigned int &Number_Of_Elements_Petsc, const unsigned int &Number_Of_TimeSteps_In_One_Period, const unsigned int &N_Petsc, const Mat &DIV)
{
    PetscScalar   sigma;
    sigma = calculate_sigma_2DEB(rho_0_Deriv, kxmode, kzmode, Fr);
    //PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);
    std::cout << "Computing Initial Condition " << std::endl;
    // Initial Condition
    // Size = sum_i Np_i, i = 1 .. Nel
    Vec VecNodesVel;
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&Initial_Condition);
    VecCreateSeq(PETSC_COMM_WORLD, 4*N_Nodes,&VecNodes);
    VecCreateSeq(PETSC_COMM_WORLD, 2*N_Nodes,&VecNodesVel);
    Vec Velocity;
    VecCreateSeq(PETSC_COMM_WORLD, 2*N_Nodes, &Velocity);

    char szFileName[255] = {0};


    std::string store_solution = "Solution/Coordinates_n"+std::to_string(Number_Of_Elements_Petsc)+"x"+std::to_string(Number_Of_Elements_Petsc)+"N"+std::to_string(N_Petsc)+"Ts"+std::to_string(Number_Of_TimeSteps_In_One_Period)+".txt";
    const char *store_solution_char = store_solution.c_str();
    //FILE *f = fopen("Solution/Coordinates.txt", "w");
    FILE *f = fopen(store_solution_char, "w");
    fprintf(f, "n \t pos \t xCoor \t zCoor \t p value \n");
    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int ID = (*k).getID();
        unsigned int pos = (*k).getPosition();
        std::vector<double> xCoor, zCoor;
        xCoor = (*k).get_node_coordinates_x();
        zCoor = (*k).get_node_coordinates_y();
        double t = 0;
        for (unsigned int n = 0; n < Np; n++)
        {
            VecSetValue(VecNodesVel, pos + n, -zCoor[n],  INSERT_VALUES);
            VecSetValue(VecNodesVel, N_Nodes + pos + n, xCoor[n], INSERT_VALUES);
            VecSetValue(VecNodes, pos + n, -zCoor[n],  INSERT_VALUES);
            VecSetValue(VecNodes, N_Nodes + pos + n, xCoor[n], INSERT_VALUES);
            VecSetValue(VecNodes, 2*N_Nodes + pos + n, 0, INSERT_VALUES);
            VecSetValue(VecNodes, 3*N_Nodes + pos + n, 0, INSERT_VALUES);

            double value = Exact_Solution_mx_2DEB(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode, Fr);
            VecSetValue(Velocity, pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, pos + n, value, INSERT_VALUES);

            value = Exact_Solution_mz_2DEB(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode, Fr);
            VecSetValue(Velocity, N_Nodes+pos + n, value, INSERT_VALUES);
            VecSetValue(Initial_Condition, N_Nodes + pos + n, value, INSERT_VALUES);

            value = Exact_Solution_r_2DEB(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode, Fr);
            VecSetValue(Initial_Condition, 2*N_Nodes + pos + n, value, INSERT_VALUES);

            value = Exact_Solution_p_2DEB(xCoor[n], zCoor[n], t, rho_0_Deriv, sigma, kxmode, kzmode, Fr);
            VecSetValue(Initial_Condition, 3*N_Nodes + pos + n, value, INSERT_VALUES);
            //std::cout << "ID = " << ID << ", n = " << n << ", pos = " << pos << ", xCoor[n] = " << xCoor[n] << ", zCoor[n] = " << zCoor[n] << ", p = " << value << std::endl;

            fprintf(f, "%1u \t %u \t %1.16e \t %1.16e \t %1.16e \n", n, pos, xCoor[n], zCoor[n], value);
        }

    }
    fclose(f);
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    VecAssemblyBegin(Velocity);
    VecAssemblyEnd(Velocity);
    VecAssemblyBegin(VecNodes);
    VecAssemblyEnd(VecNodes);
    VecAssemblyBegin(VecNodesVel);
    VecAssemblyEnd(VecNodesVel);

    compute_Divergence_Velocity(Initial_Condition, N_Nodes, DIV);
    correctInitialProjectionOfVelocity(N_Nodes, Velocity, DIV);

	PetscScalar* XTEMP;
	VecGetArray(Velocity, &XTEMP);
    PetscInt ix[2*N_Nodes];
    for (unsigned int k=0;k<2*N_Nodes; k++)
    {
        ix[k] = k;
    }
    VecSetValues(Initial_Condition, 2*N_Nodes, ix, XTEMP , INSERT_VALUES);
	VecRestoreArray(Velocity, &XTEMP);
    VecAssemblyBegin(Initial_Condition);
    VecAssemblyEnd(Initial_Condition);
    compute_Divergence_Velocity(Initial_Condition, N_Nodes, DIV);
    VecDestroy(&Velocity);

    //std::cout << "Forcing" << std::endl;
    //compute_Divergence_Velocity(VecNodes, N_Nodes, DIV);
    correctInitialProjectionOfVelocity(N_Nodes, VecNodesVel, DIV);
	PetscScalar* XTEMP2;
	VecGetArray(VecNodesVel, &XTEMP2);
    VecSetValues(VecNodes, 2*N_Nodes, ix, XTEMP2, INSERT_VALUES);
	VecRestoreArray(VecNodesVel, &XTEMP2);
    VecAssemblyBegin(VecNodes);
    VecAssemblyEnd(VecNodes);
    //compute_Divergence_Velocity(VecNodes, N_Nodes, DIV);
    VecDestroy(&VecNodesVel);
}
/*--------------------------------------------------------------------------*/
extern void correctInitialProjectionOfVelocity(const unsigned int &N_Nodes, Vec &UInit, const Mat &DIV)
{
	double reltol = 1.e-16;
	double abstol = 1.e-16;

    Vec RHS, UCorrected;
	VecCreateSeq(PETSC_COMM_SELF, N_Nodes, &RHS);
	VecDuplicate(UInit, &UCorrected);

    MatMult(DIV, UInit, RHS);
    VecScale(RHS,-1);
    PetscInt a(10);
    PetscReal max;
    VecMax(RHS, &a, &max);
    //std::cout << "Discrete divergence OLD = " << std::setprecision (20) <<max<< std::endl;

	if (max > abstol)
	{
		KSP ksp;
		PC pc;
		KSPCreate(PETSC_COMM_SELF, &ksp);
		KSPSetOperators(ksp, DIV, DIV);
		KSPSetType(ksp, KSPLSQR);
		KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
		KSPGetPC(ksp,&pc);
		KSPSetFromOptions(ksp);
		PCSetType(pc, PCNONE);

		KSPSetTolerances(ksp, reltol, abstol, PETSC_DEFAULT, PETSC_DEFAULT);
		KSPSetUp(ksp);
		KSPSolve(ksp, RHS, UCorrected);

		KSPConvergedReason reason;
		KSPGetConvergedReason(ksp, &reason);
		if(reason < 0)
		{
			PetscInt its;
			KSPGetIterationNumber(ksp, &its);
			PetscReal rnorm;
			KSPGetResidualNorm(ksp, &rnorm);
			std::cout << "\t\tPetSc: Solving the linear system has failed. Reason code: "
            << reason << std::endl << "Check KSPConvergedReason for the reason" << std::endl
            << "\t\tPetSc: Residual after " << int(its) << " iterations : ";
			std::cout.setf(std::ios::scientific, std::ios::floatfield);
			std::cout.precision(4);
			std::cout << rnorm << std::endl;
			std::cout.setf(std::ios::fixed, std::ios::floatfield);
			std::cout.precision(5);
		}

        else
		{
            /*
			if(reason > 0)
			{
				PetscInt its;
				KSPGetIterationNumber(ksp, &its);
				PetscReal rnorm;
				KSPGetResidualNorm(ksp, &rnorm);
				std::cout << "\t\tPetSc: Solving the linear system has succeeded. Reason code: "
                << reason << std::endl << "Check KSPConvergedReason for the reason" << std::endl
                << "\t\tPetsc: Convergence in " << int(its) << " iterations : ";
				std::cout.setf(std::ios::scientific, std::ios::floatfield);
				std::cout.precision(4);
				std::cout << rnorm << std::endl;
				std::cout.setf(std::ios::fixed, std::ios::floatfield);
				std::cout.precision(5);
			}
			else
			{
				std::cout << "\t\tPetSc: Solving the linear system is still under way" << std::endl;
			}
			*/
		}

        a	= 10;
        max	= 0;
        //VecMax(UCorrected, &a, &max);
        //std::cout << "Difference between UInit and UCorrected "<<max<< std::endl;
			//[maxE=max(u+velocity);]
        //VecAXPY(UCorrected, 1, UInit);
        VecAXPY(UInit, 1, UCorrected);
            //   [error=DIV*u;]
        MatMult(DIV, UInit, RHS);
        VecMax(RHS, &a, &max);
        //std::cout << "Discrete divergence NEW = " << std::setprecision (20) <<max<< std::endl;
        KSPDestroy(&ksp);
	}
	else
	{
		//VecCopy(UInit,UCorrected);
	}

    std::cout << "Initial Velocity Projected"<<std::endl;

    VecDestroy(&RHS);
    VecDestroy(&UCorrected);

}
/*--------------------------------------------------------------------------*/
extern void Simulate(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Vec &Initial_Condition, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables)
{

    double H0 = calculate_Hamiltonian2D_full(M1_small, M2_small, Initial_Condition, List_Of_Elements, N_Nodes);
    std::cout << "Initial Energy = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Start Simulations " << std::endl;

    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    KSPSetUp(ksp);
    KSPSetTolerances(ksp, 1e-12, 1e-12, 1e30, PETSC_DEFAULT);

    //KSPSetType(ksp,KSPPREONLY);
    //KSPSetType(ksp,KSPCG);
    KSPSetType(ksp,KSPGMRES);
    //KSPSetType(ksp,KSPBCGS);
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

    KSPGetPC(ksp,&pc);
    //PCSetType(pc,PCLU);
    PCSetType(pc,PCILU);
    //PCSetType(pc,PCNONE);
    //PCSetType(pc,PCSOR);

    KSPSetFromOptions(ksp);

    Vec QX;
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Variables*N_Nodes, &Sol);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Variables*N_Nodes, &QX);
    VecCopy(Initial_Condition, Sol);
    double H1 = 0.0;

    //PetscPrintf(PETSC_COMM_SELF,"Size Global Matrices %6.4e\n",(double)sigma);
    // MatView(A, viewer_info);

    char szFileName[255] = {0};
    FILE *f = fopen("Solution/Energy.txt", "w");
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

            H1 = calculate_Hamiltonian2D_full(M1_small, M2_small, Sol, List_Of_Elements, N_Nodes);

            //std::cout << "Energy Diff= " << std::setprecision(16) << H1-Hold <<std::endl;
            Hold = H1;

           // std::cout << "QX = " << std::endl;
            //VecView(QX, PETSC_VIEWER_STDOUT_SELF);
            //std::cout << "Solution = " << std::endl;
            //VecView(Sol, PETSC_VIEWER_STDOUT_SELF);

        PetscViewer viewer2;
        sprintf(szFileName, "Solution/solution%d.txt", t);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
        VecView(Sol, viewer2);
        PetscViewerDestroy(&viewer2);


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
void compute_Divergence_Velocity(const Vec &Initial_Condition, const double &N_Nodes, const Mat &DIV)
{
    Vec RHS, UVel;
	VecCreateSeq(PETSC_COMM_SELF, N_Nodes, &RHS);
    PetscScalar* XTEMP;
	VecGetArray(Initial_Condition, &XTEMP);
    VecCreateSeqWithArray(PETSC_COMM_SELF, 2.0*N_Nodes, 2.0*N_Nodes, XTEMP, &UVel);
    VecAssemblyBegin(UVel);
    VecAssemblyEnd(UVel);
	VecRestoreArray(Initial_Condition, &XTEMP);
    MatMult(DIV, UVel, RHS);
    PetscInt a(10);
    PetscReal max;
    VecMax(RHS, &a, &max);
    std::cout << "Max Divergence = " << std::setprecision (20) <<max<< std::endl;
    VecDestroy(&UVel);
    VecDestroy(&RHS);
}
/*--------------------------------------------------------------------------*/
extern void Simulate_IC(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const Vec &VecNodes, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega)
{

    double M0;
    double H0 = calculate_Hamiltonian2D_IC(M1_small, M2_small, Initial_Condition, List_Of_Elements, N_Nodes, M0);
    std::cout << "Initial Energy = " << std::setprecision(16) << H0 << std::endl;

    compute_Divergence_Velocity(Initial_Condition, N_Nodes, DIV);

    std::cout << "Start Simulations " << std::endl;

    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    KSPSetUp(ksp);
    KSPSetTolerances(ksp, 1e-14, 1e-14, 1e30, PETSC_DEFAULT);

    //KSPSetType(ksp,KSPCG);
    KSPSetType(ksp,KSPGMRES);
    //KSPSetType(ksp,KSPBCGS);
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

    KSPGetPC(ksp,&pc);
    //PCSetType(pc,PCLU);
    PCSetType(pc,PCILU);
    //PCSetType(pc,PCNONE);
    //PCSetType(pc,PCSOR);

    KSPSetFromOptions(ksp);

    Vec QX;
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Variables*N_Nodes, &Sol);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Variables*N_Nodes, &QX);
    VecCopy(Initial_Condition, Sol);
    double H1 = H0;
    double M1 = M0;

    //PetscPrintf(PETSC_COMM_SELF,"Size Global Matrices %6.4e\n",(double)sigma);
    // MatView(A, viewer_info);

    char szFileName[255] = {0};
    FILE *f = fopen("Solution/Energy.txt", "w");
    // Solve Linear System
    std::cout << "Start Time Stepping" << std::endl;
    double time = 0.0;

    double Hold = H0;
    double Mold = M0;
    for (unsigned int t = 0; t < Number_Of_TimeSteps; t++) //
    {
        // Forcing at half time Step
        PetscScalar Forcing = 0.5*F0*sin(omega*(time+DeltaT/2.0));
        fprintf(f, "%1.16e \t %1.16e \t %1.16e\n", time, H1, M1);
        time = (t+1)*DeltaT;
            MatMult(B, Sol, QX);
            VecAXPY(QX, Forcing*DeltaT, VecNodes);


            KSPSolve(ksp, QX, Sol);

            H1 = calculate_Hamiltonian2D_IC(M1_small, M2_small, Sol, List_Of_Elements, N_Nodes, M1);

            //std::cout << "Energy Diff= " << std::setprecision(16) << H1-Hold <<std::endl;
            Hold = H1;
            Mold = M1;

           // std::cout << "QX = " << std::endl;
            //VecView(QX, PETSC_VIEWER_STDOUT_SELF);
            //std::cout << "Solution = " << std::endl;
            //VecView(Sol, PETSC_VIEWER_STDOUT_SELF);

        PetscViewer viewer2;
        sprintf(szFileName, "Solution/solution%d.txt", t);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
        VecView(Sol, viewer2);
        PetscViewerDestroy(&viewer2);


    }
    fprintf(f, "%1.16e \t %1.16e \t %1.16e\n", time, H1, M1);
    fclose(f);
    compute_Divergence_Velocity(Sol, N_Nodes, DIV);
    std::cout << "Final Time = " << time << std::endl;
    std::cout << "End Time Stepping" << std::endl;
    KSPDestroy(&ksp);
    VecDestroy(&QX);
    std::cout << "Initial Energy    = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Final Energy      = " << std::setprecision(16) << H1 << std::endl;
    std::cout << "Difference Energy = " << std::setprecision(16) << H1-H0 << std::endl;
}
/*--------------------------------------------------------------------------*/
extern void Simulate_WA(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const Vec &VecNodes, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega)
{
    double M0;
    double H0 = calculate_Hamiltonian2D_IC(M1_small, M2_small, Initial_Condition, List_Of_Elements, N_Nodes, M0);
    std::cout << "Initial Energy = " << std::setprecision(16) << H0 << std::endl;




    std::cout << "Start Simulations " << std::endl;

    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    KSPSetUp(ksp);
    KSPSetTolerances(ksp, 1e-14, 1e-14, 1e30, PETSC_DEFAULT);

    //KSPSetType(ksp,KSPPREONLY);
    //KSPSetType(ksp,KSPCG);
    KSPSetType(ksp,KSPGMRES);
    //KSPSetType(ksp,KSPBCGS);
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

    KSPGetPC(ksp,&pc);
    //PCSetType(pc,PCLU);
    PCSetType(pc,PCILU);
    //PCSetType(pc,PCNONE);
    //PCSetType(pc,PCSOR);

    KSPSetFromOptions(ksp);

    Vec QX;
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Variables*N_Nodes, &Sol);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Variables*N_Nodes, &QX);
    VecCopy(Initial_Condition, Sol);
    double H1 = H0;
    double M1 = M0;

    //PetscPrintf(PETSC_COMM_SELF,"Size Global Matrices %6.4e\n",(double)sigma);
    // MatView(A, viewer_info);

    //VecView(VecNodes, PETSC_VIEWER_STDOUT_SELF);
    char szFileName[255] = {0};
    FILE *f = fopen("Solution/Energy.txt", "w");
    // Solve Linear System
    std::cout << "Start Time Stepping" << std::endl;
    double time = 0.0;
    std::cout << "Forcing Amplitude = " << F0 << std::endl;


    double Hold = H0;
    double Mold = M0;
    for (unsigned int t = 0; t < Number_Of_TimeSteps; t++) //
    {
        // Forcing at half time Step
        PetscScalar Forcing = F0*sin(omega*(time+DeltaT/2.0));
        fprintf(f, "%1.16e \t %1.16e \t %1.16e\n", time, H1, M1);
        time = (t+1)*DeltaT;
            MatMult(B, Sol, QX);
            if (t < (2.0/3.0*Number_Of_TimeSteps))
            {
                VecAXPY(QX, Forcing*DeltaT, VecNodes);
            }

            KSPSolve(ksp, QX, Sol);

            H1 = calculate_Hamiltonian2D_IC(M1_small, M2_small, Sol, List_Of_Elements, N_Nodes, M1);

            //std::cout << "Energy Diff= " << std::setprecision(16) << H1-Hold <<std::endl;
            Hold = H1;
            Mold = M1;

           // std::cout << "QX = " << std::endl;
            //VecView(QX, PETSC_VIEWER_STDOUT_SELF);
            //std::cout << "Solution = " << std::endl;
            //VecView(Sol, PETSC_VIEWER_STDOUT_SELF);

        PetscViewer viewer2;
        sprintf(szFileName, "Solution/solution%d.txt", t);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
        VecView(Sol, viewer2);
        PetscViewerDestroy(&viewer2);


    }
    fprintf(f, "%1.16e \t %1.16e \t %1.16e\n", time, H1, M1);
    fclose(f);
    std::cout << "Final Time = " << time << std::endl;
    std::cout << "End Time Stepping" << std::endl;
    KSPDestroy(&ksp);
    VecDestroy(&QX);
    std::cout << "Initial Energy    = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Final Energy      = " << std::setprecision(16) << H1 << std::endl;
    std::cout << "Difference Energy = " << std::setprecision(16) << H1-H0 << std::endl;
    std::cout << "Initial " ; compute_Divergence_Velocity(Initial_Condition, N_Nodes, DIV);
    std::cout << "Final   "  ; compute_Divergence_Velocity(Sol, N_Nodes, DIV);
}
/*--------------------------------------------------------------------------*/
extern void Simulate_WA_Forced(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega)
{
    double M0 = 0.0;
    double H0 = calculate_Hamiltonian2D_IC(M1_small, M2_small, Initial_Condition, List_Of_Elements, N_Nodes, M0);
    std::cout << "Initial Energy = " << std::setprecision(16) << H0 << std::endl;

    std::cout << "Start Simulations " << std::endl;

    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    KSPSetUp(ksp);
    KSPSetTolerances(ksp, 1e-14, 1e-14, 1e30, PETSC_DEFAULT);

    //KSPSetType(ksp,KSPPREONLY);
    //KSPSetType(ksp,KSPCG);
    KSPSetType(ksp,KSPGMRES);
    //KSPSetType(ksp,KSPBCGS);
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

    KSPGetPC(ksp,&pc);
    //PCSetType(pc,PCLU);
    PCSetType(pc,PCILU);
    //PCSetType(pc,PCNONE);
    //PCSetType(pc,PCSOR);

    KSPSetFromOptions(ksp);

    Vec QX;
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Variables*N_Nodes, &Sol);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Variables*N_Nodes, &QX);
    VecCopy(Initial_Condition, Sol);
    double H1 = H0;
    double M1 = M0;

    //PetscPrintf(PETSC_COMM_SELF,"Size Global Matrices %6.4e\n",(double)sigma);
    // MatView(A, viewer_info);

    //VecView(VecNodes, PETSC_VIEWER_STDOUT_SELF);
    char szFileName[255] = {0};
    FILE *f = fopen("Solution/Energy.txt", "w");
    //fprintf(f, "%1.16e \t %1.16e \t %1.16e\n", 0.0, H0, M0);
    // Solve Linear System
    std::cout << "Start Time Stepping" << std::endl;
    double time = 0.0;
    std::cout << "Forcing Amplitude = " << F0 << std::endl;


    double Hold = H0;
    double Mold = M0;

    PetscInt ix[N_Nodes];
    for (unsigned int k=0;k<N_Nodes; k++)
        ix[k] = 3*N_Nodes+k;


    for (unsigned int t = 0; t < Number_Of_TimeSteps; t++)
    {
        time = (t)*DeltaT;//(t+1)*DeltaT;
        fprintf(f, "%1.16e \t %1.16e \t %1.16e\n", time, H1, M1);


        // Compute the Forcing term and replace the pressure variable from the previous time step (since these are not used)
        // Forcing at half time Step
        PetscScalar Forcing = 0.0;
        if (t < (2.0/3.0*Number_Of_TimeSteps))
            { Forcing = F0*sin(omega*(time+DeltaT/2.0));}
        PetscScalar Fx[N_Nodes];
        for (unsigned int k=0;k<N_Nodes; k++)
            {Fx[k] = Forcing;}
        VecSetValues(Sol,N_Nodes,ix,Fx,INSERT_VALUES);
        VecAssemblyBegin(Sol);
        VecAssemblyEnd(Sol);




        MatMult(B, Sol, QX);
        KSPSolve(ksp, QX, Sol);

        H1 = calculate_Hamiltonian2D_IC(M1_small, M2_small, Sol, List_Of_Elements, N_Nodes, M1);

        //std::cout << "Energy Diff= " << std::setprecision(16) << H1-Hold <<std::endl;
        Hold = H1;
        Mold = M1;

       // std::cout << "QX = " << std::endl;
        //VecView(QX, PETSC_VIEWER_STDOUT_SELF);
        //std::cout << "Solution = " << std::endl;
        //VecView(Sol, PETSC_VIEWER_STDOUT_SELF);
        PetscViewer viewer2;
        sprintf(szFileName, "Solution/solution%d.txt", t);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
        VecView(Sol, viewer2);
        PetscViewerDestroy(&viewer2);


    }
    fprintf(f, "%1.16e \t %1.16e \t %1.16e\n", time, H1, M1);
    fclose(f);
    std::cout << "Final Time = " << time << std::endl;
    std::cout << "End Time Stepping" << std::endl;
    KSPDestroy(&ksp);
    VecDestroy(&QX);
    std::cout << "Initial Energy    = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Final Energy      = " << std::setprecision(16) << H1 << std::endl;
    std::cout << "Difference Energy = " << std::setprecision(16) << H1-H0 << std::endl;
    std::cout << "Initial " ; compute_Divergence_Velocity(Initial_Condition, N_Nodes, DIV);
    std::cout << "Final   "  ; compute_Divergence_Velocity(Sol, N_Nodes, DIV);
}
/*--------------------------------------------------------------------------*/
extern void Simulate_EB(const Mat &A, const Mat &B, const Mat &M1_small, const Mat &M2_small, const Mat &DIV, const Vec &Initial_Condition, const Vec &VecNodes, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, const unsigned int &Number_Of_TimeSteps, const double &DeltaT, Vec &Sol, const unsigned int &Number_Of_Variables, const double &F0, const double &omega)
{
    double M0;
    double H0 = calculate_Hamiltonian2D_IC(M1_small, M2_small, Initial_Condition, List_Of_Elements, N_Nodes, M0);
    std::cout << "Initial Energy = " << std::setprecision(16) << H0 << std::endl;




    std::cout << "Start Simulations " << std::endl;

    KSP ksp;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    KSPSetUp(ksp);
    KSPSetTolerances(ksp, 1e-14, 1e-14, 1e30, PETSC_DEFAULT);

    //KSPSetType(ksp,KSPPREONLY);
    //KSPSetType(ksp,KSPCG);
    KSPSetType(ksp,KSPGMRES);
    //KSPSetType(ksp,KSPBCGS);
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

    KSPGetPC(ksp,&pc);
    //PCSetType(pc,PCLU);
    PCSetType(pc,PCILU);
    //PCSetType(pc,PCNONE);
    //PCSetType(pc,PCSOR);

    KSPSetFromOptions(ksp);

    Vec QX;
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Variables*N_Nodes, &Sol);
    VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Variables*N_Nodes, &QX);
    VecCopy(Initial_Condition, Sol);
    double H1 = H0;
    double M1 = M0;

    //PetscPrintf(PETSC_COMM_SELF,"Size Global Matrices %6.4e\n",(double)sigma);
    // MatView(A, viewer_info);

    //VecView(VecNodes, PETSC_VIEWER_STDOUT_SELF);
    char szFileName[255] = {0};
    FILE *f = fopen("Solution/Energy.txt", "w");
    // Solve Linear System
    std::cout << "Start Time Stepping" << std::endl;
    double time = 0.0;
    std::cout << "Forcing Amplitude = " << F0 << std::endl;


    double Hold = H0;
    double Mold = M0;
    for (unsigned int t = 0; t < Number_Of_TimeSteps; t++) //
    {
        // Forcing at half time Step
        PetscScalar Forcing = F0*sin(omega*(time+DeltaT/2.0));
        fprintf(f, "%1.16e \t %1.16e \t %1.16e\n", time, H1, M1);
        time = (t+1)*DeltaT;
            MatMult(B, Sol, QX);
            VecAXPY(QX, Forcing*DeltaT, VecNodes);


            KSPSolve(ksp, QX, Sol);

            H1 = calculate_Hamiltonian2D_IC(M1_small, M2_small, Sol, List_Of_Elements, N_Nodes, M1);

            //std::cout << "Energy Diff= " << std::setprecision(16) << H1-Hold <<std::endl;
            Hold = H1;
            Mold = M1;

           // std::cout << "QX = " << std::endl;
            //VecView(QX, PETSC_VIEWER_STDOUT_SELF);
            //std::cout << "Solution = " << std::endl;
            //VecView(Sol, PETSC_VIEWER_STDOUT_SELF);

        PetscViewer viewer2;
        sprintf(szFileName, "Solution/solution%d.txt", t);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
        VecView(Sol, viewer2);
        PetscViewerDestroy(&viewer2);


    }
    fprintf(f, "%1.16e \t %1.16e \t %1.16e\n", time, H1, M1);
    fclose(f);
    std::cout << "Final Time = " << time << std::endl;
    std::cout << "End Time Stepping" << std::endl;
    KSPDestroy(&ksp);
    VecDestroy(&QX);
    std::cout << "Initial Energy    = " << std::setprecision(16) << H0 << std::endl;
    std::cout << "Final Energy      = " << std::setprecision(16) << H1 << std::endl;
    std::cout << "Difference Energy = " << std::setprecision(16) << H1-H0 << std::endl;
    std::cout << "Initial " ; compute_Divergence_Velocity(Initial_Condition, N_Nodes, DIV);
    std::cout << "Final   "  ; compute_Divergence_Velocity(Sol, N_Nodes, DIV);
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
                //H += 0.5*massij*(XTemp[pos+i]*XTemp[pos+j]+XTemp[N_Nodes+pos+i]*XTemp[N_Nodes+pos+j]+XTemp[3*N_Nodes+pos+i]*XTemp[3*N_Nodes+pos+j]);
                H += 0.5*massij*(XTemp[pos+i]*XTemp[pos+j]+XTemp[N_Nodes+pos+i]*XTemp[N_Nodes+pos+j]+XTemp[2*N_Nodes+pos+i]*XTemp[2*N_Nodes+pos+j]);

            }
        }

        MatDestroy(&M1_Elemental);
        ISDestroy(&isrow);
    }
    VecRestoreArray(Solution, &XTemp);
    return H;
}
/*--------------------------------------------------------------------------*/
extern double calculate_Hamiltonian2D_full(const Mat &M1, const Mat &M2, const Vec &Solution, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes)
{
    double H = 0.0;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int pos = (*k).getPosition();

        // Get Submatrix from M1
        Mat M1_Elemental, M2_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, pos, 1, &isrow);
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
                H += 0.5*massij1*(XTemp[pos+i]*XTemp[pos+j] + XTemp[N_Nodes+pos+i]*XTemp[N_Nodes+pos+j]) // Kinetic Energy
                    + 0.5*massij1*(XTemp[3*N_Nodes+pos+i]*XTemp[3*N_Nodes+pos+j])
                    + 0.5*massij2*(XTemp[2*N_Nodes+pos+i]*XTemp[2*N_Nodes+pos+j]-XTemp[2*N_Nodes+pos+i]*XTemp[3*N_Nodes+pos+j]-XTemp[3*N_Nodes+pos+i]*XTemp[2*N_Nodes+pos+j]+XTemp[3*N_Nodes+pos+i]*XTemp[3*N_Nodes+pos+j]);

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
extern double calculate_Hamiltonian2D_IC(const Mat &M1, const Mat &M2, const Vec &Solution, const std::vector<Squares2D> &List_Of_Elements, const unsigned int &N_Nodes, double &M)
{
    double H = 0.0;
    PetscScalar *XTemp;
    VecGetArray(Solution, &XTemp);

    for (auto k = List_Of_Elements.begin(); k < List_Of_Elements.end(); k++)
    {
        unsigned int Np = (*k).get_Number_Of_Nodes();
        unsigned int pos = (*k).getPosition();

        // Get Submatrix from M1
        Mat M1_Elemental, M2_Elemental;
        IS isrow;
        ISCreateStride(PETSC_COMM_WORLD, Np, pos, 1, &isrow);
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
                H += 0.5*massij1*(XTemp[pos+i]*XTemp[pos+j] + XTemp[N_Nodes+pos+i]*XTemp[N_Nodes+pos+j]) // Kinetic Energy
                    + 0.5*massij2*(XTemp[2*N_Nodes+pos+i]*XTemp[2*N_Nodes+pos+j]); // Potential Energy
                M += massij1*XTemp[2*N_Nodes+pos+j]; //Total Mass (when Boussinesq Approximation)
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
    error *= sqrt(DeltaX*DeltaY); //0.5* for triangles
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

    for (auto e = List_Of_Elements.begin(); e < List_Of_Elements.end(); e++)
    {
        unsigned int Np = (*e).get_Number_Of_Nodes();
        //std::cout << "ID = " << (*e).getID() << std::endl;
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

        double sum = 0.0;
        double sum_p = 0.0;
        for (unsigned int i = 0; i < (N+1); i++)
        {
            for (unsigned int j = 0; j < (N+1); j++)
            {
                double ex_u = exact_a[pos+i*(N+1)+j];
                double sol_u = sol_a[pos+i*(N+1)+j];
                double ex_w = exact_a[N_Nodes+pos+i*(N+1)+j];
                double sol_w = sol_a[N_Nodes+pos+i*(N+1)+j];
                //double ex_r = exact_a[2*N_Nodes+pos+i*(N+1)+j];
                //double sol_r = sol_a[2*N_Nodes+pos+i*(N+1)+j];
                //double ex_p = exact_a[3*N_Nodes+pos+i*(N+1)+j];
                //double sol_p = sol_a[3*N_Nodes+pos+i*(N+1)+j];
                double ex_p = exact_a[2*N_Nodes+pos+i*(N+1)+j];
                double sol_p = sol_a[2*N_Nodes+pos+i*(N+1)+j];
                //double ex_p = exact_a[3*N_Nodes+pos+i*(N+1)+j];
                //double sol_p = sol_a[3*N_Nodes+pos+i*(N+1)+j];

                //std::cout << "n = " << (*e).getID() << ", J =  " << J << ", pos = " << pos << ", xCoor[n] = " << xCoor[i*(N+1)+j] << ", zCoor[n] = " << zCoor[i*(N+1)+j] << ", ex p = " << ex_p << ", sol p = " << sol_p << std::endl;

                exact_u += w[i]*w[j]*ex_u*ex_u*J;
                diff_u += w[i]*w[j]*(ex_u-sol_u)*(ex_u-sol_u)*J;
                exact_w += w[i]*w[j]*ex_w*ex_w*J;
                diff_w += w[i]*w[j]*(ex_w-sol_w)*(ex_w-sol_w)*J;
                //exact_r += w[i]*w[j]*ex_r*ex_r*J;
                //diff_r += w[i]*w[j]*(ex_r-sol_r)*(ex_r-sol_r)*J;
                exact_p += w[i]*w[j]*ex_p*ex_p*J;
                diff_p += w[i]*w[j]*(ex_p-sol_p)*(ex_p-sol_p)*J;

                //std::cout << "w[i] = " << w[i] << ", w[j] = " << w[j] << ", w[i]*w[j] = " << w[i]*w[j] << std::endl;
                //std::cout << "ex_p = " << ex_p << ", sol_p = " << sol_p << std::endl;
                //std::cout << "(ex_p-sol_p) = " << (ex_p-sol_p) << std::endl;
                sum += w[i]*w[j];
                sum_p += abs(ex_p-sol_p);
                //std::cout << "r[i] r[j] = " << qp[i] << " " << qp[j] << std::endl;
            }
        }
        //std::cout << "sum = " << sum << std::endl;
        //std::cout << "sum_p = " << sum_p << std::endl;
        //std::cout << "diff_u = " << diff_u << std::endl;
        //std::cout << "diff_w = " << diff_w << std::endl;
        //std::cout << "diff_p = " << diff_p << std::endl;
        //std::cout << std::endl;
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
    std::cout << "Error u = " << error_u << std::endl;
    std::cout << "Error w = " << error_w << std::endl;
    std::cout << "Error r = " << error_r << std::endl;
    std::cout << "Error p = " << error_p << std::endl;

    //return sqrt(error_u+error_w+error_r+error_p);
    //return error_p;
    return error_w;

}
/*--------------------------------------------------------------------------*/
extern double calculate_Error2D_Quad(const Vec &Exact, const Vec &Solution, const unsigned int &Norm_Type, const std::vector<Squares2D> &List_Of_Elements, const std::vector<VertexCoordinates2D> &List_Of_Vertices, const unsigned int &N_Nodes)
{
    std::cout << "Start Elemental Calculations " << std::endl;
    PetscScalar *exact_a, *sol_a;
    VecGetArray(Exact, &exact_a);
    VecGetArray(Solution, &sol_a);
    double error_u = 0, error_w = 0, error_r = 0, error_p = 0;

    for (auto e = List_Of_Elements.begin(); e < List_Of_Elements.end(); e++)
    {
        unsigned int Np = (*e).get_Number_Of_Nodes();
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

        double sum = 0.0;
        double sum_p = 0.0;
        for (unsigned int i = 0; i < (N+1); i++)
        {
            for (unsigned int j = 0; j < (N+1); j++)
            {
                double J, drdx, drdy, dsdx, dsdy, x, y;
                Calculate_Jacobian_Quadrilateral((*e), List_Of_Vertices, qp[j], qp[i], J, drdx, drdy, dsdx, dsdy, x, y);

                double ex_u = exact_a[pos+i*(N+1)+j];
                double sol_u = sol_a[pos+i*(N+1)+j];
                double ex_w = exact_a[N_Nodes+pos+i*(N+1)+j];
                double sol_w = sol_a[N_Nodes+pos+i*(N+1)+j];
                //double ex_r = exact_a[2*N_Nodes+pos+i*(N+1)+j];
                //double sol_r = sol_a[2*N_Nodes+pos+i*(N+1)+j];
                //double ex_p = exact_a[3*N_Nodes+pos+i*(N+1)+j];
                //double sol_p = sol_a[3*N_Nodes+pos+i*(N+1)+j];
                double ex_p = exact_a[2*N_Nodes+pos+i*(N+1)+j];
                double sol_p = sol_a[2*N_Nodes+pos+i*(N+1)+j];
                //double ex_p = exact_a[3*N_Nodes+pos+i*(N+1)+j];
                //double sol_p = sol_a[3*N_Nodes+pos+i*(N+1)+j];

                //std::cout << "n = " << (*e).getID() << ", J =  " << J << ", pos = " << pos << ", xCoor[n] = " << xCoor[i*(N+1)+j] << ", zCoor[n] = " << zCoor[i*(N+1)+j] << ", ex p = " << ex_p << ", sol p = " << sol_p << std::endl;

                exact_u += w[i]*w[j]*ex_u*ex_u*J;
                diff_u += w[i]*w[j]*(ex_u-sol_u)*(ex_u-sol_u)*J;
                exact_w += w[i]*w[j]*ex_w*ex_w*J;
                diff_w += w[i]*w[j]*(ex_w-sol_w)*(ex_w-sol_w)*J;
                //exact_r += w[i]*w[j]*ex_r*ex_r*J;
                //diff_r += w[i]*w[j]*(ex_r-sol_r)*(ex_r-sol_r)*J;
                exact_p += w[i]*w[j]*ex_p*ex_p*J;
                diff_p += w[i]*w[j]*(ex_p-sol_p)*(ex_p-sol_p)*J;

                //std::cout << "w[i] = " << w[i] << ", w[j] = " << w[j] << ", w[i]*w[j] = " << w[i]*w[j] << std::endl;
                //std::cout << "ex_p = " << ex_p << ", sol_p = " << sol_p << std::endl;
                //std::cout << "(ex_p-sol_p) = " << (ex_p-sol_p) << std::endl;
                sum += w[i]*w[j];
                sum_p += abs(ex_p-sol_p);
                //std::cout << "r[i] r[j] = " << qp[i] << " " << qp[j] << std::endl;
            }
        }
        //std::cout << "sum = " << sum << std::endl;
        //std::cout << "sum_p = " << sum_p << std::endl;
        //std::cout << "diff_u = " << diff_u << std::endl;
        //std::cout << "diff_w = " << diff_w << std::endl;
        //std::cout << "diff_p = " << diff_p << std::endl;
        //std::cout << std::endl;
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
    std::cout << "Error u = " << error_u << std::endl;
    std::cout << "Error w = " << error_w << std::endl;
    std::cout << "Error r = " << error_r << std::endl;
    std::cout << "Error p = " << error_p << std::endl;

    //return sqrt(error_u+error_w+error_r+error_p);
    return error_p;

}
/*--------------------------------------------------------------------------*/
void Calculate_Jacobian_Quadrilateral(const Squares2D &Quad, const std::vector<VertexCoordinates2D> &List_Of_Vertices, const double &r_p, const double &s_p, double &Jacobian, double &drdx, double &drdy, double &dsdx, double &dsdy, double &x, double &y)
{

    double x1 = List_Of_Vertices[Quad.getVertex_V1()].getxCoordinate();
    double y1 = List_Of_Vertices[Quad.getVertex_V1()].getyCoordinate();
    double x2 = List_Of_Vertices[Quad.getVertex_V2()].getxCoordinate();
    double y2 = List_Of_Vertices[Quad.getVertex_V2()].getyCoordinate();
    double x3 = List_Of_Vertices[Quad.getVertex_V3()].getxCoordinate();
    double y3 = List_Of_Vertices[Quad.getVertex_V3()].getyCoordinate();
    double x4 = List_Of_Vertices[Quad.getVertex_V4()].getxCoordinate();
    double y4 = List_Of_Vertices[Quad.getVertex_V4()].getyCoordinate();

    double Area = 0.5*abs(x1*y2+x2*y3+x3*y4+x4*y1-x2*y1-x3*y2-x4*y3-x1*y4);

    //std::cout << "ID = " << Quad.getID() << ", Area = " << Area << std::endl;

    x = ((1.0-s_p)*(1.0-r_p)*x1 + (1.0-s_p)*(1.0+r_p)*x2 + (1.0+s_p)*(1.0+r_p)*x3 + (1.0+s_p)*(1.0-r_p)*x4)/4.0;
    y = ((1.0-s_p)*(1.0-r_p)*y1 + (1.0-s_p)*(1.0+r_p)*y2 + (1.0+s_p)*(1.0+r_p)*y3 + (1.0+s_p)*(1.0-r_p)*y4)/4.0;


    double dxdr = (x2+x3-x1-x4)/4.0 + (x1-x2+x3-x4)*s_p/4.0;
    double dxds = (x3+x4-x2-x1)/4.0 + (x1-x2+x3-x4)*r_p/4.0;
    double dydr = (y2+y3-y1-y4)/4.0 + (y1-y2+y3-y4)*s_p/4.0;
    double dyds = (y3+y4-y2-y1)/4.0 + (y1-y2+y3-y4)*r_p/4.0;

    Jacobian = dxdr*dyds-dxds*dydr;

    drdx = dyds/Jacobian;
    drdy = -dxds/Jacobian;
    dsdx = -dydr/Jacobian;
    dsdy = dxdr/Jacobian;
    //std::cout << "Jacobian = " << Jacobian << std::endl;
    /*
    std::cout << r_p << " " << s_p << std::endl;
    std::cout << "Area = " << Area << std::endl;
    std::cout << x1 << " " << x2 << " " << x3 << " " << x4 << std::endl;
    std::cout << y1 << " " << y2 << " " << y3 << " " << y4 << std::endl;

    std::cout << drdx << " " << drdy << " " << dsdx << " " << dsdy << std::endl;
    std::cout << std::endl;
    */
}
/*--------------------------------------------------------------------------*/
