#include "Legendre_Gauss_Lobatto.hpp"

/*--------------------------------------------------------------------------*/
extern Vec JacobiGL(const double &alpha, const double &beta, const unsigned int &N)
{
// function [x] = JacobiGL(alpha,beta,N)
// Purpose: Compute the N'th order Gauss Lobatto quadrature
//          points, x, associated with the Jacobi polynomial,
//          of type (alpha,beta) > -1 ( <> -0.5).

    Vec x;
    VecCreateSeq(PETSC_COMM_WORLD, N+1,&x);
    if (N==0)
    {
        VecSetValue(x, 0, 0.0, INSERT_VALUES);
    }
    else if (N==1)
    {
        VecSetValue(x, 0, -1.0, INSERT_VALUES);
        VecSetValue(x, 1, 1.0, INSERT_VALUES);
    }
    else
    {
        VecSetValue(x, 0, -1.0, INSERT_VALUES);
        Vec xint;
        xint = JacobiGQ(alpha+1, beta+1, N-2);
        PetscScalar    *array;
        VecGetArray(xint,&array);
        for (unsigned int i = 1; i < N; i++)
            VecSetValue(x, i, array[i-1], INSERT_VALUES);
        VecRestoreArray(xint, &array);
        VecSetValue(x, N, 1.0, INSERT_VALUES);
        VecDestroy(&xint);
    }
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    return x;
}
/*--------------------------------------------------------------------------*/
extern Vec JacobiGL_withWeights(const double &alpha, const double &beta, const unsigned int &N, Vec &Weights)
{
// function [x] = JacobiGL(alpha,beta,N)
// Purpose: Compute the N'th order Gauss Lobatto quadrature
//          points, x, associated with the Jacobi polynomial,
//          of type (alpha,beta) > -1 ( <> -0.5).

    Vec x;
    VecCreateSeq(PETSC_COMM_WORLD, N+1, &x);
    VecCreateSeq(PETSC_COMM_WORLD, N+1, &Weights);
    if (N==0)
    {
        VecSetValue(x, 0, 0.0, INSERT_VALUES);
    }
    else if (N==1)
    {
        VecSetValue(x, 0, -1.0, INSERT_VALUES);
        VecSetValue(x, 1, 1.0, INSERT_VALUES);
    }
    else
    {
        VecSetValue(x, 0, -1.0, INSERT_VALUES);
        VecSetValue(Weights, 0, 2.0/(N+1)/(N), INSERT_VALUES);
        Vec xint;
        xint = JacobiGQ(alpha+1, beta+1, N-2);
        PetscScalar    *array;
        VecGetArray(xint,&array);
        Vec ww = JacobiP(xint, 0, 0, N);
        VecScale(ww, sqrt(2.0/(2.0*(N)+1.0)));
        PetscScalar    *warray;
        VecGetArray(ww,&warray);
        for (unsigned int i = 1; i < N; i++)
        {
            VecSetValue(x, i, array[i-1], INSERT_VALUES);
            double w = warray[i-1];
            w = w*w;
            w = 2.0/N/(N+1.0)/w;
            VecSetValue(Weights, i, w, INSERT_VALUES);
        }
        VecRestoreArray(xint, &array);
        VecRestoreArray(ww,&warray);
        VecSetValue(x, N, 1.0, INSERT_VALUES);
        VecSetValue(Weights, N, 2.0/N/(N+1.0), INSERT_VALUES);
        VecDestroy(&xint);

    }
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    VecAssemblyBegin(Weights);
    VecAssemblyEnd(Weights);
    return x;
}
/*--------------------------------------------------------------------------*/
Vec JacobiGQ(const double &alpha, const double &beta, const unsigned int &N)
{
    Vec x;
    VecCreateSeq(PETSC_COMM_WORLD, N+1, &x);
    if (N==0)
    {
        VecSetValue(x, 0, -(alpha-beta)/(alpha+beta+2.0), INSERT_VALUES);
    }
    else
    {
        PetscScalar h1[N+1];
        PetscScalar a[(N+1)*(N+1)]={0};
        for (unsigned int j=0; j<=N; j++)
        {
            h1[j] = 2.0*j+alpha+beta;
        }
        for (unsigned int j = 0; j<=N; j++)
        {
            double value = -0.5*(alpha*alpha-beta*beta)/(h1[j]+2.0)/h1[j];
            a[j+(N+1)*j]= value;
            if (j < N)
            {
                double value1 = 2.0/(h1[j]+2.0)*sqrt((j+1.0)*(j+1.0+alpha+beta)*(j+1.0+alpha)*(j+1.0+beta)/(h1[j]+1.0)/(h1[j]+3.0));
                a[j+(N+1)*(j+1)] = value1;
                a[j+1+(N+1)*(j)] = value1;
            }
        }
        if ((alpha+beta) <= 10.0*std::numeric_limits<double>::epsilon())
            a[0] = 0.0;
        Mat J;
        MatCreate(PETSC_COMM_WORLD,&J);
        MatSetSizes(J, N+1, N+1, N+1, N+1);
        MatSetType(J, MATSEQDENSE);
        MatSeqDenseSetPreallocation(J,a);
        MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
        //MatView(J, PETSC_VIEWER_STDOUT_SELF);

        // Solve for eigenvalues
        EPS eps;
        SlepcInitializeNoArguments();

        EPSCreate(PETSC_COMM_WORLD,&eps);
        EPSSetOperators(eps,J,NULL);
        EPSSetProblemType(eps,EPS_HEP);
        EPSSetFromOptions(eps);
        EPSSetDimensions(eps, N+1, 2*(N+1), PETSC_DEFAULT);
        EPSSolve(eps);

        PetscReal      re;
        PetscScalar    kr,ki;
        Vec            xr,xi;
        PetscInt       i,nconv;

        MatCreateVecs(J,NULL,&xr);
        MatCreateVecs(J,NULL,&xi);

        EPSGetConverged(eps,&nconv);

        PetscScalar lam[nconv];
        if (nconv>0)
        {
            for (i=0;i<nconv;i++)
            {
                EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);
                #if defined(PETSC_USE_COMPLEX)
                    re = PetscRealPart(kr);
                #else
                    re = kr;
                #endif
                lam[i] = (double)re;
            }
        }
        EPSDestroy(&eps);
        VecDestroy(&xr);
        VecDestroy(&xi);
        MatDestroy(&J);

        SlepcFinalize();

        PetscSortReal(N+1,lam);
        PetscInt ix[N+1];
        for (unsigned int k=0;k<=N+1; k++)
            ix[k] = k;

        VecSetValues(x, N+1, ix, lam, INSERT_VALUES);
    }
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

    //VecView(x, PETSC_VIEWER_STDOUT_SELF);
    return x;
}/*--------------------------------------------------------------------------*/
Vec JacobiGQ_withWeights(const double &alpha, const double &beta, const unsigned int &N, Vec &Weights)
{
    Vec x;
    VecCreateSeq(PETSC_COMM_WORLD, N+1, &x);
    VecCreateSeq(PETSC_COMM_WORLD, N+1, &Weights);
    if (N==0)
    {
        VecSetValue(x, 0, -(alpha-beta)/(alpha+beta+2.0), INSERT_VALUES);
        VecSetValue(x, 1, 2.0, INSERT_VALUES);
    }
    else
    {
        PetscScalar h1[N+1];
        PetscScalar a[(N+1)*(N+1)]={0};
        for (unsigned int j=0; j<=N; j++)
        {
            h1[j] = 2.0*j+alpha+beta;
        }
        for (unsigned int j = 0; j<=N; j++)
        {
            double value = -0.5*(alpha*alpha-beta*beta)/(h1[j]+2.0)/h1[j];
            a[j+(N+1)*j]= value;
            if (j < N)
            {
                double value1 = 2.0/(h1[j]+2.0)*sqrt((j+1.0)*(j+1.0+alpha+beta)*(j+1.0+alpha)*(j+1.0+beta)/(h1[j]+1.0)/(h1[j]+3.0));
                a[j+(N+1)*(j+1)] = value1;
                a[j+1+(N+1)*(j)] = value1;
            }
        }
        if ((alpha+beta) <= 10.0*std::numeric_limits<double>::epsilon())
            a[0] = 0.0;
        Mat J;
        MatCreate(PETSC_COMM_WORLD,&J);
        MatSetSizes(J, N+1, N+1, N+1, N+1);
        MatSetType(J, MATSEQDENSE);
        MatSeqDenseSetPreallocation(J,a);
        MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
        //MatView(J, PETSC_VIEWER_STDOUT_SELF);

        // Solve for eigenvalues and eigenvectors
        EPS eps;
        SlepcInitializeNoArguments();

        EPSCreate(PETSC_COMM_WORLD,&eps);
        EPSSetOperators(eps,J,NULL);
        EPSSetProblemType(eps,EPS_HEP);
        EPSSetFromOptions(eps);
        EPSSetDimensions(eps, N+1, 2*(N+1), PETSC_DEFAULT);
        EPSSetTarget(eps, -1.0);
        EPSSetWhichEigenpairs(eps, EPS_TARGET_REAL);
        EPSSolve(eps);

        PetscReal      re;
        PetscScalar    kr,ki;
        Vec            xr,xi;
        PetscInt       i,nconv;

        MatCreateVecs(J,NULL,&xr);
        MatCreateVecs(J,NULL,&xi);

        EPSGetConverged(eps,&nconv);

        PetscScalar lam[nconv];
        PetscScalar w[nconv];
        if (nconv>0)
        {
            for (i=0;i<nconv;i++)
            {
                EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);
                #if defined(PETSC_USE_COMPLEX)
                    re = PetscRealPart(kr);
                #else
                    re = kr;
                #endif
                lam[i] = (double)re;
                PetscScalar *xa;
                VecGetArray(xr, &xa);
                w[i] = xa[0]*xa[0]*pow(2,(alpha+beta+1.0)/(alpha+beta+1.0)*tgamma(alpha+1.0)*tgamma(beta+1.0)/tgamma(alpha+beta+1.0));
                VecRestoreArray(xr, &xa);
            }
        }
        EPSDestroy(&eps);
        VecDestroy(&xr);
        VecDestroy(&xi);
        MatDestroy(&J);
        SlepcFinalize();

        PetscInt ix[N+1];
        for (unsigned int k=0;k<=N+1; k++)
            ix[k] = k;

        VecSetValues(x, N+1, ix, lam, INSERT_VALUES);
        //VecSetValues(x, N+1, ix, w, INSERT_VALUES);
        VecSetValues(Weights, N+1, ix, w, INSERT_VALUES);
    }
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    VecAssemblyBegin(Weights);
    VecAssemblyEnd(Weights);

    //VecView(x, PETSC_VIEWER_STDOUT_SELF);
    return x;
}
/*--------------------------------------------------------------------------*/
extern Vec JacobiP(const Vec &x, const double &alpha, const double &beta, const unsigned int &N)
{
    Vec P;
    Mat PL;
    MatCreate(PETSC_COMM_WORLD,&PL);
    PetscInt size_r;
    VecGetSize(x, &size_r);
    VecCreateSeq(PETSC_COMM_WORLD, size_r,&P);
    MatSetSizes(PL, PETSC_DECIDE, PETSC_DECIDE, N+1, size_r);
    MatSetType(PL, MATSEQAIJ);
    MatSeqAIJSetPreallocation(PL,size_r,NULL);

    // Initial values P_0(x) and P_1(x)
    double gamma0 = pow(2.0, alpha+beta+1.0)/(alpha+beta+1.0)*tgamma(alpha+1.0)*tgamma(beta+1.0)/tgamma(alpha+beta+1.0);
    PetscInt ix[size_r];
    PetscScalar u[size_r];
    PetscInt ir[1]={0};
    for (PetscInt k=0;k<=size_r-1; k++)
    {
        ix[k] = k;
        u[k] = 1.0/sqrt(gamma0);
    }
    MatSetValues(PL, 1, ir, size_r, ix, u, INSERT_VALUES);

    if (N>0)
    {
        double gamma1 = (alpha+1.0)*(beta+1.0)/(alpha+beta+3.0)*gamma0;
        PetscScalar *a;
        VecGetArray(x, &a);
        PetscScalar v[size_r];
        for (PetscInt k=0; k <=size_r-1; k++)
        {
            v[k] = ((alpha+beta+2.0)*a[k]/2.0 + (alpha-beta)/2.0)/sqrt(gamma1);
        }
        ir[0]=1;
        MatSetValues(PL, 1, ir, size_r, ix, v, INSERT_VALUES);
        if (N>1)
        {
            // Recurrence Relation
            PetscScalar aold = 2.0/(2.0+alpha+beta)*sqrt((alpha+1.0)*(beta+1.0)/(alpha+beta+3.0));
            for (unsigned int i=1; i<=N-1;i++)
            {
                PetscScalar h1 = 2.0*(double)i+alpha+beta;
                PetscScalar anew = 2.0/(h1+2.0)*sqrt( ((double)i+1.0)*((double)i+1.0+alpha+beta)*((double)i+1.0+alpha)*((double)i+1.0+beta)/(h1+1.0)/(h1+3.0));
                PetscScalar bnew = - (alpha*alpha-beta*beta)/h1/(h1+2.0);
                PetscScalar w[size_r];
                for (PetscInt k=0; k <=size_r-1; k++)
                {
                    w[k] = 1.0/anew*(-aold*u[k]+(a[k]-bnew)*v[k]);
                    u[k] = v[k];
                    v[k] = w[k];
                }
                ir[0] = i+1;
                MatSetValues(PL, 1, ir, size_r, ix, w, INSERT_VALUES);
                aold = anew;
            }
        }
    VecRestoreArray(x,&a);
    }
    MatAssemblyBegin(PL, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(PL, MAT_FINAL_ASSEMBLY);
    //MatView(PL, PETSC_VIEWER_STDOUT_SELF);

    // Get Last Row
    const PetscScalar *a;
    PetscInt ncols;
    const PetscInt    *cols;
    MatGetRow(PL, N, &ncols,&cols,&a);
    VecSetValues(P, size_r, ix, a, INSERT_VALUES);
    MatRestoreRow(PL, N, &ncols,&cols,&a);
    MatDestroy(&PL);
    VecAssemblyBegin(P);
    VecAssemblyEnd(P);
    //VecView(P, PETSC_VIEWER_STDOUT_SELF);
    return P;
}
/*--------------------------------------------------------------------------*/
extern Mat Vandermonde1D(const Vec &r, const unsigned int &N)
{
    Mat V;
    MatCreate(PETSC_COMM_WORLD,&V);
    MatSetType(V,MATSEQAIJ);
    PetscInt size_r;
    VecGetSize(r, &size_r);
    //MatSetSizes(V, N+1, N+1, PETSC_DECIDE, PETSC_DECIDE);
    MatSetSizes(V, size_r, N+1, PETSC_DECIDE, PETSC_DECIDE);
    MatSeqAIJSetPreallocation(V, size_r, NULL);

    //PetscInt ix[N+1];
    PetscInt ix[size_r];
    PetscInt ir[1]={0};
    //for (unsigned int k=0;k<=N+1; k++)
    for (unsigned int k=0;k<size_r; k++)
    {
        ix[k] = k;
    }
    for (unsigned int l = 0; l<=N; l++)
    {
        Vec P;
        P = JacobiP(r, 0, 0, l);
        //VecView(P, PETSC_VIEWER_STDOUT_SELF);
        PetscScalar *a;
        VecGetArray(P, &a);
        ir[0]=l;
        //MatSetValues(V, N+1, ix, 1, ir, a, INSERT_VALUES);
        MatSetValues(V, size_r, ix, 1, ir, a, INSERT_VALUES);
        VecRestoreArray(P, &a);
        VecDestroy(&P);
    }
    MatAssemblyBegin(V, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(V, MAT_FINAL_ASSEMBLY);

    //MatView(V, PETSC_VIEWER_STDOUT_SELF);
    return V;
}
/*--------------------------------------------------------------------------*/
extern Mat DMatrix1D(const Vec &r, const unsigned int &N, const Mat &V)
{
    Mat Dr;
    Mat Vr;
    Vr = GradVandermonde1D(r, N);
    //MatView(Vr, PETSC_VIEWER_STDOUT_SELF);

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
    info.dtcol=1.0;

    MatLUFactor(A, row, col, &info);
    MatMatSolve(A, B, X);
    //MatView(X, PETSC_VIEWER_STDOUT_SELF);
    // X is the inverse

    //Dr = Vr *inv(V)
    MatMatMult(Vr, X, MAT_INITIAL_MATRIX, 1, &Dr);

    MatDestroy(&Vr);
    MatDestroy(&A);
    MatDestroy(&B);
    MatDestroy(&X);
    ISDestroy(&row);
    ISDestroy(&col);
    return Dr;
}
/*--------------------------------------------------------------------------*/
Mat GradVandermonde1D(const Vec &r, const unsigned int &N)
{
    Mat DVr;
    MatCreate(PETSC_COMM_WORLD,&DVr);
    MatSetType(DVr,MATSEQAIJ);
    MatSetSizes(DVr, N+1, N+1, PETSC_DECIDE, PETSC_DECIDE);
    MatSeqAIJSetPreallocation(DVr, N+1, NULL);

    PetscInt ix[N+1];
    PetscInt ir[1]={0};
    for (unsigned int k=0;k<=N+1; k++)
    {
        ix[k] = k;
    }
    for (unsigned int l = 0; l<=N; l++)
    {
        Vec P;
        P = GradJacobiP(r, 0, 0, l);
        //VecView(P, PETSC_VIEWER_STDOUT_SELF);
        PetscScalar *a;
        VecGetArray(P, &a);
        ir[0]=l;
        MatSetValues(DVr, N+1, ix, 1, ir, a, INSERT_VALUES);
        VecRestoreArray(P, &a);
        VecDestroy(&P);
    }
    MatAssemblyBegin(DVr, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(DVr, MAT_FINAL_ASSEMBLY);

    //MatView(V, PETSC_VIEWER_STDOUT_SELF);

    return DVr;
}
/*--------------------------------------------------------------------------*/
Vec GradJacobiP(const Vec &r, const double &alpha, const double &beta, const unsigned int &N)
{
    PetscInt size_r;
    VecGetSize(r, &size_r);
    Vec dP;
    if (N==0)
    {
        VecCreateSeq(PETSC_COMM_WORLD, size_r, &dP);
        VecSet(dP, 0.0);
    }
    else if (N>0)
    {
        dP=JacobiP(r,alpha+1.0,beta+1.0, N-1);
        VecScale(dP, sqrt(N*(N+alpha+beta+1.0)));
    }
    VecAssemblyBegin(dP);
    VecAssemblyEnd(dP);
    //VecView(dP, PETSC_VIEWER_STDOUT_SELF);
    return dP;
}
/*--------------------------------------------------------------------------*/
extern Mat Lift1D(const unsigned int &N, const Mat &V)
{
    unsigned int Np = N + 1;
    unsigned int Nfp = 1;
    unsigned int Nfaces = 2;

    Mat Emat;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Np, Nfaces*Nfp, Nfaces*Nfp, NULL, &Emat);
    MatSetValue(Emat, 0, 0, 1, INSERT_VALUES);
    MatSetValue(Emat, Np-1, 1, 1, INSERT_VALUES);
    MatAssemblyBegin(Emat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Emat, MAT_FINAL_ASSEMBLY);

    //MatView(Emat, PETSC_VIEWER_STDOUT_SELF);
    Mat Lift;
    Mat Intermediate;
    MatTransposeMatMult(V, Emat, MAT_INITIAL_MATRIX, 1.0, &Intermediate);
    MatMatMult(V, Intermediate, MAT_INITIAL_MATRIX, 1.0, &Lift);
    MatDestroy(&Emat);
    MatDestroy(&Intermediate);
    return Lift;
}
/*--------------------------------------------------------------------------*/
extern Mat GeometricFactors1D(const Mat &x, const Mat &Dr)
{
    Mat J;
    MatMatMult(Dr, x, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &J);
    return J;
}
/*--------------------------------------------------------------------------*/
extern Mat normals1D(const unsigned int &N, const unsigned int &Number_Of_Elements)
{
    //unsigned int Np = N + 1;
    unsigned int Nfp = 1;
    unsigned int Nfaces = 2;

    Mat nx;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, Nfp*Nfaces, Number_Of_Elements, Number_Of_Elements, NULL, &nx);

    for (unsigned int j=0; j<Number_Of_Elements; j++)
    {
        MatSetValue(nx, 0, j, -1, INSERT_VALUES);
        MatSetValue(nx, 1, j, 1, INSERT_VALUES);
    }
    MatAssemblyBegin(nx, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(nx, MAT_FINAL_ASSEMBLY);

    return nx;
}
/*--------------------------------------------------------------------------*/
extern double LagrangePolynomial(const Vec &r, const double &x, const unsigned int &i)
{
    PetscInt size_r;
    VecGetSize(r, &size_r);
    PetscScalar *ra;
    VecGetArray(r, &ra);
    double returnvalue = 1;
    for (PetscInt j = 0; j < size_r; j++)
    {
        if (i!=j)
        {
            returnvalue *= (x-ra[j])/(ra[i]-ra[j]);
        }
    }
    VecRestoreArray(r, &ra);

    return returnvalue;
}
/*--------------------------------------------------------------------------*/
extern double LagrangePolynomialDeriv(const Vec &r, const double &x, const unsigned int &i)
{
    PetscInt size_r;
    VecGetSize(r, &size_r);
    PetscScalar *ra;
    VecGetArray(r, &ra);
    double returnvalue = 0;
    for (PetscInt k = 0; k < size_r; k++)
    {
        if (k!=i)
        {
            double product=1;
            for (PetscInt l = 0; l < size_r; l++)
            {
                if(l!=k && l!=i)
                {
                    product *= x-ra[l];
                }
            }
            returnvalue += product;
        }
    }
    double denom=1;
    for (PetscInt k = 0; k < size_r; k++)
    {
        if(k!=i)
        {
            denom *= ra[i]-ra[k];
        }
    }
    VecRestoreArray(r, &ra);
    returnvalue /= denom;

    return returnvalue;
}
/*--------------------------------------------------------------------------*/
// 2D //
/*--------------------------------------------------------------------------*/
extern void Nodes2D(unsigned int N, Vec &XX, Vec &YY)
{
    // function [x,y] = Nodes2D(N);
    // Purpose  : Compute (x,y) nodes in equilateral triangle for
    //             polynomial of order N


    double alpopt [] = {0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999,  1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258};

    // Set optimized parameter, alpha, depending on order N
    double alpha = 5/3;
    if (N<16)
    {
        alpha = alpopt[N];
    }
    // total number of nodes
    unsigned int Np = (N+1)*(N+2)/2;

    // Create equidistributed nodes on equilateral triangle
    double NN = double(N);


    double L1 [Np] = {};
    double L2 [Np] = {};
    double L3 [Np] = {};
    double x [Np] = {};
    double y [Np] = {};

  int sk = 0;
  for (int n=1; n<=(N+1); ++n) {
    for (int m=1; m<=(N+2-n); ++m) {
      L1[sk] = (n-1)/NN;
      L3[sk] = (m-1)/NN;
      L2[sk] = 1.0-L1[sk]-L3[sk];
      x[sk] = -L2[sk]+L3[sk];
      y[sk] = (-L2[sk]+L3[sk]+2.0*L1[sk])/sqrt(3.0);
      std::cout << L1[sk] << " " << L2[sk] << " " << L3[sk] << " " << x[sk] << " " << y[sk] << std::endl;
      ++sk;
    }
  }
  std::cout << "sk = " << sk << ", Np = " << Np << std::endl;

    Vec X, Y, L1Vec, L2Vec, L3Vec;
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,x,&X);
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,y,&Y);
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,L1,&L1Vec);
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,L2,&L2Vec);
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,L3,&L3Vec);

    VecAssemblyBegin(X);
    VecAssemblyEnd(X);
    VecAssemblyBegin(Y);
    VecAssemblyEnd(Y);
    VecAssemblyBegin(L1Vec);
    VecAssemblyEnd(L1Vec);
    VecAssemblyBegin(L2Vec);
    VecAssemblyEnd(L2Vec);
    VecAssemblyBegin(L3Vec);
    VecAssemblyEnd(L3Vec);



    // Compute blending function at each node for each edge
    // blend1 = 4*L2.*L3; blend2 = 4*L1.*L3; blend3 = 4*L1.*L2;

    Vec blend1, blend2, blend3;
    VecDuplicate(L1Vec, &blend1);
    VecDuplicate(L1Vec, &blend2);
    VecDuplicate(L1Vec, &blend3);

    VecPointwiseMult(L2Vec, L3Vec, blend1);
    VecScale(blend1, 4.0);
    VecPointwiseMult(L1Vec, L3Vec, blend2);
    VecScale(blend2, 4.0);
    VecPointwiseMult(L1Vec, L2Vec, blend3);
    VecScale(blend3, 4.0);

    //  Amount of warp for each node, for each edge
    Vec W1;
    VecDuplicate(L1Vec, &W1);
    VecWAXPY(W1,-1.0,L2Vec,L3Vec);
    Vec warpf1 = Warpfactor(N, W1);
    VecWAXPY(W1,-1.0,L3Vec,L1Vec);
    Vec warpf2 = Warpfactor(N, W1);
    VecWAXPY(W1,-1.0,L1Vec,L2Vec);
    Vec warpf3 = Warpfactor(N, W1);

    PetscInt size_warpf1, size_blend1, size_L1Vec, size_W1;
    VecGetSize(W1,&size_W1);
    VecGetSize(warpf1,&size_warpf1);
    VecGetSize(blend1,&size_blend1);
    VecGetSize(L1Vec,&size_L1Vec);
    std::cout << " Sizes : " << size_W1 << " " << size_warpf1 << " " << size_blend1 << " " << size_L1Vec << std::endl;

    // Combine blend & warp
    PetscScalar *blend_a, *warpf_a, *L_a;
    PetscInt size_rout;
    VecGetSize(blend1,&size_rout);
    PetscScalar warp[size_rout];

    VecGetArray(blend1, &blend_a);
    VecGetArray(warpf1, &warpf_a);
    VecGetArray(L1Vec, &L_a);
    for (PetscInt i = 0; i < size_rout; i++)
    {
        warp[i] =   blend_a[i]*warpf_a[i]*(1.0+pow(alpha*L_a[i],2.0));
        std::cout << i << " " << blend_a[i] << " " << warpf_a[i] << " " << L_a[i] << std::endl;
    }
    VecRestoreArray(blend1, &blend_a);
    VecRestoreArray(warpf1, &warpf_a);
    VecRestoreArray(L1Vec, &L_a);

    Vec warp1;
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,size_rout,warp,&warp1);
    VecAssemblyBegin(warp1);
    VecAssemblyEnd(warp1);

    VecGetArray(blend2, &blend_a);
    VecGetArray(warpf2, &warpf_a);
    VecGetArray(L2Vec, &L_a);
    for (PetscInt i = 0; i < size_rout; i++)
    {
        warp[i] =   blend_a[i]*warpf_a[i]*(1.0+pow(alpha*L_a[i],2.0));
    }
    VecRestoreArray(blend2, &blend_a);
    VecRestoreArray(warpf2, &warpf_a);
    VecRestoreArray(L2Vec, &L_a);

    Vec warp2;
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,size_rout,warp,&warp2);
    VecAssemblyBegin(warp2);
    VecAssemblyEnd(warp2);

    VecGetArray(blend3, &blend_a);
    VecGetArray(warpf3, &warpf_a);
    VecGetArray(L3Vec, &L_a);
    for (PetscInt i = 0; i < size_rout; i++)
    {
        warp[i] =   blend_a[i]*warpf_a[i]*(1.0+pow(alpha*L_a[i],2.0));
    }
    VecRestoreArray(blend3, &blend_a);
    VecRestoreArray(warpf3, &warpf_a);
    VecRestoreArray(L3Vec, &L_a);

    Vec warp3;
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,size_rout,warp,&warp3);
    VecAssemblyBegin(warp3);
    VecAssemblyEnd(warp3);

    // Accumulate deformations associated with each edge
    // x = x + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
    // y = y + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;

    VecAXPY(X, 1.0, warp1);
    VecAXPY(X, cos(2.0/3.0*PETSC_PI), warp2);
    VecAXPY(X, cos(4.0/3.0*PETSC_PI), warp3);

    //VecAXPY(X, 0.0, warp1);
    VecAXPY(Y, sin(2.0/3.0*PETSC_PI), warp2);
    VecAXPY(Y, sin(4.0/3.0*PETSC_PI), warp3);

    VecDestroy(&warp1);
    VecDestroy(&warp2);
    VecDestroy(&warp3);
    VecDestroy(&L1Vec);
    VecDestroy(&L2Vec);
    VecDestroy(&L3Vec);

    VecDestroy(&warpf1);
    VecDestroy(&warpf2);
    VecDestroy(&warpf3);
    VecDestroy(&W1);
    VecDestroy(&blend1);
    VecDestroy(&blend2);
    VecDestroy(&blend3);

    VecDuplicate(X,&XX);
    VecCopy(X, XX);

    VecDuplicate(Y,&YY);
    VecCopy(Y, YY);

    VecDestroy(&X);
    VecDestroy(&Y);

}
/*--------------------------------------------------------------------------*/
extern void XYtoRS(const Vec &X, const Vec &Y, Vec &R, Vec &S)
{
    //function [r,s] = xytors(x, y)
    // Purpose : From (x,y) in equilateral triangle to (r,s) coordinates in standard triangle

    Vec L1;
    VecDuplicate(Y, &L1);
    VecCopy(Y, L1);
    VecScale(L1, sqrt(3.0));
    VecShift(L1, 1.0);
    VecScale(L1, 1.0/3.0);

    Vec L2;
    VecDuplicate(X, &L2);
    VecCopy(X, L2);
    VecScale(L2, -3.0);
    VecAXPY(L2, -sqrt(3.0), Y);
    VecShift(L2, 2.0);
    VecScale(L2, 1.0/6.0);

    Vec L3;
    VecDuplicate(X, &L3);
    VecCopy(X, L3);
    VecScale(L3, 3.0);
    VecAXPY(L3, -sqrt(3.0), Y);
    VecShift(L3, 2.0);
    VecScale(L3, 1.0/6.0);

//L1 = (sqrt(3.0)*y+1.0)/3.0;
//L2 = (-3.0*x - sqrt(3.0)*y + 2.0)/6.0;
//L3 = ( 3.0*x - sqrt(3.0)*y + 2.0)/6.0;

    Vec RR;
    VecDuplicate(L3, &RR);
    VecCopy(L3, RR);
    VecAXPBYPCZ(RR, -1.0, -1.0, 1.0, L1, L2);
    Vec SS;
    VecDuplicate(L3, &SS);
    VecCopy(L3, SS);
    VecAXPBYPCZ(SS, 1.0, -1.0, -1.0, L1, L2);

    VecDuplicate(RR,&R);
    VecCopy(RR, R);
    VecDestroy(&RR);

    VecDuplicate(SS,&S);
    VecCopy(SS, S);
    VecDestroy(&SS);

    VecDestroy(&L1);
    VecDestroy(&L2);
    VecDestroy(&L3);
}
/*--------------------------------------------------------------------------*/
Vec Warpfactor(const unsigned int &N, const Vec &rout)
{
    // function warp = Warpfactor(N, rout)
    // Purpose : Compute scaled warp function at order N based on rout interpolation nodes

    //Compute LGL and equidistant node distribution
    Vec LGLr = JacobiGL(0, 0, N);
    double re [N+1] = {};
    double NN = double(N);
    for (unsigned int i = 0; i <= N; i++)
    {
        re[i] = -1.0+2.0*(double)i/NN;
    }
    // Compute V based on req
    Vec req;
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,N+1,re,&req);
    VecAssemblyBegin(req);
    VecAssemblyEnd(req);
    Mat Veq = Vandermonde1D(req, N);

    // Evaluate Lagrange polynomial at rout
    //PetscInt size_r;
    //VecGetSize(rout, &Nr);

    Mat PL;
    MatCreate(PETSC_COMM_WORLD,&PL);
    PetscInt size_r;
    VecGetSize(rout, &size_r);
    MatSetSizes(PL, PETSC_DECIDE, PETSC_DECIDE, N+1, size_r);
    MatSetType(PL, MATSEQAIJ);
    MatSeqAIJSetPreallocation(PL,size_r,NULL);

    PetscInt ix[size_r];
    PetscInt ir[1]={0};
    for (PetscInt k=0;k<=size_r-1; k++)
    {
        ix[k] = k;
    }

    for (unsigned int i = 1; i <= N+1; i++)
    {
        Vec P;
        P = JacobiP(rout, 0, 0, i-1);
        PetscScalar *a;
        VecGetArray(P, &a);
        ir[0]=i-1;
        MatSetValues(PL, 1, ir, size_r, ix, a, INSERT_VALUES);
        VecRestoreArray(P, &a);
        VecDestroy(&P);
    }

    MatAssemblyBegin(PL, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(PL, MAT_FINAL_ASSEMBLY);

    // Tranpose Vandermonde Matrix
    MatTranspose(Veq, MAT_INPLACE_MATRIX, &Veq);
    // Calculate inverse Vandermonde Matrix
    Mat A, B, X;
    MatDuplicate(Veq,MAT_COPY_VALUES,&A);
    MatDuplicate(Veq,MAT_DO_NOT_COPY_VALUES,&X);
    MatDuplicate(Veq,MAT_DO_NOT_COPY_VALUES,&B);
    MatConvert(B, MATSEQDENSE, MAT_INPLACE_MATRIX, &B);
    MatConvert(X, MATSEQDENSE, MAT_INPLACE_MATRIX, &X);

    MatShift(B, 1.0);

    MatOrderingType rtype = MATORDERINGNATURAL;
    IS row, col;
    MatGetOrdering(A, rtype, &row, &col);
    MatFactorInfo info;
    MatFactorInfoInitialize(&info);
    info.fill=1.0;
    info.dtcol=1.0;

    MatLUFactor(A, row, col, &info);
    MatMatSolve(A, B, X);
    //std::cout << "Veq = " << std::endl;
    //MatView(Veq, PETSC_VIEWER_STDOUT_SELF);
    //std::cout << "X = " << std::endl;
    //MatView(X, PETSC_VIEWER_STDOUT_SELF);
    // X is the inverse

    Mat LMat;
    MatMatMult(X, PL, MAT_INITIAL_MATRIX, 1, &LMat);


    Vec W1;
    VecDuplicate(LGLr, &W1);
    VecWAXPY(W1,-1.0,req,LGLr);
    Vec Warp;
    VecDuplicate(rout,&Warp);
    MatMultTranspose(LMat, W1, Warp);

    PetscScalar *warp_a;
    VecGetArray(Warp, &warp_a);

    PetscScalar *rout_a;
    VecGetArray(rout, &rout_a);
    PetscInt size_rout;
    VecGetSize(Warp,&size_rout);
    for (PetscInt i = 0; i < size_rout; i++)
    {
        double zerof = abs(rout_a[i]) < (1.0-1e-10);
        double sf = 1.0-pow(zerof*rout_a[i],2.0);
        warp_a[i] = warp_a[i]/sf+warp_a[i]*(zerof-1);
    }
    VecRestoreArray(rout, &rout_a);
    VecRestoreArray(Warp, &warp_a);

    VecDestroy(&W1);
    VecDestroy(&req);
    VecDestroy(&LGLr);

    MatDestroy(&Veq);
    MatDestroy(&LMat);
    MatDestroy(&A);
    MatDestroy(&B);
    MatDestroy(&X);
    ISDestroy(&row);
    ISDestroy(&col);


    MatDestroy(&PL);

    return Warp;
}
/*--------------------------------------------------------------------------*/
/*
function [V2D] = Vandermonde2D(N, r, s);
% function [V2D] = Vandermonde2D(N, r, s);
% Purpose : Initialize the 2D Vandermonde Matrix,
%
V_{ij} = phi_j(r_i, s_i);
V2D = zeros(length(r),(N+1)*(N+2)/2);
% Transfer to (a,b) coordinates
[a, b] = rstoab(r, s);
% build the Vandermonde matrix
sk = 1;
for i=0:N
for j=0:N - i
V2D(:,sk) = Simplex2DP(a,b,i,j);
sk = sk+1;
end
end
return;
*/
/*--------------------------------------------------------------------------*/
/*
function [a,b] = rstoab(r,s)
% function [a,b] = rstoab(r,s)
% Purpose : Transfer from (r,s) -> (a,b) coordinates in triangle
Np = length(r); a = zeros(Np,1);
for n=1:Np
if(s(n) ~= 1)
a(n) = 2*(1+r(n))/(1-s(n))-1;
else a(n) = -1; end
end
b = s;
return;
/*
/*--------------------------------------------------------------------------*/
/*
function [P] = Simplex2DP(a,b,i,j);
% function [P] = Simplex2DP(a,b,i,j);
% Purpose : Evaluate 2D orthonormal polynomial
%
on simplex at (a,b) of order (i,j).
h1 = JacobiP(a,0,0,i); h2 = JacobiP(b,2*i+1,0,j);
P = sqrt(2.0)*h1.*h2.*(1-b).^i;
return;
*/
/*--------------------------------------------------------------------------*/
