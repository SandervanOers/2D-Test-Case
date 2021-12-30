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
        VecSetValue(Weights, 0, 2.0, INSERT_VALUES);
    }
    else if (N==1)
    {
        VecSetValue(x, 0, -1.0, INSERT_VALUES);
        VecSetValue(x, 1, 1.0, INSERT_VALUES);
        VecSetValue(Weights, 0, 1.0, INSERT_VALUES);
        VecSetValue(Weights, 1, 1.0, INSERT_VALUES);
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
        VecDestroy(&ww);
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
Vec SecondDerivJacobiP(const Vec &r, const double &alpha, const double &beta, const unsigned int &N)
{
    PetscInt size_r;
    VecGetSize(r, &size_r);
    Vec d2P;
    if (N==0 || N == 1)
    {
        VecCreateSeq(PETSC_COMM_WORLD, size_r, &d2P);
        VecSet(d2P, 0.0);
    }
    else if (N>0)
    {
        d2P = JacobiP(r,alpha+2.0,beta+2.0, N-2);
        VecScale(d2P, sqrt(N*(N+alpha+beta+1.0)*(N-1)*(N+alpha+beta+2.0)));
    }
    VecAssemblyBegin(d2P);
    VecAssemblyEnd(d2P);
    //VecView(d2P, PETSC_VIEWER_STDOUT_SELF);
    return d2P;
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
extern double LagrangePolynomial_Test(const Vec &r, const double &x, const unsigned int &i)
{
    PetscInt size_r;
    VecGetSize(r, &size_r);
    unsigned int N = size_r-1;

    double returnvalue = 1.0;
    if (N == 0)
    {
        return 1;
    }
    else
    {

    PetscScalar *r_a;
    VecGetArray(r, &r_a);
    if (abs(x-r_a[i]) < NODETOL)
    {
        VecRestoreArray(r, &r_a);
        return 1;
    }
    else
    {
    Vec X, Ri;
    VecCreateSeq(PETSC_COMM_SELF, 1, &X);
    VecCreateSeq(PETSC_COMM_SELF, 1, &Ri);
    VecSetValue(X,0, x, INSERT_VALUES);
    VecSetValue(Ri,0, r_a[i], INSERT_VALUES);

    Vec P = JacobiP(Ri, 0, 0, N);
    Vec dP = GradJacobiP(X, 0, 0, N);


    PetscScalar *P_a, *dP_a;
    VecGetArray(P, &P_a);
    VecGetArray(dP, &dP_a);
    returnvalue = - 1.0/N/(N+1.0)*(1.0-x*x)/(x-r_a[i])*dP_a[0]/P_a[0];
    VecRestoreArray(P, &P_a);
    VecRestoreArray(dP, &dP_a);
    VecRestoreArray(r, &r_a);
    VecDestroy(&X);
    VecDestroy(&Ri);
    VecDestroy(&P);
    VecDestroy(&dP);
    }

    }
    return returnvalue;
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
            //std::cout << "(x-ra[j])/(ra[i]-ra[j]) = " << (x-ra[j])/(ra[i]-ra[j]) << std::endl;
        }
    }
    VecRestoreArray(r, &ra);

    return returnvalue;
}
/*--------------------------------------------------------------------------*/
extern double LagrangePolynomialDeriv_Test(const Vec &r, const double &x, const unsigned int &i)
{
    PetscInt size_r;
    VecGetSize(r, &size_r);
    unsigned int N = size_r-1;

    double returnvalue = 0.0;
    if (N == 0)
    {
    }
    else
    {

    PetscScalar *r_a;
    VecGetArray(r, &r_a);
    if (abs(x-r_a[i]) < NODETOL)
    {

    }
    else
    {
    Vec X, Ri;
    VecCreateSeq(PETSC_COMM_SELF, 1, &X);
    VecCreateSeq(PETSC_COMM_SELF, 1, &Ri);
    VecSetValue(X,0, x, INSERT_VALUES);
    VecSetValue(Ri,0, r_a[i], INSERT_VALUES);

    Vec P = JacobiP(Ri, 0, 0, N);
    Vec dP = GradJacobiP(X, 0, 0, N);
    Vec d2P = SecondDerivJacobiP(X, 0, 0, N);


    PetscScalar *P_a, *dP_a, *d2P_a;
    VecGetArray(P, &P_a);
    VecGetArray(dP, &dP_a);
    VecGetArray(d2P, &d2P_a);
    //std::cout << "d2P_a[0] = " << d2P_a[0] << std::endl;
    returnvalue = - 1.0/N/(N+1.0)/P_a[0]/(x-r_a[i])/(x-r_a[i])*(-2.0*x*(x-r_a[i])*dP_a[0]+(1.0-x*x)*(d2P_a[0]*(x-r_a[i])-dP_a[0]));
    VecRestoreArray(P, &P_a);
    VecRestoreArray(dP, &dP_a);
    VecRestoreArray(d2P, &d2P_a);
    VecRestoreArray(r, &r_a);
    VecDestroy(&X);
    VecDestroy(&Ri);
    VecDestroy(&P);
    VecDestroy(&dP);
    VecDestroy(&d2P);
    }

    }

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

    // total number of nodes
    unsigned int Np = (N+1)*(N+2)/2;
    double x [Np] = {};
    double y [Np] = {};


    if (N > 0)
    {

    double alpopt [] = {0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999,  1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258};

    double NN = double(N);

    // Set optimized parameter, alpha, depending on order N
    double alpha = 5/3;
    if (N<16)
    {
        alpha = alpopt[N-1];
    }

    // Create equidistributed nodes on equilateral triangle
    double L1 [Np] = {};
    double L2 [Np] = {};
    double L3 [Np] = {};

  int sk = 0;
  for (int n=1; n<=(N+1); ++n) {
    for (int m=1; m<=(N+2-n); ++m) {
      L1[sk] = (n-1)/NN;
      L3[sk] = (m-1)/NN;
      L2[sk] = 1.0-L1[sk]-L3[sk];
      ++sk;
    }
  }

    Vec L1Vec, L2Vec, L3Vec;
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,L1,&L1Vec);
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,L2,&L2Vec);
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,L3,&L3Vec);
    VecAssemblyBegin(L1Vec);
    VecAssemblyEnd(L1Vec);
    VecAssemblyBegin(L2Vec);
    VecAssemblyEnd(L2Vec);
    VecAssemblyBegin(L3Vec);
    VecAssemblyEnd(L3Vec);
    // Compute blending function at each node for each edge
    // blend1 = 4*L2.*L3; blend2 = 4*L1.*L3; blend3 = 4*L1.*L2;

    //  Amount of warp for each node, for each edge
    Vec W1;
    VecDuplicate(L1Vec, &W1);
    VecWAXPY(W1,-1.0,L2Vec,L3Vec);
    Vec warpf1 = Warpfactor(N, W1);
    VecWAXPY(W1,-1.0,L3Vec,L1Vec);
    Vec warpf2 = Warpfactor(N, W1);
    VecWAXPY(W1,-1.0,L1Vec,L2Vec);
    Vec warpf3 = Warpfactor(N, W1);

    // Combine blend & warp
    PetscScalar *warpf1_a, *warpf2_a, *warpf3_a;
    VecGetArray(warpf1, &warpf1_a);
    VecGetArray(warpf2, &warpf2_a);
    VecGetArray(warpf3, &warpf3_a);
    for (PetscInt i = 0; i < Np; i++)
    {
        double blend1 = 4.0*L2[i]*L3[i];
        double blend2 = 4.0*L1[i]*L3[i];
        double blend3 = 4.0*L1[i]*L2[i];

        double warp1 = blend1*warpf1_a[i]*(1.0+pow(alpha*L1[i],2.0));
        double warp2 = blend2*warpf2_a[i]*(1.0+pow(alpha*L2[i],2.0));
        double warp3 = blend3*warpf3_a[i]*(1.0+pow(alpha*L3[i],2.0));

        x[i] = -L2[i]+L3[i] + warp1 - 0.5*warp2 - 0.5*warp3;
        y[i] = (-L2[i]-L3[i]+2.0*L1[i])/sqrt(3.0) + 0.5*sqrt(3)*warp2 -0.5*sqrt(3)*warp3;
        //std::cout << i << " " << blend1 << " " << blend2 << " " << blend3 << " " << warp1 << " " << blend2 << " " << blend3 << " " << x[i] << " " << y[i] << std::endl;
    }
    VecRestoreArray(warpf1, &warpf1_a);
    VecRestoreArray(warpf2, &warpf2_a);
    VecRestoreArray(warpf3, &warpf3_a);

    VecDestroy(&L1Vec);
    VecDestroy(&L2Vec);
    VecDestroy(&L3Vec);

    VecDestroy(&warpf1);
    VecDestroy(&warpf2);
    VecDestroy(&warpf3);
    VecDestroy(&W1);
    }
    else if (N == 0)
    {
        x[0] = 0;
        y[0] = 0;
    }

    Vec X, Y;
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,x,&X);
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,y,&Y);


    VecAssemblyBegin(X);
    VecAssemblyEnd(X);
    VecAssemblyBegin(Y);
    VecAssemblyEnd(Y);

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
    // Purpose : From (x,y) in equilateral triangle to (r,s) coordinates in reference triangle

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

    // Assumes N > 0
    // if (N == 0) NN = 1.0;

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
extern Mat Vandermonde2D(const unsigned int &N, const Vec &r, const Vec &s)
{
    // function [V2D] = Vandermonde2D(N, r, s);
    // Purpose : Initialize the 2D Vandermonde Matrix
    // V_{ij} = phi_j(r_i, s_i);

    // Transfer to (a,b) coordinates
    Vec A, B;
    RStoAB(r, s, A, B);

    unsigned int Np = (N+1)*(N+2)/2;

    Mat V;
    MatCreate(PETSC_COMM_WORLD,&V);
    MatSetType(V,MATSEQAIJ);
    PetscInt size_r;
    VecGetSize(r, &size_r);
    MatSetSizes(V, size_r, Np, PETSC_DECIDE, PETSC_DECIDE);
    MatSeqAIJSetPreallocation(V, size_r, NULL);

    PetscInt ix[size_r];
    PetscInt ir[1]={0};
    //for (unsigned int k=0;k<=N+1; k++)
    for (unsigned int k=0;k<size_r; k++)
    {
        ix[k] = k;
    }
    unsigned int sk = 0;
    for (unsigned int i = 0; i <= N; i++)
    {
        for (unsigned int j = 0; j <= N - i; j ++)
        {
            Vec P;
            P = Simplex2DP(A, B, i, j);

            PetscScalar *p_a;
            VecGetArray(P, &p_a);
            ir[0]=sk;

            MatSetValues(V, size_r, ix, 1, ir, p_a, INSERT_VALUES);
            VecRestoreArray(P, &p_a);
            VecDestroy(&P);
            sk++;
        }
    }
    MatAssemblyBegin(V, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(V, MAT_FINAL_ASSEMBLY);

    VecDestroy(&A);
    VecDestroy(&B);
    return V;

}
/*--------------------------------------------------------------------------*/
void RStoAB(const Vec &R, const Vec &S, Vec &A, Vec &B)
{
    // function [a,b] = rstoab(r,s)
    // Purpose : Transfer from (r,s) -> (a,b) coordinates in triangle

    PetscScalar *r_a, *s_a;
    VecGetArray(R, &r_a);
    VecGetArray(S, &s_a);
    PetscInt size_r;
    VecGetSize(R,&size_r);
    double a [size_r] = {};
    for (unsigned int i = 0; i < size_r; i++)
    {
        if (s_a[i] != 1)
        {
            a[i] = 2.0*(1.0+r_a[i])/(1.0-s_a[i])-1.0;
        }
        else
        {
            a[i] = 1.0;
        }
    }
    VecRestoreArray(R, &r_a);
    VecRestoreArray(S, &s_a);

    Vec AA;
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,size_r,a,&AA);


    VecAssemblyBegin(AA);
    VecAssemblyEnd(AA);

    VecDuplicate(AA,&A);
    VecCopy(AA, A);

    VecDuplicate(S,&B);
    VecCopy(S, B);

    VecDestroy(&AA);

}
/*--------------------------------------------------------------------------*/
Vec Simplex2DP(const Vec &A, const Vec &B, const unsigned int &i, const unsigned int &j)
{
    // function [P] = Simplex2DP(a,b,i,j);
    // Purpose : Evaluate 2D orthonormal polynomial
    // on simplex at (a,b) of order (i,j).

    Vec h1, h2;
    h1 = JacobiP(A, 0, 0, i);
    h2 = JacobiP(B, 2*i+1, 0, j);

    // P = sqrt(2.0)*h1.*h2.*(1-b).^i;
    PetscScalar *h1_a, *h2_a, *b_a;
    VecGetArray(h1, &h1_a);
    VecGetArray(h2, &h2_a);
    VecGetArray(B, &b_a);
    PetscInt size_r;
    VecGetSize(B,&size_r);
    double p [size_r] = {};
    for (unsigned int k = 0; k < size_r; k++)
    {
        p[k] = sqrt(2.0)*h1_a[k]*h2_a[k]*pow(1.0-b_a[k],i);
    }
    VecRestoreArray(h1, &h1_a);
    VecRestoreArray(h2, &h2_a);
    VecRestoreArray(B, &b_a);

    VecDestroy(&h1);
    VecDestroy(&h2);

    Vec P, PP;
    VecCreateSeqWithArray(PETSC_COMM_WORLD,1,size_r,p,&PP);

    VecAssemblyBegin(PP);
    VecAssemblyEnd(PP);

    VecDuplicate(PP, &P);
    VecCopy(PP, P);
    VecDestroy(&PP);

    return P;
}
/*--------------------------------------------------------------------------*/
void set_Node_Coordinates_Uniform(std::vector<std::unique_ptr<Element>> &List_Of_Elements, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, const unsigned int &Nx, const unsigned int &Ny, const unsigned int &Nz)
{
    Vec R, S, T;
    // Compute Nodes on a Reference Element (Cuboid)
    NodesCuboid(Nx, Ny, Nz, R, S, T);

      PetscScalar *r_p, *s_p, *t_p;
      VecGetArray(R, &r_p);
      VecGetArray(S, &s_p);
      VecGetArray(T, &t_p);
      PetscInt size_r;
      VecGetSize(R, &size_r);
      // Compute the physical location of each node on each element
      for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
      {
          double x0 = List_Of_Vertices[(*i)->getVertex_V1()]->getxCoordinate();
          double y0 = List_Of_Vertices[(*i)->getVertex_V1()]->getyCoordinate();
          double z0 = List_Of_Vertices[(*i)->getVertex_V1()]->getzCoordinate();
          double x1 = List_Of_Vertices[(*i)->getVertex_V2()]->getxCoordinate();
          double y1 = List_Of_Vertices[(*i)->getVertex_V2()]->getyCoordinate();
          double z1 = List_Of_Vertices[(*i)->getVertex_V2()]->getzCoordinate();
          double x2 = List_Of_Vertices[(*i)->getVertex_V3()]->getxCoordinate();
          double y2 = List_Of_Vertices[(*i)->getVertex_V3()]->getyCoordinate();
          double z2 = List_Of_Vertices[(*i)->getVertex_V3()]->getzCoordinate();
          double x3 = List_Of_Vertices[(*i)->getVertex_V4()]->getxCoordinate();
          double y3 = List_Of_Vertices[(*i)->getVertex_V4()]->getyCoordinate();
          double z3 = List_Of_Vertices[(*i)->getVertex_V4()]->getzCoordinate();
          double x4 = List_Of_Vertices[(*i)->getVertex_V5()]->getxCoordinate();
          double y4 = List_Of_Vertices[(*i)->getVertex_V5()]->getyCoordinate();
          double z4 = List_Of_Vertices[(*i)->getVertex_V5()]->getzCoordinate();
          double x5 = List_Of_Vertices[(*i)->getVertex_V6()]->getxCoordinate();
          double y5 = List_Of_Vertices[(*i)->getVertex_V6()]->getyCoordinate();
          double z5 = List_Of_Vertices[(*i)->getVertex_V6()]->getzCoordinate();
          double x6 = List_Of_Vertices[(*i)->getVertex_V7()]->getxCoordinate();
          double y6 = List_Of_Vertices[(*i)->getVertex_V7()]->getyCoordinate();
          double z6 = List_Of_Vertices[(*i)->getVertex_V7()]->getzCoordinate();
          double x7 = List_Of_Vertices[(*i)->getVertex_V8()]->getxCoordinate();
          double y7 = List_Of_Vertices[(*i)->getVertex_V8()]->getyCoordinate();
          double z7 = List_Of_Vertices[(*i)->getVertex_V8()]->getzCoordinate();

          for(unsigned int k = 0; k < size_r; k++)
          {
              // All nodes
              double x = (1.0-r_p[k])*(1.0-s_p[k])*(1.0-t_p[k])*x0 + (1.0+r_p[k])*(1.0-s_p[k])*(1.0-t_p[k])*x1 + (1.0+r_p[k])*(1.0+s_p[k])*(1.0-t_p[k])*x2 + (1.0-r_p[k])*(1.0+s_p[k])*(1.0-t_p[k])*x3 + (1.0-r_p[k])*(1.0-s_p[k])*(1.0+t_p[k])*x4 + (1.0+r_p[k])*(1.0-s_p[k])*(1.0+t_p[k])*x5 + (1.0+r_p[k])*(1.0+s_p[k])*(1.0+t_p[k])*x6 + (1.0-r_p[k])*(1.0+s_p[k])*(1.0+t_p[k])*x7;
              double y = (1.0-r_p[k])*(1.0-s_p[k])*(1.0-t_p[k])*y0 + (1.0+r_p[k])*(1.0-s_p[k])*(1.0-t_p[k])*y1 + (1.0+r_p[k])*(1.0+s_p[k])*(1.0-t_p[k])*y2 + (1.0-r_p[k])*(1.0+s_p[k])*(1.0-t_p[k])*y3 + (1.0-r_p[k])*(1.0-s_p[k])*(1.0+t_p[k])*y4 + (1.0+r_p[k])*(1.0-s_p[k])*(1.0+t_p[k])*y5 + (1.0+r_p[k])*(1.0+s_p[k])*(1.0+t_p[k])*y6 + (1.0-r_p[k])*(1.0+s_p[k])*(1.0+t_p[k])*y7;
              double z = (1.0-r_p[k])*(1.0-s_p[k])*(1.0-t_p[k])*z0 + (1.0+r_p[k])*(1.0-s_p[k])*(1.0-t_p[k])*z1 + (1.0+r_p[k])*(1.0+s_p[k])*(1.0-t_p[k])*z2 + (1.0-r_p[k])*(1.0+s_p[k])*(1.0-t_p[k])*z3 + (1.0-r_p[k])*(1.0-s_p[k])*(1.0+t_p[k])*z4 + (1.0+r_p[k])*(1.0-s_p[k])*(1.0+t_p[k])*z5 + (1.0+r_p[k])*(1.0+s_p[k])*(1.0+t_p[k])*z6 + (1.0-r_p[k])*(1.0+s_p[k])*(1.0+t_p[k])*z7;
              x = x/8.0;
              y = y/8.0;
              z = z/8.0;

              // store node physical coordinates for each element
              (*i)->set_node_coordinates_x(x);
              (*i)->set_node_coordinates_y(y);
              (*i)->set_node_coordinates_z(z);

                // Boundary Nodes
                if (abs(t_p[k]+1.0) < NODETOL || Nz == 0)
                {
                    (*i)->set_node_on_face0(k);
                }
                if (abs(t_p[k]-1.0) < NODETOL || Nz == 0)
                {
                    (*i)->set_node_on_face1(k);
                }
                if (abs(s_p[k]+1.0) < NODETOL || Ny == 0)
                {
                    (*i)->set_node_on_face2(k);
                }
                if (abs(s_p[k]-1.0) < NODETOL || Ny == 0)
                {
                    (*i)->set_node_on_face3(k);
                }
                if (abs(r_p[k]+1.0) < NODETOL || Nx == 0)
                {
                    (*i)->set_node_on_face4(k);
                }
                if (abs(r_p[k]-1.0) < NODETOL || Nx == 0)
                {
                    (*i)->set_node_on_face5(k);
                }

            }
      }
      VecRestoreArray(R, &r_p);
      VecRestoreArray(S, &s_p);
      VecRestoreArray(T, &t_p);

      VecDestroy(&R);
      VecDestroy(&S);
      VecDestroy(&T);
}
/*--------------------------------------------------------------------------*/
void set_Node_Coordinates_Uniform_Square2D(std::vector<Squares2D> &List_Of_Elements, const std::vector<std::unique_ptr<Vertex>> &List_Of_Vertices, const unsigned int &N)
{
    Vec R, S;
    // Compute Nodes on a Reference Element (Quadrilateral)
    NodesSquares2D(N, R, S);

    PetscScalar *r_a, *s_a;
    VecGetArray(R, &r_a);
    VecGetArray(S, &s_a);
    PetscInt size_r;
    VecGetSize(R, &size_r);
    // Compute the physical location of each node on each element
    for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
    {
        double x_v1 = List_Of_Vertices[(*i).getVertex_V1()]->getxCoordinate();
        double y_v1 = List_Of_Vertices[(*i).getVertex_V1()]->getyCoordinate();
        double x_v2 = List_Of_Vertices[(*i).getVertex_V2()]->getxCoordinate();
        double y_v2 = List_Of_Vertices[(*i).getVertex_V2()]->getyCoordinate();
        double x_v3 = List_Of_Vertices[(*i).getVertex_V3()]->getxCoordinate();
        double y_v3 = List_Of_Vertices[(*i).getVertex_V3()]->getyCoordinate();
        double x_v4 = List_Of_Vertices[(*i).getVertex_V4()]->getxCoordinate();
        double y_v4 = List_Of_Vertices[(*i).getVertex_V4()]->getyCoordinate();

        for(unsigned int k = 0; k < size_r; k++)
        {
            double x = 0.25*((1.0-s_a[k]-r_a[k]+s_a[k]*r_a[k])*x_v1+(1.0-s_a[k]+r_a[k]-s_a[k]*r_a[k])*x_v2+(1.0+s_a[k]+r_a[k]+s_a[k]*r_a[k])*x_v3+(1.0+s_a[k]-r_a[k]-s_a[k]*r_a[k])*x_v4);
            double y = 0.25*((1.0-s_a[k]-r_a[k]+s_a[k]*r_a[k])*y_v1+(1.0-s_a[k]+r_a[k]-s_a[k]*r_a[k])*y_v2+(1.0+s_a[k]+r_a[k]+s_a[k]*r_a[k])*y_v3+(1.0+s_a[k]-r_a[k]-s_a[k]*r_a[k])*y_v4);

            (*i).set_node_coordinates_x(x);
            (*i).set_node_coordinates_y(y);

            unsigned int type_node = 0;

            // Boundary Elements
            if (abs(s_a[k]+1.0) < NODETOL || N == 0)
            {
                // face 1
                ///unsigned int ID_B1 = (*i).getBoundary_B1();
                ///List_Of_Boundaries2D[ID_B1].set_node_coordinates_x(x);
                ///List_Of_Boundaries2D[ID_B1].set_node_coordinates_y(y);
                type_node = 1;
                (*i).set_node_on_boundary_1(k);
            }
            if (abs(r_a[k]+1.0) < NODETOL || N == 0)
            {
                // face 4
                ///unsigned int ID_B3 = (*i).getBoundary_B3();
                ///List_Of_Boundaries2D[ID_B3].set_node_coordinates_x(x);
                ///List_Of_Boundaries2D[ID_B3].set_node_coordinates_y(y);
                type_node = 4;
                (*i).set_node_on_boundary_4(k);
            }
            if (abs(r_a[k]-1.0) < NODETOL || N == 0)
            {
                // face 2
                ///unsigned int ID_B2 = (*i).getBoundary_B2();
                ///List_Of_Boundaries2D[ID_B2].set_node_coordinates_x(x);
                ///List_Of_Boundaries2D[ID_B2].set_node_coordinates_y(y);
                type_node = 2;
                (*i).set_node_on_boundary_2(k);
            }
            if (abs(s_a[k]-1.0) < NODETOL || N == 0)
            {
                // face 3
                ///unsigned int ID_B2 = (*i).getBoundary_B2();
                ///List_Of_Boundaries2D[ID_B2].set_node_coordinates_x(x);
                ///List_Of_Boundaries2D[ID_B2].set_node_coordinates_y(y);
                type_node = 3;
                (*i).set_node_on_boundary_3(k);
            }
            (*i).set_node_on_boundary(type_node);

        }
    }
    VecRestoreArray(R, &r_a);
    VecRestoreArray(S, &s_a);

    VecDestroy(&R);
    VecDestroy(&S);
}
/*--------------------------------------------------------------------------*/
void store_Nodes_Reference_Triangle()
{

    for (unsigned int N = 0; N < 20; N ++)
    {
    auto t1 = std::chrono::high_resolution_clock::now();



        // Compute the nodes for a equilateral element
        Vec X, Y;
        Nodes2D(N, X, Y);
        // Convert nodes to a reference element
        Vec R, S;
        XYtoRS(X, Y, R, S);
        {
        std::string LocationName = "Nodes2D/";
        LocationName.append(std::to_string(N));
        LocationName.append("R.dat");
        PetscViewer    viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_WRITE,&viewer);
        VecView(R,viewer);
        PetscViewerDestroy(&viewer);
        }

        {
        std::string LocationName = "Nodes2D/";
        LocationName.append(std::to_string(N));
        LocationName.append("S.dat");
        PetscViewer    viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_WRITE,&viewer);
        VecView(S,viewer);
        PetscViewerDestroy(&viewer);
        }


        Mat V2D = Vandermonde2D(N, R, S);
            //std::cout << "V2D  = " << std::endl;
            //MatView(V2D , PETSC_VIEWER_STDOUT_SELF);
        {
        std::string LocationName = "Vandermonde2D/";
        LocationName.append(std::to_string(N));
        LocationName.append("V2D.dat");
        PetscViewer    viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_WRITE,&viewer);
        MatView(V2D,viewer);
        PetscViewerDestroy(&viewer);
        }
        Mat V2DInv = Inverse_Matrix(V2D);


            //std::cout << "V2DInv  = " << std::endl;
            //MatView(V2DInv , PETSC_VIEWER_STDOUT_SELF);
        {
        std::string LocationName = "Vandermonde2D/";
        LocationName.append(std::to_string(N));
        LocationName.append("V2DInv.dat");
        PetscViewer    viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_WRITE,&viewer);
        MatView(V2DInv,viewer);
        PetscViewerDestroy(&viewer);
        }
        MatDestroy(&V2DInv);

        MatDestroy(&V2D);
        VecDestroy(&R);
        VecDestroy(&S);
        VecDestroy(&X);
        VecDestroy(&Y);
        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << "N = " << N << " Execution took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";



    }
}
/*--------------------------------------------------------------------------*/
extern void Read_RS_Coordinates_Reference_Triangle(const unsigned int &N, Vec &R, Vec &S)
{
    Vec RR, SS;
            // Read the nodes for a Reference Element
            {
            std::string LocationName = "Nodes2D/";
            LocationName.append(std::to_string(N));
            LocationName.append("R.dat");
            PetscViewer    viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_READ,&viewer);
            VecCreate(PETSC_COMM_WORLD,&RR);
            VecLoad(RR,viewer);
            PetscViewerDestroy(&viewer);
            }

            {
            std::string LocationName = "Nodes2D/";
            LocationName.append(std::to_string(N));
            LocationName.append("S.dat");
            PetscViewer    viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_READ,&viewer);
            VecCreate(PETSC_COMM_WORLD,&SS);
            VecLoad(SS,viewer);
            PetscViewerDestroy(&viewer);
            }

        VecDuplicate(RR, &R);
        VecCopy(RR, R);
        VecDuplicate(SS, &S);
        VecCopy(SS, S);

        VecDestroy(&RR);
        VecDestroy(&SS);
}
/*--------------------------------------------------------------------------*/
void GradSimplex2D(const Vec &a, const Vec &b, const unsigned int &id, const unsigned int &jd, Vec &dmodedr, Vec &dmodeds)
{
  // function [dmodedr, dmodeds] = GradSimplex2DP(a,b,id,jd)
  // Purpose: Return the derivatives of the modal basis (id,jd)
  //          on the 2D simplex at (a,b).

  Vec fa = JacobiP(a, 0, 0, id);
  Vec gb = JacobiP(b, 2*id+1, 0, jd);

  Vec dfa = GradJacobiP(a, 0, 0, id);
  Vec dgb = GradJacobiP(b, 2*id+1, 0, jd);

  PetscScalar *fa_a, *dfa_a, *gb_a, *dgb_a;
  VecGetArray(fa, &fa_a);
  VecGetArray(dfa, &dfa_a);
  VecGetArray(gb, &gb_a);
  VecGetArray(dgb, &dgb_a);

  // r-derivative
  // d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
  VecDuplicate(fa, &dmodedr);
  VecPointwiseMult(dmodedr, dfa, gb);

  // s-derivative
  // d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
  VecDuplicate(fa, &dmodeds);
  VecPointwiseMult(dmodeds, dfa, gb);

  Vec tmp;
  VecDuplicate(a, &tmp);
  VecCopy(a, tmp);
  VecShift(tmp, 1.0);
  VecScale(tmp, 0.5);
  VecPointwiseMult(dmodeds, dmodeds, tmp);

  PetscScalar *dmodedr_a, *dmodeds_a;
  VecGetArray(dmodedr, &dmodedr_a);
  VecGetArray(dmodeds, &dmodeds_a);
  VecCopy(b, tmp);
  PetscScalar *b_a;
  VecGetArray(b, &b_a);
  PetscScalar *tmp_a;
  VecGetArray(tmp, &tmp_a);
  PetscInt size_b;
  VecGetSize(b, &size_b);
  for (unsigned int i = 0; i < size_b; i++)
  {
      tmp_a[i] = dgb_a[i]*pow(0.5*(1.0-b_a[i]), double(id));
  }
  if (id > 0)
  {
    for (unsigned int i = 0; i < size_b; i++)
    {
        dmodedr_a[i] = dmodedr_a[i]*pow(0.5*(1.0-b_a[i]), (id-1.0));
        dmodeds_a[i] = dmodeds_a[i]*pow(0.5*(1.0-b_a[i]), (id-1.0));
        tmp_a[i] = tmp_a[i] - 0.5*double(id)*gb_a[i]*pow(0.5*(1.0-b_a[i]),id-1.0);
    }
  }
  VecRestoreArray(dmodedr, &dmodedr_a);
  VecRestoreArray(dmodeds, &dmodeds_a);
  VecRestoreArray(b, &b_a);
  VecRestoreArray(fa, &fa_a);
  VecRestoreArray(dfa, &dfa_a);
  VecRestoreArray(gb, &gb_a);
  VecRestoreArray(dgb, &dgb_a);
  VecRestoreArray(tmp, &tmp_a);

  VecPointwiseMult(tmp, fa, tmp);
  VecAXPY(dmodeds,1.0,tmp); /// Check this

  // Normalize
  VecScale(dmodedr, pow(2.0, (id+0.5)));
  VecScale(dmodeds, pow(2.0, (id+0.5)));

  VecDestroy(&fa);
  VecDestroy(&gb);
  VecDestroy(&dfa);
  VecDestroy(&dgb);
  VecDestroy(&tmp);
}
/*--------------------------------------------------------------------------*/
void GradVandermonde2D(const unsigned int &N, const Vec &R, const Vec &S, Mat &V2Dr, Mat &V2Ds)
{
    unsigned int Np = (N+1)*(N+2)/2;
    PetscInt size_r;
    VecGetSize(R, &size_r);

    MatCreate(PETSC_COMM_WORLD,&V2Dr);
    MatCreate(PETSC_COMM_WORLD,&V2Ds);
    MatSetType(V2Dr,MATSEQAIJ);
    MatSetType(V2Ds,MATSEQAIJ);
    MatSetSizes(V2Dr, size_r, Np, PETSC_DECIDE, PETSC_DECIDE);
    MatSetSizes(V2Ds, size_r, Np, PETSC_DECIDE, PETSC_DECIDE);
    MatSeqAIJSetPreallocation(V2Dr, size_r, NULL);
    MatSeqAIJSetPreallocation(V2Ds, size_r, NULL);

    PetscInt ix[size_r];
    PetscInt ir[1]={0};
    for (unsigned int k=0;k<size_r; k++)
    {
        ix[k] = k;
    }

    // find tensor-product coordinates
    Vec A, B;
    RStoAB(R, S, A, B);

    // Initialize matrices
    unsigned int sk = 0;
    for (unsigned int i = 0; i <= N; i++)
    {
        for (unsigned int j = 0; j <= (N-i); j++)
        {
            Vec ddr, dds;
            GradSimplex2D(A, B, i, j, ddr, dds);

            PetscScalar *ddr_a, *dds_a;
            VecGetArray(ddr, &ddr_a);
            VecGetArray(dds, &dds_a);
            ir[0]=sk;

            MatSetValues(V2Dr, size_r, ix, 1, ir, ddr_a, INSERT_VALUES);
            MatSetValues(V2Ds, size_r, ix, 1, ir, dds_a, INSERT_VALUES);

            sk++;
            VecRestoreArray(ddr, &ddr_a);
            VecRestoreArray(dds, &dds_a);
            VecDestroy(&ddr);
            VecDestroy(&dds);
        }
    }
    VecDestroy(&A);
    VecDestroy(&B);
    MatAssemblyBegin(V2Dr, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(V2Dr, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(V2Ds, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(V2Ds, MAT_FINAL_ASSEMBLY);
}
/*--------------------------------------------------------------------------*/
extern void DMatrices2D(const unsigned int &N, const Vec &R, const Vec &S, const Mat &V, Mat &Dr, Mat &Ds)
{
  // function [Dr,Ds] = Dmatrices2D(N,r,s,V)
  // Purpose : Initialize the (r,s) differentiation matrices
  //	    on the simplex, evaluated at (r,s) at order N

    Mat Vr, Vs;
    GradVandermonde2D(N, R, S, Vr, Vs);
    MatConvert(Vr, MATSEQDENSE, MAT_INPLACE_MATRIX, &Vr);
    MatConvert(Vs, MATSEQDENSE, MAT_INPLACE_MATRIX, &Vs);

    //std::cout << "Vr = " << std::endl;
    //MatView(Vr, PETSC_VIEWER_STDOUT_SELF);

    Mat X = Inverse_Matrix(V);
    /*
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
    ISDestroy(&row);
    ISDestroy(&col);
    MatDestroy(&A);
    MatDestroy(&B);
    */
    //std::cout << "X = " << std::endl;
    //MatView(X, PETSC_VIEWER_STDOUT_SELF);
    // X is the inverse

    //Dr = Vr *inv(V)
    MatMatMult(Vr, X, MAT_INITIAL_MATRIX, 1, &Dr);
    MatMatMult(Vs, X, MAT_INITIAL_MATRIX, 1, &Ds);

    //MatDestroy(&Vr);
    MatDestroy(&X);

    MatDestroy(&Vr);
    MatDestroy(&Vs);
}
/*--------------------------------------------------------------------------*/
extern Mat InterpMatrix2D(const unsigned int &N, const Vec &R, const Vec &S)
{
  // function [IM] = InterpMatrix2D(rout, sout)
  // purpose: compute local elemental interpolation matrix

  // N = order polynomials, not number of cubature points
    Mat V2Dc = Vandermonde2D(N, R, S);

    // Build interpolation matrix (nodes->cubature nodes)
    // Need Element inverse Vandermonde matrix
    Mat VInv = load_InverseVandermondeMatrix(N);
    MatConvert(VInv, MATSEQDENSE, MAT_INPLACE_MATRIX, &VInv);
    Mat Interp;
    MatMatMult(V2Dc, VInv,  MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Interp);

    MatDestroy(&V2Dc);
    MatDestroy(&VInv);
    return Interp;
}
/*--------------------------------------------------------------------------*/
/* extern Vec LagrangePolynomial2D(const Mat &V, const Vec &P, const unsigned int &Np)
{
    Vec L;
    VecCreateSeq(PETSC_COMM_SELF, Np, &L);
    //VecDuplicate(P, &L);

    Mat VT;
    MatTranspose(V,MAT_INITIAL_MATRIX, &VT);

    std::cout << "V = " << std::endl;
    MatView(V,  PETSC_VIEWER_STDOUT_SELF);
    std::cout << "VT = " << std::endl;
    MatView(VT,  PETSC_VIEWER_STDOUT_SELF);

    KSP                solver;
    PC                 prec;
    KSPConvergedReason reason;
    PetscInt           i,j,its;
    PetscErrorCode     ierr;
    MPI_Comm comm = MPI_COMM_SELF;

    std::cout << " 1 " << std::endl;
    KSPCreate(comm,&solver);
    std::cout << " 2 " << std::endl;
    KSPSetOperators(solver,VT,VT);
    std::cout << " 3 " << std::endl;

    KSPSetType(solver, KSPPREONLY);
    std::cout << " 4 " << std::endl;
    KSPGetPC(solver,&prec);
    std::cout << " 4.5 " << std::endl;
    PCSetType(prec, PCLU);
    std::cout << " 5 " << std::endl;

    //KSPSetFromOptions(solver);
    std::cout << " 6 " << std::endl;
    KSPSetUp(solver);
    std::cout << " 7 " << std::endl;

                    PetscInt size_V1, size_V2;
                    MatGetSize(VT, &size_V1, &size_V2);
                    PetscInt size_P, size_L;
                    VecGetSize(P, &size_P);
                    VecGetSize(L, &size_L);

                    std::cout << "matrix size VT = " << size_V1 << " x " << size_V2 << std::endl;
                    MatGetSize(V, &size_V1, &size_V2);
                    std::cout << "matrix size V = " << size_V1 << " x " << size_V2 << std::endl;
                    std::cout << "vec P size = " << size_P << std::endl;
                    std::cout << "vec L size = " << size_L << std::endl;

    KSPSolve(solver,P,L);
    std::cout << " 8 " << std::endl;

    KSPDestroy(&solver);
    MatDestroy(&VT);

    return L;

}*/
/*--------------------------------------------------------------------------*/
extern Mat Inverse_Matrix(const Mat &V)
{
    // Calculate inverse Vandermonde Matrix
    Mat A, B, X;
    MatDuplicate(V,MAT_COPY_VALUES,&A);
    MatDuplicate(V,MAT_DO_NOT_COPY_VALUES,&X);
    MatDuplicate(V,MAT_DO_NOT_COPY_VALUES,&B);
    MatConvert(A, MATSEQDENSE, MAT_INPLACE_MATRIX, &A);
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
    //std::cout << "X = " << std::endl;
    //MatView(X, PETSC_VIEWER_STDOUT_SELF);
    // X is the inverse
    MatDestroy(&A);
    MatDestroy(&B);
    ISDestroy(&row);
    ISDestroy(&col);


    return X;
}
/*--------------------------------------------------------------------------*/
extern Mat load_VandermondeMatrix(const unsigned int N)
{
    Mat V;

    std::string LocationName = "Vandermonde2D/";
    LocationName.append(std::to_string(N));
    LocationName.append("V2D.dat");
    PetscViewer    viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_READ,&viewer);
    MatCreate(PETSC_COMM_WORLD,&V);
    MatLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    return V;
}
/*--------------------------------------------------------------------------*/
extern Mat load_InverseVandermondeMatrix(const unsigned int N)
{
    Mat VInv;

    std::string LocationName = "Vandermonde2D/";
    LocationName.append(std::to_string(N));
    LocationName.append("V2DInv.dat");
    PetscViewer    viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_READ,&viewer);
    MatCreate(PETSC_COMM_WORLD,&VInv);
    MatLoad(VInv,viewer);
    PetscViewerDestroy(&viewer);

    return VInv;
}
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix2D(const unsigned int &N)
{
    Mat Product;
    Mat VInv = load_InverseVandermondeMatrix(N);
    MatTransposeMatMult(VInv, VInv, MAT_INITIAL_MATRIX, 1.0, &Product);

    MatDestroy(&VInv);
    return Product;
}
/*--------------------------------------------------------------------------*/
extern Mat InverseMassMatrix2D(const unsigned int &N)
{
    Mat Product;
    Mat V = load_VandermondeMatrix(N);
    MatMatTransposeMult(V, V, MAT_INITIAL_MATRIX, 1.0, &Product);

    MatDestroy(&V);
    return Product;
}
/*--------------------------------------------------------------------------*/
extern Mat MassMatrix2D_Cubature(const unsigned int &N, const Mat &cubV, const Vec &cubW, const unsigned int &Ncub)
{
    Mat Product, Product2;

    Mat W;
    //MatCreateSeqAIJ(PETSC_COMM_SELF, Ncub, Ncub,1, NULL,&W);
    MatCreateSeqDense(PETSC_COMM_SELF, Ncub, Ncub, NULL,&W);
    MatDiagonalSet(W,cubW,INSERT_VALUES );

    //std::cout << "W = " << std::endl;
    //MatView(W, PETSC_VIEWER_STDOUT_SELF);
    //MatTransposeMatMult(Product, W, MAT_INITIAL_MATRIX, 1.0, &Product2);
    //std::cout << "Product2 = " << std::endl;
    //MatView(Product2, PETSC_VIEWER_STDOUT_SELF);
    //Mat Product3;
    //MatMatMult(Product2, Product, MAT_INITIAL_MATRIX, 1.0, &Product3);
    //std::cout << "Product3 = " << std::endl;
    //MatView(Product3, PETSC_VIEWER_STDOUT_SELF);


    MatTransposeMatMult(cubV, W, MAT_INITIAL_MATRIX, 1.0, &Product);
    MatMatMult(Product, cubV, MAT_INITIAL_MATRIX, 1.0, &Product2);

    MatDestroy(&Product);
    //MatDestroy(&Product2);
    //MatDestroy(&VInv);
    MatDestroy(&W);
    return Product2;
}
/*--------------------------------------------------------------------------*/
extern void store_LagrangePolynomial_Cubature()
{
    // Not Langrange, but Legendre
    for (unsigned int cubatureOrder = 1; cubatureOrder <= 28; cubatureOrder++)
    {
        Vec cubR, cubS, cubW;
        unsigned int Ncub;
        Cubature2D(cubatureOrder, cubR, cubS, cubW, Ncub);
        Vec cubA, cubB;
        RStoAB(cubR, cubS, cubA, cubB);

        PetscInt cubin[Ncub];
        for (unsigned int n=0;n<Ncub; n++)
        {
            cubin[n] = n;
        }

        for (unsigned int Order_Polynomials = 0; Order_Polynomials <= 20; Order_Polynomials++)
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            PetscInt sk[1]={0};
            unsigned int Np = (Order_Polynomials+1)*(Order_Polynomials+2)/2;

            Mat L;
            MatCreate(PETSC_COMM_WORLD,&L);
            MatSetType(L,MATSEQAIJ);
            MatSetSizes(L, Ncub, Np, PETSC_DECIDE, PETSC_DECIDE);
            MatSeqAIJSetPreallocation(L, Np, NULL);

            for (unsigned int i = 0; i <= Order_Polynomials; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials-i; j++)
                {
                    Vec P;
                    P = Simplex2DP(cubA, cubB, i, j);
                    PetscScalar *P_a;
                    VecGetArray(P, &P_a);
                    MatSetValues(L, Ncub, cubin, 1, sk, P_a, INSERT_VALUES);
                    VecRestoreArray(P, &P_a);
                    VecDestroy(&P);
                    sk[0] += 1;
                }
            }
            MatAssemblyBegin(L, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(L, MAT_FINAL_ASSEMBLY);

            Mat LTranspose;
            MatTranspose(L, MAT_INITIAL_MATRIX, &LTranspose);

            std::string LocationName = "CubatureLagrange/N";
            LocationName.append(std::to_string(Order_Polynomials));
            LocationName.append("CO");
            LocationName.append(std::to_string(cubatureOrder));
            LocationName.append(".dat");
            PetscViewer    viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_WRITE,&viewer);
            MatView(LTranspose,viewer);
            PetscViewerDestroy(&viewer);


            auto t2 = std::chrono::high_resolution_clock::now();
            std::cout << "Cubature Order = " << cubatureOrder << ", Polynomial Order = " << Order_Polynomials << " Execution took "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
                  << " milliseconds\n";
            MatDestroy(&L);
            MatDestroy(&LTranspose);
        }

            VecDestroy(&cubR);
            VecDestroy(&cubS);
            VecDestroy(&cubW);
            VecDestroy(&cubA);
            VecDestroy(&cubB);
    }
}
/*--------------------------------------------------------------------------*/
extern void store_DerivativeLagrangePolynomial_Cubature()
{
    // Not Langrange, but Legendre
    for (unsigned int cubatureOrder = 1; cubatureOrder <= 28; cubatureOrder++)
    {
        Vec cubR, cubS, cubW;
        unsigned int Ncub;
        Cubature2D(cubatureOrder, cubR, cubS, cubW, Ncub);
        Vec cubA, cubB;
        RStoAB(cubR, cubS, cubA, cubB);
        VecDestroy(&cubR);
        VecDestroy(&cubS);

        PetscInt cubin[Ncub];
        for (unsigned int n=0;n<Ncub; n++)
        {
            cubin[n] = n;
        }

        for (unsigned int Order_Polynomials = 0; Order_Polynomials <= 20; Order_Polynomials++)
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            PetscInt sk[1]={0};
            unsigned int Np = (Order_Polynomials+1)*(Order_Polynomials+2)/2;

            Mat dLdr;
            MatCreate(PETSC_COMM_WORLD,&dLdr);
            MatSetType(dLdr,MATSEQAIJ);
            MatSetSizes(dLdr, Ncub, Np, PETSC_DECIDE, PETSC_DECIDE);
            MatSeqAIJSetPreallocation(dLdr, Np, NULL);
            Mat dLds;
            MatCreate(PETSC_COMM_WORLD,&dLds);
            MatSetType(dLds,MATSEQAIJ);
            MatSetSizes(dLds, Ncub, Np, PETSC_DECIDE, PETSC_DECIDE);
            MatSeqAIJSetPreallocation(dLds, Np, NULL);
            for (unsigned int i = 0; i <= Order_Polynomials; i++)
            {
                for (unsigned int j = 0; j <= Order_Polynomials-i; j++)
                {
                    Vec dPdr, dPds;
                    GradSimplex2D(cubA, cubB, i, j, dPdr, dPds);

                    PetscScalar *dPdr_a;
                    VecGetArray(dPdr, &dPdr_a);
                    MatSetValues(dLdr, Ncub, cubin, 1, sk, dPdr_a, INSERT_VALUES);
                    VecRestoreArray(dPdr, &dPdr_a);
                    VecDestroy(&dPdr);

                    PetscScalar *dPds_a;
                    VecGetArray(dPds, &dPds_a);
                    MatSetValues(dLds, Ncub, cubin, 1, sk, dPds_a, INSERT_VALUES);
                    VecRestoreArray(dPds, &dPds_a);
                    VecDestroy(&dPds);

                    sk[0] += 1;
                }
            }

            MatAssemblyBegin(dLdr, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(dLdr, MAT_FINAL_ASSEMBLY);
            MatAssemblyBegin(dLds, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(dLds, MAT_FINAL_ASSEMBLY);

            Mat dLdrTranspose;
            MatTranspose(dLdr, MAT_INITIAL_MATRIX, &dLdrTranspose);
            Mat dLdsTranspose;
            MatTranspose(dLds, MAT_INITIAL_MATRIX, &dLdsTranspose);

            {
            std::string LocationName = "CubatureLagrangeDeriv/N";
            LocationName.append(std::to_string(Order_Polynomials));
            LocationName.append("CO");
            LocationName.append(std::to_string(cubatureOrder));
            LocationName.append("dr.dat");
            PetscViewer    viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_WRITE,&viewer);
            MatView(dLdrTranspose,viewer);
            PetscViewerDestroy(&viewer);
            }
            {
            std::string LocationName = "CubatureLagrangeDeriv/N";
            LocationName.append(std::to_string(Order_Polynomials));
            LocationName.append("CO");
            LocationName.append(std::to_string(cubatureOrder));
            LocationName.append("ds.dat");
            PetscViewer    viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_WRITE,&viewer);
            MatView(dLdsTranspose,viewer);
            PetscViewerDestroy(&viewer);
            }

            auto t2 = std::chrono::high_resolution_clock::now();
            std::cout << "Cubature Order = " << cubatureOrder << ", Polynomial Order = " << Order_Polynomials << " Execution took "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
                  << " milliseconds\n";

            MatDestroy(&dLdr);
            MatDestroy(&dLds);
            MatDestroy(&dLdrTranspose);
            MatDestroy(&dLdsTranspose);
        }
            VecDestroy(&cubW);
            VecDestroy(&cubA);
            VecDestroy(&cubB);
    }
}
/*--------------------------------------------------------------------------*/
extern Mat load_LagrangePolynomial_Cubature(const unsigned int &Order_Polynomials, const unsigned int &Cubature_Order)
{
    Mat L;
    // Read the nodes for a Reference Element
    std::string LocationName = "CubatureLagrange/N";
    LocationName.append(std::to_string(Order_Polynomials));
    LocationName.append("CO");
    LocationName.append(std::to_string(Cubature_Order));
    LocationName.append(".dat");
    PetscViewer    viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_READ,&viewer);
    MatCreate(PETSC_COMM_WORLD,&L);
    MatLoad(L,viewer);
    PetscViewerDestroy(&viewer);
    return L;
}
/*--------------------------------------------------------------------------*/
extern Mat load_DerivativeLagrangePolynomial_Cubature(const unsigned int &Order_Polynomials, const unsigned int &Cubature_Order, const bool &deriv)
{
    Mat dL;
    MatCreate(PETSC_COMM_WORLD,&dL);
        std::string LocationName = "CubatureLagrangeDeriv/N";
        LocationName.append(std::to_string(Order_Polynomials));
        LocationName.append("CO");
        LocationName.append(std::to_string(Cubature_Order));
    if (deriv == 1) // dr
    {
        LocationName.append("dr.dat");
        PetscViewer    viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_READ,&viewer);
        MatLoad(dL,viewer);
        PetscViewerDestroy(&viewer);
    }
    else  // ds
    {
        LocationName.append("ds.dat");
        PetscViewer    viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,LocationName.c_str(),FILE_MODE_READ,&viewer);
        MatLoad(dL,viewer);
        PetscViewerDestroy(&viewer);
    }
    return dL;
}
/*--------------------------------------------------------------------------*/
extern Mat VandermondeMultiply(const Mat &V2DInv, const Mat &Matrix)
{
    Mat Intermediate, Return;
    MatTransposeMatMult(V2DInv, Matrix, MAT_INITIAL_MATRIX, 1.0, &Intermediate);
    MatMatMult(Intermediate,V2DInv,MAT_INITIAL_MATRIX,1.0, &Return);
    MatDestroy(&Intermediate);
    return Return;
}
/*--------------------------------------------------------------------------*/
void GeometricFactors2D(const Vec &x, const Vec &y, const Mat &Dr, const Mat &Ds, const unsigned int &Np, Vec &rx, Vec &ry, Vec &sx, Vec &sy, double &J)
{
    // % function [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds)
    // Purpose : Compute the metric elements for the local mappings
    // of the elements

    Vec xr, xs, yr, ys;
    VecCreateSeq(PETSC_COMM_SELF,Np, &xr);
    VecCreateSeq(PETSC_COMM_SELF,Np, &xs);
    VecCreateSeq(PETSC_COMM_SELF,Np, &yr);
    VecCreateSeq(PETSC_COMM_SELF,Np, &ys);
    MatMult(Dr, x, xr);
    MatMult(Ds, x, xs);
    MatMult(Dr, y, yr);
    MatMult(Ds, y, ys);

    Vec Interm, Interm2;
    VecCreateSeq(PETSC_COMM_SELF,Np, &Interm);
    VecCreateSeq(PETSC_COMM_SELF,Np, &Interm2);
    //PetscScalar Interm;
    //VecDot(xs, yr, &Interm);
    VecPointwiseMult(Interm, xs, yr);
    //PetscScalar Interm2;
    //VecDot(xr, ys, &Interm2);
    VecPointwiseMult(Interm2, xr, ys);

    //J = -Interm+Interm2;
    Vec JVec;
    VecCreateSeq(PETSC_COMM_SELF,Np, &JVec);
    VecWAXPY(JVec,-1.0, Interm, Interm2);

    std::cout << " J = "  << std::endl;
    VecView(JVec, PETSC_VIEWER_STDOUT_SELF);
    //std::cout << " J = " << J << std::endl;
    std::cout << "xr = " << std::endl;
    VecView(xr, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "xs = " << std::endl;
    VecView(xs, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "yr = " << std::endl;
    VecView(yr, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "ys = " << std::endl;
    VecView(ys, PETSC_VIEWER_STDOUT_SELF);

    // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
    VecDuplicate(xr, &rx);
    VecDuplicate(xr, &ry);
    VecDuplicate(xr, &sx);
    VecDuplicate(xr, &sy);
    //VecCopy(ys, rx);
    //VecCopy(xs, ry);
    //VecCopy(yr, sx);
    //VecCopy(xr, sy);
    //VecScale(rx, 1.0/J);
    //VecScale(ry, -1.0/J);
    //VecScale(sx, -1.0/J);
    //VecScale(sy, 1.0/J);

    VecPointwiseDivide(rx, ys, JVec);
    VecPointwiseDivide(ry, xs, JVec);
    VecScale(ry, -1.0);
    VecPointwiseDivide(sx, yr, JVec);
    VecScale(ry, -1.0);
    VecPointwiseDivide(sy, xr, JVec);

    VecDestroy(&Interm);
    VecDestroy(&Interm2);
    VecDestroy(&JVec);
    VecDestroy(&xr);
    VecDestroy(&xs);
    VecDestroy(&yr);
    VecDestroy(&ys);
}
/*--------------------------------------------------------------------------*/
extern void stdVectorToPetscVec(const std::vector<double> VecIn, Vec &PetscVec)
{
    unsigned int Size_Vector = VecIn.size();
    //std::cout << "Vec Size = " << Size_Vector << std::endl;
    double x_arr[Size_Vector];
    std::copy(VecIn.begin(), VecIn.end(), x_arr);
    Vec xx;
    VecCreateSeqWithArray(PETSC_COMM_SELF, Size_Vector, Size_Vector, x_arr, &xx);
    //std::cout << "xx = " << std::endl;
    VecView(xx, PETSC_VIEWER_STDOUT_SELF);

    VecDuplicate(xx, &PetscVec);
    VecCopy(xx, PetscVec);

    VecDestroy(&xx);
}
/*--------------------------------------------------------------------------*/
extern void NodesSquares2D(unsigned int N, Vec &XX, Vec &YY)
{
    // function [x,y] = Nodes2D(N);
    // Purpose  : Compute (x,y) nodes in quadrilaterals for
    //             polynomial of order N
    Vec X, Y;
    if (N==0)
    {
        double x [1] = {0};
        double y [1] = {0};
        VecCreateSeqWithArray(PETSC_COMM_WORLD,1,1,x,&X);
        VecCreateSeqWithArray(PETSC_COMM_WORLD,1,1,y,&Y);
        VecAssemblyBegin(X);
        VecAssemblyEnd(X);
        VecAssemblyBegin(Y);
        VecAssemblyEnd(Y);
        VecDuplicate(X,&XX);
        VecCopy(X, XX);

        VecDuplicate(Y,&YY);
        VecCopy(Y, YY);

    }
    else
    {
        // total number of nodes
        unsigned int  Np = (N+1)*(N+1);
        double x [Np] = {};
        double y [Np] = {};
        Vec r = JacobiGL(0, 0, N);
        //std::cout << "r = " << std::endl;
        //VecView(r, PETSC_VIEWER_STDOUT_SELF);

        PetscScalar *r_a;
        VecGetArray(r, &r_a);
        for (unsigned int i = 0; i < (N+1); i++)
        {
            for (unsigned int j = 0; j < (N+1); j++)
            {
                x[i*(N+1)+j] = r_a[j];
                y[i*(N+1)+j] = r_a[i];
                //std::cout << x[i*(N+1)+j] << " " << y[i*(N+1)+j] << std::endl;
            }
        }
        VecRestoreArray(r, &r_a);
        VecDestroy(&r);
        VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,x,&X);
        VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,y,&Y);
        VecAssemblyBegin(X);
        VecAssemblyEnd(X);
        VecAssemblyBegin(Y);
        VecAssemblyEnd(Y);
        //std::cout << "X = " << std::endl;
        //VecView(X, PETSC_VIEWER_STDOUT_SELF);
        //std::cout << "Y = " << std::endl;
        //VecView(Y, PETSC_VIEWER_STDOUT_SELF);

        VecDuplicate(X,&XX);
        VecCopy(X, XX);

        VecDuplicate(Y,&YY);
        VecCopy(Y, YY);

    }

        VecDestroy(&X);
        VecDestroy(&Y);


}
/*--------------------------------------------------------------------------*/
void NodesCuboid(unsigned int Nx, unsigned int Ny, unsigned int Nz, Vec &XX, Vec &YY, Vec &ZZ)
{
    // function [x,y, z] = NodesCuboid(Nx, Ny, Nz, X, Y, Z);
    // Purpose  : Compute (x,y,z) nodes in cuboids for
    //             polynomial of order N

    // What if only one of Nx, Ny and Nz is zero?
    Vec X, Y, Z;
    if (Nx==0 && Ny==0 && Nz==0)
    {
        double x [1] = {0};
        double y [1] = {0};
        double z [1] = {0};
        VecCreateSeqWithArray(PETSC_COMM_WORLD,1,1,x,&X);
        VecCreateSeqWithArray(PETSC_COMM_WORLD,1,1,y,&Y);
        VecCreateSeqWithArray(PETSC_COMM_WORLD,1,1,z,&Z);
        VecAssemblyBegin(X);
        VecAssemblyEnd(X);
        VecAssemblyBegin(Y);
        VecAssemblyEnd(Y);
        VecAssemblyBegin(Z);
        VecAssemblyEnd(Z);
        VecDuplicate(X,&XX);
        VecCopy(X, XX);
        VecDuplicate(Y,&YY);
        VecCopy(Y, YY);
        VecDuplicate(Z,&ZZ);
        VecCopy(Z, ZZ);
    }
    else
    {
        // total number of nodes
        unsigned int  Np = (Nx+1)*(Ny+1)*(Nz+1);
        double x [Np] = {};
        double y [Np] = {};
        double z [Np] = {};
        Vec rx = JacobiGL(0, 0, Nx);
        Vec ry = JacobiGL(0, 0, Ny);
        Vec rz = JacobiGL(0, 0, Nz);
        //std::cout << "r = " << std::endl;
        //VecView(r, PETSC_VIEWER_STDOUT_SELF);

        PetscScalar *rx_a, *ry_a, *rz_a;
        VecGetArray(rx, &rx_a);
        VecGetArray(ry, &ry_a);
        VecGetArray(rz, &rz_a);
        for (unsigned int i = 0; i < (Nz+1); i++)
        {
            for (unsigned int j = 0; j < (Ny+1); j++)
            {
                for (unsigned int k = 0; k < (Nx+1); k++)
                {
                x[i*(Ny+1)*(Nx+1)+j*(Nx+1)+k] = rx_a[k];
                y[i*(Ny+1)*(Nx+1)+j*(Nx+1)+k] = ry_a[j];
                z[i*(Ny+1)*(Nx+1)+j*(Nx+1)+k] = rz_a[i];
                //std::cout << x[i*(N+1)+j] << " " << y[i*(N+1)+j] << std::endl;
                }
            }
        }
        VecRestoreArray(rx, &rx_a);
        VecRestoreArray(ry, &ry_a);
        VecRestoreArray(rz, &rz_a);
        VecDestroy(&rx);
        VecDestroy(&ry);
        VecDestroy(&rz);
        VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,x,&X);
        VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,y,&Y);
        VecCreateSeqWithArray(PETSC_COMM_WORLD,1,Np,z,&Z);
        VecAssemblyBegin(X);
        VecAssemblyEnd(X);
        VecAssemblyBegin(Y);
        VecAssemblyEnd(Y);
        VecAssemblyBegin(Z);
        VecAssemblyEnd(Z);

        VecDuplicate(X,&XX);
        VecCopy(X, XX);
        VecDuplicate(Y,&YY);
        VecCopy(Y, YY);
        VecDuplicate(Z,&ZZ);
        VecCopy(Z, ZZ);
    }
        VecDestroy(&X);
        VecDestroy(&Y);
        VecDestroy(&Z);
}
/*--------------------------------------------------------------------------*/
/*
void set_Node_Coordinates_Cuboid(std::vector<Cuboid> &List_Of_Elements, const std::vector<VertexCoordinates3D> &List_Of_Vertices, const unsigned int &Nx, const unsigned int &Ny, const unsigned int &Nz)
{
    Vec R, S, T;
    // Compute Nodes on a Reference Element (Quadrilateral)
    NodesCuboid(Nx, Ny, Nz, R, S, T);
    /*
    PetscViewer viewer_dense;
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer_dense);
    PetscViewerPushFormat(viewer_dense, PETSC_VIEWER_ASCII_DENSE);
    PetscViewerSetType(viewer_dense, PETSCVIEWERASCII);
    std::cout << "R = " << std::endl;
    VecView(R, viewer_dense);
    std::cout << "S = " << std::endl;
    VecView(S, viewer_dense);
    std::cout << "T = " << std::endl;
    VecView(T, viewer_dense);
    PetscViewerDestroy(&viewer_dense);
    * /

    PetscScalar *r_a, *s_a, *t_a;
    VecGetArray(R, &r_a);
    VecGetArray(S, &s_a);
    VecGetArray(T, &t_a);

    PetscInt size_r;
    VecGetSize(R, &size_r); // Vec R S T are same size
    for(auto i = List_Of_Elements.begin(); i < List_Of_Elements.end(); i++)
    {
        /*
        double x_v1 = List_Of_Vertices[(*i).getVertex_V1()].getxCoordinate();
        double y_v1 = List_Of_Vertices[(*i).getVertex_V1()].getyCoordinate();
        double z_v1 = List_Of_Vertices[(*i).getVertex_V1()].getzCoordinate();
        double x_v2 = List_Of_Vertices[(*i).getVertex_V2()].getxCoordinate();
        double y_v2 = List_Of_Vertices[(*i).getVertex_V2()].getyCoordinate();
        double z_v2 = List_Of_Vertices[(*i).getVertex_V2()].getzCoordinate();
        double x_v3 = List_Of_Vertices[(*i).getVertex_V3()].getxCoordinate();
        double y_v3 = List_Of_Vertices[(*i).getVertex_V3()].getyCoordinate();
        double z_v3 = List_Of_Vertices[(*i).getVertex_V3()].getzCoordinate();
        double x_v4 = List_Of_Vertices[(*i).getVertex_V4()].getxCoordinate();
        double y_v4 = List_Of_Vertices[(*i).getVertex_V4()].getyCoordinate();
        double z_v4 = List_Of_Vertices[(*i).getVertex_V4()].getzCoordinate();
        double x_v5 = List_Of_Vertices[(*i).getVertex_V5()].getxCoordinate();
        double y_v5 = List_Of_Vertices[(*i).getVertex_V5()].getyCoordinate();
        double z_v5 = List_Of_Vertices[(*i).getVertex_V5()].getzCoordinate();
        double x_v6 = List_Of_Vertices[(*i).getVertex_V6()].getxCoordinate();
        double y_v6 = List_Of_Vertices[(*i).getVertex_V6()].getyCoordinate();
        double z_v6 = List_Of_Vertices[(*i).getVertex_V6()].getzCoordinate();
        double x_v7 = List_Of_Vertices[(*i).getVertex_V7()].getxCoordinate();
        double y_v7 = List_Of_Vertices[(*i).getVertex_V7()].getyCoordinate();
        double z_v7 = List_Of_Vertices[(*i).getVertex_V7()].getzCoordinate();
        double x_v8 = List_Of_Vertices[(*i).getVertex_V8()].getxCoordinate();
        double y_v8 = List_Of_Vertices[(*i).getVertex_V8()].getyCoordinate();
        double z_v8 = List_Of_Vertices[(*i).getVertex_V8()].getzCoordinate();
        * /

        for(unsigned int k = 0; k < size_r; k++)
        {
            unsigned int type_node = 0;
            // Boundary of Elements
            if (abs(t_a[k]+1.0) < NODETOL || Nz == 0)
            {
                // face 0
                type_node = 0;
                (*i).set_node_on_face0(k);
            }
            if (abs(t_a[k]-1.0) < NODETOL || Nz == 0)
            {
                // face 1
                type_node = 1;
                (*i).set_node_on_face1(k);
            }
            if (abs(s_a[k]+1.0) < NODETOL || Ny == 0)
            {
                // face 2
                type_node = 2;
                (*i).set_node_on_face2(k);
            }
            if (abs(s_a[k]-1.0) < NODETOL || Ny == 0)
            {
                // face 3
                type_node = 3;
                (*i).set_node_on_face3(k);
            }
            if (abs(r_a[k]+1.0) < NODETOL || Nx == 0)
            {
                // face 4
                type_node = 4;
                (*i).set_node_on_face4(k);
            }
            if (abs(r_a[k]-1.0) < NODETOL || Nx == 0)
            {
                // face 5
                type_node = 5;
                (*i).set_node_on_face5(k);
            }
        }
    }
    VecRestoreArray(R, &r_a);
    VecRestoreArray(S, &s_a);
    VecRestoreArray(R, &t_a);

    VecDestroy(&R);
    VecDestroy(&S);
    VecDestroy(&T);
}*/
/*--------------------------------------------------------------------------*/
