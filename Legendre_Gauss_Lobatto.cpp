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

    VecView(x, 	PETSC_VIEWER_STDOUT_SELF);

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
