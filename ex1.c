static char help[] = "Solves a 1D test case for a stratified fluid.\n\n";

#include <petscksp.h>
static double calculate_sigma(double N2, unsigned int k);
static double Exact_Solution_m(double x, double t, double N2, double sigma, unsigned int k);
static double Exact_Solution_p(double x, double t, double N2, double sigma, unsigned int k);
static double Inital_Condition_m(double x1, double x2, double t, double N2, double sigma, unsigned int k);
static double Inital_Condition_p(double x1, double x2, double t, double N2, double sigma, unsigned int k);
static double Background_Density(double x, double N2);
static double Background_Density_Derivative(double x, double N2);
static double Backgound_Density_Finite_Volume(double x1, double x2, double N2);
static double Hamiltonian_FV(Vec Sol, double DeltaX, Vec R0);
static double Hamiltonian_continuous(Vec Sol, double DeltaX, Vec One_Over_R0);
static double calculate_Energy(Vec Sol, double DeltaX, Vec rho0_approximation, int Method);
static double One_Over_Background_Density(double x, double N2);
static double Hamiltonian_FV_SV(Vec VecU, Vec VecP, double DeltaX, Vec R0);

static double calculate_sigma(double N2, unsigned int k)
{
    return sqrt(1/4.0*N2*N2+2.0*PETSC_PI*2.0*PETSC_PI*k*k)/2.0/PETSC_PI;
}
static double Exact_Solution_m(double x, double t, double N2, double sigma, unsigned int k)
{
    return exp(-0.5*N2*x)*sin(2.0*PETSC_PI*k*x)*sin(2.0*PETSC_PI*sigma*(t+0.0));
}
static double Exact_Solution_p(double x, double t, double N2, double sigma, unsigned int k)
{
    return exp(-0.5*N2*x)*(N2/4.0/sigma/PETSC_PI*sin(2.0*PETSC_PI*k*x)+k/sigma*cos(2.0*PETSC_PI*k*x))*cos(2.0*PETSC_PI*sigma*(t+0.0));
}
static double Inital_Condition_m(double x1, double x2, double t, double N2, double sigma, unsigned int k)
{
    // int m dx = -(exp(-(N2*x)/2)*sin(2*pi*sigma*(t + 1/10))*((N2*sin(2*pi*k*x))/2 + 2*k*pi*cos(2*pi*k*x)))/(N2^2/4 + 4*pi^2*k^2)
    // int m dx = -(sin(2*pi*sigma*(t + 1/10))*exp(-(N2*x)/2)*((N2*sin(2*pi*k*x))/2 + 2*k*pi*cos(2*pi*k*x)))/(N2^2/4 + 4*pi^2*k^2)
    return -(exp(-(N2*x2)/2.0)*sin(2.0*PETSC_PI*sigma*(t+0.0))*((N2*sin(2.0*PETSC_PI*k*x2))/2.0 + 2.0*k*PETSC_PI*cos(2.0*PETSC_PI*k*x2)))/sigma/sigma+(exp(-(N2*x1)/2.0)*sin(2.0*PETSC_PI*sigma*(t+0.0))*((N2*sin(2.0*PETSC_PI*k*x1))/2.0 + 2.0*k*PETSC_PI*cos(2.0*PETSC_PI*k*x1)))/sigma/sigma;
}
static double Inital_Condition_p(double x1, double x2, double t, double N2, double sigma, unsigned int k)
{
    // int p dx = -(exp(-(N2*x)/2)*cos(2*pi*sigma*t)*(N2^2*sin(2*pi*k*x) - 16*k^2*pi^2*sin(2*pi*k*x) + 8*N2*k*pi*cos(2*pi*k*x)))/(2*sigma*pi*(N2^2 + 16*pi^2*k^2))
    return -(exp(-(N2*x2)/2.0)*cos(2.0*PETSC_PI*sigma*(t+0.0))*(N2*N2*sin(2.0*PETSC_PI*k*x2) - 16.0*k*k*PETSC_PI*PETSC_PI*sin(2*PETSC_PI*k*x2) + 8.0*N2*k*PETSC_PI*cos(2.0*PETSC_PI*k*x2)))/(2.0*sigma*PETSC_PI*(N2*N2 + 16.0*PETSC_PI*PETSC_PI*k*k)) + (exp(-(N2*x1)/2.0)*cos(2.0*PETSC_PI*sigma*(t+0.0))*(N2*N2*sin(2.0*PETSC_PI*k*x1) - 16.0*k*k*PETSC_PI*PETSC_PI*sin(2.0*PETSC_PI*k*x1) + 8.0*N2*k*PETSC_PI*cos(2.0*PETSC_PI*k*x1)))/(2.0*sigma*PETSC_PI*(N2*N2 + 16.0*PETSC_PI*PETSC_PI*k*k));
}
static double Background_Density(double x, double N2)
{
    double b00 = 1.0;
    return b00*exp(-N2*x);
}
static double One_Over_Background_Density(double x, double N2)
{
    double b00 = 1.0;
    if (N2!=0)
        return b00*exp(N2*x)/N2;
    else
        return 0.0;
}
static double Background_Density_Derivative(double x, double N2)
{
    double b00 = 1.0;
    return -N2*b00*exp(-N2*x);
}
static double Backgound_Density_Finite_Volume(double x1, double x2, double N2)
{
    double b00 = 1.0;
    if (N2!=0)
        return b00/N2*(exp(-N2*x1)-exp(-N2*x2));
    else
        return 1.0;
}
static double Hamiltonian_FV(Vec Sol, double DeltaX, Vec R0)
{
    unsigned int i;
    double result=0.0;
    for (i=0; i<1/DeltaX; i++)
    {
      PetscScalar rho0[1];
      PetscInt ix[1] = {i};
      VecGetValues(R0,1,ix,rho0);
      PetscScalar M[1];
      PetscScalar P[1];
      PetscInt ixx[1] = {i+1/DeltaX};
      VecGetValues(Sol,1,ix,M);
      VecGetValues(Sol,1,ixx,P);

      result+= 0.5*DeltaX/rho0[0]*(M[0]*M[0]+P[0]*P[0]);
    }
    return result;
}
static double Hamiltonian_FV_SV(Vec VecU, Vec VecP, double DeltaX, Vec R0)
{
    unsigned int i;
    double result=0.0;
    for (i=0; i<1/DeltaX; i++)
    {
      PetscScalar rho0[1];
      PetscInt ix[1] = {i};
      VecGetValues(R0,1,ix,rho0);
      PetscScalar M[1];
      PetscScalar P[1];
      VecGetValues(VecU,1,ix,M);
      VecGetValues(VecP,1,ix,P);

      result+= 0.5*DeltaX/rho0[0]*(M[0]*M[0]+P[0]*P[0]);
    }
    return result;
}
static double Hamiltonian_continuous(Vec Sol, double DeltaX, Vec One_Over_R0)
{
    unsigned int i;
    double result=0.0;
    for (i=0; i<1/DeltaX; i++)
    {
      PetscScalar oneoverrho0[1];
      PetscInt ix[1] = {i};
      VecGetValues(One_Over_R0,1,ix,oneoverrho0);
      PetscScalar M[1];
      PetscScalar P[1];
      PetscInt ixx[1] = {i+1/DeltaX};
      VecGetValues(Sol,1,ix,M);
      VecGetValues(Sol,1,ixx,P);

      result+= 0.5*oneoverrho0[0]*(M[0]*M[0]+P[0]*P[0]);
    }
    return result;
}
static double calculate_Energy(Vec Sol, double DeltaX, Vec rho0_approximation, int Method)
{
    switch (Method)
    {
        case 0:
         break;
        case 1:
         return Hamiltonian_FV(Sol, DeltaX, rho0_approximation);
         break;
        case 2:
         return Hamiltonian_continuous(Sol, DeltaX, rho0_approximation);
         break;
        case 3:
         break;
        default:
         break;
    }
    return -1;
}


int main(int argc,char **args)
{
  KSP            ksp;          /* linear solver context */
  PC             pc;           /* preconditioner context */
  PetscReal      norm;         /* norm of solution error */
  PetscInt       i,its;

  PetscInitialize(&argc,&args,(char*)0,help);
  Mat           P,Q;
  Vec           Initial_Condition, Exact_Solution;
  PetscInt      Number_Of_Elements=10, Number_Of_TimeSteps_In_One_Period=10, Method=1;

  PetscInt   Number_Of_Periods=1;
  PetscScalar N2 = 1.0;
  PetscScalar   theta = 0.5;
  PetscOptionsGetInt(NULL, NULL, "-n", &Number_Of_Elements, NULL);
  PetscOptionsGetInt(NULL, NULL, "-t", &Number_Of_TimeSteps_In_One_Period, NULL);
  PetscOptionsGetInt(NULL, NULL, "-P", &Number_Of_Periods, NULL);
  PetscOptionsGetInt(NULL, NULL, "-Method", &Method, NULL);
  PetscOptionsGetScalar(NULL, NULL, "-N2", &N2, NULL);
  PetscOptionsGetScalar(NULL, NULL, "-theta", &theta, NULL);
  if (!( Method == 1 || Method == 2 || Method == 3 || Method == 4 || Method == 5 || Method == 6))
  {
    PetscPrintf(PETSC_COMM_SELF, "Incorrect Method. Please use:\n Method == 1 => Implicit Midpoint, Finite Volume Approximation rho_0\n Method == 2 => Implicit Midpoint, continuous rho_0, quadrature\n Method == 3 => Implicit Midpoint, continuous rho_0, exact integration\n                                    \n Method == 4 => Stormer-Verlet, Finite Volume Approximation rho_0\n Method == 5 => Stormer-Verlet, continuous rho_0, quadrature\n Method == 6 => Stormer-Verlet, continuous rho_0, exact integration\n");
    return -1;
  }

  PetscScalar   sigma;
  PetscInt      k = 1;
  sigma = calculate_sigma(N2, k);
  ///PetscPrintf(PETSC_COMM_SELF,"Frequency %6.4e\n",(double)sigma);
  PetscScalar   DeltaX = 1.0/(double)Number_Of_Elements, DeltaT=1.0/(double)Number_Of_TimeSteps_In_One_Period/sigma;

  // Viewer Full Digits
  PetscViewer viewer;
  PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  PetscViewerSetType(viewer, PETSCVIEWERASCII);

  VecCreateSeq(PETSC_COMM_WORLD, 2*Number_Of_Elements,&Exact_Solution);
  for (i=0; i<Number_Of_Elements; i++)
  {
        VecSetValue(Exact_Solution,i,(double)Exact_Solution_m(((double)i+0.5)/Number_Of_Elements,0,N2,sigma,1), INSERT_VALUES);
        VecSetValue(Exact_Solution,Number_Of_Elements+i,(double)Exact_Solution_p(((double)i+0.5)/Number_Of_Elements,0,N2,sigma,1), INSERT_VALUES);
  }
  VecAssemblyBegin(Exact_Solution);
  VecAssemblyEnd(Exact_Solution);
  ///VecView(Exact_Solution,viewer);

  // Initial Condition
  VecCreateSeq(PETSC_COMM_WORLD, 2*Number_Of_Elements,&Initial_Condition);
  Vec VecU, VecP;
  VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements, &VecU);
  VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements, &VecP);
  for (i=0; i<Number_Of_Elements; i++)
  {
        double value = (double)Inital_Condition_m((double)i*DeltaX, (double)(i+1.0)*DeltaX, 0, N2, sigma, k)/DeltaX;
        VecSetValue(Initial_Condition,i,value, INSERT_VALUES);
        VecSetValue(VecU,i,value, INSERT_VALUES);
        value = (double)Inital_Condition_p((double)i*DeltaX, (double)(i+1.0)*DeltaX, 0, N2, sigma, k)/DeltaX;
        VecSetValue(Initial_Condition,Number_Of_Elements+i,value, INSERT_VALUES);
        VecSetValue(VecP,i,value, INSERT_VALUES);
  }
  VecAssemblyBegin(Initial_Condition);
  VecAssemblyEnd(Initial_Condition);
  VecAssemblyBegin(VecU);
  VecAssemblyEnd(VecU);
  VecAssemblyBegin(VecP);
  VecAssemblyEnd(VecP);
  ///VecView(Initial_Condition,viewer);

  Vec R0;
  VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements, &R0);
  Vec One_Over_rho0;
  VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements, &One_Over_rho0);
  Vec drho_0dz;
  VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements, &drho_0dz);

  Vec Sol, RHS;
  VecDuplicate(Initial_Condition, &Sol);
  VecDuplicate(Initial_Condition, &RHS);
  VecCopy(Initial_Condition, Sol);

  // Implicit Midpoint
  if (Method==1)
  {
  // Finite Volume Approximation of rho_0
  for (i=0;i<Number_Of_Elements;i++)
  {
    VecSetValue(R0, i, (double) Backgound_Density_Finite_Volume((double)i*DeltaX, (double)(i+1.0)*DeltaX, N2)/DeltaX, INSERT_VALUES);
  }
  VecAssemblyBegin(R0);
  VecAssemblyEnd(R0);
  //VecView(R0,viewer);

  MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*Number_Of_Elements, 2*Number_Of_Elements, 4,NULL,&P);
  // Interior Points
  for (i=1; i<Number_Of_Elements-1; i++)
  {
    // U Equations
    MatSetValue(P, i, i, 1.0, INSERT_VALUES);
    MatSetValue(P, i, Number_Of_Elements+i-1, theta*DeltaT/2.0/DeltaX, INSERT_VALUES);
    MatSetValue(P, i, Number_Of_Elements+i, (1.0-2.0*theta)*DeltaT/2.0/DeltaX, INSERT_VALUES);
    MatSetValue(P, i, Number_Of_Elements+i+1, -(1.0-theta)*DeltaT/2.0/DeltaX, INSERT_VALUES);
    // P Equations
    MatSetValue(P, Number_Of_Elements+i, Number_Of_Elements+i, 1.0, INSERT_VALUES);
    PetscScalar rho0[3];
    PetscInt ix[3] = {i-1,i,i+1};
    VecGetValues(R0,3,ix,rho0);
    MatSetValue(P, Number_Of_Elements+i, i-1, (1.0-theta)*DeltaT/2.0/DeltaX*rho0[1]/rho0[0], INSERT_VALUES);
    MatSetValue(P, Number_Of_Elements+i, i, -(1.0-2.0*theta)*DeltaT/2.0/DeltaX, INSERT_VALUES);
    MatSetValue(P, Number_Of_Elements+i, i+1, -theta*DeltaT/2.0/DeltaX*rho0[1]/rho0[2], INSERT_VALUES);
  }
  // Left Boundary
  {
  MatSetValue(P, 0, 0, 1.0, INSERT_VALUES);
  MatSetValue(P, 0, Number_Of_Elements, (1.0-theta)*DeltaT/2.0/DeltaX, INSERT_VALUES);
  MatSetValue(P, 0, Number_Of_Elements+1, -(1.0-theta)*DeltaT/2.0/DeltaX, INSERT_VALUES);
  MatSetValue(P, Number_Of_Elements,Number_Of_Elements, 1.0, INSERT_VALUES);
  MatSetValue(P, Number_Of_Elements, 0, -(1.0-theta)*DeltaT/2.0/DeltaX, INSERT_VALUES);
  PetscScalar rho0[2];
  PetscInt ix[2] = {0,1};
  VecGetValues(R0,2,ix,rho0);
  MatSetValue(P, Number_Of_Elements, 1, -theta*DeltaT/2.0/DeltaX*rho0[0]/rho0[1], INSERT_VALUES);
  }
  //Right Boundary
  {
  MatSetValue(P,Number_Of_Elements-1,Number_Of_Elements-1, 1.0, INSERT_VALUES);
  MatSetValue(P,Number_Of_Elements-1,2*Number_Of_Elements-2, theta*DeltaT/2.0/DeltaX, INSERT_VALUES);
  MatSetValue(P,Number_Of_Elements-1,2*Number_Of_Elements-1, -theta*DeltaT/2.0/DeltaX, INSERT_VALUES);
  MatSetValue(P,2*Number_Of_Elements-1, 2*Number_Of_Elements-1, 1.0, INSERT_VALUES);
  PetscScalar rho02[2];
  PetscInt ixx[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
  VecGetValues(R0,2,ixx,rho02);
  MatSetValue(P,2*Number_Of_Elements-1, Number_Of_Elements-2, (1.0-theta)*DeltaT/2.0/DeltaX*rho02[1]/rho02[0], INSERT_VALUES);
  MatSetValue(P,2*Number_Of_Elements-1, Number_Of_Elements-1, theta*DeltaT/2.0/DeltaX, INSERT_VALUES);
  }
  MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);
  //MatView(P, viewer);
  MatDuplicate(P,MAT_COPY_VALUES,&Q);
  MatScale(Q, -1.0);
  MatShift(Q, 2.0*1.0);
  MatAssemblyBegin(Q,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Q,MAT_FINAL_ASSEMBLY);
  //MatView(Q, viewer);
  }
  // Continuous rho_0 Simpson's rule
  if (Method==2)
  {
  // Simpson's rule of continuous 1/rho_0
  for (i=0;i<Number_Of_Elements;i++)
  {
    double value = 0.0;
    value += 1.0/Background_Density((double)i*DeltaX, N2)+4.0/Background_Density(((double)(i)+0.5)*DeltaX, N2)+1.0/Background_Density((double)(i+1.0)*DeltaX, N2);
    value *= DeltaX/6.0;
    VecSetValue(One_Over_rho0, i, (double)value, INSERT_VALUES);

    value = 0.0;
    value += 1.0*Background_Density_Derivative((double)i*DeltaX, N2)+4.0*Background_Density_Derivative(((double)(i)+0.5)*DeltaX, N2)+1.0*Background_Density_Derivative((double)(i+1.0)*DeltaX, N2);
    value *= DeltaX/6.0;
    //value += Background_Density((double)(i+1.0)*DeltaX, N2) - Background_Density((double)(i)*DeltaX, N2);
    //value *= DeltaX;
    VecSetValue(drho_0dz, i, (double)value, INSERT_VALUES);

    value = 0.0;
    value = Background_Density((double)(i)*DeltaX, N2);
    VecSetValue(R0, i, (double)value, INSERT_VALUES);
  }
  VecAssemblyBegin(One_Over_rho0);
  VecAssemblyEnd(One_Over_rho0);
  VecAssemblyBegin(drho_0dz);
  VecAssemblyEnd(drho_0dz);
  VecAssemblyBegin(R0);
  VecAssemblyEnd(R0);
  //VecView(R0, viewer);
  //VecView(One_Over_rho0, viewer);
  //PetscPrintf(PETSC_COMM_SELF, "drho/dz Begin\n");
  //VecView(drho_0dz, viewer);
  //PetscPrintf(PETSC_COMM_SELF, "drho/dz End\n");

  MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*Number_Of_Elements, 2*Number_Of_Elements, 4,NULL,&P);
  // Interior Points
  int a0 = -1, a1 = 1, a2 = 1, a3 = 1, a4 = 1;
  //DeltaX=sqrt(DeltaX);
  for (i=1; i<Number_Of_Elements-1; i++)
  {
    MatSetValue(P, i, i, 1.0, INSERT_VALUES);
    PetscScalar oneoverrho0[3];
    PetscInt ix[3] = {i-1,i,i+1};
    VecGetValues(One_Over_rho0,3,ix,oneoverrho0);
    PetscScalar rho0[2];
    PetscInt iix[2] = {i,i+1};
    VecGetValues(R0,2,iix,rho0);
    PetscScalar drhodz[1];
    PetscInt iiix[1] = {i};
    VecGetValues(drho_0dz, 1, iiix, drhodz);
    // U Equations
    MatSetValue(P, i, Number_Of_Elements+i-1, a2*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0], INSERT_VALUES);
    MatSetValue(P, i, Number_Of_Elements+i, -(a4*theta*rho0[0]-a1*(1.0-theta)*rho0[1])*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1] + a0*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*drhodz[0], INSERT_VALUES);
    MatSetValue(P, i, Number_Of_Elements+i+1, -a3*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[2]*rho0[1], INSERT_VALUES);
    // P Equations
    MatSetValue(P, Number_Of_Elements+i, Number_Of_Elements+i, 1.0, INSERT_VALUES);
    MatSetValue(P, Number_Of_Elements+i, i-1, a3*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0], INSERT_VALUES);
    MatSetValue(P, Number_Of_Elements+i, i, -(a1*(1.0-theta)*rho0[1]-a4*theta*rho0[0])*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1] - a0*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*drhodz[0], INSERT_VALUES);
    MatSetValue(P, Number_Of_Elements+i, i+1, -a2*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[2]*rho0[1], INSERT_VALUES);
  }
  // Left Boundary
  {
  PetscScalar oneoverrho0[2];
  PetscInt ix[2] = {0,1};
  VecGetValues(One_Over_rho0,2,ix,oneoverrho0);
  PetscScalar rho0[1];
  PetscInt iix[1] = {1};
  VecGetValues(R0,1,iix,rho0);
  PetscScalar drhodz[1];
  PetscInt iiix[1] = {0};
  VecGetValues(drho_0dz, 1, iiix, drhodz);

  MatSetValue(P, 0, 0, 1.0, INSERT_VALUES);
  MatSetValue(P, 0, Number_Of_Elements, a1*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0] + a0*oneoverrho0[0]*drhodz[0]*DeltaT/2.0/DeltaX/DeltaX, INSERT_VALUES);
  MatSetValue(P, 0, Number_Of_Elements+1, -a3*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*rho0[0], INSERT_VALUES);

  MatSetValue(P, Number_Of_Elements, Number_Of_Elements, 1.0, INSERT_VALUES);
  MatSetValue(P, Number_Of_Elements, 0, -a1*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0] - a0*oneoverrho0[0]*drhodz[0]*DeltaT/2.0/DeltaX/DeltaX, INSERT_VALUES);
  MatSetValue(P, Number_Of_Elements, 1, -a2*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*rho0[0], INSERT_VALUES);
  }
  //Right Boundary
  {
  PetscScalar oneoverrho0[2];
  PetscInt ix[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
  VecGetValues(One_Over_rho0,2,ix,oneoverrho0);
  PetscScalar rho0[1];
  PetscInt iix[1] = {Number_Of_Elements-1};
  VecGetValues(R0,1,iix,rho0);
  PetscScalar drhodz[1];
  PetscInt iiix[1] = {Number_Of_Elements-1};
  VecGetValues(drho_0dz, 1, iiix, drhodz);

  MatSetValue(P, Number_Of_Elements-1,Number_Of_Elements-1, 1.0, INSERT_VALUES);
  MatSetValue(P, Number_Of_Elements-1,2*Number_Of_Elements-2, a2*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0], INSERT_VALUES);
  MatSetValue(P, Number_Of_Elements-1,2*Number_Of_Elements-1, -a4*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*rho0[0] + a0*oneoverrho0[1]*drhodz[0]*DeltaT/2.0/DeltaX/DeltaX, INSERT_VALUES);

  MatSetValue(P, 2*Number_Of_Elements-1, 2*Number_Of_Elements-1, 1.0, INSERT_VALUES);
  MatSetValue(P, 2*Number_Of_Elements-1, Number_Of_Elements-2, a3*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0], INSERT_VALUES);
  MatSetValue(P, 2*Number_Of_Elements-1, Number_Of_Elements-1, a4*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*rho0[0] - a0*oneoverrho0[1]*drhodz[0]*DeltaT/2.0/DeltaX/DeltaX, INSERT_VALUES);
  }
  MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);
  //MatView(P, viewer);
  MatDuplicate(P,MAT_COPY_VALUES,&Q);
  MatScale(Q, -1.0);
  MatShift(Q, 2.0);
  MatAssemblyBegin(Q,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Q,MAT_FINAL_ASSEMBLY);
  //MatView(Q, viewer);
  }

  // Continuous rho_0 exact integration
  if (Method==3)
  {
  // Simpson's rule of continuous 1/rho_0
  for (i=0;i<Number_Of_Elements;i++)
  {
    double value = 0.0;
    if (N2!=0)
        value = One_Over_Background_Density((double)(i+1.0)*DeltaX,N2)-One_Over_Background_Density((double)(i)*DeltaX,N2);
    else
        value = 1.0;
    VecSetValue(One_Over_rho0, i, (double)value, INSERT_VALUES);

    value = 0.0;
    if (N2!=0)
        value = Background_Density((double)(i+1.0)*DeltaX,N2)-Background_Density((double)(i)*DeltaX,N2);
    else
        value = 0.0;
    VecSetValue(drho_0dz, i, (double)value, INSERT_VALUES);

    value = 0.0;
    if (N2!=0)
        value = Background_Density((double)(i)*DeltaX, N2);
    else
        value = 1.0;
    VecSetValue(R0, i, (double)value, INSERT_VALUES);
  }
  VecAssemblyBegin(One_Over_rho0);
  VecAssemblyEnd(One_Over_rho0);
  VecAssemblyBegin(drho_0dz);
  VecAssemblyEnd(drho_0dz);
  VecAssemblyBegin(R0);
  VecAssemblyEnd(R0);
  //VecView(R0, viewer);
  //VecView(One_Over_rho0, viewer);
  //PetscPrintf(PETSC_COMM_SELF, "drho/dz Begin\n");
  //VecView(drho_0dz, viewer);
  //PetscPrintf(PETSC_COMM_SELF, "drho/dz End\n");

  MatCreateSeqAIJ(PETSC_COMM_WORLD, 2*Number_Of_Elements, 2*Number_Of_Elements, 4,NULL,&P);
  // Interior Points
  int a0 = -1, a1 = 1, a2 = 1, a3 = 1, a4 = 1;
  //DeltaX=sqrt(DeltaX);
  for (i=1; i<Number_Of_Elements-1; i++)
  {
    MatSetValue(P, i, i, 1.0, INSERT_VALUES);
    PetscScalar oneoverrho0[3];
    PetscInt ix[3] = {i-1,i,i+1};
    VecGetValues(One_Over_rho0,3,ix,oneoverrho0);
    PetscScalar rho0[2];
    PetscInt iix[2] = {i,i+1};
    VecGetValues(R0,2,iix,rho0);
    PetscScalar drhodz[1];
    PetscInt iiix[1] = {i};
    VecGetValues(drho_0dz, 1, iiix, drhodz);
    // U Equations
    MatSetValue(P, i, Number_Of_Elements+i-1, a2*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0], INSERT_VALUES);
    MatSetValue(P, i, Number_Of_Elements+i, -(a4*theta*rho0[0]-a1*(1.0-theta)*rho0[1])*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1] + a0*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*drhodz[0], INSERT_VALUES);
    MatSetValue(P, i, Number_Of_Elements+i+1, -a3*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[2]*rho0[1], INSERT_VALUES);
    // P Equations
    MatSetValue(P, Number_Of_Elements+i, Number_Of_Elements+i, 1.0, INSERT_VALUES);
    MatSetValue(P, Number_Of_Elements+i, i-1, a3*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0], INSERT_VALUES);
    MatSetValue(P, Number_Of_Elements+i, i, -(a1*(1.0-theta)*rho0[1]-a4*theta*rho0[0])*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1] - a0*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*drhodz[0], INSERT_VALUES);
    MatSetValue(P, Number_Of_Elements+i, i+1, -a2*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[2]*rho0[1], INSERT_VALUES);
  }
  // Left Boundary
  {
  PetscScalar oneoverrho0[2];
  PetscInt ix[2] = {0,1};
  VecGetValues(One_Over_rho0,2,ix,oneoverrho0);
  PetscScalar rho0[1];
  PetscInt iix[1] = {1};
  VecGetValues(R0,1,iix,rho0);
  PetscScalar drhodz[1];
  PetscInt iiix[1] = {0};
  VecGetValues(drho_0dz, 1, iiix, drhodz);

  MatSetValue(P, 0, 0, 1.0, INSERT_VALUES);
  MatSetValue(P, 0, Number_Of_Elements, a1*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0] + a0*oneoverrho0[0]*drhodz[0]*DeltaT/2.0/DeltaX/DeltaX, INSERT_VALUES);
  MatSetValue(P, 0, Number_Of_Elements+1, -a3*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*rho0[0], INSERT_VALUES);

  MatSetValue(P, Number_Of_Elements, Number_Of_Elements, 1.0, INSERT_VALUES);
  MatSetValue(P, Number_Of_Elements, 0, -a1*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0] - a0*oneoverrho0[0]*drhodz[0]*DeltaT/2.0/DeltaX/DeltaX, INSERT_VALUES);
  MatSetValue(P, Number_Of_Elements, 1, -a2*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*rho0[0], INSERT_VALUES);
  }
  //Right Boundary
  {
  PetscScalar oneoverrho0[2];
  PetscInt ix[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
  VecGetValues(One_Over_rho0,2,ix,oneoverrho0);
  PetscScalar rho0[1];
  PetscInt iix[1] = {Number_Of_Elements-1};
  VecGetValues(R0,1,iix,rho0);
  PetscScalar drhodz[1];
  PetscInt iiix[1] = {Number_Of_Elements-1};
  VecGetValues(drho_0dz, 1, iiix, drhodz);

  MatSetValue(P, Number_Of_Elements-1,Number_Of_Elements-1, 1.0, INSERT_VALUES);
  MatSetValue(P, Number_Of_Elements-1,2*Number_Of_Elements-2, a2*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0], INSERT_VALUES);
  MatSetValue(P, Number_Of_Elements-1,2*Number_Of_Elements-1, -a4*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*rho0[0] + a0*oneoverrho0[1]*drhodz[0]*DeltaT/2.0/DeltaX/DeltaX, INSERT_VALUES);

  MatSetValue(P, 2*Number_Of_Elements-1, 2*Number_Of_Elements-1, 1.0, INSERT_VALUES);
  MatSetValue(P, 2*Number_Of_Elements-1, Number_Of_Elements-2, a3*(1.0-theta)*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[0]*rho0[0], INSERT_VALUES);
  MatSetValue(P, 2*Number_Of_Elements-1, Number_Of_Elements-1, a4*theta*DeltaT/2.0/DeltaX/DeltaX*oneoverrho0[1]*rho0[0] - a0*oneoverrho0[1]*drhodz[0]*DeltaT/2.0/DeltaX/DeltaX, INSERT_VALUES);
  }
  MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);
  //MatView(P, viewer);
  MatDuplicate(P,MAT_COPY_VALUES,&Q);
  MatScale(Q, -1.0);
  MatShift(Q, 2.0);
  MatAssemblyBegin(Q,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Q,MAT_FINAL_ASSEMBLY);
  //MatView(Q, viewer);
  }

  double Initial_Energy = 0.0;
  if (Method == 1 || Method == 4)
  {
    Initial_Energy=Hamiltonian_FV(Initial_Condition, DeltaX, R0);
  }
  if (Method == 2 || Method == 3)
  {
    Initial_Energy=Hamiltonian_continuous(Initial_Condition, DeltaX, One_Over_rho0);
  }
  //PetscPrintf(PETSC_COMM_SELF,"Initial Energy %6.4e\n",(double)Initial_Energy);
  double Energy = 0.0;


  if (Method == 1 || Method == 2 || Method == 3)
  {
  // Time Integration
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp,P,P);
  //KSPGetPC(ksp,&pc);
  //PCSetType(pc,PCJACOBI);

  KSPSetType(ksp, KSPPREONLY);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCLU);

  KSPSetTolerances(ksp,1e-12, 1e-12,PETSC_DEFAULT,PETSC_DEFAULT);
  KSPSetFromOptions(ksp);

  FILE *f = fopen("Energy.txt", "w");
  PetscViewer viewer2;
  //PetscViewerPushFormat(viewer2, PETSC_VIEWER_ASCII_MATLAB);
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "solution0.txt", &viewer2);

  //PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  //PetscViewerSetType(viewer, PETSCVIEWERASCII);
  VecView(Initial_Condition, viewer2);
  PetscViewerDestroy(&viewer2);


  fprintf(f, "%1.16e\n", Initial_Energy-Initial_Energy);
  int t;
  char szFileName[255] = {0};
  double T = 0.0;
  for (t=1;t<=Number_Of_TimeSteps_In_One_Period*Number_Of_Periods;t++)
  {
    T += DeltaT;
    ///PetscPrintf(PETSC_COMM_SELF,"Time = %6.4e\n",t*DeltaT);
    MatMult(Q, Sol, RHS);
    KSPSolve(ksp,RHS,Sol);
    PetscViewer viewer2;
  //PetscViewerPushFormat(viewer2, PETSC_VIEWER_ASCII_MATLAB);
    sprintf(szFileName, "solution%d.txt", t);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);
    VecView(Sol, viewer2);
    PetscViewerDestroy(&viewer2);
    //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    //VecCopy(Sol, RHS);
      if (Method==1)
      {
        Energy=Hamiltonian_FV(Sol, DeltaX, R0);
      }
      if (Method == 2 || Method == 3)
      {
        Energy=Hamiltonian_continuous(Sol, DeltaX, One_Over_rho0);
      }
    fprintf(f, "%1.16e\n", Energy-Initial_Energy);
  }
  fclose(f);
  PetscPrintf(PETSC_COMM_SELF,"Energy %6.4e\n",(double)Energy);
  ///PetscPrintf(PETSC_COMM_SELF,"Sigma %6.4e\n",(double)sigma );
  ///PetscPrintf(PETSC_COMM_SELF,"Time %6.4e\n",(double)T );
  PetscPrintf(PETSC_COMM_SELF,"Difference in Energy %6.4e\n",(double)Energy-Initial_Energy);
  MatDestroy(&P);
  MatDestroy(&Q);
  KSPDestroy(&ksp);
    }
  // Stormer-Verlet

  // Finite Volume Approximation of rho_0
  if (Method==4)
  {
  // Finite Volume Approximation of rho_0
  for (i=0;i<Number_Of_Elements;i++)
  {
    VecSetValue(R0, i, (double) Backgound_Density_Finite_Volume((double)i*DeltaX, (double)(i+1.0)*DeltaX, N2)/DeltaX, INSERT_VALUES);
  }
  VecAssemblyBegin(R0);
  VecAssemblyEnd(R0);
  //VecView(R0,viewer);
  Vec UIntermediate;
  //VecDuplicate(VecU, &UIntermediate);

  VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements, &UIntermediate);

  double Initial_Energy = 0.0;
  if (Method == 4)
  {
    Initial_Energy=Hamiltonian_FV_SV(VecU, VecP, DeltaX, R0);
  }

  //PetscPrintf(PETSC_COMM_SELF,"Initial Energy %6.4e\n",(double)Initial_Energy);
  double Energy = 0.0;


  FILE *f = fopen("Energy.txt", "w");
  PetscViewer viewer2;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "solution0.txt", &viewer2);
  VecView(Initial_Condition, viewer2);
  PetscViewerDestroy(&viewer2);

  fprintf(f, "%1.16e\n", Initial_Energy-Initial_Energy);
  int t;
  char szFileName[255] = {0};
  double T = 0.0;
  for (t=1;t<=Number_Of_TimeSteps_In_One_Period*Number_Of_Periods;t++)
  //t = 1;
  {
    T += DeltaT;
  // Half Time Step
  // Interior Points
  for (i=1; i<Number_Of_Elements-1; i++)
  {
    PetscScalar pold[3];
    PetscInt ix[3] = {i-1,i,i+1};
    VecGetValues(VecP,3,ix,pold);
    PetscScalar uold[1];
    PetscInt iu[1] = {i};
    VecGetValues(VecU, 1, iu, uold);
    VecSetValue(UIntermediate, i, uold[0]+DeltaT/DeltaX/2.0*(pold[2]*(1.0-theta)-pold[1]*(1.0-2.0*theta)-pold[0]*theta), INSERT_VALUES);
  }
  {
  // left boundary
  PetscScalar pold[2];
  PetscInt ix[2] = {0, 1};
  VecGetValues(VecP, 2, ix, pold);
  PetscScalar uold[1];
  PetscInt iu[1] = {0};
  VecGetValues(VecU,1,iu, uold);
  VecSetValue(UIntermediate, 0, uold[0]+DeltaT/DeltaX/2.0*(pold[1]*(1.0-theta)-pold[0]*(1.0-theta)), INSERT_VALUES);
  }
  {
  // right boundary
  PetscScalar pold[2];
  PetscInt ix[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
  VecGetValues(VecP, 2, ix, pold);
  PetscScalar uold[1];
  PetscInt iu[1] = {Number_Of_Elements-1};
  VecGetValues(VecU,1,iu, uold);
  VecSetValue(UIntermediate, Number_Of_Elements-1, uold[0]+DeltaT/DeltaX/2.0*(pold[1]*theta-pold[0]*theta), INSERT_VALUES);
  }
  VecAssemblyBegin(UIntermediate);
  VecAssemblyEnd(UIntermediate);
  // Full Time Step
  for (i=1; i<Number_Of_Elements-1; i++)
  {
    PetscScalar uold[3];
    PetscInt ix[3] = {i-1,i,i+1};
    VecGetValues(UIntermediate,3,ix,uold);
    PetscScalar pold[1];
    PetscInt ip[1] = {i};
    VecGetValues(VecP, 1, ip, pold);
    PetscScalar rho0[3];
    PetscInt irho0[3] = {i-1,i,i+1};
    VecGetValues(R0,3,irho0,rho0);
    VecSetValue(VecP, i, pold[0]+DeltaT/DeltaX*(uold[2]*theta*rho0[1]/rho0[2]+uold[1]*(1.0-2.0*theta)-uold[0]*(1.0-theta)*rho0[1]/rho0[0]), INSERT_VALUES);
  }
  {
  // left boundary
    PetscScalar uold[2];
    PetscInt ix[2] = {0,1};
    VecGetValues(UIntermediate,2,ix,uold);
    PetscScalar pold[1];
    PetscInt ip[1] = {0};
    VecGetValues(VecP, 1, ip, pold);
    PetscScalar rho0[2];
    PetscInt irho0[2] = {0,1};
    VecGetValues(R0,2,irho0,rho0);
    VecSetValue(VecP, 0, pold[0]+DeltaT/DeltaX*(uold[1]*theta*rho0[0]/rho0[1]+uold[0]*(1.0-theta)), INSERT_VALUES);
  }
  {
  // right boundary
    PetscScalar uold[2];
    PetscInt ix[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
    VecGetValues(UIntermediate,2,ix,uold);
    PetscScalar pold[1];
    PetscInt ip[1] = {Number_Of_Elements-1};
    VecGetValues(VecP, 1, ip, pold);
    PetscScalar rho0[2];
    PetscInt irho0[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
    VecGetValues(R0,2,irho0,rho0);
    VecSetValue(VecP, Number_Of_Elements-1, pold[0]+DeltaT/DeltaX*(-uold[1]*theta-uold[0]*(1.0-theta)*rho0[1]/rho0[0]), INSERT_VALUES);
  }
  VecAssemblyBegin(VecP);
  VecAssemblyEnd(VecP);
  // Half Time Step
  // Interior Points
  for (i=1; i<Number_Of_Elements-1; i++)
  {
    PetscScalar pold[3];
    PetscInt ix[3] = {i-1,i,i+1};
    VecGetValues(VecP,3,ix,pold);
    PetscScalar uold[1];
    PetscInt iu[1] = {i};
    VecGetValues(UIntermediate, 1, iu, uold);
    VecSetValue(VecU, i, uold[0]+DeltaT/DeltaX/2.0*(pold[2]*(1.0-theta)-pold[1]*(1.0-2.0*theta)-pold[0]*theta), INSERT_VALUES);
  }
  {
  // left boundary
  PetscScalar pold[2];
  PetscInt ix[2] = {0, 1};
  VecGetValues(VecP, 2, ix, pold);
  PetscScalar uold[1];
  PetscInt iu[1] = {0};
  VecGetValues(UIntermediate,1,iu, uold);
  VecSetValue(VecU, 0, uold[0]+DeltaT/DeltaX/2.0*(pold[1]*(1.0-theta)-pold[0]*(1.0-theta)), INSERT_VALUES);
  }
  {
  // right boundary
  PetscScalar pold[2];
  PetscInt ix[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
  VecGetValues(VecP, 2, ix, pold);
  PetscScalar uold[1];
  PetscInt iu[1] = {Number_Of_Elements-1};
  VecGetValues(UIntermediate,1,iu, uold);
  VecSetValue(VecU, Number_Of_Elements-1, uold[0]+DeltaT/DeltaX/2.0*(pold[1]*theta-pold[0]*theta), INSERT_VALUES);
  }
  VecAssemblyBegin(VecU);
  VecAssemblyEnd(VecU);

    PetscViewer viewer2;
    sprintf(szFileName, "solution%d.txt", t);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);

      if (Method==4)
      {
        Energy=Hamiltonian_FV_SV(VecU, VecP, DeltaX, R0);
      }

    fprintf(f, "%1.16e\n", Energy-Initial_Energy);
  for (i=0;i<Number_Of_Elements;i++)
  {
    PetscScalar u[1],p[1];
    PetscInt ix[1] = {i};
    VecGetValues(VecU, 1, ix, u);
    VecGetValues(VecP, 1, ix, p);
    VecSetValue(Sol, i, u[0], INSERT_VALUES);
    VecSetValue(Sol, Number_Of_Elements+i, p[0], INSERT_VALUES);
  }


    VecAssemblyBegin(Sol);
    VecAssemblyEnd(Sol);
    VecView(Sol, viewer2);
    PetscViewerDestroy(&viewer2);

    //PetscPrintf(PETSC_COMM_SELF,"Delta t %6.4e\n",(double)DeltaT );
    //PetscPrintf(PETSC_COMM_SELF,"Delta x %6.4e\n",(double)DeltaX );
    //VecView(Initial_Condition, viewer);
    //  VecView(UIntermediate, viewer);
  //VecView(Sol, viewer);
  }
  fclose(f);
  PetscPrintf(PETSC_COMM_SELF,"Energy %6.4e\n",(double)Energy);
  ///PetscPrintf(PETSC_COMM_SELF,"Sigma %6.4e\n",(double)sigma );
  ///PetscPrintf(PETSC_COMM_SELF,"Time %6.4e\n",(double)T );
  PetscPrintf(PETSC_COMM_SELF,"Difference in Energy %6.4e\n",(double)Energy-Initial_Energy);
  VecDestroy(&UIntermediate);


  }
  // Continuous rho_0 Simpson's rule
  if (Method==5)
  {
  // Simpson's rule of continuous 1/rho_0
  for (i=0;i<Number_Of_Elements;i++)
  {
    double value = 0.0;
    value += 1.0/Background_Density((double)i*DeltaX, N2)+4.0/Background_Density(((double)(i)+0.5)*DeltaX, N2)+1.0/Background_Density((double)(i+1.0)*DeltaX, N2);
    value *= DeltaX/6.0;
    VecSetValue(One_Over_rho0, i, (double)value, INSERT_VALUES);

    value = 0.0;
    value += 1.0*Background_Density_Derivative((double)i*DeltaX, N2)+4.0*Background_Density_Derivative(((double)(i)+0.5)*DeltaX, N2)+1.0*Background_Density_Derivative((double)(i+1.0)*DeltaX, N2);
    value *= DeltaX/6.0;
    //value += Background_Density((double)(i+1.0)*DeltaX, N2) - Background_Density((double)(i)*DeltaX, N2);
    //value *= DeltaX;
    VecSetValue(drho_0dz, i, (double)value, INSERT_VALUES);

    value = 0.0;
    value = Background_Density((double)(i)*DeltaX, N2);
    VecSetValue(R0, i, (double)value, INSERT_VALUES);
  }
  VecAssemblyBegin(One_Over_rho0);
  VecAssemblyEnd(One_Over_rho0);
  VecAssemblyBegin(drho_0dz);
  VecAssemblyEnd(drho_0dz);
  VecAssemblyBegin(R0);
  VecAssemblyEnd(R0);
  Vec UIntermediate;
  //VecDuplicate(VecU, &UIntermediate);

  VecCreateSeq(PETSC_COMM_WORLD, Number_Of_Elements, &UIntermediate);

  double Initial_Energy = 0.0;
  if (Method == 4)
  {
    Initial_Energy=Hamiltonian_FV_SV(VecU, VecP, DeltaX, R0);
  }

  //PetscPrintf(PETSC_COMM_SELF,"Initial Energy %6.4e\n",(double)Initial_Energy);
  double Energy = 0.0;


  FILE *f = fopen("Energy.txt", "w");
  PetscViewer viewer2;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "solution0.txt", &viewer2);
  VecView(Initial_Condition, viewer2);
  PetscViewerDestroy(&viewer2);

  fprintf(f, "%1.16e\n", Initial_Energy-Initial_Energy);
  int t;
  char szFileName[255] = {0};
  double T = 0.0;
  for (t=1;t<=Number_Of_TimeSteps_In_One_Period*Number_Of_Periods;t++)
  //t = 1;
  {
    T += DeltaT;
  // Half Time Step
  // Interior Points
  for (i=1; i<Number_Of_Elements-1; i++)
  {
    PetscScalar pold[3];
    PetscInt ix[3] = {i-1,i,i+1};
    VecGetValues(VecP,3,ix,pold);
    PetscScalar uold[1];
    PetscInt iu[1] = {i};
    VecGetValues(VecU, 1, iu, uold);
    VecSetValue(UIntermediate, i, uold[0]+DeltaT/DeltaX/2.0*(pold[2]*(1.0-theta)-pold[1]*(1.0-2.0*theta)-pold[0]*theta), INSERT_VALUES);
  }
  {
  // left boundary
  PetscScalar pold[2];
  PetscInt ix[2] = {0, 1};
  VecGetValues(VecP, 2, ix, pold);
  PetscScalar uold[1];
  PetscInt iu[1] = {0};
  VecGetValues(VecU,1,iu, uold);
  VecSetValue(UIntermediate, 0, uold[0]+DeltaT/DeltaX/2.0*(pold[1]*(1.0-theta)-pold[0]*(1.0-theta)), INSERT_VALUES);
  }
  {
  // right boundary
  PetscScalar pold[2];
  PetscInt ix[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
  VecGetValues(VecP, 2, ix, pold);
  PetscScalar uold[1];
  PetscInt iu[1] = {Number_Of_Elements-1};
  VecGetValues(VecU,1,iu, uold);
  VecSetValue(UIntermediate, Number_Of_Elements-1, uold[0]+DeltaT/DeltaX/2.0*(pold[1]*theta-pold[0]*theta), INSERT_VALUES);
  }
  VecAssemblyBegin(UIntermediate);
  VecAssemblyEnd(UIntermediate);
  // Full Time Step
  for (i=1; i<Number_Of_Elements-1; i++)
  {
    PetscScalar uold[3];
    PetscInt ix[3] = {i-1,i,i+1};
    VecGetValues(UIntermediate,3,ix,uold);
    PetscScalar pold[1];
    PetscInt ip[1] = {i};
    VecGetValues(VecP, 1, ip, pold);
    PetscScalar rho0[3];
    PetscInt irho0[3] = {i-1,i,i+1};
    VecGetValues(R0,3,irho0,rho0);
    VecSetValue(VecP, i, pold[0]+DeltaT/DeltaX*(uold[2]*theta*rho0[1]/rho0[2]+uold[1]*(1.0-2.0*theta)-uold[0]*(1.0-theta)*rho0[1]/rho0[0]), INSERT_VALUES);
  }
  {
  // left boundary
    PetscScalar uold[2];
    PetscInt ix[2] = {0,1};
    VecGetValues(UIntermediate,2,ix,uold);
    PetscScalar pold[1];
    PetscInt ip[1] = {0};
    VecGetValues(VecP, 1, ip, pold);
    PetscScalar rho0[2];
    PetscInt irho0[2] = {0,1};
    VecGetValues(R0,2,irho0,rho0);
    VecSetValue(VecP, 0, pold[0]+DeltaT/DeltaX*(uold[1]*theta*rho0[0]/rho0[1]+uold[0]*(1.0-theta)), INSERT_VALUES);
  }
  {
  // right boundary
    PetscScalar uold[2];
    PetscInt ix[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
    VecGetValues(UIntermediate,2,ix,uold);
    PetscScalar pold[1];
    PetscInt ip[1] = {Number_Of_Elements-1};
    VecGetValues(VecP, 1, ip, pold);
    PetscScalar rho0[2];
    PetscInt irho0[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
    VecGetValues(R0,2,irho0,rho0);
    VecSetValue(VecP, Number_Of_Elements-1, pold[0]+DeltaT/DeltaX*(-uold[1]*theta-uold[0]*(1.0-theta)*rho0[1]/rho0[0]), INSERT_VALUES);
  }
  VecAssemblyBegin(VecP);
  VecAssemblyEnd(VecP);
  // Half Time Step
  // Interior Points
  for (i=1; i<Number_Of_Elements-1; i++)
  {
    PetscScalar pold[3];
    PetscInt ix[3] = {i-1,i,i+1};
    VecGetValues(VecP,3,ix,pold);
    PetscScalar uold[1];
    PetscInt iu[1] = {i};
    VecGetValues(UIntermediate, 1, iu, uold);
    VecSetValue(VecU, i, uold[0]+DeltaT/DeltaX/2.0*(pold[2]*(1.0-theta)-pold[1]*(1.0-2.0*theta)-pold[0]*theta), INSERT_VALUES);
  }
  {
  // left boundary
  PetscScalar pold[2];
  PetscInt ix[2] = {0, 1};
  VecGetValues(VecP, 2, ix, pold);
  PetscScalar uold[1];
  PetscInt iu[1] = {0};
  VecGetValues(UIntermediate,1,iu, uold);
  VecSetValue(VecU, 0, uold[0]+DeltaT/DeltaX/2.0*(pold[1]*(1.0-theta)-pold[0]*(1.0-theta)), INSERT_VALUES);
  }
  {
  // right boundary
  PetscScalar pold[2];
  PetscInt ix[2] = {Number_Of_Elements-2, Number_Of_Elements-1};
  VecGetValues(VecP, 2, ix, pold);
  PetscScalar uold[1];
  PetscInt iu[1] = {Number_Of_Elements-1};
  VecGetValues(UIntermediate,1,iu, uold);
  VecSetValue(VecU, Number_Of_Elements-1, uold[0]+DeltaT/DeltaX/2.0*(pold[1]*theta-pold[0]*theta), INSERT_VALUES);
  }
  VecAssemblyBegin(VecU);
  VecAssemblyEnd(VecU);

    PetscViewer viewer2;
    sprintf(szFileName, "solution%d.txt", t);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, szFileName, &viewer2);

      if (Method==4)
      {
        Energy=Hamiltonian_FV_SV(VecU, VecP, DeltaX, R0);
      }

    fprintf(f, "%1.16e\n", Energy-Initial_Energy);


  for (i=0;i<Number_Of_Elements;i++)
  {
    PetscScalar u[1],p[1];
    PetscInt ix[1] = {i};
    VecGetValues(VecU, 1, ix, u);
    VecGetValues(VecP, 1, ix, p);
    VecSetValue(Sol, i, u[0], INSERT_VALUES);
    VecSetValue(Sol, Number_Of_Elements+i, p[0], INSERT_VALUES);
  }
  VecAssemblyBegin(Sol);
  VecAssemblyEnd(Sol);

    VecView(Sol, viewer2);
    PetscViewerDestroy(&viewer2);

    //PetscPrintf(PETSC_COMM_SELF,"Delta t %6.4e\n",(double)DeltaT );
    //PetscPrintf(PETSC_COMM_SELF,"Delta x %6.4e\n",(double)DeltaX );
    //VecView(Initial_Condition, viewer);
    //  VecView(UIntermediate, viewer);
  //VecView(Sol, viewer);
  }
  fclose(f);
  PetscPrintf(PETSC_COMM_SELF,"Energy %6.4e\n",(double)Energy);
  ///PetscPrintf(PETSC_COMM_SELF,"Sigma %6.4e\n",(double)sigma );
  ///PetscPrintf(PETSC_COMM_SELF,"Time %6.4e\n",(double)T );
  PetscPrintf(PETSC_COMM_SELF,"Difference in Energy %6.4e\n",(double)Energy-Initial_Energy);
  VecDestroy(&UIntermediate);


  }

  // Continuous rho_0 exact integration


  //VecView(Initial_Condition, viewer);
  //VecView(Exact_Solution, viewer);
  //VecView(Sol, viewer);

  //VecAXPY(Sol,-1.0,Initial_Condition);
  VecAXPY(Sol,-1.0,Exact_Solution);
  VecNorm(Sol,NORM_2,&norm);
  norm *= sqrt(DeltaX);
  PetscPrintf(PETSC_COMM_WORLD,"L2-Norm of error %1.3e\n",(double)norm);

  VecDestroy(&Sol);
  VecDestroy(&RHS);
  VecDestroy(&Initial_Condition);
  VecDestroy(&Exact_Solution);
  VecDestroy(&R0);
  VecDestroy(&One_Over_rho0);
  VecDestroy(&drho_0dz);


  PetscViewerDestroy(&viewer);
  PetscFinalize();
  return 1;
}
