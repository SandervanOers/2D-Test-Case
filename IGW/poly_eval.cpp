#include "poly_eval.hpp"

double Maths::Algebra::Polynomial::polyEval(double x, const double coef[], int N)
{
  const double *p = coef;
  double ans = *p++;
  int i = N;

  do
    ans = ans * x  +  *p++;
  while( --i );

  return( ans );
}

//! Evaluates a polynomial of degree N, with $C_N=1.0$.

double Maths::Algebra::Polynomial::polyEval1(double x, const double coef[], int N)
{
  const double *p = coef;
  double ans = x + *p++;
  int i = N-1;

  do
    ans = ans * x  + *p++;
  while( --i );

  return( ans );
}