#include "tepper-garcia.h"
#include <stdio.h>
#include <math.h>
double tau_i(double N, double lam, double lam_i, double b, double f_i, double Gamma_i)
{
  //equation 1 of tepper-garcia
  double a  = voigt_a(lam_i, b, Gamma_i);
  double Ca = voigt_Ca(lam_i, b, f_i);
  double x  = voigt_x(lam, lam_i, b);
  double H  = voigt_H(x,&a);

  //printf("lam %e x %e\n",lam,x);

  // tau(lam) = C_i a N H(a,x)
  //printf("Ca %e N %e H %e Ca*N*H %e\n",Ca,N,H,Ca*N*H);
  return Ca*N*H;
}
double voigt_x(double lam, double lam_i, double b)
{
  //see section 2.1 of Tepper Garcia 2006
  //note b is in km/sec
  //lam and lam_i in Angstrom
  double c = 2.99792458e5; //km/sec
  double Delta_lambda_D = (b/c)*lam_i;
  return (lam-lam_i)/Delta_lambda_D;
}
double voigt_a(double lam_i, double b, double Gamma_i)
{
  //see section 2.1 of Tepper Garcia 2006
  //note b is in km/sec
  //lam and lam_i in Angstrom
  //Gamma_i is Einstein A?
  double c = 2.99792458e5; //km/sec
  double Delta_lambda_D = (b/c)*lam_i;
  return lam_i*(lam_i*1e-13*Gamma_i)/(4*M_PI*c*Delta_lambda_D);
}
double voigt_Ca(double lam_i, double b, double f_i)
{
  //see section 2.1 of Tepper Garcia 2006
  //and equation 10 of Inoue and Iwata 2008
  //note b is in km/sec
  //lam_i in Angstrom
  /*double c = 2.99792458e10; //cm/sec
  double nu_D = (b/(c*1e-5))*(c*1.0e8/lam_i);
  double m_e = 9.1095e-28;//grams
  double e_e = 4.8032e-10;//esu
  double A = sqrt(M_PI)*e_e*e_e*f_i;
  double B = m_e*c*c*nu_D;*/
  //printf("A %e B %e A/B %e\n",A,B,A/B);
  double li = lam_i*1.0e-8; //to cm
  double c = 2.99792458e10; //cm/sec
  double m_e = 9.1095e-28;//grams
  double e_e = 4.8032e-10;//esu
  double A = sqrt(M_PI)*e_e*e_e*f_i;
  double Delta_lambda_D = (b*1.0e5/c)*li;
  double B = m_e*c*c*Delta_lambda_D;
  return A*li*li/B;
}
double voigt_H(double x, void *params)
{
  //see footnote 4 of Tepper Garcia 2006
  //lim x->0 H1(x) -> -2/sqrt(pi)
  //lim x->\infty H1(x) -> 0
  //
  double x2 = x*x;
  double H0 = exp(-1.*x2);
  double Q  = 1.5*x2;
  double *pa = (double *) params;
  double a = *pa;
  //double A = a/sqrt(M_PI)/x2 * (H0*H0*(4*x2*x2 + 7*x2 + 4 + Q) - Q -1.);
  double ans;
  double sh = sinh(x2);

  double Kx = 0.5/x2*((4*x2+3)*(x2+1)*H0 - (1./x2)*(2*x2 +3)*sh);

  ans = H0*(1-2*a/sqrt(M_PI)*Kx);
  printf("H0 %e x2 %e sh %e Kx %e\n",H0,x2,sh,Kx);
  //printf("H0 %e a %e x2 %e A %e H0-A %e\n",H0, a, x2,A, H0-A);
  //ans = H0 - A;
  //printf("H0 %e a %e x2 %e ans %e\n",H0, a, x2, ans);
  return ans;

}