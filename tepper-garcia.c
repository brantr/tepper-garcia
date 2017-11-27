#include "tepper-garcia.h"
#include <math.h>
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
  double c = 2.99792458e10; //cm/sec
  double nu_D = (b/(c*1e-5))*(c*1.0e8/lam_i);
  double m_e = 9.1095e-28;//grams
  double e_e = 4.8032e-10;//esu
  double A = sqrt(M_PI)*e_e*e_e*f_i;
  double B = m_e*c*nu_D;
  return A/B;
}
double voigt_H(double x, void *params)
{
  //see footnote 4 of Tepper Garcia 2006
  double x2 = x*x;
  double H0 = exp(-1.*x2);
  double Q  = 1.5*x2;
  double a = (double) params;

  return H0 - a/sqrt(M_PI)/x2 * (H0*H0*(4*x2*x2 + 7*x2 + 4 + Q) - Q -1.);

}