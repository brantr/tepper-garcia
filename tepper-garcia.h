#ifndef TEPPER_GARCIA_H
#define TEPPER_GARCIA_H
double voigt_x(double lam, double lam_i, double b);
double voigt_a(double lam_i, double b, double Gamma_i);
double voigt_Ca(double lam_i, double b, double f_i);
double voigt_H(double x, void *params);
#endif //TEPPER_GARCIA_H
