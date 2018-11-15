#ifndef _STATISTIC_H
#define _STATISTIC_H

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <cmath>
#include <vector>
#include <algorithm>


#define M_EPS 1e-15



void Mean_SD(std::vector<double> &vec, double &ret_mean, double &ret_stdev);

std::vector<double> OutlierElimination_IQR(std::vector<double> vecInp, double coef = 1.5);

int gsl_polynomialfit_robust(size_t n, size_t deg,
                             std::vector<double> &dx, std::vector<double> &dy,
                             double ret_coef[], double ret_cov[] = NULL,
                             const gsl_multifit_robust_type *robust_type = gsl_multifit_robust_bisquare);

#endif
