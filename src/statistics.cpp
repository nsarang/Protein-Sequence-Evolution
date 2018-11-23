#include "statistics.h"


std::tuple<double, double> Mean_SD(std::vector<double> &vec) {
    double ret_mean, ret_stdev;
    /*
    double sum = 0, sq_sum = 0;
    int n = vec.size();

    for (auto &e : vec) {
        sum += e;
        sq_sum += e * e;
    }
    ret_mean = sum / n;
    double var = sq_sum / n - ret_mean * ret_mean;
    var = (var < M_EPS ? 0 : var);
    ret_stdev = std::sqrt(var);
    */

    ret_mean = gsl_stats_mean(vec.data(), 1, vec.size());
    ret_stdev = gsl_stats_sd_m(vec.data(), 1, vec.size(), ret_mean);
    return std::tuple(ret_mean, ret_stdev);
}

#include <iostream>
#include <cassert>
std::vector<double> OutlierElimination_IQR(std::vector<double> vecInp, double coef) {
    assert(vecInp.size());
    int Q1, Q2, Q3;
    Q1 = vecInp.size() / 4;
    Q2 = vecInp.size() / 2;
    Q3 = Q1 + Q2;
    std::cerr << "IQR " << Q1 << " " << Q2 << " " << Q3 << "\n";
    std::sort(vecInp.begin(), vecInp.end());
    double IQR = vecInp[Q3] - vecInp[Q1];

    auto low = std::lower_bound(vecInp.begin(), vecInp.end(), vecInp[Q1] - coef * IQR);
    auto up = std::upper_bound(vecInp.begin(), vecInp.end(), vecInp[Q3] + coef * IQR);
    std::cerr << IQR << " " << coef * IQR << " " << vecInp[Q1] << " " << vecInp[Q3] << "\n" << low -  vecInp.begin() << " " << up - vecInp.begin() << "   Passing\n";
    return std::vector(low, up);
}



int gsl_polynomialfit_robust(size_t n, size_t deg,
                             std::vector<double> &dx, std::vector<double> &dy,
                             double ret_coef[], double ret_cov[],
                             const gsl_multifit_robust_type *robust_type)
{
    gsl_matrix *X, *cov;
    gsl_vector *x, *y, *coef;

    X = gsl_matrix_alloc (n, deg);
    x = gsl_vector_alloc (n);
    y = gsl_vector_alloc (n);
    coef = gsl_vector_alloc (deg);
    cov = gsl_matrix_alloc (deg, deg);


    /* construct design matrix X for linear fit */
    for (size_t i = 0; i < n; ++i)
    {
        gsl_vector_set (x, i, dx[i]);
        gsl_vector_set (y, i, dy[i]);
        gsl_matrix_set (X, i, 0, 1.0);
        gsl_matrix_set (X, i, 1, dx[i]);
    }

    /* perform robust fit */
    gsl_multifit_robust_workspace *ws = gsl_multifit_robust_alloc(robust_type, X->size1, X->size2);
    int s = gsl_multifit_robust (X, y, coef, cov, ws);


#define C(i) (gsl_vector_get(coef,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

    ret_coef[0] = C(0);
    ret_coef[1] = C(1);
    if (ret_cov != NULL) {
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                ret_cov[i + j] = COV(i, j);
    }

    /* output data and model */
//   for (size_t i = 0; i < n; ++i)
//   {
//     double xi = gsl_vector_get(x, i);
//     double yi = gsl_vector_get(y, i);
//     gsl_vector_view v = gsl_matrix_row(X, i);
//     double y_ols, y_rob, y_err;

//     gsl_multifit_robust_est(&v.vector, coef, cov, &y_rob, &y_err);

//     printf("%g %g %g %g\n", xi, yi, y_rob, y_ols);
//   }

//   {
//     printf ("# best fit: Y = %g + %g X\n",
//             C(0), C(1));

//     printf ("# covariance matrix:\n");
//     printf ("# [ %+.5e, %+.5e\n",
//             COV(0, 0), COV(0, 1));
//     printf ("#   %+.5e, %+.5e\n",
//             COV(1, 0), COV(1, 1));
//   }

    gsl_matrix_free(X);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_vector_free(coef);
    gsl_matrix_free(cov);
    gsl_multifit_robust_free(ws);

    return s;
}
