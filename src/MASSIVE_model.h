#ifndef MASSIVE_MODEL_H
#define MASSIVE__MODEL_H

#include <RcppArmadillo.h>
#include <boost/math/distributions/normal.hpp>
#include <vector>

#define LOG_PDF_NORM_VEC(p, sd) (- 0.5 * log(2 * M_PI) - log(sd) - (p % p) / (2 * sd % sd))
#define LOG_PDF_NORM(p, sd) (- 0.5 * log(2 * M_PI) - log(sd) - (p * p) / (2 * sd * sd))

struct Parameters {
  arma::vec sgamma;
  arma::vec salpha;
  double sbeta;
  double skappa_X;
  double skappa_Y;
  double sigma_X;
  double sigma_Y;
};

#endif // MASSIVE_MODEL_H