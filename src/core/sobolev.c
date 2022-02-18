/**
 * @file core/sobolev.c
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "core/sobolev.h"

#include "rxi_common.h"

inline double
rxi_f_gauss (double freq, double fwhm)
{
  double sigma = fwhm * log (2) / 2;
  return gsl_ran_gaussian_pdf (freq, sigma);
}

inline double
rxi_f_velocity_distr (double r, double rc)
{
  return sqrt (1 - rc / r) * RXI_INF_VELOCITY;
}

inline double
rxi_f_optical_depth_0 (double r, double rc, double freq_0, double lstat,
                       double ustat)
{
  double velocity = rxi_f_velocity_distr (r, rc);
  return
          (M_PI * gsl_pow_2 (RXI_E) * (lstat - ustat) * r)
      / //------------------------------------------------
                  (RXI_EM * freq_0 * velocity);
}
