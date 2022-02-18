/**
 * @file core/sobolev.h
 * @brief Defines functions for more precise LVG method calculations.
 */

#ifndef RXI_CORE_SOBOLEV_H
#define RXI_CORE_SOBOLEV_H

#define RXI_INF_VELOCITY 10000

extern inline double rxi_f_gauss (double freq, double fwhm);

extern inline double rxi_f_velocity_distr (double r, double rc);

extern inline double rxi_f_optical_depth_0 (double r, double rc, double freq_0,
                                            double lstat, double ustat);

extern inline double rxi_f_optical_depth_inf (double r, double rc,
    double freq_0, double lstat, double ustat);

#endif  // RXI_CORE_SOBOLEV_H
