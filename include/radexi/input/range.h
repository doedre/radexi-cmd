#pragma once

#include <stdlib.h>

#include <radexi/status.h>

typedef struct {
	double current;
	double begin;
	double end;
	size_t n;
}
rxi_range_t;

rxi_range_t rxi_range(const double beg, const double end, const size_t n);
RXI_STAT rxi_range_get(rxi_range_t* range, const size_t i);
double rxi_range_get_cur(const rxi_range_t* const range);

