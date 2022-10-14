#include "radexi/input/range.h"

rxi_range_t
rxi_range(const double beg, const double end, const size_t n)
{
	rxi_range_t r;
	r.begin = beg;
	r.end = end;
	r.n = n;
	r.current = beg;

	return r;
}

RXI_STAT
rxi_range_get(rxi_range_t* range, const size_t i)
{
	if ((range->n != 0) && (i <= range->n)) {
		double step = (range->end - range->begin) / range->n;
		range->current = range->begin + i * step;
		return RXI_OK;
	} else if (range->n == 0) {
		return RXI_RANGE_INVALID;
	} else {
		return RXI_RANGE_END;
	}
}

double
rxi_range_get_cur(const rxi_range_t* const range)
{
	return range->current;
}
