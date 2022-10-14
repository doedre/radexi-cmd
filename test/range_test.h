#pragma once

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>

#include <cmocka.h>

#include "radexi/input/range.h"

void
valid_range_test(void** state)
{
	(void)state;

	rxi_range_t range = rxi_range(0, 10, 10);
	assert_double_equal(range.begin, 0., 1e-5);
	assert_double_equal(range.end, 10., 1e-5);
	assert_double_equal(range.current, range.begin, 1e-5);
	assert_uint_equal(range.n, 10);

	RXI_STAT r = RXI_OK;
	for (size_t i = 0; (r = rxi_range_get(&range, i)) == RXI_OK;
			r = rxi_range_get(&range, ++i)) {
		assert_double_equal(rxi_range_get_cur(&range), (double)i, 1e-5);
	}

	assert_true(r == RXI_RANGE_END);
}

void
invalid_range_test(void** state)
{
	(void)state;

	rxi_range_t range = rxi_range(0, 10, 0);
	assert_double_equal(range.begin, 0., 1e-5);
	assert_double_equal(range.end, 10., 1e-5);
	assert_double_equal(range.current, range.begin, 1e-5);
	assert_uint_equal(range.n, 0);

	RXI_STAT r = RXI_OK;
	for (size_t i = 0; (r = rxi_range_get(&range, i)) == RXI_OK;
			r = rxi_range_get(&range, ++i)) {
		fail_msg("Range is invalid and can't be used in this loop");
	}

	assert_true(r == RXI_RANGE_INVALID);
}

