#pragma once

#include <stdint.h>
#include <stdlib.h>

typedef struct {
	int32_t* buf;
	size_t size;
}
vec_i32_t;

typedef struct {
	uint32_t* buf;
	size_t size;
}
vec_u32_t;

typedef struct {
	float* buf;
	size_t size;
}
vec_f32_t;

typedef struct {
	double* buf;
	size_t size;
}
vec_f64_t;

