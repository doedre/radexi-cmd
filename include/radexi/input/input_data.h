#pragma once

#include <radexi/defines.h>
#include <radexi/input/range.h>

enum GEOMETRY {
	SPHERE = 0,
	SLAB,
	LVG
};

enum COLLISION_PARTNER {
	H2 = 0,
	PARA_H2,
	ORTHO_H2,
	ELECTRONS,
	H_I,
	HELIUM,
	H_II,
};

struct rxi_collision_partner {
	enum COLLISION_PARTNER name;
	rxi_range_t density;
};

struct rxi_input_data {
	char molecule_name[RXI_STRING_LEN];
	rxi_range_t freq;
	rxi_range_t temp_kin;
	rxi_range_t temp_bg;
	rxi_range_t col_dens;
	rxi_range_t line_width;
	enum GEOMETRY geom;
	struct rxi_collision_partner coll_part[RXI_MAX_COLLISION_PARTNERS];
};

