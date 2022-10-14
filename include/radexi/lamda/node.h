#pragma once

#include <stdbool.h>
#include <stdlib.h>

#include <radexi/defines.h>
#include <radexi/status.h>

enum RXI_LAMDA_NODE_TYPE {
	RXI_LAMDA_NODE_UNKNOWN = 0,
	RXI_LAMDA_NODE_MOLECULE_NAME,
	RXI_LAMDA_NODE_MOLECULE_WEIGHT,
	RXI_LAMDA_NODE_NUMOF_ENERGY_LEVELS,
	RXI_LAMDA_NODE_NUMOF_RADIATIVE_TRANSITIONS,
	RXI_LAMDA_NODE_NUMOF_COLLISION_PARTNERS,
	RXI_LAMDA_NODE_COLLISION_PARTNER_NAME,
	RXI_LAMDA_NODE_COLLISION_PARTNER_NUMOF_TRANSIITONS,
	RXI_LAMDA_NODE_COLLISION_PARTNER_NUMOF_TEMPERATURES,
	RXI_LAMDA_NODE_COLLISION_PARTNER_TEMPERATURES,
	RXI_LAMDA_NODE_ENERGY_LEVELS,
	RXI_LAMDA_NODE_RADIATIVE_TRANSITIONS,
	RXI_LAMDA_NODE_COLLISION_PARTNER_TRANSITIONS,
	RXI_LAMDA_NODE_ADDITIONAL_INFO,
};

typedef struct {
	enum RXI_LAMDA_NODE_TYPE type;
	size_t capacity;
	char (*lines)[RXI_LAMDA_LINE_LEN];
}
rxi_lamda_node_t;

RXI_STAT rxi_lamda_node_init(rxi_lamda_node_t* node, const size_t size);
bool rxi_lamda_node_inited(rxi_lamda_node_t* node);
size_t rxi_lamda_node_size(rxi_lamda_node_t* node);
RXI_STAT rxi_lamda_node_change_capacity(rxi_lamda_node_t* node,
		const size_t size);
RXI_STAT rxi_lamda_node_append(rxi_lamda_node_t* node, const char* str);
void rxi_lamda_node_free(rxi_lamda_node_t* node);

