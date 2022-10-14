#pragma once

#include <stdlib.h>

#include <radexi/status.h>
#include <radexi/lamda/node.h>

typedef struct {
	rxi_lamda_node_t* nodes;
	size_t size;
}
rxi_lamda_doc_t;

RXI_STAT rxi_lamda_doc_init(rxi_lamda_doc_t* doc);

RXI_STAT rxi_lamda_doc_append(rxi_lamda_doc_t* doc,
		rxi_lamda_node_t* node);

rxi_lamda_node_t* rxi_lamda_doc_get_last_node(rxi_lamda_doc_t* doc);

void rxi_lamda_doc_print(const rxi_lamda_doc_t* const doc);

void rxi_lamda_doc_free(rxi_lamda_doc_t* doc);
	
