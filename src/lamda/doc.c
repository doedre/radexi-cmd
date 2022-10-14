#include <stdio.h>
#include <string.h>

#include <radexi/status.h>
#include <radexi/lamda/doc.h>

RXI_STAT
rxi_lamda_doc_init(rxi_lamda_doc_t* doc)
{
	doc->nodes = NULL;
	doc->size = 0;
	return RXI_OK;
}

RXI_STAT
rxi_lamda_doc_append(rxi_lamda_doc_t* doc, rxi_lamda_node_t* node)
{
	rxi_lamda_node_t* new_nodes;
	new_nodes = calloc(doc->size + 1, sizeof(*new_nodes));
	if (!new_nodes) {
		perror("Error allocating memory for new LAMDA node");
		return RXI_ERR_MALLOC;
	}

	if (doc->size != 0) {
		memcpy(new_nodes, doc->nodes, doc->size * sizeof(*doc->nodes));
		free(doc->nodes);
	}

	new_nodes[doc->size] = *node;
	++doc->size;
	doc->nodes = new_nodes;

	return RXI_OK;
}

rxi_lamda_node_t*
rxi_lamda_doc_get_last_node(rxi_lamda_doc_t* doc)
{
	return &doc->nodes[doc->size - 1];
}

void
rxi_lamda_doc_print(const rxi_lamda_doc_t* const doc)
{
	for (size_t i = 0; i < doc->size; ++i) {
		size_t ns = rxi_lamda_node_size(&doc->nodes[i]);
		printf("Node type %d with size %lu\n", doc->nodes[i].type, ns);

		for (size_t j = 0; j < ns; ++j) {
			if (j < 3)
				printf("\tLine %lu: %s\n", j, doc->nodes[i].lines[j]);
			else if (j == 3)
				printf("...\n");
			else if (j == ns - 1)
				printf("\tLine %lu: %s\n", j, doc->nodes[i].lines[j]);
			else
				continue;
		}
	}
}

void
rxi_lamda_doc_free(rxi_lamda_doc_t* doc)
{
	free(doc->nodes);
}
