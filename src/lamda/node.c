#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include <radexi/lamda/node.h>

RXI_STAT
rxi_lamda_node_init(rxi_lamda_node_t* node, const size_t size)
{
	if ((size + 1) > RXI_LAMDA_MAX_LINES_IN_NODE) {
		return RXI_ERR_WRONG_ARGUMENT;
	}

	errno = 0;
	node->capacity = size;
	node->lines = calloc(size + 1, sizeof(*node->lines));
	if (errno) {
		perror("Cannot allocate memory for LAMDA parser node");
		return RXI_ERR_MALLOC;
	}

	return RXI_OK;
}

bool
rxi_lamda_node_inited(rxi_lamda_node_t* node)
{
	return node->capacity;
}

size_t
rxi_lamda_node_size(rxi_lamda_node_t* node)
{
	for (size_t i = 0; i < node->capacity; ++i) {
		if (strlen(node->lines[i]) == 0)
			return i;
	}
	return node->capacity;
}

RXI_STAT
rxi_lamda_node_change_capacity(rxi_lamda_node_t* node, const size_t size)
{
	if (size > RXI_LAMDA_MAX_LINES_IN_NODE) {
		return RXI_ERR_WRONG_ARGUMENT;
	}

	errno = 0;
	char (*buf)[RXI_LAMDA_LINE_LEN];
	buf = calloc(size, sizeof(*buf));
	if (errno) {
		perror("Cannot allocate memory on LAMDA node resize");
		return RXI_ERR_MALLOC;
	}

	size_t cur_node_size = rxi_lamda_node_size(node);
	if (cur_node_size > size)
		cur_node_size = size;

	for (size_t i = 0; i < cur_node_size; ++i)
		strncpy(buf[i], node->lines[i], RXI_LAMDA_LINE_LEN);

	free(node->lines);

	node->capacity = size;
	node->lines = buf;

	return RXI_OK;
}

RXI_STAT
rxi_lamda_node_append(rxi_lamda_node_t* node, const char* str)
{
	const size_t size = rxi_lamda_node_size(node);
	if (size >= node->capacity)
		return RXI_ERR_NOT_ENOUGH_SPACE;

	snprintf(node->lines[size], RXI_LAMDA_LINE_LEN, "%s", str);
	return RXI_OK;
}

void
rxi_lamda_node_free(rxi_lamda_node_t* node)
{
	free(node->lines);
}
