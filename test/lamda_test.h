#pragma once

#include <setjmp.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

#include <radexi/lamda/parse.h>

void
node_initialize_and_free_test(void** state)
{
	(void)state;

	RXI_STAT stat;
	const size_t size = 10;
	rxi_lamda_node_t node;

	stat = rxi_lamda_node_init(&node, size);
	assert_int_equal(stat, RXI_OK);
	assert_int_equal(node.type, RXI_LAMDA_NODE_UNKNOWN);
	for (size_t i = 0; i < size; ++i) {
		assert_non_null(node.lines[i]);
		sprintf(node.lines[i], "%lu", i);
		assert_int_equal(atoi(node.lines[i]), i);
	}

	assert_uint_equal(rxi_lamda_node_size(&node), size);

	rxi_lamda_node_free(&node);
}

void
node_change_capacity_test(void** state)
{
	(void)state;

	RXI_STAT stat;
	const size_t size = 10;
	rxi_lamda_node_t node;

	stat = rxi_lamda_node_init(&node, size);
	assert_int_equal(stat, RXI_OK);

	for (size_t i = 0; i < size; ++i)
		sprintf(node.lines[i], "%lu", i);

	assert_uint_equal(rxi_lamda_node_size(&node), size);

	// Increase capacity
	const size_t new_size = 50;
	stat = rxi_lamda_node_change_capacity(&node, new_size);
	assert_int_equal(stat, RXI_OK);

	assert_uint_equal(rxi_lamda_node_size(&node), size);
	sprintf(node.lines[10], "%d", 10);
	assert_uint_equal(rxi_lamda_node_size(&node), size + 1);

	for (size_t i = 0; i < size + 1; ++i)
		assert_uint_equal(atoi(node.lines[i]), i);

	// Decrease capacity
	const size_t new_size2 = 5;
	stat = rxi_lamda_node_change_capacity(&node, new_size2);
	assert_int_equal(stat, RXI_OK);

	assert_uint_equal(rxi_lamda_node_size(&node), new_size2);
	for (size_t i = 0; i < new_size2; ++i)
		assert_uint_equal(atoi(node.lines[i]), i);

	rxi_lamda_node_free(&node);
}

void
node_append_test(void** status)
{
	(void)status;

	RXI_STAT stat;
	const size_t size = 1;
	rxi_lamda_node_t node;

	stat = rxi_lamda_node_init(&node, size);
	assert_int_equal(stat, RXI_OK);

	const char* str = "append";
	stat = rxi_lamda_node_append(&node, str);
	assert_int_equal(stat, RXI_OK);
	assert_string_equal(node.lines[0], str);

	stat = rxi_lamda_node_append(&node, str);
	assert_int_equal(stat, RXI_ERR_NOT_ENOUGH_SPACE);
}

void
doc_append_test(void** status)
{
	(void)status;

	RXI_STAT stat;
	const size_t size = 2;
	const char* str = "append";
	const char* str2 = "append2";
	rxi_lamda_node_t node;

	stat = rxi_lamda_node_init(&node, size);
	assert_int_equal(stat, RXI_OK);
	stat = rxi_lamda_node_append(&node, str);
	assert_int_equal(stat, RXI_OK);

	rxi_lamda_doc_t doc;

	stat = rxi_lamda_doc_init(&doc);
	assert_int_equal(stat, RXI_OK);
	stat = rxi_lamda_doc_append(&doc, &node);
	assert_int_equal(stat, RXI_OK);

	rxi_lamda_node_t* last = rxi_lamda_doc_get_last_node(&doc);
	assert_non_null(last);
	assert_string_equal(last->lines[0], str);

	rxi_lamda_node_append(last, str2);
	assert_string_equal(last->lines[1], str2);
	stat = rxi_lamda_doc_append(&doc, last);
	assert_int_equal(stat, RXI_OK);

	rxi_lamda_node_t* last2 = rxi_lamda_doc_get_last_node(&doc);
	assert_non_null(last2);
	assert_string_equal(last2->lines[1], str2);

	rxi_lamda_doc_free(&doc);
	rxi_lamda_node_free(&node);
}

void
parse_test(void** status)
{
	(void)status;

	RXI_STAT stat;
	rxi_lamda_doc_t doc;

	stat = rxi_lamda_doc_init(&doc);
	assert_int_equal(stat, RXI_OK);
	stat = rxi_lamda_parse("../../data/hco+.dat", &doc);
	if (stat != RXI_OK)
		fail_msg("%s", rxi_lamda_parse_error());

//	rxi_lamda_doc_print(&doc);

	assert_int_equal(doc.nodes[0].type, RXI_LAMDA_NODE_MOLECULE_NAME);
	assert_string_equal(doc.nodes[0].lines[0], "HCO+");
	assert_int_equal(doc.nodes[1].type, RXI_LAMDA_NODE_MOLECULE_WEIGHT);
	assert_string_equal(doc.nodes[1].lines[0], "29.0");
	assert_int_equal(doc.nodes[2].type, RXI_LAMDA_NODE_NUMOF_ENERGY_LEVELS);
	assert_uint_equal(rxi_lamda_node_size(&doc.nodes[3]), atoi(doc.nodes[2].lines[0]));
	assert_int_equal(doc.nodes[4].type, RXI_LAMDA_NODE_NUMOF_RADIATIVE_TRANSITIONS);
	assert_uint_equal(rxi_lamda_node_size(&doc.nodes[5]), atoi(doc.nodes[4].lines[0]));
}
