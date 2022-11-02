#include "radexi/lamda/node.h"
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include <radexi/lamda/parse.h>

char rxi_lamda_parse_error_str[RXI_LAMDA_LINE_LEN];

static void
normalise_lamda_comment(char* comment)
{
	if (!comment)
		return;

	if (comment[0] != '!')
		return;

	for (int i = 0; comment[i] != '\0'; ++i) {
		if (comment[i] == '\n')
			comment[i] = ' ';

		if (isupper(comment[i]))
			comment[i] = tolower(comment[i]);

		if (comment[i] == ' ') {
			memmove(comment + i, comment + i + 1,
					RXI_LAMDA_LINE_LEN - i - 1);
			--i;
		}
	}
}

static void
remove_newline_symbol(char* line)
{
	size_t last = strlen(line) - 1;
	if ((line[last]) == '\n') {
		memmove(line + last, line + last + 1,
				RXI_LAMDA_LINE_LEN - last - 1);
	}
}

static enum RXI_LAMDA_NODE_TYPE
define_node_type(const char* norm_comment)
{
	if (!strcmp(norm_comment, "!molecule"))
		return RXI_LAMDA_NODE_MOLECULE_NAME;
	else if (!strcmp(norm_comment, "!molecularweight"))
		return RXI_LAMDA_NODE_MOLECULE_WEIGHT;
	else if (!strcmp(norm_comment, "!numberofenergylevels"))
		return RXI_LAMDA_NODE_NUMOF_ENERGY_LEVELS;
	else if (!strncmp(norm_comment, "!level+energ", 11))
		return RXI_LAMDA_NODE_ENERGY_LEVELS;
	else if (!strcmp(norm_comment, "!numberofradiativetransitions"))
		return RXI_LAMDA_NODE_NUMOF_RADIATIVE_TRANSITIONS;
	else if (!strncmp(norm_comment, "!trans+up+low+einsteina", 22))
		return RXI_LAMDA_NODE_RADIATIVE_TRANSITIONS;
	else if (!strcmp(norm_comment, "!numberofcollpartners")
			|| !strcmp(norm_comment, "!numberofcollisionpartners")
			)
		return RXI_LAMDA_NODE_NUMOF_COLLISION_PARTNERS;
	else if (!strncmp(norm_comment, "!partner", 7)
			|| !strcmp(norm_comment, "!collisionpartner")
			|| !strcmp(norm_comment, "!collisionbetween")
			)
		return RXI_LAMDA_NODE_COLLISION_PARTNER_NAME;
	else if (!strcmp(norm_comment, "!numberofcolltrans")
			|| !strcmp(norm_comment, "!numberofcollisionaltransitions")
			)
		return RXI_LAMDA_NODE_COLLISION_PARTNER_NUMOF_TRANSIITONS;
	else if (!strcmp(norm_comment, "!numberofcolltemps")
			|| !strcmp(norm_comment, "!numberofcollisiontemperatures")
			|| !strcmp(norm_comment, "!numberofcollisionaltemperatures")
			)
		return RXI_LAMDA_NODE_COLLISION_PARTNER_NUMOF_TEMPERATURES;
	else if (!strcmp(norm_comment, "!colltemps")
			|| !strcmp(norm_comment, "!collisionaltemperatures")
			|| !strcmp(norm_comment, "!collisiontemperatures")
			)
		return RXI_LAMDA_NODE_COLLISION_PARTNER_TEMPERATURES;
	else if (!strncmp(norm_comment, "!trans+up+low+rate", 17)
			|| !strncmp(norm_comment, "!trans+up+low+coll", 17)
			)
		return RXI_LAMDA_NODE_COLLISION_PARTNER_TRANSITIONS;
	else
		return RXI_LAMDA_NODE_UNKNOWN;
}

static bool
is_one_line_node(enum RXI_LAMDA_NODE_TYPE nt)
{
	if (nt > RXI_LAMDA_NODE_COLLISION_PARTNER_TEMPERATURES)
		return false;
	else
		return true;
}

RXI_STAT
rxi_lamda_parse(const char* path, rxi_lamda_doc_t* doc)
{
	FILE* file = fopen(path, "r");
	if (!file) {
		perror("Error opening file");
		return RXI_ERR_FILE;
	}

	RXI_STAT status;
	int nline = 0;
	size_t node_sizes = 1;
	char line[RXI_LAMDA_LINE_LEN];
	char orig_line[RXI_LAMDA_LINE_LEN];
	while (fgets(line, sizeof(line), file) != NULL) {
		++nline;

		remove_newline_symbol(line);
		if (line[0] == '!') {
			strncpy(orig_line, line, RXI_LAMDA_LINE_LEN);
			normalise_lamda_comment(line);
			enum RXI_LAMDA_NODE_TYPE nt = define_node_type(line);
			if (nt == RXI_LAMDA_NODE_UNKNOWN) {
				snprintf(rxi_lamda_parse_error_str,
						RXI_LAMDA_LINE_LEN,
						"In file %s, line %d: unknown node type `%s`",
						path, nline, orig_line);
				status = RXI_LAMDA_BAD_NODE;
				goto ANY_ERROR;
			}
			rxi_lamda_node_t new_node;
			new_node.type = nt;
			new_node.capacity = 0;
			new_node.lines = NULL;

			status = rxi_lamda_doc_append(doc, &new_node);
			if (status == RXI_ERR_MALLOC) {
				snprintf(rxi_lamda_parse_error_str,
						RXI_LAMDA_LINE_LEN,
						"In file %s, line %d: Failed to allocate memory to insert new node `%s` in LAMDA document",
						path, nline, orig_line);
				goto ANY_ERROR;
			}

			continue;
		}

		rxi_lamda_node_t* node = rxi_lamda_doc_get_last_node(doc);
		if (!node)
			continue;

		if ((node->type == RXI_LAMDA_NODE_NUMOF_COLLISION_PARTNERS)
				|| (node->type == RXI_LAMDA_NODE_NUMOF_RADIATIVE_TRANSITIONS)
				|| (node->type == RXI_LAMDA_NODE_NUMOF_ENERGY_LEVELS)
				|| (node->type == RXI_LAMDA_NODE_COLLISION_PARTNER_NUMOF_TRANSIITONS)
				) {
			node_sizes = atoi(line);
		}

		if (is_one_line_node(node->type)) {
			if (!rxi_lamda_node_inited(node)) {
				status = rxi_lamda_node_init(node, 1);
				if (status != RXI_OK) {
					snprintf(rxi_lamda_parse_error_str,
							RXI_LAMDA_LINE_LEN,
							"In file %s, line %d: Failed to allocate memory for one line node `%s`",
							path, nline, orig_line);
					goto ANY_ERROR;
				}
			}

			rxi_lamda_node_append(node, line);
			continue;
		}

		if (!rxi_lamda_node_inited(node)) {
			status = rxi_lamda_node_init(node, node_sizes);
			if (status == RXI_ERR_MALLOC) {
				snprintf(rxi_lamda_parse_error_str,
						RXI_LAMDA_LINE_LEN,
						"In file %s, line %d: Failed to allocate memory for multiline node of size %lu",
						path, nline, node_sizes);
				goto ANY_ERROR;
			} else if (status == RXI_ERR_WRONG_ARGUMENT) {
				status = rxi_lamda_node_init(node, RXI_LAMDA_MAX_LINES_IN_NODE - 1);
				if (status != RXI_OK) {
					snprintf(rxi_lamda_parse_error_str,
							RXI_LAMDA_LINE_LEN,
							"In file %s, line %d: Failed to allocate memory for multiline node of size %d",
							path, nline, RXI_LAMDA_MAX_LINES_IN_NODE - 1);
					goto ANY_ERROR;
				}
			}
		}

		status = rxi_lamda_node_append(node, line);
		if (status == RXI_ERR_NOT_ENOUGH_SPACE) {
			status = rxi_lamda_node_change_capacity(node,
					++node_sizes);
			if (status == RXI_ERR_MALLOC) {
				snprintf(rxi_lamda_parse_error_str,
						RXI_LAMDA_LINE_LEN,
						"In file %s, line %d: Failed to allocate memory when increasing capacity of document's node to %lu",
						path, nline, node_sizes);
				goto ANY_ERROR;
			}
		}
	}

	if (ferror(file)) {
		perror("Error reading from file: ");
		fclose(file);
		return RXI_ERR_FILE;
	}

	fclose(file);

	return RXI_OK;

ANY_ERROR:
	if (file)
		fclose(file);

	if (doc->size != 0)
		rxi_lamda_doc_free(doc);

	return status;
}

const char*
rxi_lamda_parse_error()
{
	return rxi_lamda_parse_error_str;
}

