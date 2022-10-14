#pragma once

#include <radexi/lamda/doc.h>

RXI_STAT rxi_lamda_parse(const char* path, rxi_lamda_doc_t* doc);

const char* rxi_lamda_parse_error();
