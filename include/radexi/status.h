#pragma once

typedef enum {
	RXI_OK = 0,
	RXI_ERR,
	RXI_ERR_MALLOC,
	RXI_ERR_WRONG_ARGUMENT,
	RXI_ERR_NOT_ENOUGH_SPACE,
	RXI_ERR_FILE,
	RXI_RANGE_END = 100,
	RXI_RANGE_INVALID,
	RXI_LAMDA_BAD_NODE,
}
RXI_STAT;
