INCLUDE := -I./src -I./3rdparty
LDFLAGS := -lgsl -lm -lgslcblas
VPATH := 3rdparty/linenoise 3rdparty/minIni src src/core src/utils

CC := clang
CFLAGS := -std=gnu11 -Wall -Wextra --pedantic ${INCLUDE} -DNDEBUG

BUILD_DIR := bin
OBJ_DIR := .obj
SRC_DIR := src 3rdparty

SRC := \
	3rdparty/linenoise/linenoise.c \
	3rdparty/minIni/minIni.c \
	src/core/background.c \
	src/core/calculation.c \
	src/core/dialogue.c \
	src/core/output.c \
	src/utils/cli_tools.c \
	src/utils/csv.c \
	src/utils/database.c \
	src/utils/options.c \
	src/main.c \
	src/rxi_common.c

OBJ := ${SRC:.c=.o}

.PHONY: all options debug clean install uninstall

all: options radexi

options:
	@echo dwm build options:
	@echo "CFLAGS   = ${CFLAGS}"
	@echo "LDFLAGS  = ${LDFLAGS}"
	@echo "CC       = ${CC}"
	@echo

.c.o:
	@echo [CC] $<
	@${CC} -c ${CFLAGS} -o ${OBJ_DIR}/${notdir $@} $<

${OBJ_DIR}:
	mkdir -p $@

${BUILD_DIR}:
	mkdir -p $@

radexi: ${BUILD_DIR} ${OBJ_DIR} ${OBJ}
	@echo [LD] ${BUILD_DIR}/$@
	@${CC} ${CFLAGS} ${LDFLAGS} ${addprefix ${OBJ_DIR}/,${notdir ${OBJ}}} -o ${BUILD_DIR}/$@

debug: CFLAGS := $(filter-out -DNDEBUG,$(CFLAGS))
debug: radexi 

install: all
	cp -f bin/radexi /usr/local/bin/
	chmod 755 /usr/local/bin/radexi

uninstall:
	rm -f /usr/local/bin/radexi

clean:
	rm -rf ${OBJ_DIR}
	rm bin/radexi
