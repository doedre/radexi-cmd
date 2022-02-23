CFLAGS := -Wall -Wextra --pedantic -DNDEBUG
CC := clang
INCLUDE := -I./src -I./3rdparty
LIBS := -lgsl -lm -lgslcblas

TARGET_NAME := radexi
BUILD_DIR := ./bin
OBJ_DIR := ./.obj
SRC_DIR := ./src ./3rdparty

SOURCES := \
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

OBJECTS := $(patsubst %.c,$(OBJ_DIR)/%.o,$(SOURCES))

TARGET := $(BUILD_DIR)/$(TARGET_NAME)
	
$(OBJ_DIR)/%.o: %.c
	@echo [CC] $<
	@mkdir -p $(dir $@)
	@$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBS) $(OBJECTS) -o $@

debug: CFLAGS := $(filter-out -DNDEBUG,$(CFLAGS))
debug: $(TARGET)

clean:
	rm -rf $(OBJ_DIR)
	rm $(TARGET)