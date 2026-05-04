# Makefile - packing functions (DP, GTP, BP)
# -------------------------------------------
#
# Targets principales:
#   make              - compila biblioteca y tests
#   make test         - corre todos los tests
#   make clean        - limpia objetos y binarios
#   make debug        - compila con flags de depuracion (-g -O0 -fsanitize)
#   make run-arith    - corre solo el test de aritmetica

CC      ?= cc
CFLAGS  ?= -std=c11 -Wall -Wextra -Wpedantic -Wshadow -Wconversion \
           -Wno-sign-conversion -O2
CPPFLAGS = -Iinclude
LDFLAGS  =

ifeq ($(MODE),debug)
    CFLAGS  += -g -O0 -fsanitize=address,undefined
    LDFLAGS += -fsanitize=address,undefined
endif

SRC_DIR   = src
TEST_DIR  = tests
BUILD_DIR = build

LIB_SRC   = $(wildcard $(SRC_DIR)/*.c)
LIB_OBJ   = $(LIB_SRC:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

TEST_SRC  = $(wildcard $(TEST_DIR)/*.c)
TEST_BIN  = $(TEST_SRC:$(TEST_DIR)/%.c=$(BUILD_DIR)/%)

# ---------- targets ----------

.PHONY: all test clean debug run-arith help

all: $(TEST_BIN)

# Construir cada test enlazando con todos los objetos de la libreria.
$(BUILD_DIR)/%: $(TEST_DIR)/%.c $(LIB_OBJ) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LIB_OBJ) -o $@ $(LDFLAGS)

# Compilar fuentes de la libreria.
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Correr todos los tests en secuencia.  Falla en el primer error.
test: $(TEST_BIN)
	@set -e; for t in $(TEST_BIN); do \
		echo "=== ejecutando $$t ==="; \
		$$t; \
	done
	@echo "\n*** todos los tests pasaron ***"

run-arith: $(BUILD_DIR)/test_arith
	@$(BUILD_DIR)/test_arith

debug:
	$(MAKE) MODE=debug

clean:
	rm -rf $(BUILD_DIR)

help:
	@echo "Targets disponibles:"
	@echo "  make            - compila biblioteca y tests"
	@echo "  make test       - corre todos los tests"
	@echo "  make run-arith  - corre solo test_arith"
	@echo "  make debug      - compila con sanitizers"
	@echo "  make clean      - elimina build/"
