CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -Wpedantic -O2 -g
CFLAGS  += -I libsinter

SRC_DIR  = libsinter
TEST_DIR = $(SRC_DIR)/test
BUILD    = build

# Core source files (added as implemented)
CORE_SRC = \
	$(SRC_DIR)/core/transform/lapping.c \
	$(SRC_DIR)/core/transform/dct.c \
	$(SRC_DIR)/core/quantize/pvq.c \
	$(SRC_DIR)/core/quantize/qp_table.c \
	$(SRC_DIR)/core/entropy/rans.c \
	$(SRC_DIR)/core/bitstream/bitstream.c \
	$(SRC_DIR)/core/prediction/intra.c \
	$(SRC_DIR)/core/prediction/obmc.c \
	$(SRC_DIR)/core/grain/grain_synth.c \
	$(SRC_DIR)/perceptual/csf_weights.c \
	$(SRC_DIR)/perceptual/activity_mask.c \
	$(SRC_DIR)/perceptual/ptd.c \
	$(SRC_DIR)/perceptual/grain_analysis.c \
	$(SRC_DIR)/codec/encoder.c \
	$(SRC_DIR)/codec/decoder.c \
	$(SRC_DIR)/codec/rdo.c

CORE_OBJ = $(patsubst $(SRC_DIR)/%.c, $(BUILD)/%.o, $(CORE_SRC))

# Test executables
TESTS = \
	$(BUILD)/test_transform \
	$(BUILD)/test_pvq \
	$(BUILD)/test_entropy \
	$(BUILD)/test_roundtrip

.PHONY: all clean test

all: $(BUILD)/libsinter.a

$(BUILD)/libsinter.a: $(CORE_OBJ)
	@mkdir -p $(dir $@)
	ar rcs $@ $^

$(BUILD)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

# Test targets — each links against the library
$(BUILD)/test_%: $(TEST_DIR)/test_%.c $(BUILD)/libsinter.a
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $< -L$(BUILD) -lsinter -lm -o $@

test: $(TESTS)
	@for t in $(TESTS); do \
		echo "=== Running $$t ==="; \
		$$t || exit 1; \
	done
	@echo "All tests passed."

clean:
	rm -rf $(BUILD)
