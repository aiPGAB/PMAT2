CC := gcc
CFLAGS := -g -O3 -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-result
LIBS := -lm -lz -lpthread

SOURCES := PMAT.c log.c misc.c autoMito.c graphBuild.c hitseeds.c DFSseed.c \
           graphtools.c break_long_reads.c fastq2fa.c runassembly.c path2fa.c\
           get_subsample.c correct_sequences.c yak-count.c kthread.c \
		   graphPath.c
TARGET := PMAT

EXCLUDE_MAINS := -DHITSEEDS_MAIN -DDFSSEED_MAIN -DSUBSAMPLE_MAIN -DFQ2FA_MAIN -DRUNASSEMBLY_MAIN -DYAK_MAIN

BLUE := \033[1;34m
GREEN := \033[1;32m
RED := \033[1;31m
YELLOW := \033[1;33m
RESET := \033[0m

.PHONY: all clean

all: info $(TARGET)
	@echo "$(GREEN)Build complete: $(TARGET)$(RESET)"

$(TARGET): $(SOURCES)
	@echo "$(BLUE)Compiling $(TARGET)...$(RESET)"
	@$(CC) $(CFLAGS) $(EXCLUDE_MAINS) -o $@ $^ $(LIBS)

clean:
	@rm -f $(TARGET)

info:
	@echo "$(BLUE)Build Information:$(RESET)"
	@echo "  Compiler: $(CC)"
	@echo "  Libraries: $(LIBS)"
#	@echo "  Source files: $(SOURCES)"
