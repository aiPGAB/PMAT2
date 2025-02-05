CC := gcc
CFLAGS := -g -O3 -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-result
LIBS := -lm -lz -lpthread

SOURCES := PMAT.c log.c misc.c autoMito.c graphBuild.c hitseeds.c BFSseed.c \
           graphtools.c break_long_reads.c fastq2fa.c runassembly.c path2fa.c\
           get_subsample.c correct_sequences.c yak-count.c kthread.c \
		   graphPath.c orgAss.c
TARGET := PMAT

EXCLUDE_MAINS := -DHITSEEDS_MAIN -DBFSSEED_MAIN -DSUBSAMPLE_MAIN -DFQ2FA_MAIN -DRUNASSEMBLY_MAIN -DYAK_MAIN


.PHONY: all clean

all: info $(TARGET)
	@echo "Build complete: $(TARGET)"

$(TARGET): $(SOURCES)
	@echo "Compiling $(TARGET)..."
	@$(CC) $(CFLAGS) $(EXCLUDE_MAINS) -o $@ $^ $(LIBS)

clean:
	@rm -f $(TARGET)

info:
	@echo "Build Information:"
	@echo "  Compiler: $(CC)"
	@echo "  Libraries: $(LIBS)"
#	@echo "  Source files: $(SOURCES)"
