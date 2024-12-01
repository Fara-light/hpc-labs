
CC=mpicc
C_OPTS=-Iinclude

BIN_FILE=parallel_prim.out

CPU_COUNT?=4
RUN=mpirun
RUN_OPTS=-np $(CPU_COUNT)

LIB_FILES = \
	lib/matrix2d.h

SOURCE_FILES= \
	lib/matrix2d.c

.PHONY: build run

$(BIN_FILE): main.c $(SOURCE_FILES)
	$(CC) $(C_OPTS) $^ -o $(BIN_FILE)

IN_FILE?=input_graph.txt
OUT_FILE?=output_graph.txt

run: $(BIN_FILE)
	$(RUN) $(RUN_OPTS) ./$< $(IN_FILE) $(OUT_FILE)

clean:
	rm -f $(BIN_FILE) output_graph.txt
