
CC=mpicc
C_OPTS=-Iinclude

OUT_FILE=parallel_prim.out

CPU_COUNT?=4
RUN=mpirun
RUN_OPTS=-np $(CPU_COUNT)

LIB_FILES = \
	lib/matrix2d.h

SOURCE_FILES= \
	lib/matrix2d.c

.PHONY: build run

$(OUT_FILE): main.c $(SOURCE_FILES)
	$(CC) $(C_OPTS) $^ -o $(OUT_FILE)

run: $(OUT_FILE)
	$(RUN) $(RUN_OPTS) ./$<

clean:
	rm -f $(OUT_FILE) output_graph.txt
