#ifndef MATRIX2D_H
#define MATRIX2D_H

#include <string.h>
#include <stdlib.h>

void** matrix2d_new(int el_size, int i_dim, int j_dim);
void** matrix2d_new_val(int el_size, int i_dim, int j_dim, void* val);
void matrix2d_free(void** matrix);

#endif // MATRIX2D_H
