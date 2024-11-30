#include <lib/matrix2d.h>

void** matrix2d_new(int el_size, int i_dim, int j_dim)
{
    void* data_pool = malloc(el_size * i_dim * j_dim);
    void** matrix = malloc(sizeof(void*) * i_dim);
    for (int i = 0; i < i_dim; ++i)
    {
        matrix[i] = (data_pool + i * j_dim * el_size);
    }
    return matrix;
}

void** matrix2d_new_val(int el_size, int i_dim, int j_dim, void* val)
{
    void** matrix = matrix2d_new(el_size, i_dim, j_dim);
    for (int i = 0; i < i_dim; ++i)
    {
        for (int j = 0; j < j_dim; ++j)
        {
            memcpy(matrix[i] + j * el_size, val, el_size);
        }
    }
    return matrix;
}

void matrix2d_free(void** matrix)
{
    free(matrix[0]);
    free(matrix);
}
