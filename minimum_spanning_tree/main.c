/*
 * MPI realization of Prim's algorithm
 *
 * (c) Fara-light
 */

#include <lib/matrix2d.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define INF 10000000

int get_real_index(int n, int rank, int size, int im_index)
{
    int vertices_to_work = n / size;
    int offset = vertices_to_work * rank;
    return offset + im_index;
}

int** MPI_minimumSpanningTree(int** adjacency_matrix, int n, int rank, int size)
{
    int** result;
    if (rank == 0)
    {
        int tmp = 0;
        result = (int**)matrix2d_new_val(sizeof(int), n, n, &tmp);
    }

    int vertices_to_work = n / size;
    if (rank + 1 == size)
    {
        vertices_to_work += n % size;
    }

    char* used = (char*)calloc(n, sizeof(char));

    int* thread_cheapest_connection_cost = (int*)malloc(sizeof(int) * vertices_to_work); // C_i
    for (size_t i = 0; i < vertices_to_work; ++i)
    {
        thread_cheapest_connection_cost[i] = INF;
    }

    int* thread_cheapest_connection_edge = (int*)malloc(sizeof(int) * vertices_to_work); // E_i
    for (size_t i = 0; i < vertices_to_work; ++i)
    {
        thread_cheapest_connection_edge[i] = -1;
    }

    for (int node = 0; node < n; ++node)
    {
        struct {
            int cost;
            int index;
        } thread_iteration_cheapest = { .cost = INF, .index = rank }, iteration_cheapest;

        int thread_iteration_cheapest_edge = -1;
        int thread_prev_cheapest_edge = -1;

        for (size_t i = 0; i < vertices_to_work; ++i)
        {
            int real_index = get_real_index(n, rank, size, i);
            if (!used[real_index] && thread_cheapest_connection_cost[i] < thread_iteration_cheapest.cost)
            {
                thread_iteration_cheapest.cost = thread_cheapest_connection_cost[i];
                thread_prev_cheapest_edge = thread_cheapest_connection_edge[i];
                thread_iteration_cheapest_edge = i;
            }
        }

        int thread_previous_edge_index = thread_prev_cheapest_edge == -1 ? 0 : thread_prev_cheapest_edge;
        thread_iteration_cheapest.index = thread_previous_edge_index * n + get_real_index(n, rank, size, thread_iteration_cheapest_edge);

        MPI_Reduce(
            &thread_iteration_cheapest,
            &iteration_cheapest,
            1,
            MPI_2INT,
            MPI_MINLOC,
            0,
            MPI_COMM_WORLD
        );

        int previous_index, cheapest_index;
        if (rank == 0)
        {
            if (node == 0)
            {
                previous_index = 0;
                cheapest_index = 0;
            } else
            {
                previous_index = iteration_cheapest.index / n;
                cheapest_index = iteration_cheapest.index % n;
            }
            if (previous_index != cheapest_index)
            {
                result[previous_index][cheapest_index] = adjacency_matrix[previous_index][cheapest_index];
            }
        }

        MPI_Bcast(
            &cheapest_index,
            1,
            MPI_INT,
            0,
            MPI_COMM_WORLD
        );

        used[cheapest_index] = 1;

        for (size_t i = 0; i < vertices_to_work; ++i)
        {
            int real_index = get_real_index(n, rank, size, i);
            if (adjacency_matrix[cheapest_index][real_index] < thread_cheapest_connection_cost[i])
            {
                thread_cheapest_connection_cost[i] = adjacency_matrix[cheapest_index][real_index];
                thread_cheapest_connection_edge[i] = cheapest_index;
            }
        }
    }


    MPI_Barrier(MPI_COMM_WORLD); 

    free(used);
    free(thread_cheapest_connection_cost);
    free(thread_cheapest_connection_edge);

    return result;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    FILE *in_fp = fopen("input_graph.txt", "r");
    if (in_fp == NULL)
    {
        printf("Error: can't open input_graph.txt file\n");
        exit(1);
    }

    int n;
    fscanf(in_fp, "%d", &n);

    int tmp_value = INF;
    int** adjacency_matrix = (int**)matrix2d_new_val(sizeof(int), n, n, &tmp_value);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            fscanf(in_fp, "%d", &tmp_value);
            adjacency_matrix[i][j] = tmp_value == 0 ? INF : tmp_value;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int** res = MPI_minimumSpanningTree(adjacency_matrix, n, rank, size);

    if (rank == 0)
    {
        FILE *out_fp = fopen("output_graph.txt", "w");
        if (out_fp == NULL)
        {
            printf("Error: can't open output_graph.txt file\n");
            exit(2);
        }
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                fprintf(out_fp, "%d\t", res[i][j]);
            }
            fprintf(out_fp, "\n");
        }
        fclose(out_fp);
        matrix2d_free((void**)res);
    }

    matrix2d_free((void**)adjacency_matrix);

    MPI_Finalize();
    return 0;
}
