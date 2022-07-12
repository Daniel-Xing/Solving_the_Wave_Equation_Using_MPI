//
// Created by 邢政 on 2022/2/16.
//

#ifndef TESTS_MPI_COURSEWORK_H
#define TESTS_MPI_COURSEWORK_H

#include <iostream>
#include <mpi.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <string>

namespace parallel
{

    // location info stored the important location info of the processor
    // All the data points are distributed into a square space, divided into many
    // small squares and assigned to different processors. For example, the size of
    // data points is 10 x 10, and the number of processors is 4. So, every processor
    // is responsible for 10 x 10 / 4 = 5 x 5 data points.
    // For example,
    // y = [0, 10)
    // ^ ___________________
    // | |5x5     |5x5     |
    // | |p2 1,0  |p3 1,1  |
    // | ------------
    // | |5x5     |5x5     |
    // | |p0 0,0  |p1 0,1  |
    // | ___________________
    // | -- -- -- -- -- > x = [0, 10)
    // the location info of p0 is:
    // id = 0 (processor id) p = 4 (total processors)
    // _local_row = 0, _local_col = 0, represent the two dimension idnex of chunks p0 is responsible for.
    // _chunk_col_size = 5, _chunk_row_size = 5
    // _col_low_index = 0, _col_high_index = 5 , indicate the low and high col indexs of chunk, the range is [_col_low_index, _col_high_index).
    // _row_low_index = 0, _row_high_index = 5 , indicate the low and high row indexs of chunk, the range is [_row_low_index, _row_high_index).
    //
    // Ghost chunk is designed to be bigger than chunck, which add one row on the top and bottom and add one column on the left and right of chunk.
    // So that when updating the boundary data, there is no need to communication to the other processor and just read the data from local memory.
    //
    // ghost_chunl_col_size = 6, ghost_chuck_row_size = 6, For those chunks not located in the bounday, it is should be 7, but cause p0 is in the corner,
    // so, no need to add on the left and bottom.
    // ghost_col_low_index = 0, ghost_col_low_index = 6, range is [0, 6).
    // ghost_row_low_index = 0, ghost_row_high_index = 6, range is [0, 6).
    // total_size = (ghost_col_high_index - ghost_col_low_index) * (ghost_row_high_index - ghost_row_low_index) = 6 * 6 = 36
    struct _location_info
    {
        int id, p;

        int _local_row; // the index of chunk taken by local processor
        int _local_col;

        int _chunk_col_size;
        int _chunk_row_size;

        int _col_low_index;
        int _col_high_index;

        int _row_low_index;
        int _row_high_index;

        int ghost_chunk_col_size;
        int ghost_chunk_row_size;

        int ghost_col_low_index;
        int ghost_col_high_index;

        int ghost_row_low_index;
        int ghost_row_high_index;

        int total_size;
    };

    // boundary_data is designed to recieve the date form neighbour processors.
    // if there is no neighbour in some direction, the data related willl be nullptr.
    struct boundary_data
    {
        double *left_boundary;
        double *right_boundary;
        double *top_boundary;
        double *bottom_boundary;
    };

    // comm_data is designed to storage the single data point and will be send to other
    // processor in comm_data array.
    class comm_data
    {
    public:
        double data;

        static void buildMPIType();
        static MPI_Datatype MPI_type;
    };

    // calculate the wave equation with:
    // boundary condition:
    int parallel_wave_equation_basic(int argc, char **argv);
}

#endif // TESTS_MPI_COURSEWORK_H
