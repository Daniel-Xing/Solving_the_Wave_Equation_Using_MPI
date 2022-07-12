#ifndef MPI_COURSE_WORK_CPP
#define MPI_COURSE_WORK_CPP

#define collect_time

#include "Mpi_CourseWork.h"

namespace parallel
{

    /*** run-configuration defined ***/
    int imax = 100;
    int jmax = 100;      // the size of grid points
    double t_max = 30.0; // upper bound of t
    double t, t_out = 0.0, dt_out = 1, dt;
    double y_max = 10.0, x_max = 10.0, dx, dy;
    double c = 1;

    int id; // processor id
    int p;  // the number of processors

    int rows = 0; // how many data chuncks in direction of rows.
    int cols = 0; // how many data chuncks in direction of cols.

    _location_info locationInfo;        // global variable: location info
    boundary_data boundary;             // define boundary datas for storaging the boundary revieved.
    double *new_grid, *grid, *old_grid; // grid datas in the iteration of n -1, n , n +1

    MPI_Datatype comm_data::MPI_type;
    /*** run-configuration defined ***/

    /** Implemment the member function of comm_data -
     * buildMPIType will be used to create a mpi struct
     */
    void comm_data::buildMPIType()
    {
        int block_lengths[1];
        MPI_Aint displacements[1];
        MPI_Aint addresses[1], add_start;
        MPI_Datatype typelist[1];

        comm_data temp;

        typelist[0] = MPI_DOUBLE;
        block_lengths[0] = 1;
        MPI_Get_address(&temp.data, &addresses[0]);

        MPI_Get_address(&temp, &add_start);
        for (int i = 0; i < 1; i++)
            displacements[i] = addresses[i] - add_start;

        MPI_Type_create_struct(1, block_lengths, displacements, typelist, &MPI_type);
        MPI_Type_commit(&MPI_type);
    }

    /**
     * @brief find_dimensions decompose p into a product of two numbers that are as
     * close as possible to each other. P = m x n
     *
     * @param processor : the number of processors
     * @param m : first number
     * @param n : second number
     */
    void find_dimensions(const int &processor, int &m, int &n)
    {
        // init the values
        m = 1;
        n = processor;

        // try to find the divisible number from large to small
        // in this case, if found, the pair of numbers will the closest to each other.
        for (int i = round(sqrt(processor)); i > 1; i--)
        {
            if (processor % i == 0)
            {
                m = i;
                n = processor / m;
                break;
            }
        }
    }

    /**
     * @brief init_location_info init the global variable locationInfo.
     * At the beginning, init the processor id and number of processor.
     * the determinate the position of the row and col of chunk, then use
     * this information to call the length of row and col of chunks, which
     * is _chunk_col_size and _chunk_row_size. And the index of chunk will
     * be determinated.
     * all the index is Left open and right closed interval -> [low_index, high_index)
     *
     */
    void init_location_info()
    {
        locationInfo.id = id;
        locationInfo.p = p;

        // calculate the processor location
        locationInfo._local_row = id / cols;
        locationInfo._local_col = id % cols;

        locationInfo._chunk_col_size = locationInfo._local_col != cols - 1 ? imax / cols : imax / cols + imax % cols;
        locationInfo._chunk_row_size = locationInfo._local_row != rows - 1 ? jmax / rows : jmax / rows + jmax % rows;

        locationInfo._col_low_index = (imax / cols) * locationInfo._local_col;
        locationInfo._col_high_index = locationInfo._local_col != cols - 1 ? (imax / cols) * (locationInfo._local_col + 1) : imax;

        locationInfo._row_low_index = (jmax / rows) * locationInfo._local_row;
        locationInfo._row_high_index = locationInfo._local_row != rows - 1 ? (jmax / rows) * (locationInfo._local_row + 1) : jmax;

        // calculate the chunk low_index and high_index
        locationInfo.ghost_col_low_index = locationInfo._local_col != 0 ? locationInfo._local_col * imax / cols - 1 : 0;
        locationInfo.ghost_col_high_index = locationInfo._local_col != cols - 1 ? (locationInfo._local_col + 1) * imax / cols + 1 : imax;

        locationInfo.ghost_row_low_index = locationInfo._local_row != 0 ? locationInfo._local_row * jmax / rows - 1 : 0;
        locationInfo.ghost_row_high_index = locationInfo._local_row != rows - 1 ? (locationInfo._local_row + 1) * jmax / rows + 1 : jmax;

        locationInfo.ghost_chunk_col_size = locationInfo.ghost_col_high_index - locationInfo.ghost_col_low_index;
        locationInfo.ghost_chunk_row_size = locationInfo.ghost_row_high_index - locationInfo.ghost_row_low_index;

        locationInfo.total_size = (locationInfo.ghost_row_high_index - locationInfo.ghost_row_low_index) *
                                  (locationInfo.ghost_col_high_index - locationInfo.ghost_col_low_index);
    }

    /**
     * @brief print_location_info print the value of lcoationInfo for checking,
     *
     */
    void print_location_info()
    {
        std::cout << "processor " << id << " location info:" << std::endl;
        std::cout << " local_row: " << locationInfo._local_row << "   local_col: " << locationInfo._local_col << std::endl;
        std::cout << " _chunk_row_size: " << locationInfo._chunk_row_size << " chunk_col_size: " << locationInfo._chunk_col_size << std::endl;
        std::cout << " _col_low_index: " << locationInfo._col_low_index << " _col_high_index: " << locationInfo._col_high_index << std::endl;
        std::cout << " _row_low_index: " << locationInfo._row_low_index << " _row_high_index: " << locationInfo._row_high_index << std::endl;
        std::cout << " ghost_chunk_row_size: " << locationInfo.ghost_chunk_row_size << " ghost_chunk_col_size: " << locationInfo.ghost_chunk_col_size << std::endl;
        std::cout << " ghost_col_low_index: " << locationInfo.ghost_col_low_index << " ghost_col_high_index: " << locationInfo.ghost_col_high_index << std::endl;
        std::cout << " ghost_row_low_index: " << locationInfo.ghost_row_low_index << " ghost_row_high_index: " << locationInfo.ghost_row_high_index << std::endl;
        std::cout << " total_size: " << locationInfo.total_size << std::endl;
        std::cout << std::endl;
        std::cout.flush();
    }

    /**
     * @brief index_mapper transfer the 2D index (i, j) into actual 1D index of grid, which is 1D array.
     *
     * @param i row
     * @param j col
     * @return int value, if return -1, means not found
     */
    int index_mapper(int i, int j)
    {
        if (i < locationInfo.ghost_row_low_index || i >= locationInfo.ghost_row_high_index ||
            j < locationInfo.ghost_col_low_index || j >= locationInfo.ghost_col_high_index)
            return -1;

        return (i - locationInfo.ghost_row_low_index) * locationInfo.ghost_chunk_col_size + (j - locationInfo.ghost_col_low_index);
    }

    /**
     * @brief id_from_index return the processor id using id_row and id_column
     *
     * @param id_row : rows
     * @param id_column : cols
     * @return int : processor id
     */
    int id_from_index(int id_row, int id_column)
    {
        if (id_row >= rows || id_row < 0)
            return -1;
        if (id_column >= cols || id_column < 0)
            return -1;

        return id_row * cols + id_column;
    }

    /**
     * @brief grid_to_file wirte the data of grid intp the file.
     *
     * @param out the number of iteration.
     */
    void grid_to_file(int out)
    {
        // Write the output for a single time step to file
        std::stringstream fname;
        std::fstream f1;
        // fname << "../out/output" <<  "_" << out << "_" << id << ".dat";
        fname << "/rds/general/user/dx121/home/mpi-coursework-acse-dx121"
              << "/out/output"
              << "_" << out << "_" << id << ".dat";
        f1.open(fname.str().c_str(), std::ios_base::out);
        for (int i = locationInfo._row_low_index; i < locationInfo._row_high_index; i++)
        {
            for (int j = locationInfo._col_low_index; j < locationInfo._col_high_index; j++)
                f1 << grid[index_mapper(i, j)] << "\t";
            f1 << std::endl;
        }
        f1.close();
    }

    /**
     * @brief grid_to_file_with_flag write the data of new grid into file
     * using during processing to check if all thing work corcertly
     *
     * @param out the number of iteration
     * @param flag a string to identify
     */
    void grid_to_file_with_flag(int out, std::string flag)
    {
        // Write the output for a single time step to file
        std::stringstream fname;
        std::fstream f1;
        fname << "../out/output" << flag << "_" << out << "_" << id << ".dat";
        f1.open(fname.str().c_str(), std::ios_base::out);
        for (int i = locationInfo._row_low_index; i < locationInfo._row_high_index; i++)
        {
            for (int j = locationInfo._col_low_index; j < locationInfo._col_high_index; j++)
                f1 << new_grid[index_mapper(i, j)] << "\t";
            f1 << std::endl;
        }
        f1.close();
    }

    /**
     * @brief half_sinusoidal_intitial_disturbance init the disturbance at the begnining
     *
     */
    void half_sinusoidal_intitial_disturbance(double r_splash, double x_splash, double y_splash)
    {
        // double r_splash = 1.0;
        // double x_splash = 3.0;
        // double y_splash = 3.0;
        for (int i = locationInfo.ghost_row_low_index; i < locationInfo.ghost_row_high_index; i++)
        {
            if (i == 0 || i == jmax - 1)
                continue;
            for (int j = locationInfo.ghost_col_low_index; j < locationInfo.ghost_col_high_index; j++)
            {
                if (j == 0 || j == imax - 1)
                    continue;

                double x = dx * i;
                double y = dy * j;
                double dist = sqrt(pow(x - x_splash, 2.0) + pow(y - y_splash, 2.0));

                int id = index_mapper(i, j);
                if (dist < r_splash)
                {
                    double h = 5.0 * (cos(dist / r_splash * M_PI) + 1.0);
                    grid[id] += h;
                    old_grid[id] += h;
                }
            }
        }
    }

    /**
     * @brief fake_initial_disturbance init the disturbance using the const value id + 10.
     * using to check if every thing works, easy to calculate the boundary values.
     *
     */
    void fake_initial_disturbance()
    {
        for (int i = locationInfo.ghost_row_low_index; i < locationInfo.ghost_row_high_index; i++)
        {
            if (i == 0 || i == jmax - 1)
                continue;
            for (int j = locationInfo.ghost_col_low_index; j < locationInfo.ghost_col_high_index; j++)
            {
                if (j == 0 || j == imax - 1)
                    continue;

                grid[index_mapper(i, j)] = id + 10;
                old_grid[index_mapper(i, j)] = id + 10;
            }
        }
    }

    /**
     * @brief do_calculation calculate the new value using the data from grid and update new_grid
     * no need to worry about the out of range, the checker has been done when using this.
     *
     * @param i
     * @param j
     * @return double
     */
    double do_calculation(int i, int j)
    {
        // if(id == 0 && i == 1 && j == 3) {
        //             std::cout << "upper: " << grid[index_mapper(i + 1, j)] << " ";
        //             std::cout << "down: " << grid[index_mapper(i - 1 , j)] << " ";
        //             std::cout << "left: " << grid[index_mapper(i , j -1)] << " ";
        //             std::cout << "right: " << grid[index_mapper(i , j + 1)] << " ";

        //             std::cout << "self: " << grid[index_mapper(i , j)] << " ";
        //             std::cout << "dx: " << dx << " ";
        //             std::cout << "dy: " << dy << " ";
        //             std::cout<< std::endl;
        //             std::cout.flush();
        //         }
        double res = pow(dt * c, 2.0) * ((grid[index_mapper(i + 1, j)] - 2.0 * grid[index_mapper(i, j)] + grid[index_mapper(i - 1, j)]) / pow(dx, 2.0) + (grid[index_mapper(i, j + 1)] - 2.0 * grid[index_mapper(i, j)] + grid[index_mapper(i, j - 1)]) / pow(dy, 2.0)) + 2.0 * grid[index_mapper(i, j)] - old_grid[index_mapper(i, j)];
        return res;
    }

    /**
     * @brief do_iteration_no_boundary calcualte the new wave value expect for the edges
     * used during communication to improve the runtime
     *
     */
    void do_iteration_no_boundary()
    {
        // TODO: take care of case of out of range.
        for (int i = locationInfo._row_low_index + 1; i < locationInfo._row_high_index - 1; i++)
        {
            for (int j = locationInfo._col_low_index + 1; j < locationInfo._col_high_index - 1; j++)
            {
                // if(id == 0 && i == 2 && j == 1) {
                //     std::cout << "upper: " << index_mapper(i + 1, j) << " ";
                //     std::cout << "down: " << index_mapper(i - 1 , j) << " ";
                //     std::cout << "left: " << index_mapper(i , j -1) << " ";
                //     std::cout << "right: " << index_mapper(i , j + 1) << std::endl;
                //     std::cout.flush();
                // }
                new_grid[index_mapper(i, j)] = do_calculation(i, j);
            }
        }
    }

    /**
     * @brief do_iteration_boundary calculte the new value of edges and apply the boundary condition
     * and then swap the new grid, grid and old_grid. The update order is important. First we need to
     * try to update the Nonboundary data and then to apply the boudary condition.
     *
     * @param out
     */
    void do_iteration_boundary(int out)
    {

        // update the left
        if (locationInfo._col_low_index != 0)
        {
            for (int i = locationInfo._row_low_index; i < locationInfo._row_high_index; i++)
            {
                if (i == 0 || i == jmax - 1)
                    continue;
                int j = locationInfo._col_low_index;
                new_grid[index_mapper(i, j)] = do_calculation(i, j);
            }
        }

        // update the right
        if (locationInfo._col_high_index != imax)
        {
            for (int i = locationInfo._row_low_index; i < locationInfo._row_high_index; i++)
            {
                if (i == 0 || i == jmax - 1)
                    continue;
                int j = locationInfo._col_high_index - 1;
                new_grid[index_mapper(i, j)] = do_calculation(i, j);
            }
        }

        // update the bottom
        if (locationInfo._row_low_index != 0)
        {
            for (int j = locationInfo._col_low_index; j < locationInfo._col_high_index; j++)
            {
                if (j == 0 || j == imax - 1)
                    continue;
                int i = locationInfo._row_low_index;
                new_grid[index_mapper(i, j)] = do_calculation(i, j);
            }
        }

        // update the top
        if (locationInfo._row_high_index != jmax)
        {
            for (int j = locationInfo._col_low_index; j < locationInfo._col_high_index; j++)
            {
                if (j == 0 || j == imax - 1)
                    continue;
                int i = locationInfo._row_high_index - 1;
                new_grid[index_mapper(i, j)] = do_calculation(i, j);
            }
        }

        //// start to apply the boundary condition
        // update the left boundary
        if (locationInfo._col_low_index == 0)
        {
            for (int i = locationInfo._row_low_index; i < locationInfo._row_high_index; i++)
            {
                new_grid[index_mapper(i, locationInfo._col_low_index)] =
                    new_grid[index_mapper(i, locationInfo._col_low_index + 1)];
            }
        }

        // update the right boundary
        if (locationInfo._col_high_index == imax)
        {
            for (int i = locationInfo._row_low_index; i < locationInfo._row_high_index; i++)
            {
                new_grid[index_mapper(i, locationInfo._col_high_index - 1)] =
                    new_grid[index_mapper(i, locationInfo._col_high_index - 2)];
            }
        }

        // update the bottom boundary
        if (locationInfo._row_low_index == 0)
        {
            for (int j = locationInfo._col_low_index; j < locationInfo._col_high_index; j++)
            {
                new_grid[index_mapper(locationInfo.ghost_row_low_index, j)] =
                    new_grid[index_mapper(locationInfo._row_low_index + 1, j)];
            }
        }

        // update the top boundary
        if (locationInfo._row_high_index == jmax)
        {
            for (int j = locationInfo._col_low_index; j < locationInfo._col_high_index; j++)
            {
                new_grid[index_mapper(locationInfo._row_high_index - 1, j)] =
                    new_grid[index_mapper(locationInfo._row_high_index - 2, j)];
            }
        }

        t += dt;

        // swap the grid data
        double *temp = old_grid;
        old_grid = grid;
        grid = new_grid;
        new_grid = temp;
    }

    /**
     * @brief do_munication send the edge data to other processors and revieve the data from them to
     * storage it into global variable boundary
     *
     * @param requests
     * @param cnt
     * @param tag_num
     */
    void do_communication(MPI_Request *requests, int &cnt, int tag_num)
    {
        comm_data **send_data = new comm_data *[4];
        send_data[0] = locationInfo.ghost_row_high_index == jmax ? nullptr : new comm_data[locationInfo._chunk_col_size]();
        send_data[2] = locationInfo.ghost_row_low_index == 0 ? nullptr : new comm_data[locationInfo._chunk_col_size]();
        send_data[1] = locationInfo.ghost_col_high_index == imax ? nullptr : new comm_data[locationInfo._chunk_row_size]();
        send_data[3] = locationInfo.ghost_col_low_index == 0 ? nullptr : new comm_data[locationInfo._chunk_row_size]();

        // send and recv data to the right
        if (send_data[1] != nullptr)
        { // not the left boundary
            // find the comm processor id
            int com_id = id_from_index(locationInfo._local_row, locationInfo._local_col + 1);

            // collect the send data
            // std::vector<int> index(0);
            for (int i = locationInfo._row_low_index, j = 0; i < locationInfo._row_high_index; i++, j++)
            {
                // index.push_back(index_mapper(i, locationInfo._col_high_index -1));
                send_data[1][j].data = grid[index_mapper(i, locationInfo._col_high_index - 1)];
            }

            // send data
            MPI_Isend(send_data[1], locationInfo._chunk_row_size, comm_data::MPI_type, com_id, tag_num, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;

            // print the send data:
            // std::cout << "processor " << id << " send the right data to the " << com_id << "    ";
            // std::cout << "Data:";
            // for(int i = 0; i < locationInfo._chunk_row_size; i ++) {
            // std::cout << send_data[1][i].data << ":" << index[i] << " ";
            // }
            // // std::cout << std::endl; std::cout.flush();

            boundary.right_boundary = new double[locationInfo._chunk_row_size]();
            // recv the data from comm_id
            MPI_Irecv(boundary.right_boundary, locationInfo._chunk_row_size, comm_data::MPI_type, com_id, tag_num, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;
        }

        // send the left
        if (send_data[3] != nullptr)
        {
            // find the comm processor id
            int com_id = id_from_index(locationInfo._local_row, locationInfo._local_col - 1);

            // collect the send data
            for (int i = locationInfo._row_low_index, j = 0; i < locationInfo._row_high_index; i++, j++)
            {
                send_data[3][j].data = grid[index_mapper(i, locationInfo._col_low_index)];
            }

            // send data
            MPI_Isend(send_data[3], locationInfo._chunk_row_size, comm_data::MPI_type, com_id, tag_num, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;

            // print the send data:
            // std::cout << "processor " << id << " send the left data to the " << com_id << "    ";
            // std::cout << "Data:";
            // for(int i = 0; i < locationInfo._chunk_row_size; i ++) {
            // std::cout << send_data[3][i].data << " ";
            // }
            // // std::cout << std::endl; std::cout.flush();

            boundary.left_boundary = new double[locationInfo._chunk_row_size]();
            // recv the data from comm_id
            MPI_Irecv(boundary.left_boundary, locationInfo._chunk_row_size, comm_data::MPI_type, com_id, tag_num, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;
        }

        // send to the top
        if (send_data[0] != nullptr)
        {
            // find the comm processor id
            int com_id = id_from_index(locationInfo._local_row + 1, locationInfo._local_col);

            // collect the send data
            for (int i = locationInfo._col_low_index, j = 0; i < locationInfo._col_high_index; i++, j++)
            {
                send_data[0][j].data = grid[index_mapper(locationInfo._row_high_index - 1, i)];
            }

            // send data
            MPI_Isend(send_data[0], locationInfo._chunk_col_size, comm_data::MPI_type, com_id, tag_num, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;

            // print the send data:
            // std::cout << "processor " << id << " send the top data to the " << com_id << "    ";
            // std::cout << "Data:";
            // for(int i = 0; i < locationInfo._chunk_row_size; i ++) {
            // std::cout << send_data[0][i].data << " ";
            // }
            // // std::cout << std::endl; std::cout.flush();

            boundary.top_boundary = new double[locationInfo._chunk_col_size]();
            // recv the data from comm_id
            MPI_Irecv(boundary.top_boundary, locationInfo._chunk_col_size, comm_data::MPI_type, com_id, tag_num, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;
        }

        // send to the bottom
        if (send_data[2] != nullptr)
        {
            // find the comm processor id
            int com_id = id_from_index(locationInfo._local_row - 1, locationInfo._local_col);

            // collect the send data
            for (int i = locationInfo._col_low_index, j = 0; i < locationInfo._col_high_index; i++, j++)
            {

                send_data[2][j].data = grid[index_mapper(locationInfo._row_low_index, i)];
            }

            // send data
            MPI_Isend(send_data[2], locationInfo._chunk_col_size, comm_data::MPI_type, com_id, tag_num, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;

            // print the send data:
            // std::cout << "processor " << id << " send the bottom data to the " << com_id << "    ";
            // std::cout << "Data:";
            // for(int i = 0; i < locationInfo._chunk_row_size; i ++) {
            // std::cout << send_data[2][i].data << " ";
            // }
            // // std::cout << std::endl; std::cout.flush();

            boundary.bottom_boundary = new double[locationInfo._chunk_col_size]();
            // recv the data from comm_id
            MPI_Irecv(boundary.bottom_boundary, locationInfo._chunk_col_size, comm_data::MPI_type, com_id, tag_num, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;
        }
    }

    /**
     * @brief update_boundary use the global variable boudary
     *
     */
    void update_boundary()
    {
        // print the ghost state
        if (boundary.top_boundary != nullptr)
        {
            // std::cout << "processor " << id << "has the top data  ";
            // for(int i = 0; i < locationInfo._chunk_col_size; i ++) {
            // std::cout << boundary.top_boundary[i] << " ";
            // }
            // // std::cout << std::endl; std::cout.flush();
            std::memcpy(&grid[index_mapper(locationInfo.ghost_row_high_index - 1, locationInfo._col_low_index)],
                        boundary.top_boundary, sizeof(double) * locationInfo._chunk_col_size);
        }

        // print the ghost state
        if (boundary.right_boundary != nullptr)
        {
            // std::cout << "processor " << id << "has the right data  ";
            for (int i = 0; i < locationInfo._chunk_row_size; i++)
            {
                // std::cout << boundary.right_boundary[i] << " ";
                grid[index_mapper(i + locationInfo._row_low_index, locationInfo.ghost_col_high_index - 1)] = boundary.right_boundary[i];
            }
            // // std::cout << std::endl; std::cout.flush();
        }

        // print the ghost state
        if (boundary.bottom_boundary != nullptr)
        {
            // std::cout << "processor " << id << "has the bottom data  ";
            for (int i = 0; i < locationInfo._chunk_col_size; i++)
            {
                // std::cout << boundary.bottom_boundary[i] << " ";
            }
            // // std::cout << std::endl; std::cout.flush();
            std::memcpy(&grid[index_mapper(locationInfo.ghost_row_low_index, locationInfo._col_low_index)],
                        boundary.bottom_boundary, sizeof(double) * locationInfo._chunk_col_size);
        }

        // print the ghost state
        if (boundary.left_boundary != nullptr)
        {
            // std::cout << "processor " << id << "has the left data  ";
            for (int i = 0; i < locationInfo._chunk_row_size; i++)
            {
                // std::cout << boundary.left_boundary[i] << " ";
                grid[index_mapper(i + locationInfo._row_low_index, locationInfo.ghost_col_low_index)] = boundary.left_boundary[i];
            }
            // // std::cout << std::endl; std::cout.flush();
        }
    }

    /**
     * @brief do_iteration first communicate with the processor to exchange the edge values,
     * then while doing communication, do_iteration update the wave value expect for the edgas using
     * do_iteration_no_boundary, after that using wait_all to make a barrier untill all the communication
     * finish. After that update the edge data and apply the boundary condition.
     *
     * @param out
     */
    void do_iteration(int out)
    {
        MPI_Request *requests = new MPI_Request[8];
        int cnt = 0;
        int tag_num = 1;

        // communication will send and recv the boundary data
        //    std::cout << "processor " << id << " starts to do_communication "  << std::endl; std::cout.flush();
        do_communication(requests, cnt, tag_num);
        //    std::cout << "processor " << id << " ends to do_communication "  << std::endl; std::cout.flush();

        // update the data expect the boundary
        //    std::cout << "processor " << id << " starts to do_iteration_no_boundary "  << std::endl; std::cout.flush();
        do_iteration_no_boundary();
        //    std::cout << "processor " << id << " ends to do_iteration_no_boundary "  << std::endl; std::cout.flush();
        // grid_to_file_with_flag(out, "_no_boundart_");
        // wait for communication finished
        //    std::cout << "processor " << id << " starts to wait for communication "  << std::endl; std::cout.flush();
        MPI_Waitall(cnt, requests, MPI_STATUS_IGNORE);
        //    std::cout << "processor " << id << " ends to wait for communication "  << std::endl; std::cout.flush();

        update_boundary();
        // grid_temp_to_file(out);
        // grid_to_file_with_flag(out, "_update_ghost_state_");

        // update the boundary data
        //    std::cout << "processor " << id << " starts to do_iteration_boundary "  << std::endl; std::cout.flush();
        do_iteration_boundary(out);
        //    std::cout << "processor " << id << " starts to do_iteration_boundary "  << std::endl; std::cout.flush();
    }

    void parse_paramater(int argc, char **argv)
    {
        if (argc == 2)
        {
            imax = atoi(argv[1]);
            jmax = imax;
        }
    }

    int parallel_wave_equation_basic(int argc, char **argv)
    {
        // Basic MPI Init
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &id);
        MPI_Comm_size(MPI_COMM_WORLD, &p);

        double elapsed_time = -MPI_Wtime();

        comm_data::buildMPIType();

        // user pass imax, jmax
        parse_paramater(argc, argv);

        std::cout << "processor " << id << " has finished the init " << std::endl;
        std::cout.flush();

        // calculate the grid size
        std::cout << "processor " << id << " has started the find " << std::endl;
        std::cout.flush();
        find_dimensions(p, rows, cols);
        std::cout << "processor " << id << " has finished the find " << std::endl;
        std::cout.flush();

        // grid size configuration -
        dx = x_max / ((double)imax - 1);
        dy = y_max / ((double)jmax - 1);
        t = 0.0;
        dt = 0.1 * std::min(dx, dy) / c;

        // calculate the processor location
        std::cout << "processor " << id << " starts to init location info " << std::endl;
        std::cout.flush();

        init_location_info();
        if (id == 1)
            print_location_info();
        boundary.left_boundary = nullptr;
        boundary.right_boundary = nullptr;
        boundary.top_boundary = nullptr;
        boundary.bottom_boundary = nullptr;
        std::cout << "processor " << id << " ends to init location info " << std::endl;
        std::cout.flush();

        // create the grid
        grid = new double[locationInfo.total_size]();
        new_grid = new double[locationInfo.total_size]();
        old_grid = new double[locationInfo.total_size]();

        // init the grid
        std::cout << "processor " << id << " starts to half_sinusoidal_intitial_disturbance " << std::endl;
        std::cout.flush();
        half_sinusoidal_intitial_disturbance(1, 3, 3);
        // half_sinusoidal_intitial_disturbance(2, 5, 5);
        // fake_initial_disturbance();
        std::cout << "processor " << id << " starts to ihalf_sinusoidal_intitial_disturbance " << std::endl;
        std::cout.flush();

        std::cout << "processor " << id << " starts to grid_to_file " << std::endl;
        std::cout.flush();
        int out_cnt = 0, it = 0;
        grid_to_file(out_cnt);
        out_cnt++;
        t_out += dt_out;
        std::cout << "processor " << id << " ends to grid_to_file " << std::endl;
        std::cout.flush();

        while (t < t_max)
        {
            {

                // while (it < 2) {
                do_iteration(out_cnt);

                // Note that I am outputing at a fixed time interval rather than after a fixed number of time steps.
                // This means that the output time interval will be independent of the time step (and thus the resolution)
#ifdef output
                if (t_out <= t)
                {
                    std::cout << "Processor: " << id << " generate the output: " << out_cnt << "\tt: " << t << "\titeration: " << it << std::endl;
                    std::cout.flush();
                    grid_to_file(out_cnt);
                    out_cnt++;
                    t_out += dt_out;
                }
#endif
                it++;
            }
        }

#ifdef collect_time
        elapsed_time += MPI_Wtime();
        if (!id)
        {
            printf("wave_equation (%d) %10.6f\n", p, elapsed_time);
        }
#endif

        // delete []grid;
        // delete []new_grid;
        // delete []old_grid;

        MPI_Finalize();
        return 0;
    }
}

#endif
