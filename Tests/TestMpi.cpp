#include <iostream>
#include <mpi.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <string>

int id, p;
int rows, cols;

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

int main(int argc, char **argv)
{
    // Basic MPI Init
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    find_dimensions(p, rows, cols);
    //    std::cout << "processor " << id << " has finished the find "  << std::endl; std::cout.flush();


    MPI_Finalize();
    return 0;
}
