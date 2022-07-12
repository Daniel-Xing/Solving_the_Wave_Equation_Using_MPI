#include "../src/Mpi_CourseWork.h"
#include "iostream"
#include <gtest/gtest.h>

TEST(generate_test, BasicAssertions) {
   int m, n;
   find_dimensions(10, m, n);
   EXPECT_TRUE(m == 2 && n == 5);

   find_dimensions(110, m, n);
   EXPECT_TRUE(m == 10 && n == 11);
}

TEST(location_info_test, BasicAssertions) {
   _location_info locationInfo;
   int p = 15;
   int rows, cols;
   find_dimensions(p, rows, cols);

   int imax = 10, jmax = 10;
   for(int i = 0 ; i < 15; i ++) {
       init_location_info(i, p, imax, jmax, rows, cols, locationInfo);
       std::cout << "location info:" << std::endl;
       std::cout << locationInfo._local_row << " " << locationInfo._local_col << std::endl;
       std::cout << locationInfo._chunk_row_size << " " << locationInfo._chunk_col_size << std::endl;
       std::cout << locationInfo.col_low_index << " " << locationInfo.col_high_index << std::endl;
       std::cout << locationInfo.row_low_index << " " << locationInfo.row_high_index << std::endl;
       std::cout << locationInfo.total_size << std::endl;
       std::cout << std::endl;
   }

   for(int i = locationInfo.row_low_index ; i < locationInfo.row_high_index ; i ++) {
//        if (i == 0 || i == jmax -1) continue;
       for(int j = locationInfo.col_low_index; j < locationInfo.col_high_index ; j ++) {
//            if (j == 0 || j == jmax -1) continue;
           std::cout << i << " " <<
           j << " " <<  index_mapper(i, j) << std::endl;
       }
   }
}