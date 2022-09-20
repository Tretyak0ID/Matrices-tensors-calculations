#include <iostream>
#include <omp.h>
#include "LU-decomposition.cpp"

int main(){
    int SIZE = 256 * 4;
    Matrix<double> Mat(SIZE,SIZE), L(SIZE,SIZE), U(SIZE,SIZE);
    Mat.fill_rand(1, 50);

    double start_time = omp_get_wtime();

    CholeskyLU(Mat, L, U);

    double end_time = omp_get_wtime();

    std::cout << end_time - start_time << std::endl;

    return 0;
}