#include <mpi.h>
#include <cassert>
#include <iostream>
#include <random>
#include <ctime>
#include <vector>
#include <cmath>
#include <stdexcept>

static int offset = 0;

std::vector <double> getRandomMatrix(int rows, int cols, double min_value, double max_value) {
    std::mt19937 gen;
    gen.seed((unsigned)time(0) + ++offset);
    std::uniform_real_distribution<> dis(min_value, max_value);
    std::vector <double> a(rows * cols);
    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            a[index++] = dis(gen);
        }
    }
    return a;
}

std::vector <double> solveSequential(const std::vector <double> &a, size_t rows, size_t cols) {
    if (rows * cols != a.size()) {
        throw std::runtime_error("Matrix sizes does not match");
    }
    if (rows + 1 != cols) {
        throw std::runtime_error("Incorrect amount of rows and cols");
    }

    std::vector <double> result(rows);
    std::vector <double> b(a);
    for (size_t k = 0; k < rows; ++k) {
        for (size_t i = k; i < rows; ++i) {
            double temp = b[i * cols + k];
            for (size_t j = 0; j < cols; ++j)
                b[i * cols + j] /= temp;
            if (i != k) {
                for (size_t j = 0; j < cols; ++j) {
                    b[i * cols + j] -= b[k * cols + j];
                }
            }
        }
    }

    for (int k = static_cast<int>(rows) - 1; k >= 0; --k) {
        result[k] = b[k * cols + cols - 1];
        for (int i = 0; i < k; i++) {
            b[i * cols + cols - 1] -= b[i * cols + k] * result[k];
        }
    }
    return result;
}

std::vector <double> solveParallel(const std::vector <double> &a, size_t rows, size_t cols) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = cols / size;
    const int rem = cols % size;

    int code = 0;

    if (rows * cols != a.size()) {
        code = 1;
    }
    MPI_Bcast(&code, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (code != 0) {
        throw std::runtime_error("Matrix sizes does not match");
    }

    if (rows + 1 != cols) {
        code = 2;
    }
    MPI_Bcast(&code, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (code != 0) {
        throw std::runtime_error("Incorrect amount of rows and cols");
    }

    std::vector <double> v((delta + (rank < rem ? 1 : 0)) * rows);

    if (rank == 0) {
        for (int proc = size - 1; proc >= 0; --proc) {
            int index = 0;
            for (size_t j = proc; j < cols; j += size) {
                for (size_t i = 0; i < rows; ++i) {
                    v[index++] = a[i * cols + j];
                }
            }
            if (proc > 0) {
                MPI_Send(v.data(), index, MPI_DOUBLE, proc, 1, MPI_COMM_WORLD);
            }
        }
    } else {
        MPI_Status stat;
        MPI_Recv(v.data(), v.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
    }

    std::vector <double> pivotCol(rows);
    for (size_t row = 0; row < rows; ++row) {
        if (static_cast<int>(row) % size == rank) {
            int index = 0;
            for (size_t i = rows * (row / size); i < rows * (row / size + 1); ++i) {
                pivotCol[index++] = v[i];
            }
            assert(index == rows);
        }
        MPI_Bcast(pivotCol.data(), rows, MPI_DOUBLE, row % size, MPI_COMM_WORLD);
        double pivotRow = pivotCol[row];
        for (int j = row / size; j < (delta + (rank < rem ? 1 : 0)); ++j) {
            double pivotC = v[j * rows + row];
            for (size_t k = 0; k < rows; ++k) {
                if (k == row) {
                    v[j * rows + k] /= pivotRow;
                } else {
                    v[j * rows + k] -= pivotC * pivotCol[k] / pivotRow;
                }
            }
        }
    }

    if ((cols - 1) % size == (size_t)rank) {
        MPI_Request rq;
        MPI_Isend(v.data() + ((cols - 1) / size) * rows, rows, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &rq);
    }
    if (rank == 0) {
        v.resize(rows);
        MPI_Status stat;
        MPI_Recv(v.data(), rows, MPI_DOUBLE, (cols - 1) % size, 2, MPI_COMM_WORLD, &stat);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    return v;
}

static const double EPS = 1e-9;

static bool checkEqual(const std::vector <double> &a, const std::vector <double> &b) {
    if (a.size() != b.size()) {
        return 0;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::abs(a[i] - b[i]) > EPS) {
            return false;
        }
    }
    return true;
}

static bool checkSolution(const std::vector <double> &a, size_t rows, size_t cols,
                          const std::vector <double> &x) {
    if (rows * cols != a.size()) {
        throw std::runtime_error("Matrix sizes does not match");
    }
    if (rows + 1 != cols) {
        throw std::runtime_error("Incorrect amount of rows and cols");
    }

    std::vector <double> result(rows, 0);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < rows; ++j) {
            result[i] += a[i * cols + j] * x[j];
        }
    }
    for (size_t i = 0; i < rows; ++i) {
        if (result[i] - a[i * cols + rows] > EPS) {
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> a;
    const int rows = 5;
    const int cols = 6;

    if (rank == 0) {
        a = {
            3, -2, 1, 5, 6, 0,
            5, 1, 3, 4, 5, 0,
            1, 2, 3, 9, 8, 0,
            7, 5, -5, 9, 1, 0,
            1, -8, -2, 9, 1, 0
        };
    }

    std::vector <double> answer = solveParallel(a, rows, cols);
    if (rank == 0) {
        std::vector <double> seqAnswer = solveSequential(a, rows, cols);
        if (checkEqual(seqAnswer, answer) == true) { printf("true");} else { puts("false (error)"); };
        if (checkSolution(a, rows, cols, answer) == true) { printf("true");} else { puts("false (error)"); };
    }
    
    MPI_Finalize();
    return 0;
}
