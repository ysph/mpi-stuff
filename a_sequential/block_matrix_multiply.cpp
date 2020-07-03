#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define UPPER_BOUND	10000
#define NUMBER_OF_ELEMENTS 500
#define RANGE 1000.0f

int squareRoot(int n){
  int lo = 0, hi = n, mid;
  for(int i = 0 ; i < 1000; i++){
      mid = (lo+hi) >> 1;
      if(mid * mid == n) return mid;
      if(mid * mid > n) hi = mid;
      else lo = mid;
  }
  return mid;
}

void printMatrix(double* A, int N) {
	for (int i = 0; i < N * N; i++) {
		printf("%f ", A[i]);
	}
}

void generateRandomMatrix(double *matrix, int N) {
	int i = 0;
	
	while(1) {
		if (i == (N * N)) return; else matrix[i] = (rand() % UPPER_BOUND) / RANGE;
		i += 1;
	}
}
void dumbMultiply(double* matrix1, double* matrix2, double* result, int block_size, int n) {
    for (int i = block_size - 1; i >= 0; i--)
        for (int j = block_size - 1; j >= 0; j--)
            for (int k = block_size - 1; k >= 0; k--) {
                result[i * n + j] += matrix1[i * n + k] * matrix2[k * n + j];
            }
}
void blockMultiplication(double *A, double *B, double* C, int n) {
	int divider = squareRoot(n) >> 1;
    int block_size = n / divider;
    int Aindex, Bindex, Cindex;
    
    for (int i = 0; i < divider; ++i) {
        for (int j = 0; j < divider; ++j) {
            for (int k = 0; k < divider; ++k) {
            	Aindex = (i * n + (j + i + k) % divider) * block_size;
            	Bindex = (((i + j + k) % divider) * n + j) * block_size;
            	Cindex = (i * n + j) * block_size;
            	
                dumbMultiply(&A[Aindex], &B[Bindex], &C[Cindex], block_size, n);
            }
        }
    }
}

void checkEquality(double* A, double* C, int N) {
	double temp;
    for (int i = 0; i < N * N; i++) {
    	temp = (A[i] > C[i]) ? A[i] - C[i] : C[i] - A[i];
   		if (temp > 0.000001) {
   			puts("Not equal matrices.");
   			return;
   		}
    }
    puts("Equal matrices.");
}

int main(int argc, char** argv) {
	int N;
	
	srand(time(NULL));
	
	(argc == 2) ? N = atoi(argv[1]) : N = NUMBER_OF_ELEMENTS;
    
    double* matrix1 = new double[N * N]; 
    double* matrix2 = new double[N * N];
    double* result1 = new double[N * N];
    double* result2 = new double[N * N];
    
    generateRandomMatrix(matrix1, N);
    generateRandomMatrix(matrix2, N);

    clock_t begin = clock();
    blockMultiplication(matrix1, matrix2, result1, N);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Sequential first time: %f\n", time_spent);
    
    begin = clock();
    blockMultiplication(matrix1, matrix2, result2, N);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    
    printf("Sequential second time: %f\n", time_spent);
    
    checkEquality(result1, result2, N);
    
    delete[] matrix1;
    delete[] matrix2;
    delete[] result1;
    delete[] result2;
    
    return 0;
}
