#include <mpi.h> 
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>

typedef uint8_t BYTE;
typedef uint32_t DWORD;
typedef int32_t LONG;
typedef int64_t LONGLONG;

typedef union _LARGE_INTEGER {
  struct {
    DWORD LowPart;
    LONG  HighPart;
  };
  struct {
    DWORD LowPart;
    LONG  HighPart;
  } u;
  LONGLONG QuadPart;
} LARGE_INTEGER, *PLARGE_INTEGER;

int main(int argc, char* argv[]) {	
	int TotalSum, ProcSum = 0;
	srand(time(NULL));
	int ProcRank, ProcNum, N = 1000000 + rand() % 20000000;
	int* x = new int[N];
	double t1, t2;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	// подготовка данных
	int part = N/ProcNum;
	int *partOfX = new int[part];
	if (ProcRank == 0) {
		printf("Number of element =  %d", N);
		srand(time(NULL));
		for (int i = 0; i < N; i++) {
			x[i] = 10 + rand() % 1000;
		};
		t1 = MPI_Wtime();
		for (int i = 0; i < ProcNum - 1; i++)
		{
			partOfX = x + i * part;
			MPI_Send(partOfX, part, MPI_INT, i+1, 0, MPI_COMM_WORLD);
		}
	}
	else { // все процессы отсылают свои частичные суммы
		MPI_Recv(partOfX, part, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
		for (int i = 0; i < part; i++)
			ProcSum += partOfX[i];
		MPI_Send(&ProcSum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	// сборка частичных сумм на процессе с рангом 0
	if (ProcRank == 0) {
		TotalSum = ProcSum;
		for (int i = (ProcNum - 1)*part; i < N; i++) {
			TotalSum += x[i];
		}
		for (int i = 1; i < ProcNum; i++) {
			MPI_Recv(&ProcSum, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &Status);
			TotalSum += ProcSum;
		}
	}
	
	// вывод результата
	if (ProcRank == 0) {
		t2 = MPI_Wtime();
		printf("\nTotal Average Sum = %i", (TotalSum/N));
		printf("\nTotal time = %1.10f", t2-t1);
		LARGE_INTEGER frequency;        
		LARGE_INTEGER t1Consistent, t2Consistent;           
		/*
		windows' high-resolution timer
		QueryPerformanceFrequency(&frequency); 
		QueryPerformanceCounter(&t1Consistent); 
		*/
		int sumConsistent = 0;
		for (int i = 0; i < N; i++) {
			sumConsistent += x[i];
		}
		//QueryPerformanceCounter(&t2Consistent);
		//printf("\nTotal Average Sum Consistent = %i", (sumConsistent/N));
		//printf("\nTotal time Consistent = %1.10f", (t2Consistent.QuadPart - t1Consistent.QuadPart)*1.0/frequency.QuadPart);
	}
	MPI_Finalize();
	return 0;
}
