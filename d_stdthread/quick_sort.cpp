#include <cstdlib>
#include <time.h>
#include <stdio.h>
#include <limits.h>
#include <thread>

template<typename T>
void swp(T& a, T& b)
{ T temp = a;
  a = b;
  b = temp;
}

void quickSort(double* array, int left, int right)
{
    while (right > left)
	{
        int iterate_left = left;
        int iterate_right = right;
        double pivot = array[(left + right) >> 1];

        while (iterate_left <= iterate_right)
		{
            while (array[iterate_left] < pivot)
			{
                iterate_left += 1;
            }
            while (array[iterate_right] > pivot)
			{
                iterate_right -= 1;
            }
            if (iterate_left <= iterate_right)
			{
                double temp = array[iterate_left];
                array[iterate_left] = array[iterate_right];
                array[iterate_right] = temp;

                iterate_left += 1;
                iterate_right -= 1;
            }
        }

        if ((iterate_left << 1) > left + right)
		{
            quickSort(array, iterate_left, right);
            right = iterate_left - 1;
        }
		else
		{
            quickSort(array, left, iterate_left - 1);
            left = iterate_left;
        }
    }
}

void even(double* array, double* tmp, int left_part, int right_part)
{
	int i = 0, j;
    	while(i < left_part)
		{
    		tmp[i] = array[i];
    		i += (1 << 1);
    	}

    	double* array2 = array + left_part;
    	int a = 0, b = 0;
    	i = a;

    	for (i; (a < left_part) && (b < right_part); i += 2)
		{
    		array[i] = tmp[a];

    		if (tmp[a] <= array2[b])
			{
            	a += 2;
        	}
			else
			{
            	array[i] = array2[b];
            	b += 2;
        	}
    	}

    	j = b;
    	if (a == left_part)
		{
    		while(j < right_part)
			{
    			array[i] = array2[j];
    			j += 2;
    			i += 2;
    		}
    	}
		else
		{
    		j = a;
    		while(j < left_part)
			{
    			array[i] = tmp[j];
    			j += 2;
    			i += 2;
    		}
    	}
}

void odd(double* array, double* tmp, int left_part, int right_part)
{
	int i = 1, j;
    	while(i < left_part)
		{
    		tmp[i] = array[i];
    		i += 2;
    	}

    	double* array2 = array + left_part;
    	int a = 1, b = 1;
    	i = a;

    	for (i; (a < left_part) && (b < right_part); i += 2)
		{
    		array[i] = tmp[a];

    		if (tmp[a] <= array2[b])
			{
            	a += 2;
        	}
			else
			{
            	array[i] = array2[b];
            	b += 2;
        	}
    	}

    	j = b;
    	if (a == left_part)
		{
    		while(j < right_part)
			{
    			array[i] = array2[j];
    			j += 2;
    			i += 2;
    		}
    	}
		else
		{
    		j = a;
    		while(j < left_part)
			{
    			array[i] = tmp[j];
    			j += 2;
    			i += 2;
    		}
    	}
}

void quick(double* array, double* tmp, int size, int part) {
	if (size <= part) {
        	quickSort(array, 0, size - 1);
    } else {
    	int divide = size >> 1;
        int partial = divide + divide % 2;

		std::thread sortL(quick, array, tmp, partial, part);
		std::thread sortR(quick, array + partial, tmp + partial, size - partial, part);
		sortL.join();
		sortR.join();

		std::thread batcherE(even, array, tmp, partial, size - partial);
		std::thread batcherO(odd, array, tmp, partial, size - partial);
		batcherE.join();
		batcherO.join();

		auto lambda = [&]() {
			int i = 1;
            while (i < (size + 1) >> 1) {
            	if (array[i << 1] < array[(i << 1) - 1]) {
            		swp(array[(i << 1) - 1], array[i << 1]);
				}
        		i += 1;
            }
		};
		std::thread prll(lambda);
		prll.join();
    }
}
void quickSort__STD(double* array, int threads, int size)
{
    double* temporary = new double[size];

    int portion = size / threads;
    if (size % threads)
        portion += 1;

    std::thread start(quick, array, temporary, size, portion);
    start.join();
    delete[]temporary;
}

void getRandomArray(double* arr, int size)
{
	int i = 0;
	double number;

	while(i < size) {
		number = rand() / (RAND_MAX + 1.0);
		arr[i] = number;
		i += 1;
	}
}

bool isSorted(double* ar, int size) {
    const double *previous_value = ar;

    while (size) {
       if (*ar < *previous_value)
             return false;
       previous_value = ar;

       ++ar;
       --size;
     }
     return true;
}

int main(void) {
	srand(time(NULL));

	int size = 20000;
	int threads = 4;

    double* std = new double[size];
    double* seq = new double[size];
    getRandomArray(std, size);

    for (int i = 0; i < size; i++) {
		seq[i] = std[i];
	}

	clock_t start, end;
	float seconds;

    start = clock();
	quickSort__STD(std, threads, size);
    end = clock();
    seconds = (float)(end - start) / CLOCKS_PER_SEC;
    printf("(STD) time for quicksort = %f \n", seconds);

    start = clock();
    quickSort(seq, 0, size - 1);
    end = clock();
	seconds = (float)(end - start) / CLOCKS_PER_SEC;
    printf("(Sequential) time for quicksort = %f \n", seconds);

    if (isSorted(std, size) )
        printf("Correctly sorted\n");
    else
        printf("Incorretly sorted\n");

    if (isSorted(seq, size) )
        printf("Correctly sorted\n");
    else
        printf("Incorretly sorted\n");

    for (int i = 0; i < size; i++) {
    	if (std[i] != seq[i]) {
    		puts("not equal"); break;
    	}
    }

    delete[]std;
    delete[]seq;

	return 0;
}
