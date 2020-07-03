#include <cstdlib>
#include <time.h>
#include <stdio.h>
#include <limits.h>
#include <tbb/tbb.h>

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
    	
    	for (i; (a < left_part) && (b < right_part); i += (1 << 1)) 
		{
    		array[i] = tmp[a];
    		
    		if (tmp[a] <= array2[b]) 
			{
            	a += (1 << 1);
        	} 
			else
			{
            	array[i] = array2[b];
            	b += (1 << 1);
        	}
    	}
    	
    	j = b;
    	if (a == left_part) 
		{
    		while(j < right_part) 
			{
    			array[i] = array2[j];
    			j += (1 << 1);
    			i += (1 << 1);
    		}
    	} 
		else 
		{
    		j = a;
    		while(j < left_part)
			{
    			array[i] = tmp[j];
    			j += (1 << 1);
    			i += (1 << 1);
    		}
    	}
}

void odd(double* array, double* tmp, int left_part, int right_part)
{
	int i = 1, j;	
    	while(i < left_part) 
		{
    		tmp[i] = array[i];
    		i += (1 << 1);
    	}
        	
    	double* array2 = array + left_part;
    	int a = 1, b = 1;
    	i = a;
    	
    	for (i; (a < left_part) && (b < right_part); i += (1 << 1)) 
		{
    		array[i] = tmp[a];
    		
    		if (tmp[a] <= array2[b]) 
			{
            	a += (1 << 1);
        	} 
			else
			{
            	array[i] = array2[b];
            	b += (1 << 1);
        	}
    	}
    	
    	j = b;
    	if (a == left_part) 
		{
    		while(j < right_part) 
			{
    			array[i] = array2[j];
    			j += (1 << 1);
    			i += (1 << 1);
    		}
    	} 
		else
		{
    		j = a;
    		while(j < left_part)
			{
    			array[i] = tmp[j];
    			j += (1 << 1);
    			i += (1 << 1);
    		}
    	}
}

void quick(double* array, double* tmp, int size, int part) 
{
	if (size <= part) 
	{
        	quickSort(array, 0, size - 1);
    } 
	else 
	{
    		int divide = size >> 1;
        	int partial = divide + divide % (1 << 1);
        	
			tbb::task_group sort;
			tbb::task_group batcher;
			
			sort.run([&]{quick(array, tmp, partial, part);});
			sort.run_and_wait([&]{quick(array + partial, tmp + partial, size - partial, part);});
			
			batcher.run([&]{even(array, tmp, partial, size - partial);});
			batcher.run_and_wait([&]{odd(array, tmp, partial, size - partial);});
			
        	tbb::parallel_for(tbb::blocked_range<int>(1, (size + 1) >> 1),[&](const tbb::blocked_range<int>& r) 
			{ 
            		int i = r.begin();
            		while (i < r.end()) 
					{
            			if (array[i << 1] < array[(i << 1) - 1]) 
						{
            				std::swap(array[(i << 1) - 1], array[i << 1]);
        				}
        				i += 1;
            		}
            });
    	}
}
void quickSort__TBB(double* array, int threads, int size) 
{
    double* temporary = new double[size];

    int portion = size / threads;
    if (size % threads)
        portion += 1;

    tbb::task_group g;
    g.run_and_wait([&]{quick(array, temporary, size, portion);});

    delete[]temporary;
}

void getRandomArray(double* arr, int size) 
{
	int i = 0;
	double number;
	
	while(i < size) 
	{
		number = rand() / (RAND_MAX + 1.0);
		arr[i] = number;
		i += 1;
	}
}

void isSorted(double* arr, int size) 
{
	int i = 0;
	
	while(i < size - 1) 
	{
		if (arr[i] > arr[i + 1]) { printf("\n Incorretly sorted %f > %f \n", arr[i], arr[i + 1]); return; }
		i += 1;
	}
	printf("\nCorrectly sorted\n");
}

int main(void) {
	srand(time(NULL));
	
	int size = 200;
	int threads = 4;
	
    double* tbb = new double[size];
    double* seq = new double[size];
    getRandomArray(tbb, size);
    
    for (int i = 0; i < size; i++) 
	{
		seq[i] = tbb[i];
	}
    
    tbb::tick_count t0 = tbb::tick_count::now();
	quickSort__TBB(tbb, threads, size);
    tbb::tick_count t1 = tbb::tick_count::now();
    printf("(TBB) time for quicksort = %g seconds\n", (t1 - t0).seconds() );
    
    clock_t start = clock();
    quickSort(seq, 0, size - 1);
    clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
    printf("(Sequential) time for quicksort = %f \n", seconds);
    
    isSorted(tbb, size);
    isSorted(seq, size);

    for (int i = 0; i < size; i++) 
	{
    	if (tbb[i] != seq[i]) 
		{
    		puts("not equal"); break;
    	}
    }
    
    delete[]tbb;
    delete[]seq;
    
	return 0;
}
