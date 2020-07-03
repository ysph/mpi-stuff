#include <mpi.h>
#include <iostream>

bool is_sorted(double* source, int size) {
    for (int i = 0; i < size - 1; i++) {
        if (source[i + 1] < source[i]) {
            return false;
        }
    }
    return true;
}

void pass(double* source, double* dest, int size, int offset) {
    unsigned char* pmem = (unsigned char*)source;
    int counters[256];
    int sum = 0;

    for (int i = 0; i < 256; i++) counters[i] = 0;

    for (int i = 0; i < size; i++) counters[pmem[8 * i + offset]]++;

    for (int i = 0; i < 256; i++) {
        int tmp = counters[i];
        counters[i] = sum;
        sum += tmp;
    }

    for (int i = 0; i < size; i++) {
        int index = 8 * i + offset;
        dest[counters[pmem[index]]] = source[i];
        counters[pmem[index]]++;
    }
}

void last_pass(double* source, double* dest, int size, int offset) {
    unsigned char* pmem = (unsigned char*)source;
    int counters[256];
    int sum = 0;

    for (int i = 0; i < 256; i++) counters[i] = 0;

    for (int i = 255; i > 127; i--) {
        counters[i] += sum;
        sum = counters[i];
    }

    for (int i = 0; i < 128; i++) {
        int tmp = counters[i];
        counters[i] = sum;
        sum += tmp;
    }

    for (int i = 0; i < size; i++) {
        int index = 8 * i + offset;
        if (pmem[index] < 128) {
            dest[counters[pmem[index]]] = source[i];
            counters[pmem[index]]++;
        }
        else {
            counters[pmem[index]]--;
            dest[counters[pmem[index]]] = source[i];
        }
    }
}

void ordered_merge(double* source1, int size1, double* source2, int size2, double* dest) {
    int i = 0, j = 0, k = 0;

    while ((i < size1) && (j < size2)) {
        if (source1[i] < source2[j]) {
            dest[k] = source1[i];
            i++;
        }
        else {
            dest[k] = source2[j];
            j++;
        }
        k++;
    }
    while (i < size1) {
        dest[k] = source1[i];
        i++;
        k++;
    }
    while (j < size2) {
        dest[k] = source2[j];
        j++;
        k++;
    }
}

void seq_sort(double* source, int size) {
    double* temp = NULL;
    double* dest = new double[size];
    for (int i = 0; i < 8; i++) {
        pass(source, dest, size, i);
        temp = source;
        source = dest;
        dest = temp;
    }
    last_pass(source, dest, size, 7);
    delete[] dest;
}

void par_sort(double** source, int size) {
    int proc_rank = 0, proc_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    int lenght = size / proc_size, addition = size % proc_size;

    double* dest = NULL;
    double* temp = NULL;
    if (0 == proc_rank) {
        dest = new double[lenght + addition];
        temp = new double[size];
    }
    else {
        dest = new double[lenght];
    }
    int* displs = new int[proc_size];
    int* scounts = new int[proc_size];

    displs[0] = 0;
    scounts[0] = lenght + addition;
    for (int i = 1; i < proc_size; i++) {
        displs[i] = addition + lenght * i;
        scounts[i] = lenght;
    }

    MPI_Scatterv(*source, scounts, displs, MPI_DOUBLE, dest, scounts[proc_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    seq_sort(dest, scounts[proc_rank]);
    MPI_Gatherv(dest, scounts[proc_rank], MPI_DOUBLE, temp, scounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (0 == proc_rank) {
        double* tmp = NULL;

        for (int i = 0; i < proc_size - 1; i++) {
            ordered_merge(temp, displs[i + 1], &temp[displs[i + 1]], scounts[i + 1], *source);
            tmp = *source;
            *source = temp;
            temp = tmp;
        }
        tmp = *source;
        *source = temp;
        temp = tmp;
    }
    delete[] dest;
    delete[] temp;
    delete[] displs;
    delete[] scounts;
}

int main(void) {
    int initialized, finalized;

    MPI_Initialized(&initialized);
    if (!initialized)
        MPI_Init(NULL, NULL);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (0 == rank) {
        int size = 5;
        double* vector = new double[size] { 2.0, 8.66667, 1.69, 99.147, 57.3579 };
        bool result = false;

        seq_sort(vector, size);
        result = is_sorted(vector, size);

        if (result == true) {
            std::cout << "Alright" << std::endl;
        }
    }

    rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (0 == rank) {
        int size = 1000;
        double* vector = new double[size];
        for (int i = 0; i < size; i++) vector[i] = (size - i) * 0.001 + (i + size / 2);
        bool result = false;

        seq_sort(vector, size);
        result = is_sorted(vector, size);

        if (result == true) {
            std::cout << "Alright" << std::endl;
        }
    }

    rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int size = 1000;
    double* vector = new double[size];
    for (int i = 0; i < size; i++) vector[i] = (size - i) * 0.001 + (i + size / 2);

    par_sort(&vector, size);

    if (0 == rank) {
        bool result = false;
        result = is_sorted(vector, size);

        if (result == true) {
            std::cout << "Alright" << std::endl;
        }
    }
    delete[] vector;

    rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    size = 5;
    vector = new double[size] { 99.147, 8.66667, 1.69, 2.0, 57.3579 };

    par_sort(&vector, size);

    if (0 == rank) {
        bool result = false;
        result = is_sorted(vector, size);

        if (result == true) {
            std::cout << "Alright" << std::endl;
        }
    }
    delete[] vector;

    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();

    return 0;
}
