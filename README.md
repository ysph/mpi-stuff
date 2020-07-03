# parallel-stuff
Messing up with different technologies that capable of making code run in parallel.

## What projects are there?

### MPI
1. Calculating the average value of the sum of the elements of the vector
2. Solving the system of linear equations using Gaussian elimination
3. Radix sort

### OpenMP, TBB, std::thread
1. convex_hull.cpp - Building the smallest convex set of a shape  
2. quick_sort.cpp - Divide and Conquer using quicksort
3. block_matrix_multiply.cpp - Multiply matrixes of a different size

## Install
### MPI
Ubuntu:
```bash
sudo apt-get install mpic++
```

Arch:
```bash
pacman -S mpic++
```

RedHat and CentOS:
```bash
yum install mpic++
```

### OpenMP (omp)
Ubuntu:
```bash
sudo apt-get install libomp-dev
```

Arch:
```bash
pacman -S libomp-dev
```

RedHat and CentOS:
```bash
yum install libomp-dev
```

### TBB
Ubuntu:
```bash
sudo apt-get install libtbb-dev
```

Arch:
```bash
pacman -S libtbb-dev
```

RedHat and CentOS:
```bash
yum install libtbb-dev
```

## Complie

In \<folder\>
```bash
g++ <file>.cpp -pthread -ltbb -mpic++ -fopenmp
```

## Usage
Launch file you'd like to test:
```bash
./<name_of_the_project>
```

