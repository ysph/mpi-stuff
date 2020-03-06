# mpi-attempts
Messing up with MPI and making some basic stuff

1. Calculating the average value of the sum of the elements of the vector
2. Solving the system of linear equations using Gaussian elimination 
3. Radix sort

## Install MPI

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

## Complie

```bash
mpic++ first.cpp -o first && mpic++ second.cpp -o second && mpic++ third.cpp -o third
```

## Usage

Launch file you'd like to test:
```bash
./first
```

