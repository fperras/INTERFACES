# INTERFACES
Program to determine the 3D structure of a surface species using REDOR NMR data

There are two options for using the program:

1) If you are planning on submitting jobs directly on that machine using cmd or a terminal
then use the main.cpp file. Compile it using OpenMP and O3 optimizations. For instance:
g++ main.cpp -o INTERFACES -fopenmp -lm -O3 -march=native

2) If you are planning on submitting jobs to a cluster the main_cluster.cpp file takes instead
an argument for the input file. Compile it using OpenMP and O3 optimizations. For instance:
g++ main_cluster.cpp -o INTERFACES -fopenmp -lm -O3 -march=native

To submit a job to a cluster, use the following command in your job submission script:
./INTERFACES input_file > log.txt

Precompiled binaries are available in the binaries directory for windows, linux, and a linux cluster

Note: All needed REDOR library files (see REDOR libraries directory) must be included in the running 
directory of the program, together with the input file, any REDOR data files, and the starting structure, 
given as a mol2 file. Please reference the publication and the examples for more details.
