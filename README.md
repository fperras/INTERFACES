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

Precompiled binaries are available in the binaries directory for windows, linux, and a linux cluster (all x86)

Note: All needed REDOR library files (see REDOR libraries directory) must be included in the running 
directory of the program, together with the input file, any REDOR data files, and the starting structure, 
given as a mol2 file. Please reference the publication and the examples for more details. 

If a REDOR library does not already exist for the support and/or nuclide of interest one can be generated by using 
the "library_maker.c" program, or by contacting Frederic Perras. A few parameters in the library_maker.c file need 
to be edited for each case. The program can be compiled with the following command. More details are given in the 
comments in the source code.
mpicc library_maker.c -o library_maker -lm -O3 -march=native

Copyright 2022.  Iowa State University.  This material was produced under U.S. Government contract DE-AC02-07CH11358 
for the Ames National Laboratory, which is operated by Iowa State University for the U.S. Department of Energy.  The 
U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT, AMES NATIONAL 
LABORATORY, NOR IOWA STATE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF 
THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, 
so as not to confuse it with the version available from The Ames Laboratory.

Additionally, this program is free software; you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your 
option) any later version.  Accordingly, this program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
General Public License for more details. 
