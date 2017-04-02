# mpitumor
I've writen this program as an assignment of a parallel programming course.
This software simulates a glioma tumor growth inside a 3d brain. The brain is represented as a 3d matrix of difussion coeficients. 
Then an initial small tumor is placed somewhere inside the matrix and it's difussion and growth is simulated in several interations.
Difussion equation is simulated using the finite differences method.  The codes has been paralelized an can be compiled in 3 versions:
Serial (no additional libraries requiered), OpenMPI and MPI.
The program runs for 1000 iterations and outputs the total concentration of tumor cells.
Every 100 iterations, a VTK file is generating containg the 3d points occupied by tumor cells.
In order to compile, you need a working C++ compiler (I've tested it under Linux and MacOS)

- comp_serial.sh -> compiles the serial version (main_serial) 
- comp_openmp -> compiles the OpenMP version (main_openmp)
- comp_mpi -> compiles the MPI version (main_mpi)

After compiling, run directly the generated executable. In order to run the programs properly you must first decompress 
Cerebro.csv.tar.gz and place it along the desired executable.

I've edited main_serial_ijk.cpp to translate comments to english

