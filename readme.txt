Everything in the program is instructed in the console log.

The program used a special library called "Eigen". It is a libraty for vector and matrix operation.
The library file will be include in the folder called "Eigen".

The following header MUST be include:
	#include <iostream>
	#include <vector>
	#include <cstdlib>
	#include <time.h>
	#include <fstream>
	#include <cmath>
	#include <numeric>
	#include <stdio.h>
	#include <iomanip>
	#include <../package/Eigen/Eigen/Dense>
	#include <../package/Eigen/Eigen/LU>


The url of the library: https://eigen.tuxfamily.org/

The class Sim_prim refers to lattices points contructed by arbirtary primitives vector
the class Atom_Other is a inheritance class from Sim_prim.It refers to atom structure by input basis for lattice structre in Sim_prim

Other function or class object is function as its name