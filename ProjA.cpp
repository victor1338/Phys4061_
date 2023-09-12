#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <fstream>
#include <cmath>
#include <numeric>
struct double3 {
	double x, y, z;
};
struct int3 {
	int x, y, z;
};;
using namespace std;
double const pi = 3.14159265358979323846;
double const c = 2997924.58; //in armstrong per ps
double const kb = 8.6173324 * pow(10, -5); //

class Cubic {
private:
	double3 vect;
	double siz;

public:
	vector<double3> atom;
	Cubic(int x, int y, int z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		for (int i = 0; i < vect.x; i++) {
			for (int j = 0; j <vect.y ; j++) {
				for (int k = 0; k< vect.z; k++) {
					double3 point = {double(i) * siz,double(j) * siz,double(k) * siz};
					atom.push_back(point);

				}
			}
		}
	};
};

class BCC {
private:
	int3 vect;
	double siz;

public:
	vector<double3> atom;
	BCC(int x, int y, int z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		for (int i = 0; i < vect.x; i++) {
			for (int j = 0; j < vect.y; j++) {
				for (int k = 0; k < vect.z; k++) {
					double3 point = { double(i) * siz,double(j) * siz,double(k) * siz };
					double3 point2 = { (double(i) + 0.5) * siz,(double(j) + 0.5) * siz,(double(k) + 0.5) * siz };
					atom.push_back(point);
					atom.push_back(point2);

				}
			}
		}

	}
}
;

class FCC {
private:
	double3 vect;
	double siz;

public:
	vector<double3> atom;
	FCC(double x, double y, double z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		for (int i = 0; i < vect.x; i++) {
			for (int j = 0; j < vect.y; j++) {
				for (int k = 0; k < vect.z; k++) {
					double3 point = { double(i) * siz,double(j) * siz,double(k) * siz };
					double3 point2 = { (double(i) ) * siz,(double(j) + 0.5) * siz,(double(k) + 0.5) * siz };
					double3 point3 = { (double(i) + 0.5) * siz,(double(j) ) * siz,(double(k) + 0.5) * siz };
					double3 point4 = { (double(i) + 0.5) * siz,(double(j) + 0.5) * siz, g(double(k)) * siz };
					atom.push_back(point);
					atom.push_back(point2);
					atom.push_back(point3);
					atom.push_back(point4);

				}
			}
		}
	};
};
class Diamond {
private:
	double3 vect;
	double siz;

public:
	vector<double3> atom;
	Diamond(double x, double y, double z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		for (double i = 0; i < vect.x ; i++) {
			for (double j = 0; j < vect.y ; j++) {
				for (double k = 0; k < vect.z ; k++) {
					double3 point = { i * siz,j * siz,k * siz };
					atom.push_back(point);

				}
			}
		}
		for (double i = 0; i < vect.x; i++) {
			for (double j = 0.5; j < vect.y; j++) {
				for (double k = 0.5; k < vect.z; k++) {
					double3 point = { i * siz,j * siz,k * siz };
					atom.push_back(point);

				}
			}
		}
		for (double i = 0.25; i < vect.x ; i++) {
			for (double j = 0.75; j < vect.y; j++) {
				for (double k = 0.75; k < vect.z; k++) {
					double3 point = { i * siz,j * siz,k * siz };
					atom.push_back(point);

				}
			}
		}
		for (double i = 0.5; i < vect.x; i++) {
			for (double j = 0; j < vect.y ; j++) {
				for (double k = 0.5; k < vect.z; k++) {
					double3 point = { i * siz,j * siz,k * siz };
					atom.push_back(point);

				}
			}
		}
		for (double i = 0.5; i < vect.x; i++) {
			for (double j = 0.5; j < vect.y; j++) {
				for (double k = 0; k < vect.z ; k++) {
					double3 point = { i * siz,j * siz,k * siz };
					atom.push_back(point);

				}
			}
		}
		for (double i = 0.75; i < vect.x; i++) {
			for (double j = 0.25; j < vect.y ; j++) {
				for (double k = 0.5; k < vect.z; k++) {
					double3 point = { i * siz,j * siz,k * siz };
					atom.push_back(point);

				}
			}
		}
		for (double i = 0.75; i < vect.x; i++) {
			for (double j = 0.75; j < vect.y; j++) {
				for (double k = 0.25; k < vect.z ; k++) {
					double3 point = { i * siz,j * siz,k * siz };
					atom.push_back(point);

				}
			}
		}
		for (double i = 0.25; i < vect.x ; i++) {
			for (double j = 0.25; j < vect.y ; j++) {
				for (double k = 0.25; k < vect.z ; k++) {
					double3 point = { i * siz,j * siz,k * siz };
					atom.push_back(point);

				}
			}
		}

	};
};

void print(vector<double3> cell) {
	ofstream Outttfile("particle.xyz", ofstream::trunc);
	Outttfile << cell.size() << "\nframe " << 0 << "\n";
	for (double3 i : cell) {
		Outttfile << "Particle" << " " << i.x << " " << i.y << " " << i.z << "\n";
	}
};

int main() {
	double Latconst;
	int x, y, z;
	cout << "input lattice constant" << endl;
	cin >> Latconst;
	cout << "input number of period in x y z direction" << endl;
	cin >> x >> y >> z;
	BCC cell1(x, y, z,Latconst);
	print(cell1.atom);

	return 0;

}