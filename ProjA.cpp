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
	string name="Cubic";
	vector<double3> point;
	Cubic(int x, int y, int z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		for (int i = 0; i < vect.x; i++) {
			for (int j = 0; j <vect.y ; j++) {
				for (int k = 0; k< vect.z; k++) {
					point.push_back({ double(i) * siz,double(j) * siz,double(k) * siz });
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
	vector<double3> point;
	string name = "BCC";
	BCC(int x, int y, int z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		for (int i = 0; i < vect.x; i++) {
			for (int j = 0; j < vect.y; j++) {
				for (int k = 0; k < vect.z; k++) {
					point.push_back({ double(i) * siz,double(j) * siz,double(k) * siz });
					point.push_back({ (double(i) + 0.5) * siz,(double(j) + 0.5) * siz,(double(k) + 0.5) * siz });

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
	vector<double3> point;
	string name = "FCC";
	FCC(double x, double y, double z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		for (int i = 0; i < vect.x; i++) {
			for (int j = 0; j < vect.y; j++) {
				for (int k = 0; k < vect.z; k++) {
					point.push_back({ double(i) * siz,double(j) * siz,double(k) * siz });
					point.push_back({ (double(i)) * siz,(double(j) + 0.5) * siz,(double(k) + 0.5) * siz });
					point.push_back({ (double(i) + 0.5) * siz,(double(j)) * siz,(double(k) + 0.5) * siz });
					point.push_back({ (double(i) + 0.5) * siz,(double(j) + 0.5) * siz, (double(k)) * siz });

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
	string name = "Diamond";
	Diamond(double x, double y, double z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		for (int i = 0; i < vect.x ; i++) {
			for (int j = 0; j < vect.y ; j++) {
				for (int k = 0; k < vect.z ; k++) {
					atom.push_back({ double(i) * siz,double(j) * siz,double(k) * siz });
					atom.push_back({ double(i) * siz,(double(j)+0.5) * siz,(double(k)+0.5) * siz });
					atom.push_back({ (double(i)+0.25) * siz,(double(j)+0.75) * siz,(double(k)+0.75) * siz });
					atom.push_back({ (double(i)+0.5) * siz,double(j) * siz,(double(k)+0.5) * siz });
					atom.push_back({ (double(i)+0.5) * siz,(double(j)+0.5) * siz,double(k) * siz });
					atom.push_back({ (double(i)+0.75) * siz,(double(j)+0.25) * siz,(double(k)+0.75) * siz });
					atom.push_back({ (double(i)+0.75) * siz,(double(j)+0.75) * siz,(double(k)+0.25) * siz });
					atom.push_back({ (double(i)+0.25) * siz,(double(j)+0.25) * siz,(double(k)+0.25) * siz });

				}
			}
		}

	};
};

void print(vector<double3> cell, string filename) {
	ofstream Outttfile(filename+".xyz", ofstream::trunc);
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
	Cubic cell_1(x, y, z,Latconst);
	FCC cell_2(x, y, z, Latconst);
	BCC cell_3(x, y, z, Latconst);
	Diamond cell_4(x, y, z, Latconst);
	print(cell_1.point, cell_1.name);
	print(cell_2.point, cell_2.name);
	print(cell_3.point, cell_3.name);
	print(cell_4.atom, cell_4.name);
	return 0;

}