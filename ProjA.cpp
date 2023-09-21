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

double dot(double3 a,double3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
};

double3 cross(double3 a, double3 b,double c=1) {
	return { (a.y * b.z - a.z * b.y)*c , (a.z * b.x - a.x * b.z)*c , (a.x * b.y - a.y * b.x)*c };
};

double Volume(vector<double3> a) {
	return dot(a.at(0), cross(a.at(1), a.at(2)));
};

class Cubic {
private:
	double3 vect;
	double siz;

public:
	string name="Cubic";
	vector<double3> point;
	vector<double3> prim_vec;
	double prim_Vol ;
	vector<double3> recip_vec;
	double recip_Vol;
	Cubic(int x, int y, int z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;

		prim_vec.push_back({ siz,0,0 });
		prim_vec.push_back({ 0,siz,0 });
		prim_vec.push_back({ 0,0,siz });
		prim_Vol = Volume(prim_vec);
		recip_vec.push_back({ cross(prim_vec.at(1),prim_vec.at(2),2.0 * pi / prim_Vol) });
		recip_vec.push_back({ cross(prim_vec.at(2),prim_vec.at(0),2.0 * pi / prim_Vol) });
		recip_vec.push_back({ cross(prim_vec.at(0),prim_vec.at(1),2.0 * pi / prim_Vol) });
		recip_Vol = Volume(recip_vec);

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
	vector<double3> prim_vec;
	double prim_Vol;
	vector<double3> recip_vec;
	double recip_Vol;
	BCC(int x, int y, int z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;

		prim_vec.push_back({ -siz/2,siz / 2,siz / 2 });
		prim_vec.push_back({ siz / 2,-siz / 2,siz / 2 });
		prim_vec.push_back({ siz / 2,siz / 2,-siz / 2 });
		prim_Vol = Volume(prim_vec);
		recip_vec.push_back({ cross(prim_vec.at(1),prim_vec.at(2),2.0 * pi / prim_Vol) });
		recip_vec.push_back({ cross(prim_vec.at(2),prim_vec.at(0),2.0 * pi / prim_Vol) });
		recip_vec.push_back({ cross(prim_vec.at(0),prim_vec.at(1),2.0 * pi / prim_Vol) });
		recip_Vol = Volume(recip_vec);

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
	vector<double3> prim_vec;
	double prim_Vol;
	vector<double3> recip_vec;
	double recip_Vol;

	FCC(double x, double y, double z, double size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;

		prim_vec.push_back({ 0,siz / 2,siz / 2 });
		prim_vec.push_back({ siz / 2,0,siz / 2 });
		prim_vec.push_back({ siz / 2,siz / 2,0 });
		prim_Vol = Volume(prim_vec);
		recip_vec.push_back({ cross(prim_vec.at(1),prim_vec.at(2),2.0 * pi / prim_Vol) });
		recip_vec.push_back({ cross(prim_vec.at(2),prim_vec.at(0),2.0 * pi / prim_Vol) });
		recip_vec.push_back({ cross(prim_vec.at(0),prim_vec.at(1),2.0 * pi / prim_Vol) });
		recip_Vol = Volume(recip_vec);

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

void print_vect(double3 a) {
	cout <<"("<< a.x << " , " << a.y << " , " << a.z<<")";
}

void print_vect_set(vector<double3> a) {
	for (auto& i : a) {
		cout << endl;
		print_vect(i);
		cout << endl;
	}
};

void print(vector<double3> cell, string filename) {
	ofstream Outttfile("./xyz/"+filename + ".xyz", ofstream::trunc);
	Outttfile << cell.size() << "\nframe " << 0 << "\n";
	for (double3 i : cell) {
		Outttfile << "Particle" << " " << i.x << " " << i.y << " " << i.z << "\n";
	}
};

void print_volume(double a, double b, double c, double d, double e, double f) {
	cout << "The volume of primitive cell of Cubic, BCC, FCC are:" << endl;
	cout << a << " " << b << " " << c << endl;
	cout << "The volume of reciprocal primitive cell of Cubic, BCC, FCC are:" << endl;
	cout << d << " " << e << " " << f << endl;
};

void print_vector(Cubic &a, BCC &b, FCC &c) {
	cout << "The primitive vector for each lattice cell are:" << endl;
	cout << endl << a.name << endl;
	print_vect_set(a.prim_vec);

	cout << endl << b.name << endl;
	print_vect_set(b.prim_vec);

	cout << endl << c.name << endl;
	print_vect_set(c.prim_vec);

	cout << "The reciprocal vector for each lattice cell are:" << endl;

	cout << endl << a.name << endl;
	print_vect_set(a.recip_vec);

	cout << endl << b.name << endl;
	print_vect_set(b.recip_vec);

	cout << endl << c.name << endl;
	print_vect_set(c.recip_vec);
};

int main() {
	double Latconst;
	int x, y, z;
	cout << "input lattice constant" << endl;
	cin >> Latconst;
	cout << "input number of period in x y z direction" << endl;
	cin >> x >> y >> z;
	Cubic cell_1(x, y, z,Latconst);
	BCC cell_2(x, y, z, Latconst);
	FCC cell_3(x, y, z, Latconst);
	Diamond cell_4(x, y, z, Latconst);
	
	print_volume(cell_1.prim_Vol, cell_2.prim_Vol, cell_3.prim_Vol, cell_1.recip_Vol, cell_2.recip_Vol, cell_3.recip_Vol);
	print_vector(cell_1, cell_2,cell_3);


	print(cell_1.point, cell_1.name);
	print(cell_2.point, cell_2.name);
	print(cell_3.point, cell_3.name);
	print(cell_4.atom, cell_4.name);
	return 0;

}