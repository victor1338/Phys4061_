#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <fstream>
#include <cmath>
#include <numeric>
#include <stdio.h>

typedef struct double3 {
	double x, y, z;
}double3;

double3 difference(double3 b, double3 a) {
	return { a.x - b.x,a.y - b.y,a.z - b.z };
};
double3 sum(double3 a, double3 b) {
	return { a.x + b.x,a.y + b.y,a.z + b.z };
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

double norm(double3 a) {
	return sqrt(pow(a.x, 2)+pow(a.y,2)+pow(a.z,2));
};

class Lattice {
protected:

	double3 vect;
	double siz;
	string name;
	vector<double3> point;
	vector<double3> prim_vec;
	double prim_Vol;
	vector<double3> recip_vec;
	double recip_Vol;
	vector<double3> conv_vec;
	vector<double3> atom;

public:

	Lattice(double siz) {
		conv_vec.push_back({ siz,0,0 });
		conv_vec.push_back({ 0,siz,0 });
		conv_vec.push_back({ 0,0,siz });
	};

	void print_file() {
		ofstream Outttfile("./xyz/" + name + ".xyz", ofstream::trunc);
		Outttfile << point.size() << "\nframe " << 0 << "\n";
		for (auto &i : point) {
			Outttfile << "Particle" << " " << i.x << " " << i.y << " " << i.z << "\n";
		}
	};

	void print_file_atom() {
		ofstream Outttfile("./xyz/" + name + ".xyz", ofstream::trunc);
		Outttfile << atom.size() << "\nframe " << 0 << "\n";
		for (auto &i : atom) {
			Outttfile << "Particle" << " " << i.x << " " << i.y << " " << i.z << "\n";
		}
	};

	void PBC(double3 &a, double3 b = {0,0,0}) {
		double3 r= difference(b, a);
		if (r.x > siz/2) {
			r.x = fmod(r.x + siz / 2, siz) - siz / 2;
		};
		if (r.x <= -siz/2||r.x==-siz/2) {
			r.x = fmod(r.x - siz / 2, siz) + siz / 2;
		};
		if (r.y >  siz/2) {
			r.y = fmod(r.y + siz / 2, siz) - siz / 2;
		};
		if (r.y <= - siz/2||r.y==-siz/2) {
			r.y = fmod(r.y - siz / 2, siz) + siz / 2;
		};
		if (r.z > siz/2) {
			r.z = fmod(r.z + siz / 2, siz) - siz / 2;
		};
		if (r.z <= -siz/2||r.z==-siz/2) {
			r.z = fmod(r.z - siz / 2, siz) + siz / 2;
		};
		a = sum(r, b);

	}


	void print_prim_Vol() {
		cout << name << ":" << prim_Vol << endl;
	};

	void print_recip_Vol() {
		cout << name << ":" << recip_Vol << endl;
	};

	void print_prim_vect() {
		cout << endl << name << endl;
		print_vect_set(prim_vec);
	};

	void print_recip_vect() {
		cout << endl << name << endl;
		print_vect_set(recip_vec);
	};

	void print_vect(double3 a) {
		cout << "(" << a.x << " , " << a.y << " , " << a.z << ")";
	}

	void print_vect_set(vector<double3> a) {
		for (auto& i : a) {
			cout << endl;
			print_vect(i);
			cout << endl;
		}
	};
};

class Cubic : public Lattice{
public:
	 
	Cubic(int x, int y, int z, double size):Lattice(size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		name = "Cubic";
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

class BCC:public Lattice {
private:


public:
	BCC(int x, int y, int z, double size):Lattice (size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		name = "BCC";
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

class FCC:public Lattice{

public:

	FCC(double x, double y, double z, double size) :Lattice(size) {
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		name = "FCC";
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

class Diamond:public FCC{

public:
	Diamond(double x, double y, double z, double size) :FCC(x,y,z,size) {
		name = "Diamond";
		for (const double3 &i : point) {
			atom.push_back({i.x,i.y,i.z});
			atom.push_back({i.x+0.25*siz,i.y+0.25*siz,i.z+0.25*siz });
		}

	};
};


void print_volume(Cubic& a, BCC& b, FCC& c) {
	cout << "The volume of primitive cell of Cubic, BCC, FCC are:" << endl;
	a.print_prim_Vol(); b.print_prim_Vol(); c.print_prim_Vol();
	cout << "The volume of reciprocal primitive cell of Cubic, BCC, FCC are:" << endl;
	a.print_recip_Vol(); b.print_recip_Vol(); c.print_recip_Vol();
};

void print_vector(Cubic &a, BCC &b, FCC &c) {
	cout << "The primitive vector for each lattice cell are:" << endl;
	a.print_prim_vect();
	b.print_prim_vect();
	c.print_prim_vect();

	cout << "The reciprocal vector for each lattice cell are:" << endl;
	a.print_recip_vect();
	b.print_recip_vect();
	c.print_recip_vect();
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
	Diamond diamond(x, y, z, Latconst);
	double3 vector;
	cout << "input an coordinate in \( x,y,z \)" << endl;
	cin >> vector.x >> vector.y >> vector.z;
	cout << "it's cooredinate after applying periodic boundary condition in \( x,y,z \) is: (range is ["<<-Latconst*0.5<<","<<Latconst*0.5 << "])" << endl<<endl;
	cell_1.PBC(vector, {1,1,1});
	cout << vector.x << " " << vector.y << " " << vector.z << endl;


	print_vector(cell_1, cell_2,cell_3);
	print_volume(cell_1, cell_2, cell_3);

	cell_1.print_file();
	cell_2.print_file();
	cell_3.print_file();
	diamond.print_file_atom();
	return 0;

}