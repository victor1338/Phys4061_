#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <fstream>
#include <cmath>
#include <numeric>
#include <stdio.h>
#include <iomanip>
#include <map>

typedef struct double3 {
	double x, y, z;
}double3;
typedef struct lattice_vector {
	double3 a1, a2, a3;
};

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
	double3 sizz;
	string name;
	vector<double3> point;
	vector<double3> prim_vec;
	double prim_Vol;
	vector<double3> recip_vec;
	double recip_Vol;
	vector<double3> conv_vec;
	vector<double3> atom;

public:



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

	void evaluate(int n) {
		for (auto& i : point) {
			PBC(i, point.at(n));
		};
	}
	void evaluate_atom(int n) {
		for (auto& i : atom) {
			PBC(i, atom.at(n));
		};
	}

	double3 frac_coor(double3& a) {
		return { a.x / siz,a.y / siz,a.z / siz };
	};

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
	};

	

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
	 
	Cubic(int x, int y, int z, double size):Lattice() {
		conv_vec.push_back({ siz,0,0 });
		conv_vec.push_back({ 0,siz,0 });
		conv_vec.push_back({ 0,0,siz });
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

class BCC :public Lattice {
private:


public:
	BCC(int x, int y, int z, double size) :Lattice() {
		conv_vec.push_back({ siz,0,0 });
		conv_vec.push_back({ 0,siz,0 });
		conv_vec.push_back({ 0,0,siz });
		vect.x = x;
		vect.y = y;
		vect.z = z;
		siz = size;
		name = "BCC";
		prim_vec.push_back({ -siz / 2,siz / 2,siz / 2 });
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
};

class Sim_Prim: public Lattice {
private:

public:
	Sim_Prim(int x, int y , int z, double3 size):Lattice() {
		sizz = size;
		conv_vec.push_back({ sizz.x,0,0 });
		conv_vec.push_back({ 0,sizz.y,0 });
		conv_vec.push_back({ 0,0,sizz.z });
		vect.x = x;
		vect.y = y;
		vect.z = z;
		name = "Simple Orthorhombic";
		prim_vec.push_back({ sizz.x,0,0 });
		prim_vec.push_back({ 0,sizz.y,0 });
		prim_vec.push_back({ 0,0,sizz.z });
		prim_Vol = Volume(prim_vec);
		recip_vec.push_back({ cross(prim_vec.at(1),prim_vec.at(2),2.0 * pi / prim_Vol) });
		recip_vec.push_back({ cross(prim_vec.at(2),prim_vec.at(0),2.0 * pi / prim_Vol) });
		recip_vec.push_back({ cross(prim_vec.at(0),prim_vec.at(1),2.0 * pi / prim_Vol) });
		recip_Vol = Volume(recip_vec);

		for (int i = 0; i < vect.x; i++) {
			for (int j = 0; j < vect.y; j++) {
				for (int k = 0; k < vect.z; k++) {
					point.push_back({ double(i) * sizz.x,double(j) * sizz.y,double(k) * sizz.z });
				}
			}
		}
	}
};

class Atom_Othor : public Sim_Prim {
public:
	Atom_Othor(int x, int y, int z, double3 size, string atom_name, map<int,double3> Basis) : Sim_Prim(x,y,z,size) {
		name = atom_name;
		for (const double3& i : point) {
			for (auto const& j : Basis) {
				atom.push_back({ i.x + j.second.x*sizz.x,i.y + j.second.y*sizz.y  ,i.z + j.second.z*sizz.z });
			}
		}

	}
};

class FCC:public Lattice{

public:

	FCC(int x, int y, int z, double size) :Lattice() {
		conv_vec.push_back({ siz,0,0 });
		conv_vec.push_back({ 0,siz,0 });
		conv_vec.push_back({ 0,0,siz });
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
	Diamond(int x, int y, int z, double size) :FCC(x,y,z,size) {
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
	int n = 0;
	int Num = 0;
	double3 vector;
	double Latconst;
	int x, y, z;
	Cubic* cell_1 ;
	BCC* cell_2 ;
	FCC* cell_3 ;
	Diamond* diamond ;
	Atom_Othor* cell_4;
	map<int, double3> basis;
	string atom_name;
	lattice_vector latticee;
home:

	cout << "This is Project A in PHYS_4061" << endl<<endl;
	cout << setprecision(3);
	cout << "Select the Lab number to start the simulation (e.g. input 1 will go to Lab_1 ), or 0 to end the program" << endl;
	cout << "1: Lab_1" << endl << "2: Lab_2"<<endl<<"0: exit"<<endl;
	cin >> n;
	cout << "-----------------------------------------------------------------------------------------------"<<endl<<endl;
	switch (n) {
	case 1:
		goto Lab_1;
		break;
	case 2:
		goto Lab_2;
		break;
	case 0:
		goto exit;
		break;
	};
Lab_1:
	cout << "Lab 1:" << endl;
	cout << "input lattice constant" << endl;
	cin >> Latconst;
	cout << "input number of period in x y z direction" << endl;
	cin >> x >> y >> z;

	cell_1 = new Cubic(x, y, z, Latconst);
	cell_2 = new BCC(x, y, z, Latconst);
	cell_3 = new FCC(x, y, z, Latconst);
	diamond = new Diamond(x, y, z, Latconst);

	cell_1->print_file();
	cell_2->print_file();
	cell_3->print_file();
	diamond->print_file_atom();
	delete cell_1;
	delete cell_2;
	delete cell_3;
	delete diamond;
	cout << "finished printing xyz file for the Simple Cubic, FCC, BCC and diamond structures"<<endl;
	cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
	goto home;

Lab_2:

	cout << "Lab 2: " << endl;
	cout << "input lattice constant" << endl;
	cin >> Latconst;
	cout << "input number of period in x y z direction" << endl;
	cin >> x >> y >> z;
	cell_1 = new Cubic(x, y, z, Latconst);
	cell_2 = new BCC(x, y, z, Latconst);
	cell_3 = new FCC(x, y, z, Latconst);
	diamond = new Diamond(x, y, z, Latconst);

Task_2_1:

	// task 1 
	print_volume(*cell_1, *cell_2, *cell_3);
	print_vector(*cell_1, *cell_2, *cell_3);

Task_2_2:

	// task 2
	cell_1 = new Cubic(x, y, z, Latconst);
	cout << "input an coordinate in \( x,y,z \)" << endl;
	cin >> vector.x >> vector.y >> vector.z;
	cout << "it's cooredinate after applying periodic boundary condition in \( x,y,z \) is: (range in absolute coordinate is ["<<-Latconst*0.5<<","<<Latconst*0.5 << "])" << endl<<endl;
	cell_1->PBC(vector);
	cout << "(" << vector.x << " , " << vector.y << " , " << vector.z << ")" << endl;
	vector = cell_1->frac_coor(vector);
	cout << "The respective frational coordinate is: ( in range of (0.5,0.5] )";
	cout << "(" << vector.x << " , " << vector.y << " , " << vector.z << ")" << endl;

Task_2_3:
	Num = 0;
	
	cout << "input a_1 in (x,y,z)" << endl;
	cin >> vector.x >> vector.y >> vector.z;
	latticee.a1 = vector;
	cout << "input a_2 in (x,y,z)" << endl;
	cin >> vector.x >> vector.y >> vector.z;
	latticee.a2 = vector;
	cout << "input a_3 in (x,y,z)" << endl;
	cin >> vector.x >> vector.y >> vector.z;
	latticee.a3 = vector;
	cout << "input number of atom" << endl;
	cin >> Num;
	for (int i = 0; i < Num; i++) {
		cout << "Input the fractional coordinate of " << endl;
	}



	basis.insert(pair<int, double3>(Num, vector));
	cell_4 = new Atom_Othor(x, y, z, {}, "asdasd", basis);
	

	// input a_1 a_2 a_3 to form orthorhombic/Triclinic/
	//PBC on point n in the different structure





	goto home;

exit:

	return 0;

}