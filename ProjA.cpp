#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <fstream>
#include <cmath>
#include <numeric>
#include <stdio.h>
#include <iomanip>
#include <Eigen/Eigen/Dense>
#include <Eigen/Eigen/LU>
using Eigen::Matrix3d;
using Eigen::Vector3d;
typedef struct double3 {
	double x, y, z;
}double3;
typedef struct lattice_vector {
	double3 a1, a2, a3;
};
Vector3d ToVector(double3 a) {
	return { a.x,a.y,a.z };
};
double3 Todouble3(Vector3d a) {
	return { a(0),a(1),a(2) };
}

double3 difference(double3 b, double3 a) {
	return { a.x - b.x,a.y - b.y,a.z - b.z };
};
double3 sum(double3 a, double3 b) {
	return { a.x + b.x,a.y + b.y,a.z + b.z };
};
double3 operator+(double3 a,double3 b) {
	return { a.x + b.x,a.y + b.y,a.z + b.z };
}

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
	Matrix3d con_vec_mat;
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

	void PBC_Cubic(double3 &a, double3 b = {0,0,0}) {
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

	void PBC_prim(double3& a, double3 b = { 0,0,0 }, double3 vector = {0,0,0}) {
		double3 r = difference(b, a);

		if (vector.x != 0) {
			if (r.x > vector.x / 2) {
				r.x = fmod(r.x + vector.x / 2, vector.x) - vector.x / 2;
			};
			if (r.x <= -vector.x / 2 || r.x == -vector.x / 2) {
				r.x = fmod(r.x - vector.x / 2, vector.x) + vector.x / 2;
			};
		};

		if (vector.y != 0) {
			if (r.y > vector.y / 2) {
				r.y = fmod(r.y + vector.y / 2, vector.y) - vector.y / 2;
			};
			if (r.y <= -vector.y / 2 || r.y == -vector.y / 2) {
				r.y = fmod(r.y - vector.y / 2, vector.y) + vector.y / 2;
			};
		};
		if (vector.z != 0) {
			if (r.z > vector.z / 2) {
				r.z = fmod(r.z + vector.z / 2, vector.z) - vector.z / 2;
			};
			if (r.z <= -vector.z / 2 || r.z == -vector.z / 2) {
				r.z = fmod(r.z - vector.z / 2, vector.z) + vector.z / 2;
			};
		};

		a = sum(r, b);
	}

	void PBC(double3& a, double3 b = { 0,0,0 }) {
		for (auto const& i : conv_vec) {
			PBC_prim( a, b, i);
		}
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

	double3 frac_coor_othor(double3& a) {
		return { a.x / siz,a.y / siz,a.z / siz };
	};

	void frac_coor(double3& a) {
		Vector3d b(a.x,a.y,a.z);
		b = con_vec_mat.inverse() * b;
		a = {b(0),b(1),b(2)};

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
		
		for (int i = 0; i < 3; i++) {
			con_vec_mat(0, i) = conv_vec.at(i).x;
			con_vec_mat(1, i) = conv_vec.at(i).y;
			con_vec_mat(2, i) = conv_vec.at(i).z;
		}

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

		for (int i = 0; i < 3; i++) {
			con_vec_mat(0, i) = conv_vec.at(i).x;
			con_vec_mat(1, i) = conv_vec.at(i).y;
			con_vec_mat(2, i) = conv_vec.at(i).z;
		}

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
	Sim_Prim(int x, int y , int z, lattice_vector la_vec):Lattice() {
		conv_vec.push_back(la_vec.a1);
		conv_vec.push_back(la_vec.a2);
		conv_vec.push_back(la_vec.a3);
		vect.x = x;
		vect.y = y;
		vect.z = z;
		name = "Simple_Primitive_cell";
		prim_vec.push_back(la_vec.a1);
		prim_vec.push_back(la_vec.a2);
		prim_vec.push_back(la_vec.a3);
		prim_Vol = Volume(prim_vec);
		recip_vec.push_back({ cross(prim_vec.at(1),prim_vec.at(2),2.0 * pi / prim_Vol) });
		recip_vec.push_back({ cross(prim_vec.at(2),prim_vec.at(0),2.0 * pi / prim_Vol) });
		recip_vec.push_back({ cross(prim_vec.at(0),prim_vec.at(1),2.0 * pi / prim_Vol) });
		recip_Vol = Volume(recip_vec);

		for (int i = 0; i < 3; i++) {
			con_vec_mat(0, i) = conv_vec.at(i).x;
			con_vec_mat(1, i) = conv_vec.at(i).y;
			con_vec_mat(2, i) = conv_vec.at(i).z;
		}

		for (int i = 0; i < vect.x; i++) {
			for (int j = 0; j < vect.y; j++) {
				for (int k = 0; k < vect.z; k++) {
					point.push_back({	double(i) * la_vec.a1.x + double(j) * la_vec.a2.x + double(k) * la_vec.a3.x , 
										double(i) * la_vec.a1.y + double(j) * la_vec.a2.y + double(k) * la_vec.a3.y ,
										double(i) * la_vec.a1.z + double(j) * la_vec.a2.z + double(k) * la_vec.a3.z });
				}
			}
		}
	}
};

class Atom_Othor : public Sim_Prim {
public:
	Atom_Othor(int x, int y, int z, lattice_vector la_vec, string atom_name, vector<double3> Basis) : Sim_Prim(x,y,z,la_vec) {
		name = atom_name;
		for (const double3& i : point) {
			for (auto const& j : Basis) {
				atom.push_back(i+(Todouble3(con_vec_mat*ToVector(j))));
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

		for (int i = 0; i < 3; i++) {
			con_vec_mat(0, i) = conv_vec.at(i).x;
			con_vec_mat(1, i) = conv_vec.at(i).y;
			con_vec_mat(2, i) = conv_vec.at(i).z;
		}

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
	int task = 0;
	int Num = 0;
	double3 vectorr;
	vector<double3> basis;
	double Latconst;
	int x, y, z;
	Cubic* cell_1 ;
	BCC* cell_2 ;
	FCC* cell_3 ;
	Diamond* diamond ;
	Sim_Prim* cell_4;
	Atom_Othor* cell_5;
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
	cout << "Select the task number" << endl;
	cout << "1: Volume and vector of primitive cell and reciprocal cell " << endl << "2: Return fractional coordinate of input cooredinate with respect to origin" << endl << "3: Neighbor list" << endl << "0: Back to main menu" << endl;
	cin >> task;
	cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
	switch (task) {
	case 1:
		goto Task_2_1;
		break;
	case 2:
		goto Task_2_2;
		break;
	case 3:
		goto Task_2_3;
		break;
	case 0:
		goto home;
		break;
	};
	diamond = new Diamond(x, y, z, Latconst);

Task_2_1:
	// task 1
	cout << "Volume and vector of primitive cell and reciprocal cell" << endl;
	cout << "input lattice constant" << endl;
	cin >> Latconst;

	cell_1 = new Cubic(2, 2, 2, Latconst);
	cell_2 = new BCC(2, 2, 2, Latconst);
	cell_3 = new FCC(2, 2, 2, Latconst);

	print_volume(*cell_1, *cell_2, *cell_3);
	print_vector(*cell_1, *cell_2, *cell_3);
	delete cell_1;
	delete cell_2;
	delete cell_3;
	cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
	goto Lab_2;


Task_2_2:
	// task 2
	cout << "Return fractional coordinate of input cooredinate with respect to origin" << endl;
Re:

	cout << "input a_1 in (x,y,z)" << endl;
	cin >> latticee.a1.x >> latticee.a1.y >> latticee.a1.z;
	cout << "input a_2 in (x,y,z)" << endl;
	cin >> latticee.a2.x >> latticee.a2.y >> latticee.a2.z;
	cout << "input a_3 in (x,y,z)" << endl;
	cin >> latticee.a3.x >> latticee.a3.y >> latticee.a3.z;
	if (dot(latticee.a1, cross(latticee.a2, latticee.a3)) == 0) {
		cout << "No coplanar sets of vectors! Plaese eneter new sets of vectors" << endl;
		goto Re;
	}
	cell_4 = new Sim_Prim(2, 2, 2, latticee);


	cout << "input an coordinate in \( x,y,z \)" << endl;
	cin >> vectorr.x >> vectorr.y >> vectorr.z;
	cout << "it's cooredinate after applying periodic boundary condition in \( x,y,z \) is: "<<endl;
	cell_4->PBC(vectorr);
	cout << "(" << vectorr.x << " , " << vectorr.y << " , " << vectorr.z << ")" << endl;
	cell_4->frac_coor(vectorr);
	cout << "The respective frational coordinate is: ( in range of (0.5,0.5] )"<<endl;
	cout << "(" << vectorr.x << " , " << vectorr.y << " , " << vectorr.z << ")" << endl;

	delete cell_4;
	cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
	goto Lab_2;

Task_2_3:

	Num = 0;	
	cout << "input a_1 in (x,y,z)" << endl;
	cin >> latticee.a1.x >> latticee.a1.y >> latticee.a1.z;
	cout << "input a_2 in (x,y,z)" << endl;
	cin >> latticee.a2.x >> latticee.a2.y >> latticee.a2.z;
	cout << "input a_3 in (x,y,z)" << endl;
	cin >> latticee.a3.x >> latticee.a3.y >> latticee.a3.z;
	cout << "Basis setup: " << endl;
	cout << "input number of atom in the basis" << endl;
	cin >> Num;

	for (int i = 0; i < Num; i++) {
		cout << "Input the fractional coordinate of atom "<<i+1 << endl;
		basis.push_back({});
		cin >> basis.at(i).x >> basis.at(i).y >> basis.at(i).z;
	}

	cout << "input atom name" << endl;
	cin >> atom_name;

	cell_5 = new Atom_Othor(2, 2, 2, latticee, atom_name, {});
	
	delete cell_5;
	goto Lab_2;

exit:

	return 0;

}