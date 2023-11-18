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
#include <map>
double const Well = 1;

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

	

	int atom_number() {
		return atom.size();
	}
	void print_atom_position(int n) {
		cout << "The position of atom" << n << "is:";
		cout << "(" << atom.at(n).x << " " << atom.at(n).y << " " << atom.at(n).z << ")" << endl;
	}

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

	void PBC_prim(double3& a, double3 b = { 0,0,0 }, double3 vector = {0,0,0}, double n=1) {
		double3 r = difference(b, a);
		if (vector.x != 0) {
			if (r.x > vector.x*n / 2) {
				r.x = fmod(r.x + vector.x*n / 2, vector.x*n) - vector.x*n / 2;
			};
			if (r.x < -vector.x*n / 2) {
				r.x = fmod(r.x - vector.x*n / 2, vector.x*n) + vector.x*n / 2;
			};
		};

		if (vector.y != 0) {
			if (r.y > vector.y*n/ 2) {
				r.y = fmod(r.y + vector.y*n / 2, vector.y*n) - vector.y*n / 2;
			};
			if (r.y < -vector.y*n/ 2) {
				r.y = fmod(r.y - vector.y*n / 2, vector.y*n) + vector.y*n / 2;
			};
		};
		if (vector.z != 0) {
			if (r.z > vector.z*n/ 2) {
				r.z = fmod(r.z + vector.z*n / 2, vector.z*n) - vector.z*n / 2;
			};
			if (r.z < -vector.z*n / 2 ) {
				r.z = fmod(r.z - vector.z*n / 2, vector.z*n) + vector.z*n / 2;
			};
		};

		a = sum(r, b);
	}



	void PBC(double3& a, double3 b = { 0,0,0 }) {
		PBC_prim(a, b, conv_vec.at(0), vect.x);
		PBC_prim(a, b, conv_vec.at(1), vect.y);
		PBC_prim(a, b, conv_vec.at(2), vect.z);
	}

	void evaluate(int n, vector<double3>& image ) {
		for (auto& i : image) {
			PBC(i, point.at(n));
		};
	}

	void evaluate_atom(int n, vector<double3>& image ) {
		for (auto& i : image) {
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

	void evaluate_image(vector<double3> &image, int i) {
		if (atom.size() == 0) {
			image = point;
			evaluate(i, image);

		}
		else {
			image = atom;
			evaluate_atom(i,  image);
		}
	}


	void find_neighbor_list(int i) {
		vector<double3> a;
		double3 control;
		
		evaluate_image(a,i);

		vector<double> distance;
		for (auto& k : a  ) {
			distance.push_back(norm(difference(k, a.at(i))));
		}
		distance.at(i) = {};
		int n = 0; double min_r = min_value(distance);
		vector<int> list;
		for (auto& k : distance) {
			if (min_r - 2048.0 * DBL_MIN <= k && k <= min_r + 2048.0 * DBL_MIN) {
				list.push_back(n);
			}
			n++;
		}
		if (min_r==norm(conv_vec.at(0))|| min_r == norm(conv_vec.at(1))|| min_r == norm(conv_vec.at(2))) {
			cout << "The first - neibouring is the image of the particle itself" << endl;
			cout << "The distance is "<< min_r<<endl;
			cout << "The coordinates of the self image particles of the nearest neibouring are:" << endl;

			for (auto& k : min_Vector(min_r)) {
				cout << "(" <<a.at(i).x + conv_vec.at(k).x << " " << a.at(i).y + conv_vec.at(k).y << " " << a.at(i).z + conv_vec.at(k).z << ")" << endl;
				cout << "(" << -1 * (a.at(i).x + conv_vec.at(k).x) << " " << -1 * (a.at(i).y + conv_vec.at(k).y) << " " << -1 * (a.at(i).z + conv_vec.at(k).z) << ")" << endl;
			}

			return;
		};
		
		cout << "Distance for first - neibouring particle for atom"<<i<<" is : " << min_r<<endl;
		cout << "input distance cutoff" << endl; double cutoff; cin >> cutoff;
		cout << "The coordinate are" << endl<<"atom"<<i << " (" << a.at(i).x << " " << a.at(i).y << " " << a.at(i).z << ")" << endl;
		for (auto& k : list) {
			control = a.at(k);
			cout << "atom" << k << "(" << a.at(k).x << " " << a.at(k).y << " " << a.at(k).z << ")" << endl;
			frac_coor(control);
			if (fmod(control.x,0.5*vect.x)==0&& fmod(control.y, 0.5*vect.y) == 0&& fmod(control.z,vect.z* 0.5) == 0) {
				cout << "atom" << k << "(" << -1*a.at(k).x << " " << -1*a.at(k).y << " " << -1*a.at(k).z << ")" << endl;
			}
			if (fmod(control.x, 0.5 * vect.x)==0 ) {
				cout << "atom" << k << "(" <<   a.at(k).x- conv_vec.at(0).x << " " << a.at(k).y - conv_vec.at(0).y << " " << a.at(k).z - conv_vec.at(0).z << ")" << endl;
			}
			if (fmod(control.x, 0.5 * vect.x) == 0 && fmod(control.y, 0.5 * vect.y)==0) {
				cout << "atom" << k << "(" << a.at(k).x - conv_vec.at(0).x - conv_vec.at(1).x << " " <<  a.at(k).y - conv_vec.at(0).y - conv_vec.at(1).y << " " << a.at(k).z - conv_vec.at(0).z - conv_vec.at(1).z << ")" << endl;
			}
			if (fmod(control.x, 0.5 * vect.x) == 0 && fmod(control.z , vect.z* 0.5) == 0) {
				cout << "atom" << k << "(" << a.at(k).x - conv_vec.at(0).x - conv_vec.at(2).x << " " <<  a.at(k).y - conv_vec.at(0).y - conv_vec.at(2).y << " " << a.at(k).z - conv_vec.at(0).z - conv_vec.at(2).z << ")" << endl;
			}
			if (fmod(control.y, 0.5 * vect.y) == 0 && fmod(control.z * vect.z, 0.5) == 0) {
				cout << "atom" << k << "(" << a.at(k).x - conv_vec.at(1).x - conv_vec.at(2).x << " " <<  a.at(k).y - conv_vec.at(1).y - conv_vec.at(2).y << " " << a.at(k).z - conv_vec.at(1).z - conv_vec.at(2).z << ")" << endl;
			}
			if (fmod(control.y, 0.5 * vect.y) == 0) {
				cout << "atom" << k << "(" << a.at(k).x - conv_vec.at(1).x << " " << a.at(k).y - conv_vec.at(1).y << " " << a.at(k).z - conv_vec.at(1).z << ")" << endl;
			}
			if (fmod(control.z , vect.z* 0.5) == 0) {
				cout << "atom" << k << "(" << a.at(k).x - conv_vec.at(2).x << " " << a.at(k).y - conv_vec.at(2).y << " " << a.at(k).z - conv_vec.at(2).z << ")" << endl;
			}

		}
		
	}

	double LJ_potential() {

	}

	double min_value(vector<double>& a) {
		double min = *max_element(a.begin(), a.end());
		for (auto& i : a) {
			if (min > i && i != 0) {
				min = i;
			}
		}
		return min;
	};


	vector<int> min_Vector(double &r) {
		vector<int> list;
		if (r == norm(conv_vec.at(0))) {
			list.push_back(0);
		}
		if (r == norm(conv_vec.at(1))) {
			list.push_back(1);
		}
		if (r == norm(conv_vec.at(2))) {
			list.push_back(2);
		}
		return list;
	}
	
};

class Cubic : public Lattice{
public:
	 
	Cubic(int x, int y, int z, double size):Lattice() {
		siz = size;
		conv_vec.push_back({ siz,0,0 });
		conv_vec.push_back({ 0,siz,0 });
		conv_vec.push_back({ 0,0,siz });
		vect.x = x;
		vect.y = y;
		vect.z = z;

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
		siz = size;
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
	Atom_Othor(int x, int y, int z, lattice_vector la_vec, string atom_name, vector<double3> Basis) : Sim_Prim(x, y, z, la_vec) {
		name = atom_name;
		for (const double3& i : point) {
			for (auto const& j : Basis) {
				atom.push_back(i + (Todouble3(con_vec_mat * ToVector(j))));
			}
		}
	}
};

class FCC:public Lattice{

public:

	FCC(int x, int y, int z, double size) :Lattice() {
		siz = size;
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
	cout << "input number of period in x y z direction" << endl;
	cin >> x >> y >> z;
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
	cell_4 = new Sim_Prim(x, y, z, latticee);


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

	cout << "neibouring list" << endl;
	cout << "input number of period in x y z direction" << endl;
	cin >> x >> y >> z;
	Num = 0;
Re2:
	cout << "input a_1 in (x,y,z)" << endl;
	cin >> latticee.a1.x >> latticee.a1.y >> latticee.a1.z;
	cout << "input a_2 in (x,y,z)" << endl;
	cin >> latticee.a2.x >> latticee.a2.y >> latticee.a2.z;
	cout << "input a_3 in (x,y,z)" << endl;
	cin >> latticee.a3.x >> latticee.a3.y >> latticee.a3.z;
	if (dot(latticee.a1, cross(latticee.a2, latticee.a3)) == 0) {
		cout << "No coplanar sets of vectors! Plaese eneter new sets of vectors" << endl;
		goto Re2;
	}
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

	cell_5 = new Atom_Othor(x, y, z, latticee, atom_name, basis);
	cout << "Which atom number you would like to evaluate? There are total " << cell_5->atom_number()<<" atoms"<<endl;
	cin >> n;
	cout << "You are evaluating atom"<<n<<endl;
	cell_5->print_atom_position(n);
	cell_5->find_neighbor_list(n);
	delete cell_5;
	basis = {};
	cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
	goto Lab_2;

exit:

	return 0;

}