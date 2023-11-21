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
};
bool operator==(double3 a,double3 b) {
	return (a.x == b.x && a.y == b.y && a.z == b.z);
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
double3 operator-(double3 a, double3 b) {
	return { a.x - b.x,a.y - b.y,a.z - b.z };
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

double dist(double3 a,double3 b) {
	return norm(difference(a, b));
}



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

	double size() {
		return atom.size();
	}

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

	void PBC_2(double3& a, double3 b = { 0,0,0 }) {
		double3 re = difference(b,a);
		frac_coor(re);
		if ( re.x > vect.x/2 ) {
			re.x = fmod(re.x + vect.x / 2, vect.x) - vect.x / 2;
		}
		if (re.x < -vect.x / 2) {
			re.x = fmod(re.x - vect.x / 2, vect.x) + vect.x / 2;
		}
		if (re.y > vect.y / 2) {
			re.y = fmod(re.y + vect.y / 2, vect.y) - vect.y / 2;
		}
		if (re.y < -vect.y / 2) {
			re.y = fmod(re.y - vect.y / 2, vect.y) + vect.y / 2;
		}
		if ( re.z > vect.z / 2) {
			re.z = fmod(re.z + vect.z / 2, vect.z) - vect.z / 2;
		}
		if (re.z < -vect.z / 2) {
			re.z = fmod(re.z - vect.z / 2, vect.z) + vect.z / 2;
		}
		Vector3d c(re.x,re.y,re.z);
		c = con_vec_mat * c;
		a = {c(0)+b.x,c(1)+b.y,c(2)+b.z};

	}

	void evaluate(int n, vector<double3>& image ) {
		for (auto& i : image) {
			PBC(i, point.at(n));
		};
	}

	void evaluate_atom(int n, vector<double3>& image ) {
		for (auto& i : image) {
			PBC_2(i, atom.at(n));
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
		double3 bb;
		
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

		if (vect.x == 1 && vect.y == 1 && vect.z == 1 && a.size() == 1) {
			cout << "The first - neibouring is the image of the particle itself" << endl;
			cout << "The distance is " << norm(change_cor(a.at(0),1,0,0)) << endl;
			cout << "The coordinates of the self image particles of the nearest neibouring are:" << endl;
			if (norm(conv_vec.at(0)) <= norm(conv_vec.at(1))&& norm(conv_vec.at(0)) <= norm(conv_vec.at(2))) {
				for (int i = -1; i <= 1; i += 2) {
					bb = change_cor(a.at(0), i, 0, 0);
						cout << "(" << bb.x << " " << bb.y << " " << bb.z << ")" << endl;
				}
			}
			if (norm(conv_vec.at(1)) <= norm(conv_vec.at(0)) && norm(conv_vec.at(1)) <= norm(conv_vec.at(2))) {
				for (int i = -1; i <= 1; i += 2) {
					bb = change_cor(a.at(0), 0, i, 0);
					cout << "(" << bb.x << " " << bb.y << " " << bb.z << ")" << endl;
				}
			}
			if (norm(conv_vec.at(2)) <= norm(conv_vec.at(0)) && norm(conv_vec.at(2)) <= norm(conv_vec.at(1))) {
				for (int i = -1; i <= 1; i += 2) {
					bb = change_cor(a.at(0), 0, 0, i);
					cout << "(" << bb.x << " " << bb.y << " " << bb.z << ")" << endl;
				}
			}


			return;
		}
		
		cout << "Distance for first - neibouring particle for atom"<<i<<" is : " << min_r<<endl;
		cout << "input distance cutoff" << endl; double cutoff; cin >> cutoff;
		cout << "The coordinate are" << endl<<"atom"<<i << " (" << a.at(i).x << " " << a.at(i).y << " " << a.at(i).z << ")" << endl;
		for (auto& k : list) {
			control = a.at(k);
			cout << "atom" << k << "(" << a.at(k).x << " " << a.at(k).y << " " << a.at(k).z << ")" << endl;
			frac_coor(control);
			if (vect.x == 1 && vect.y == 1 && vect.z == 1 && dist(a.at(k), a.at(i)) == dist(change_cor(a.at(k), -1, -1, -1),a.at(i))) {
				cout << "atom" << k << "(" << -1*a.at(k).x << " " << -1*a.at(k).y << " " << -1*a.at(k).z << ")" << endl;
			}
			if (vect.x == 1 && dist(a.at(k), a.at(i)) == dist(change_cor(a.at(k),-1,0,0),a.at(i))) {
				cout << "atom" << k << "(" <<   a.at(k).x- conv_vec.at(0).x << " " << a.at(k).y - conv_vec.at(0).y << " " << a.at(k).z - conv_vec.at(0).z << ")" << endl;
			}
			if (vect.x == 1 &&vect.y==1&& dist(a.at(k), a.at(i)) == dist(change_cor(a.at(k), -1, -1, 0), a.at(i))) {
				cout << "atom" << k << "(" << a.at(k).x - conv_vec.at(0).x - conv_vec.at(1).x << " " <<  a.at(k).y - conv_vec.at(0).y - conv_vec.at(1).y << " " << a.at(k).z - conv_vec.at(0).z - conv_vec.at(1).z << ")" << endl;
			}
			if (vect.x == 1&&vect.z && dist(a.at(k), a.at(i)) == dist(change_cor(a.at(k), -1, 0, -1), a.at(i))) {
				cout << "atom" << k << "(" << a.at(k).x - conv_vec.at(0).x - conv_vec.at(2).x << " " <<  a.at(k).y - conv_vec.at(0).y - conv_vec.at(2).y << " " << a.at(k).z - conv_vec.at(0).z - conv_vec.at(2).z << ")" << endl;
			}
			if (vect.y == 1 && vect.z==1&& dist(a.at(k), a.at(i)) == dist(change_cor(a.at(k), 0, -1, -1), a.at(i))) {
				cout << "atom" << k << "(" << a.at(k).x - conv_vec.at(1).x - conv_vec.at(2).x << " " <<  a.at(k).y - conv_vec.at(1).y - conv_vec.at(2).y << " " << a.at(k).z - conv_vec.at(1).z - conv_vec.at(2).z << ")" << endl;
			}
			if (vect.y == 1 && dist(a.at(k), a.at(i)) == dist(change_cor(a.at(k), 0, -1, 0), a.at(i))) {
				cout << "atom" << k << "(" << a.at(k).x - conv_vec.at(1).x << " " << a.at(k).y - conv_vec.at(1).y << " " << a.at(k).z - conv_vec.at(1).z << ")" << endl;
			}
			if (vect.z == 1 && dist(a.at(k), a.at(i)) == dist(change_cor(a.at(k), 0, 0, -1), a.at(i))) {
				cout << "atom" << k << "(" << a.at(k).x - conv_vec.at(2).x << " " << a.at(k).y - conv_vec.at(2).y << " " << a.at(k).z - conv_vec.at(2).z << ")" << endl;
			}

		}
		
	}
	

	double3 change_cor(double3 a,int x=0, int y=0 , int z=0) {
		frac_coor(a);
		a = { a.x + x,a.y + y,a.z + z };
		return Todouble3(con_vec_mat * ToVector(a));
	}

	double min_dist() {
		vector<double> min;
		for (auto& i : conv_vec) {
			min.push_back(norm(i));
		}
		return  *min_element(min.begin(), min.end());
	}

	double LJ_potential(double energy, double r) {
		vector<double3> image;
		double potential = 0;
		copy_image(image);
		for (auto& i : atom) {
			for (auto& j : image) {
				potential += 0.5*potetential(i, j,energy,r);
			}
		}
		return potential;
	}

	double potetential(double3 a, double3 b,double energy, double r) {
		if (a==b) {
			return 0;
		}
		return (4 * energy * ((pow((r / dist(a, b)), 12)) - pow((r / dist(a, b)), 6) ) );
	};

	void copy_image(vector<double3> &image){
		for (auto& a : atom) {
			for (int i = -2; i <= 2; i++) {
				for (int j = -2; j <= 2; j++) {
					for (int k = -2; k <= 2; k++) {
						image.push_back({a.x + i * vect.x * conv_vec.at(0).x + j * vect.y * conv_vec.at(1).x + k * vect.z * conv_vec.at(2).x ,
										 a.y + i * vect.x * conv_vec.at(0).y + j * vect.y * conv_vec.at(1).y + k * vect.z * conv_vec.at(2).y ,
										 a.z + i * vect.x * conv_vec.at(0).z + j * vect.y * conv_vec.at(1).z + k * vect.z * conv_vec.at(2).z });
					}
				}
			}
		}
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
	
	double Tersoff(double R=3.0,double D=0.2,double A=3264.7,double B=95.373,double lamda_1=3.2394,double lamda_2=1.3258,double lamda_3=1.3258, double h=-0.0, double c=4.8381, double d=2.0417,double gamma=0.33675, double n=22.956) {
		vector<double3> image;
		double energy = 0;
		int num_atom = atom.size();

		for (int i = 0; i < atom.size(); i++) {
			cout << "calculateing energy of atom " << i << endl;
			evaluate_image(image, i);
			for (int z = 0; z < num_atom; z++) {
				if (vect.x == 1) {
					image.push_back(change_cor(image.at(z) , 1,0,0));
					image.push_back(change_cor(image.at(z), -1, 0, 0));
				}
				if (vect.y == 1) {
					image.push_back(change_cor(image.at(z), 0, 1, 0));
					image.push_back(change_cor(image.at(z), 0, -1, 0));
				}
				if (vect.z == 1) {
					image.push_back(change_cor(image.at(z), 0, 0, 1));
					image.push_back(change_cor(image.at(z), 0, 0, -1));
				}
				if (vect.x == 1&&vect.y==1) {
					image.push_back(change_cor(image.at(z), 1, 1, 0));
					image.push_back(change_cor(image.at(z), -1, -1, 0));
					image.push_back(change_cor(image.at(z), -1, 1, 0));
					image.push_back(change_cor(image.at(z), 1, -1, 0));
				}
				if (vect.x == 1 && vect.z == 1) {
					image.push_back(change_cor(image.at(z), 1, 0, 1));
					image.push_back(change_cor(image.at(z), -1, 0, -1));
					image.push_back(change_cor(image.at(z), -1, 0, 1));
					image.push_back(change_cor(image.at(z), 1, 0, -1));
				}
				if (vect.y == 1 && vect.z == 1) {
					image.push_back(change_cor(image.at(z), 0, 1, 1));
					image.push_back(change_cor(image.at(z), 0, -1, -1));
					image.push_back(change_cor(image.at(z), 0, 1, -1));
					image.push_back(change_cor(image.at(z), 0, -1, 1));
				}
				if (vect.x == 1 && vect.y == 1 && vect.z==1) {
					for (int aa = -1; aa <= 1; aa+=2) {
						for (int bb = -1; bb <= 1; bb+=2) {
							for (int cc = -1; cc <= 1; cc+=2) {
								image.push_back(change_cor(image.at(z), aa, bb, cc));
							}
						}
					}
				}
			}
				for (auto& k : image) {
					if (atom.at(i) == k) {
						continue;
					}
					else if (dist(atom.at(i), k) > 2.4) {
						continue;
					}
					else {

						energy += 0.5 * f_c(dist(atom.at(i), k), R, D) * (V_R(dist(atom.at(i), k), A, lamda_1) - beta(R, D, c, d, h, gamma, n, lamda_3, atom.at(i), k, image) * V_A(dist(k, atom.at(i)), B, lamda_2));
					}
				}
				image = {};
		}
		return energy;
		
	}
	double f_c(double r,double R,double D) {
		if (r <= R - D) {
			return 1;
		}else if (r > R + D) {
			return 0;
		}
		else{
			return (0.5 - 0.5 * (sin(pi * (r - R) / (2 * D))));
		}

	}

	double V_R(double r,double A,double lamda) {
		return (A * exp(-1 * lamda * r));
	}
	double V_A(double r, double B, double lamda) {
		return (B * exp(-1 * lamda * r));
	}

	double beta(double R, double D, double c, double d, double h,double gamma,double n,double lamda,double3 a,double3 b,  vector<double3> image) {
		return pow((1 + pow(gamma, n) * pow(eta(R,D,c,d,h,lamda, a, b, image), n)),-0.5/n);
	}
	double eta(double R, double D, double c, double d, double h, double lamda,double3 a, double3 b, vector<double3> image) {
		double eeta = 0;
		for (auto& i : image) {
			if (b == i||a==i) {
				continue;
			}else if (dist(i, a) > R + D) {
				continue;
			}
			else {
				eeta += f_c(dist(i, a), R, D) * g(c,d,h,a,b,i)*E(lamda, dist(a,b),dist(a,i));
			}
		}
		return lamda;
	}
	double g(double c, double d,double h,double3 o,double3 i,double3 j) {
		return (1 + pow((c / d), 2) - (pow(c, 2) / (pow(d, 2) + pow(h - cos_vector(o, i, j), 2))));
	}

	double cos_vector(double3 o, double3 i, double3 j) {
		return (dot(i - o, j - o)/(norm(i-o)*norm(j-o)));

	};
	double E(double lamda, double a, double b) {
		return exp(pow(lamda*(a-b),3));
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
	Cubic* cell_1;
	BCC* cell_2;
	FCC* cell_3;
	Diamond* diamond;
	Sim_Prim* cell_4;
	Atom_Othor* cell_5;
	string atom_name;
	lattice_vector latticee;
	double well;
	double min_r;

home:
	cout << "This is Project A in PHYS_4061" << endl << endl;
	cout << setprecision(5);
	cout << "Select the Lab number to start the simulation (e.g. input 1 will go to Lab_1 ), or 0 to end the program" << endl;
	cout << "1: Lab_1" << endl << "2: Lab_2" << endl<<"3: Lab_3" <<endl<< "0: exit" << endl;
	cin >> n;
	cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
	switch (n) {
	case 1:
		goto Lab_1;
		break;
	case 2:
		goto Lab_2;
		break;
	case 3:
		goto Lab_3;
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
	cout << "finished printing xyz file for the Simple Cubic, FCC, BCC and diamond structures" << endl;
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
	cout << "it's cooredinate after applying periodic boundary condition in \( x,y,z \) is: " << endl;
	cell_4->PBC(vectorr);
	cout << "(" << vectorr.x << " , " << vectorr.y << " , " << vectorr.z << ")" << endl;
	cell_4->frac_coor(vectorr);
	cout << "The respective frational coordinate is: ( in range of (0.5,0.5] )" << endl;
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
		cout << "Input the fractional coordinate of atom " << i + 1 << endl;
		basis.push_back({});
		cin >> basis.at(i).x >> basis.at(i).y >> basis.at(i).z;
	}

	cout << "input atom name" << endl;
	cin >> atom_name;

	cell_5 = new Atom_Othor(x, y, z, latticee, atom_name, basis);
	cout << "Which atom number you would like to evaluate? There are total " << cell_5->atom_number() << " atoms" << endl;
	cin >> n;
	cout << "You are evaluating atom" << n << endl;
	cell_5->print_atom_position(n);
	cell_5->find_neighbor_list(n);
	delete cell_5;
	basis = {};
	cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
	goto Lab_2;

Lab_3:

	cout << "Lab 3: " << endl;
	cout << "Select the task number" << endl;
	cout << "1: LJ potential" << endl;
	cout << "2: Tersoff energy" << endl;
	cout << "0: home" << endl;
	cin >> task;
	cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
	switch (task) {
	case 1:
		goto Task_3_1;
		break;
	case 2:
		goto Task_3_2;
		break;
	case 0:
		goto home;
		break;
	};
Task_3_1:
	cout << "Lj potential" << endl;
	cout << "Input number of period(in xyz direction), primitive lattice vector , basis of atom and parameters of lj potential to obtain potential" << endl<<endl;
	cout << "input number of period in x y z direction" << endl;
	cin >> x >> y >> z;
	cout << "input a_1 in (x,y,z) input (2.7346 2.7346 0) for silicon" << endl;
	cin >> latticee.a1.x >> latticee.a1.y >> latticee.a1.z;
	cout << "input a_2 in (x,y,z) input (2.7346 0 2.7346) for silicon" << endl;
	cin >> latticee.a2.x >> latticee.a2.y >> latticee.a2.z;
	cout << "input a_3 in (x,y,z) input (0 2.7346 2.7346) for silicon" << endl;
	cin >> latticee.a3.x >> latticee.a3.y >> latticee.a3.z;
	if (dot(latticee.a1, cross(latticee.a2, latticee.a3)) == 0) {
		cout << "No coplanar sets of vectors! Plaese eneter new sets of vectors" << endl;
		goto Re2;
	}
	cout << "Basis setup: " << endl;
	cout << "input number of atom in the basis (input 2 for silicon; basis set= {(0 0 0), (0.25 0.25 0.25)})" << endl;
	cin >> Num;

	for (int i = 0; i < Num; i++) {
		cout << "Input the fractional coordinate of atom " << i + 1 << endl;
		basis.push_back({});
		cin >> basis.at(i).x >> basis.at(i).y >> basis.at(i).z;
	}

	cout << "input atom name" << endl;
	cin >> atom_name;

Ree:
	cout << "input the potential well for LJ potential (For Si, input 1.51)" << endl;
	cin >> well;

	cout << "distance at zero energy (For Si, input 2.1099)" << endl;
	cin >> min_r;


	cell_5 = new Atom_Othor(x, y, z, latticee, atom_name, basis);

	cout << "The energy given by LJ potential is"<<endl;
	cout<<cell_5->LJ_potential(well, min_r)/(cell_5->size()) <<"eV per atom"<<endl;
	delete cell_5;
	goto Ree;

	cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
	goto Lab_3;


Task_3_2:
	cout << "tesorff energy of Si" << endl;
	cout << "input lattice constant (For Si, input 5.423) " << endl;
	cin >> Latconst;
	cout << "input number of period in x y z direction" << endl;
	cin >> x >> y >> z;
	diamond = new Diamond(x, y, z, Latconst);

	cout << "The Tersoff energy is : " << endl;
	cout << diamond->Tersoff()/(diamond->size()+pow(2,x-1)*pow(2,y-1)*pow(2,z-1)) << "eV per atom" << endl; //considering bonding, therefore per atom energy should add some numnber to the atom number

	delete diamond;

	cout << "-----------------------------------------------------------------------------------------------" << endl << endl;
	goto Lab_3;

exit:

	return 0;

};