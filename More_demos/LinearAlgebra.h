#pragma once
#include <iostream>
using namespace std;
struct V3d
{
	double x, y, z;

	V3d operator +(V3d V) {
		return V3d{ this->x + V.x,this->y + V.y,this->z + V.z };
	}
	V3d operator -(V3d V) {
		return V3d{ this->x - V.x,this->y - V.y,this->z - V.z };
	}
	V3d operator *(double k) {
		return V3d{ k*this->x,k*this->y ,k*this->z };
	}
	V3d operator /(double k) {
		return V3d{ this->x/k,this->y/k ,this->z/k };
	}
	double operator *(V3d V) {
		return this->x*V.x + this->y*V.y + this->z*V.z;
	}
};

V3d V3dcross(V3d& v1, V3d& v2);
void printV3d(V3d v);
double determinant3d(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9);
double norm(V3d v);
double distanceV3d(V3d& v1, V3d& v2);
V3d matrix3d_solve(V3d v1, V3d v2, V3d v3, V3d vp);
double point_to_plane(V3d p0, V3d p1, V3d p2, V3d p3);