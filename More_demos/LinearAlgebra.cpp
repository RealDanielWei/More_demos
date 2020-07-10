#include "LinearAlgebra.h"

V3d V3dcross(V3d& v1, V3d& v2) {
	return V3d{ v1.y*v2.z - v1.z * v2.y,v1.z*v2.x - v1.x*v2.z,v1.x*v2.y - v1.y*v2.x };
}

void printV3d(V3d v) {
	cout << v.x << "\t\t" << v.y << "\t\t" << v.z << endl;
}

double determinant3d(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9) {
	return a1*a5*a9 + a3*a4*a8 + a2*a6*a7 - a3*a5*a7 - a2*a4*a9 - a1*a6*a8;
}

V3d matrix3d_solve(V3d v1, V3d v2, V3d v3, V3d vp) {
	double x, y, z;
	double determine0 = determinant3d(v1.x, v2.x, v3.x, v1.y, v2.y, v3.y, v1.z, v2.z, v3.z);
	x = determinant3d(vp.x, v2.x, v3.x, vp.y, v2.y, v3.y, vp.z, v2.z, v3.z) / determine0;
	y = determinant3d(v1.x, vp.x, v3.x, v1.y, vp.y, v3.y, v1.z, vp.z, v3.z) / determine0;
	z = determinant3d(v1.x, v2.x, vp.x, v1.y, v2.y, vp.y, v1.z, v2.z, vp.z) / determine0;
	return V3d{ x,y,z };
}

double point_to_plane(V3d p0, V3d p1, V3d p2, V3d p3) {
	V3d v1 = p2 - p1, v2 = p3 - p1, v0 = p0 - p1;
	V3d n = V3dcross(v1, v2);
	n = n / norm(n);
	return abs(v0*n);
}

double norm(V3d v) {
	return sqrt(v*v);
}

double distanceV3d(V3d& v1, V3d& v2) {
	return norm(v1 - v2);
}