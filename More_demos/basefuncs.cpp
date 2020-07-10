#include "basefuncs.h"

vector<indices> get_base_table() {
	vector<indices> list;
	list.push_back(indices{ 2,2,1 });//N1
	list.push_back(indices{ 1,2,1 });//N2
	list.push_back(indices{ 1,1,3 });//N3
	list.push_back(indices{ 3,1,3 });//N4
	list.push_back(indices{ 4,4,1 });//N5
	list.push_back(indices{ 1,4,1 });//N6
	list.push_back(indices{ 3,3,2 });//N7
	list.push_back(indices{ 2,3,2 });//N8
	list.push_back(indices{ 2,2,4 });//N9
	list.push_back(indices{ 4,2,4 });//N10
	list.push_back(indices{ 4,4,3 });//N11
	list.push_back(indices{ 3,4,3 });//N12

	list.push_back(indices{ 2,4,3 });//N13
	list.push_back(indices{ 3,2,4 });//N14
	list.push_back(indices{ 3,4,1 });//N15
	list.push_back(indices{ 4,1,3 });//N16
	list.push_back(indices{ 4,2,1 });//N17
	list.push_back(indices{ 1,2,4 });//N18
	list.push_back(indices{ 1,3,2 });//N19
	list.push_back(indices{ 2,1,3 });//N20

	return list;
}

double tetraVolume(V3d& p1, V3d& p2, V3d& p3, V3d& p4) {
	double sum = 0.0;
	sum += determinant3d(p2.x, p3.x, p4.x, p2.y, p3.y, p4.y, p2.z, p3.z, p4.z);
	sum += -1.0 * p1.x*determinant3d(1.0, 1.0, 1.0, p2.y, p3.y, p4.y, p2.z, p3.z, p4.z);
	sum += p1.y*determinant3d(1.0, 1.0, 1.0, p2.x, p3.x, p4.x, p2.z, p3.z, p4.z);
	sum += -1.0*p1.z*determinant3d(1.0, 1.0, 1.0, p2.x, p3.x, p4.x, p2.y, p3.y, p4.y);
	return abs(sum) / 6.0;
}

double L1(V3d& p, V3d& p1, V3d& p2, V3d& p3, V3d& p4) {
	return tetraVolume(p, p2, p3, p4) / tetraVolume(p1, p2, p3, p4);
}

double L2(V3d& p, V3d& p1, V3d& p2, V3d& p3, V3d& p4) {
	return tetraVolume(p1, p, p3, p4) / tetraVolume(p1, p2, p3, p4);
}

double L3(V3d& p, V3d& p1, V3d& p2, V3d& p3, V3d& p4) {
	return tetraVolume(p1, p2, p, p4) / tetraVolume(p1, p2, p3, p4);
}

double L4(V3d& p, V3d& p1, V3d& p2, V3d& p3, V3d& p4) {
	return tetraVolume(p1, p2, p3, p) / tetraVolume(p1, p2, p3, p4);
}

V3d gradL1(V3d& p1, V3d& p2, V3d& p3, V3d& p4) {
	double gradx = -determinant3d(1.0, 1.0, 1.0, p2.y, p3.y, p4.y, p2.z, p3.z, p4.z) / (6.0*tetraVolume(p1, p2, p3, p4));
	double grady = determinant3d(1.0, 1.0, 1.0, p2.x, p3.x, p4.x, p2.z, p3.z, p4.z) / (6.0*tetraVolume(p1, p2, p3, p4));
	double gradz = -determinant3d(1.0, 1.0, 1.0, p2.x, p3.x, p4.x, p2.y, p3.y, p4.y) / (6.0*tetraVolume(p1, p2, p3, p4));
	V3d grad{ gradx,grady,gradz };
	if (grad*(p1 - p2) < 0) {
		grad = grad*(-1.0);
	}
	return grad;
}

V3d gradL2(V3d& p01, V3d& p02, V3d& p03, V3d& p04) {
	V3d p1 = p02, p2 = p03, p3 = p04, p4 = p01;
	double gradx = -determinant3d(1.0, 1.0, 1.0, p2.y, p3.y, p4.y, p2.z, p3.z, p4.z) / (6.0*tetraVolume(p1, p2, p3, p4));
	double grady = determinant3d(1.0, 1.0, 1.0, p2.x, p3.x, p4.x, p2.z, p3.z, p4.z) / (6.0*tetraVolume(p1, p2, p3, p4));
	double gradz = -determinant3d(1.0, 1.0, 1.0, p2.x, p3.x, p4.x, p2.y, p3.y, p4.y) / (6.0*tetraVolume(p1, p2, p3, p4));
	V3d grad{ gradx,grady,gradz };
	if (grad*(p1 - p2) < 0) {
		grad = grad*(-1.0);
	}
	return grad;
}

V3d gradL3(V3d& p01, V3d& p02, V3d& p03, V3d& p04) {
	V3d p1 = p03, p2 = p04, p3 = p01, p4 = p02;
	double gradx = -determinant3d(1.0, 1.0, 1.0, p2.y, p3.y, p4.y, p2.z, p3.z, p4.z) / (6.0*tetraVolume(p1, p2, p3, p4));
	double grady = determinant3d(1.0, 1.0, 1.0, p2.x, p3.x, p4.x, p2.z, p3.z, p4.z) / (6.0*tetraVolume(p1, p2, p3, p4));
	double gradz = -determinant3d(1.0, 1.0, 1.0, p2.x, p3.x, p4.x, p2.y, p3.y, p4.y) / (6.0*tetraVolume(p1, p2, p3, p4));
	V3d grad{ gradx,grady,gradz };
	if (grad*(p1 - p2) < 0) {
		grad = grad*(-1.0);
	}
	return grad;
}

V3d gradL4(V3d& p01, V3d& p02, V3d& p03, V3d& p04) {
	V3d p1 = p04, p2 = p01, p3 = p02, p4 = p03;
	double gradx = -determinant3d(1.0, 1.0, 1.0, p2.y, p3.y, p4.y, p2.z, p3.z, p4.z) / (6.0*tetraVolume(p1, p2, p3, p4));
	double grady = determinant3d(1.0, 1.0, 1.0, p2.x, p3.x, p4.x, p2.z, p3.z, p4.z) / (6.0*tetraVolume(p1, p2, p3, p4));
	double gradz = -determinant3d(1.0, 1.0, 1.0, p2.x, p3.x, p4.x, p2.y, p3.y, p4.y) / (6.0*tetraVolume(p1, p2, p3, p4));
	V3d grad{ gradx,grady,gradz };
	if (grad*(p1 - p2) < 0) {
		grad = grad*(-1.0);
	}
	return grad;
}

V3d Ni(int base_id,V3d p, V3d p1, V3d p2, V3d p3, V3d p4,vector<indices>& index_list) {
	int i = index_list[base_id - 1].i;
	int j = index_list[base_id - 1].j;
	int k = index_list[base_id - 1].k;
	double Li[5] = { 0.0,L1(p, p1, p2, p3, p4),L2(p, p1, p2, p3, p4),L3(p, p1, p2, p3, p4),L4(p, p1, p2, p3, p4) };
	V3d gradLi[5] = { V3d{ 0.0,0.0,0.0 },gradL1(p1, p2, p3, p4),gradL2(p1, p2, p3, p4),gradL3(p1, p2, p3, p4),gradL4(p1, p2, p3, p4) };
	V3d pi[5] = { p,p1,p2,p3,p4 };

	if (base_id-1 < EDGEBASENUM) {
		return (gradLi[k] * Li[j] - gradLi[j] * Li[k])*distanceV3d(pi[j], pi[k])*(3 * Li[i] - 1.0);
	}
	else {
		return (gradLi[k] * Li[j] - gradLi[j] * Li[k])*distanceV3d(pi[j], pi[k])*4.5*Li[i];
	}
}

V3d curlNi(int base_id, V3d p, V3d p1, V3d p2, V3d p3, V3d p4,vector<indices>& index_list) {
	int i = index_list[base_id - 1].i;
	int j = index_list[base_id - 1].j;
	int k = index_list[base_id - 1].k;
	double Li[5] = { 0.0,L1(p, p1, p2, p3, p4),L2(p, p1, p2, p3, p4),L3(p, p1, p2, p3, p4),L4(p, p1, p2, p3, p4) };
	V3d gradLi[5] = { V3d{ 0.0,0.0,0.0 },gradL1(p1, p2, p3, p4),gradL2(p1, p2, p3, p4),gradL3(p1, p2, p3, p4),gradL4(p1, p2, p3, p4) };
	V3d pi[5] = { p,p1,p2,p3,p4 };

	if (base_id - 1 < EDGEBASENUM) {
		return V3dcross(gradLi[j], gradLi[k])*distanceV3d(pi[j], pi[k])*(9 * Li[i] - 2.0);
	}
	else {
		return (V3dcross(gradLi[i], gradLi[k])*Li[j] + V3dcross(gradLi[j], gradLi[i])*Li[k] + V3dcross(gradLi[j], gradLi[k])*2.0*Li[i])*4.5*distanceV3d(pi[j], pi[k]);
	}
}

