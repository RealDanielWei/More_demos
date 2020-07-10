#pragma once
#define BASENUM 20
#define EDGEBASENUM 12
#include "LinearAlgebra.h"
#include <vector>
struct indices {
	int i, j, k;
};

double tetraVolume(V3d& p1, V3d& p2, V3d& p3, V3d& p4);
double L1(V3d& p, V3d& p1, V3d& p2, V3d& p3, V3d& p4);
double L2(V3d& p, V3d& p1, V3d& p2, V3d& p3, V3d& p4);
double L3(V3d& p, V3d& p1, V3d& p2, V3d& p3, V3d& p4);
double L4(V3d& p, V3d& p1, V3d& p2, V3d& p3, V3d& p4);

V3d gradL1(V3d& p1, V3d& p2, V3d& p3, V3d& p4);
V3d gradL2(V3d& p1, V3d& p2, V3d& p3, V3d& p4);
V3d gradL3(V3d& p1, V3d& p2, V3d& p3, V3d& p4);
V3d gradL4(V3d& p1, V3d& p2, V3d& p3, V3d& p4);

V3d Ni(int base_id, V3d p, V3d p1, V3d p2, V3d p3, V3d p4, vector<indices>& index_list);

V3d curlNi(int base_id, V3d p, V3d p1, V3d p2, V3d p3, V3d p4, vector<indices>& index_list);

vector<indices> get_base_table();