#include "Tetra_system.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#define BAR -1e-6
#define LOOP_SIZE 1.5

using namespace std;

namespace Tetra_system {

	//////////////////////////     Geometry related code     ////////////////////////////////////

	bool two_sort(long& n1, long& n2) {
		long temp;
		if (n1 > n2) {
			temp = n1;
			n1 = n2;
			n2 = temp;
		}
		return true;
	}

	bool three_sort(long& n1, long& n2, long& n3) {
		two_sort(n1, n2);
		two_sort(n2, n3);
		two_sort(n1, n2);
		return true;
	}

	bool load_nodes(string nodefile, vector<Point>& plist) {
		ifstream nodestream(nodefile);
		string str;
		while (getline(nodestream, str)) {
			stringstream ss(str);
			double x, y, z;
			ss >> x;
			ss >> y;
			ss >> z;
			plist.push_back(Point{ x,y,z });
		}
		nodestream.close();
		return true;
	}

	long lookup_edge(long n1, long n2, long tetra_id, vector<Edge>& edgelist, map<point_pair, long>& dic) {
		two_sort(n1, n2);
		long index = -1;
		point_pair pair = make_pair(n1, n2);
		if (dic.count(pair)) {
			index = dic[pair];
			edgelist[index].tetral_list.push_back(tetra_id);
			return index;
		}
		else {
			Edge ed; ed.p1 = n1; ed.p2 = n2; ed.tetral_list.push_back(tetra_id);
			edgelist.push_back(ed);
			index = edgelist.size() - 1;
			dic[pair] = index;
			return index;
		}
	}

	long lookup_facet(long n1, long n2, long n3, long tetra_id, vector<Triangle>& trianglelist, map<triangle_point_pair, long>& dic_facet) {
		three_sort(n1, n2, n3);
		long index = -1;
		point_pair pair1 = make_pair(n1, n2), pair2 = make_pair(n2, n3);
		triangle_point_pair pair = make_pair(pair1, pair2);
		if (dic_facet.count(pair)) {
			index = dic_facet[pair];
			trianglelist[index].tetra_list.push_back(tetra_id);
			return index;
		}
		else {
			Triangle tri; tri.p1 = n1; tri.p2 = n2; tri.p3 = n3;  tri.tetra_list.push_back(tetra_id);
			trianglelist.push_back(tri);
			index = trianglelist.size() - 1;
			dic_facet[pair] = index;
			return index;
		}
	}

	bool assign_edge_facet_tetra(long p1, long p2, long p3, long p4, long tetraid, vector<Edge>& edgelist, vector<Triangle>& trianglelist, map<point_pair, long>& dic_edge, map<triangle_point_pair, long>& dic_facet, vector<Tetrahedral>& tetrahedral_list) {
		Tetrahedral tetra = Tetrahedral();
		tetra.vertices[0] = p1; tetra.vertices[1] = p2; tetra.vertices[2] = p3; tetra.vertices[3] = p4;
		//lookup edges
		tetra.edges[0] = lookup_edge(p1, p2, tetraid, edgelist, dic_edge);
		tetra.edges[1] = lookup_edge(p1, p3, tetraid, edgelist, dic_edge);
		tetra.edges[2] = lookup_edge(p1, p4, tetraid, edgelist, dic_edge);
		tetra.edges[3] = lookup_edge(p2, p3, tetraid, edgelist, dic_edge);
		tetra.edges[4] = lookup_edge(p2, p4, tetraid, edgelist, dic_edge);
		tetra.edges[5] = lookup_edge(p3, p4, tetraid, edgelist, dic_edge);
		//lookup facets
		tetra.triangles[0] = lookup_facet(p2, p3, p4, tetraid, trianglelist, dic_facet);
		tetra.triangles[1] = lookup_facet(p1, p3, p4, tetraid, trianglelist, dic_facet);
		tetra.triangles[2] = lookup_facet(p1, p2, p4, tetraid, trianglelist, dic_facet);
		tetra.triangles[3] = lookup_facet(p1, p2, p3, tetraid, trianglelist, dic_facet);
		//finish tetra
		tetrahedral_list.push_back(tetra);
		return true;
	}

	bool load_edge_facet_tetra(string tetrafile, vector<Edge>& edgelist, vector<Triangle>& trianglelist, vector<Tetrahedral>& tetrahedral_list) {
		ifstream filestream(tetrafile);
		string str;
		map<point_pair, long> dic_edge;
		map<triangle_point_pair, long> dic_facet;
		while (getline(filestream, str)) {
			stringstream ss(str);
			long p1, p2, p3, p4, tetra_id;
			ss >> p1; ss >> p2; ss >> p3; ss >> p4; ss >> tetra_id;
			p1--; p2--; p3--; p4--; tetra_id--;
			assign_edge_facet_tetra(p1, p2, p3, p4, tetra_id, edgelist, trianglelist, dic_edge, dic_facet, tetrahedral_list);
		}
		filestream.close();
		return true;
	}

	Tetra_Geometry build_geo_tetra(string nodefile, string tetrafile) {
		Tetra_Geometry geo;
		load_nodes(nodefile, geo.plist);
		load_edge_facet_tetra(tetrafile, geo.edgelist, geo.trianglelist, geo.tetrahedral_list);
		cout << " Tetra geometry finished: ";
		cout << "N_edge=" << geo.edgelist.size() << " N_tetrahedron=" << geo.tetrahedral_list.size() << endl;
		return geo;
	}
}



namespace Tetra_system {
	///////////////////////////      DataSpace related code      /////////////////////////////////

	V3d get_standard_direction_table(int i, V3d p1, V3d p2, V3d p3, V3d p4) {
		if (i == 1 || i == 2 || i == 17) {
			return (p2 - p1) / norm(p2 - p1)*(-1);
		}
		if (i == 3 || i == 4 || i == 16 || i == 20) {
			return (p1 - p3) / norm(p1 - p3)*(-1);
		}
		if (i == 5 || i == 6 || i == 15) {
			return (p4 - p1) / norm(p4 - p1)*(-1);
		}
		if (i == 7 || i == 8 || i == 19) {
			return (p3 - p2) / norm(p3 - p2)*(-1);
		}
		if (i == 9 || i == 10 || i == 14 || i == 18) {
			return (p2 - p4) / norm(p2 - p4)*(-1);
		}
		if (i == 11 || i == 12 || i == 13) {
			return (p4 - p3) / norm(p4 - p3)*(-1);
		}
		cout << "wrong direction id!" << endl;
	}

	void redirect(double& e1, double& e2, V3d v1_old, V3d v2_old, V3d v1, V3d v2) {
		V3d v3old = V3dcross(v1_old, v2_old); v3old = v3old / norm(v3old);
		V3d v01{ v1_old.x,v2_old.x,v3old.x }, v02{ v1_old.y,v2_old.y,v3old.y }, v03{ v1_old.z,v2_old.z,v3old.z }, vp{ e1,e2,0.0 };
		V3d E = matrix3d_solve(v01, v02, v03, vp);
		e1 = v1*E; e2 = v2*E;
	}

	double tetra_min_distance(V3d p0, long id, Tetra_Geometry& geo) {
		V3d p1 = { geo.plist[geo.tetrahedral_list[id].vertices[0]].x,geo.plist[geo.tetrahedral_list[id].vertices[0]].y,geo.plist[geo.tetrahedral_list[id].vertices[0]].z };
		V3d p2 = { geo.plist[geo.tetrahedral_list[id].vertices[1]].x,geo.plist[geo.tetrahedral_list[id].vertices[1]].y,geo.plist[geo.tetrahedral_list[id].vertices[1]].z };
		V3d p3 = { geo.plist[geo.tetrahedral_list[id].vertices[2]].x,geo.plist[geo.tetrahedral_list[id].vertices[2]].y,geo.plist[geo.tetrahedral_list[id].vertices[2]].z };
		V3d p4 = { geo.plist[geo.tetrahedral_list[id].vertices[3]].x,geo.plist[geo.tetrahedral_list[id].vertices[3]].y,geo.plist[geo.tetrahedral_list[id].vertices[3]].z };
		double average_length = (distanceV3d(p1, p2) + distanceV3d(p1, p3) + distanceV3d(p1, p4) + distanceV3d(p2, p3) + distanceV3d(p2, p4) + distanceV3d(p3, p4)) / 6.0;
		vector<double> raw_values, values;
		raw_values.push_back(point_to_plane(p0, p1, p2, p3));
		raw_values.push_back(point_to_plane(p0, p1, p2, p4));
		raw_values.push_back(point_to_plane(p0, p1, p3, p4));
		raw_values.push_back(point_to_plane(p0, p2, p3, p4));

		for (int i = 0; i < raw_values.size(); i++) {
			if (raw_values[i] > average_length*1e-6) {
				values.push_back(raw_values[i]);
			}
		}
		return *min_element(values.begin(), values.end());
	}

	double get_length(long id, V3d position, Tetra_Geometry& geo) {
		vector<double> values;
		if (id < geo.edgelist.size()) {
			for (int i = 0; i < geo.edgelist[id].tetral_list.size(); i++) {
				long m = geo.edgelist[id].tetral_list[i];
				values.push_back(tetra_min_distance(position, m, geo));
			}
		}
		else {
			id = id - geo.edgelist.size();
			for (int i = 0; i<geo.trianglelist[id].tetra_list.size(); i++) {
				long m = geo.trianglelist[id].tetra_list[i];
				values.push_back(tetra_min_distance(position, m, geo));
			}
		}
		return *min_element(values.begin(), values.end());
	}

	bool check_if_point_on_segment(V3d p0, V3d p1, V3d p2) {
		double bar = 1e-15;
		double dif = (norm(p0 - p1) + norm(p0 - p2)) - norm(p1 - p2);
		if (abs(dif) < bar * norm(p1 - p2)) {
			return true;
		}
		else {
			return false;
		}
	}

	bool check_if_point_on_triangle(V3d p0, V3d p1, V3d p2, V3d p3) {
		V3d v01 = p2 - p1;
		V3d v02 = p3 - p1;
		V3d v1 = p1 - p0;
		V3d v2 = p2 - p0;
		V3d v3 = p3 - p0;
		double area = 0.5 * norm(V3dcross(v01, v02));
		double area1 = 0.5 * norm(V3dcross(v1, v2));
		double area2 = 0.5 * norm(V3dcross(v1, v3));
		double area3 = 0.5 * norm(V3dcross(v2, v3));
		double dif = (area1 + area2 + area3) - area;
		double bar = 1e-15;
		if (abs(dif) < bar * area) {
			return true;
		}
		else {
			return false;
		}
	}

	bool check_if_p_on_tetra_surface(V3d p0, V3d p1, V3d p2, V3d p3, V3d p4) {
		bool cond1 = check_if_point_on_segment(p0, p1, p2);
		bool cond2 = check_if_point_on_segment(p0, p1, p3);
		bool cond3 = check_if_point_on_segment(p0, p1, p4);
		bool cond4 = check_if_point_on_segment(p0, p2, p3);
		bool cond5 = check_if_point_on_segment(p0, p2, p4);
		bool cond6 = check_if_point_on_segment(p0, p3, p4);
		bool cond7 = check_if_point_on_triangle(p0, p2, p3, p4);
		bool cond8 = check_if_point_on_triangle(p0, p1, p3, p4);
		bool cond9 = check_if_point_on_triangle(p0, p1, p2, p4);
		bool cond10 = check_if_point_on_triangle(p0, p1, p2, p3);
		bool on_segments = cond1 || cond2 || cond3 || cond4 || cond5 || cond6;
		bool on_surface = cond7 || cond8 || cond9 || cond10;
		if (on_segments || on_surface) {
			return true;
		}
		else {
			return false;
		}
	}

	bool check_if_p_falls_into_tetra(V3d p0, long tetra_id, Tetra_system::Tetra_Geometry& geo) {
		V3d p1 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[0]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[0]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[0]].z };
		V3d p2 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[1]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[1]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[1]].z };
		V3d p3 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[2]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[2]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[2]].z };
		V3d p4 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[3]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[3]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[3]].z };
		if (check_if_p_on_tetra_surface(p0, p1, p2, p3, p4)) {
			return true;
		}
		double bar = 1e-15;
		V3d v1 = p2 - p1, v2 = p3 - p1, v3 = p4 - p1, vp = p0 - p1;
		V3d vx = matrix3d_solve(v1, v2, v3, vp);
		double sum = vx.x + vx.y + vx.z;
		if (vx.x >= -bar && vx.y >= -bar && vx.z >= -bar && vx.x + vx.y + vx.z <= 1.0 + bar) {
			return true;
		}
		else {
			return false;
		}
	}

	long get_tetra_id_for_h(V3d p0, vector<long>& tetralist, Tetra_system::Tetra_Geometry& geo) {
		for (int i = 0; i < tetralist.size(); i++) {
			if (check_if_p_falls_into_tetra(p0, tetralist[i], geo)) {
				return tetralist[i];
			}
		}
		return -1;
	}
	/*
	long get_tetra_id_for_h(V3d p0, vector<long>& tetralist, Tetra_Geometry& geo) {
		for (int i = 0; i < tetralist.size(); i++) {
			V3d p1 = { geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[0]].x,geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[0]].y,geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[0]].z };
			V3d p2 = { geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[1]].x,geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[1]].y,geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[1]].z };
			V3d p3 = { geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[2]].x,geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[2]].y,geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[2]].z };
			V3d p4 = { geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[3]].x,geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[3]].y,geo.plist[geo.tetrahedral_list[tetralist[i]].vertices[3]].z };

			V3d v1 = p2 - p1, v2 = p3 - p1, v3 = p4 - p1, vp = p0 - p1;
			V3d vx = matrix3d_solve(v1, v2, v3, vp);
			double sum = vx.x + vx.y + vx.z;
			if (vx.x >= BAR && vx.y >= BAR && vx.z >= BAR && vx.x + vx.y + vx.z <= 1.0 - BAR) {
				//cout << "Find one tetrahedron matches h position!" << endl;
				return tetralist[i];
			}
		}
		if (p0.x >= 0.0 && p0.x <= 1.0 && p0.y >= 0.0 && p0.y <= 0.5 && p0.z >= 0.0 && p0.z <= 0.75) {
			//cout << "Error! No tetrahedron matches h position!" << endl;
			//cout << p0.x << "," << p0.y << "," << p0.z << endl;
		}
		else {
			//cout << "PASS" << endl;
		}
		//cout << p0.x << "," << p0.y << "," << p0.z << endl;
		return -1;
	}
	*/
	e_coff_set set_h_for_e(V3d position, V3d direction, double L, vector<long>& tetralist, Tetra_DataSpace& ds, Tetra_Geometry& geo,bool flag=false, V3d default_dir= { 0,0,1.0 }) {
		V3d vz = { 0,0,1.0 };
		if (!(abs(norm(V3dcross(direction, vz))) > 5e-16)) {
			vz.x += 1.0;
		}
		vz = vz / norm(vz);
		//V3d v1 = direction - vz;
		V3d v1 = vz;
		if (flag) {
			v1 = default_dir;
		}
		v1 = v1 / norm(v1);
		V3d v2 = V3dcross(direction, v1);
		v2 = v2 / norm(v2);
		v1 = V3dcross(v2, direction);

		e_coff_set set;
		set.h_ids = new long[4];
		set.h_ids[0] = ds.Htable.size();
		set.h_ids[1] = ds.Htable.size() + 1;
		set.h_ids[2] = ds.Htable.size() + 2;
		set.h_ids[3] = ds.Htable.size() + 3;
		double k = LOOP_SIZE;  //control loop size
		set.h_coff = 1 / (k*L);
		//fill H parameters
		Hparameters hpara1, hpara2, hpara3, hpara4;
		hpara1.hposition = position + v1*k*0.5*L; hpara1.hdirection = v2;
		hpara1.tetra_id = get_tetra_id_for_h(hpara1.hposition, tetralist, geo);
		hpara2.hposition = position + v2*k*0.5*L; hpara2.hdirection = v1*(-1.0);
		hpara2.tetra_id = get_tetra_id_for_h(hpara2.hposition, tetralist, geo);
		hpara3.hposition = position - v1*k*0.5*L; hpara3.hdirection = v2*(-1.0);
		hpara3.tetra_id = get_tetra_id_for_h(hpara3.hposition, tetralist, geo);
		hpara4.hposition = position - v2*k*0.5*L; hpara4.hdirection = v1;
		hpara4.tetra_id = get_tetra_id_for_h(hpara4.hposition, tetralist, geo);
		ds.Htable.push_back(hpara1);
		ds.Htable.push_back(hpara2);
		ds.Htable.push_back(hpara3);
		ds.Htable.push_back(hpara4);

		return set;
	}

	long get_e_id(long base_id, long tetra_id, Tetra_Geometry& geo) {
		long standard_edge_data[12] = { 2,1,1,3,4,1,3,2,2,4,4,3 };
		long tetra_plist[5];
		tetra_plist[0] = -1;
		for (int i = 1; i < 5; i++) {
			tetra_plist[i] = geo.tetrahedral_list[tetra_id].vertices[i - 1];
		}
		if (base_id < EDGEBASENUM) {
			long edge_id = floor(base_id / 2.0);
			if (tetra_plist[standard_edge_data[2 * edge_id]] > tetra_plist[standard_edge_data[2 * edge_id + 1]]) {
				return 2 * geo.tetrahedral_list[tetra_id].edges[edge_id] + base_id - 2 * edge_id;
			}
			else {
				if (base_id - 2 * edge_id > 0) {
					return 2 * geo.tetrahedral_list[tetra_id].edges[edge_id];
				}
				else {
					return 2 * geo.tetrahedral_list[tetra_id].edges[edge_id] + 1;
				}
			}
		}
		else {
			long face_id = floor((base_id - EDGEBASENUM) / 2.0);
			return 2 * geo.edgelist.size() + 2 * geo.tetrahedral_list[tetra_id].triangles[face_id] + base_id - EDGEBASENUM - 2 * face_id;
		}
	}

	Tetra_DataSpace geo_to_dataspace_tetra(Tetra_Geometry& geo) {
		Tetra_DataSpace dataspace;
		dataspace.base_index_table = get_base_table();
		for (long i = 0; i < geo.edgelist.size(); i++) {
			V3d p1 = { geo.plist[geo.edgelist[i].p1].x,geo.plist[geo.edgelist[i].p1].y,geo.plist[geo.edgelist[i].p1].z };
			V3d p2 = { geo.plist[geo.edgelist[i].p2].x,geo.plist[geo.edgelist[i].p2].y,geo.plist[geo.edgelist[i].p2].z };
			Eparameters epara1, epara2;
			epara1.eposition = p2 + (p1 - p2) / 3.0;
			epara2.eposition = p2 + (p1 - p2)*2.0 / 3.0;
			epara1.edirection = (p1 - p2) / norm(p1 - p2);
			epara2.edirection = (p1 - p2) / norm(p1 - p2);
			epara1.ecoffs = set_h_for_e(epara1.eposition, epara1.edirection, get_length(i, epara1.eposition, geo), geo.edgelist[i].tetral_list, dataspace, geo);
			epara2.ecoffs = set_h_for_e(epara2.eposition, epara2.edirection, get_length(i, epara2.eposition, geo), geo.edgelist[i].tetral_list, dataspace, geo);
			epara1.edge_facet_id = i;
			epara2.edge_facet_id = i;
			dataspace.Etable.push_back(epara1);
			dataspace.Etable.push_back(epara2);
		}
		long edge_number = geo.edgelist.size();
		for (long i = 0; i < geo.trianglelist.size(); i++) {
			V3d p1 = { geo.plist[geo.trianglelist[i].p1].x,geo.plist[geo.trianglelist[i].p1].y,geo.plist[geo.trianglelist[i].p1].z };
			V3d p2 = { geo.plist[geo.trianglelist[i].p2].x,geo.plist[geo.trianglelist[i].p2].y,geo.plist[geo.trianglelist[i].p2].z };
			V3d p3 = { geo.plist[geo.trianglelist[i].p3].x,geo.plist[geo.trianglelist[i].p3].y,geo.plist[geo.trianglelist[i].p3].z };
			Eparameters epara1, epara2;
			epara1.eposition = (p1 + p2 + p3) / 3.0;
			epara2.eposition = (p1 + p2 + p3) / 3.0;
			epara1.edirection = (p3 - p2) / norm(p3 - p2);
			epara2.edirection = (p1 - p3) / norm(p1 - p3);
			//epara1.ecoffs = set_h_for_e(epara1.eposition, epara1.edirection, get_length(i + edge_number, epara1.eposition, geo), geo.trianglelist[i].tetra_list, dataspace, geo, true, epara2.edirection);
			//epara2.ecoffs = set_h_for_e(epara2.eposition, epara2.edirection, get_length(i + edge_number, epara2.eposition, geo), geo.trianglelist[i].tetra_list, dataspace, geo, true, epara1.edirection);
			epara1.ecoffs = set_h_for_e(epara1.eposition, epara1.edirection, get_length(i + edge_number, epara1.eposition, geo), geo.trianglelist[i].tetra_list, dataspace, geo);
			epara2.ecoffs = set_h_for_e(epara2.eposition, epara2.edirection, get_length(i + edge_number, epara2.eposition, geo), geo.trianglelist[i].tetra_list, dataspace, geo);
			epara1.edge_facet_id = edge_number + i;
			epara2.edge_facet_id = edge_number + i;
			dataspace.Etable.push_back(epara1);
			dataspace.Etable.push_back(epara2);
		}
		for (long i = 0; i < geo.tetrahedral_list.size(); i++) {
			tetra_e_list elist = new long[BASENUM];
			for (int j = 0; j < BASENUM; j++) {
				elist[j] = get_e_id(j, i, geo);
			}
			dataspace.tetra_e_table.push_back(elist);
		}

		dataspace.Ne = dataspace.Etable.size();
		dataspace.Nh = dataspace.Htable.size();
		cout << "Tetra ComputationSpace finished: ";
		cout << "Ne=" << dataspace.Etable.size() << " " << dataspace.Etable.size() << " Nh=" << dataspace.Htable.size() << " " << dataspace.Htable.size() << endl;
		return dataspace;
	}

}