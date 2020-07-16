#pragma once
#include "basefuncs.h"
#include "SparseVector.h"
#include "BasicDomain.h"
#include <algorithm>
#include<vector>
using namespace std;

namespace Tetra_system {
	struct Point
	{
		double x, y, z;
	};

	typedef pair<long, long> point_pair;

	typedef pair <point_pair, point_pair> triangle_point_pair;

	struct Edge {
		long p1, p2;
		vector<long> tetral_list;
	};

	struct Triangle {
		long p1, p2, p3;
		vector<long> tetra_list;
	};

	struct Tetrahedral {
		long *vertices;
		long *edges;
		long *triangles;
		Tetrahedral() {
			this->vertices = new long[4];
			this->edges = new long[6];
			this->triangles = new long[4];
		}
	};

	class Tetra_Geometry {
	public:
		vector<Point> plist;
		vector<Edge> edgelist;
		vector<Triangle> trianglelist;
		vector<Tetrahedral> tetrahedral_list;
		
	};

	Tetra_Geometry build_geo_tetra(string nodefile, string tetrafile);

	typedef long* tetra_e_list;

	struct e_coff_set {
		long *h_ids;
		double h_coff;
	};

	struct Eparameters {
		V3d eposition;
		V3d edirection;
		e_coff_set ecoffs;
		long edge_facet_id;
	};

	struct Hparameters {
		V3d hposition;
		V3d hdirection;
		long tetra_id;
	};

	class Tetra_DataSpace {
	public:
		long Ne, Nh;
		vector<indices> base_index_table;
		vector<Eparameters> Etable;
		vector<Hparameters> Htable;
		vector<tetra_e_list> tetra_e_table;
	};

	Tetra_DataSpace geo_to_dataspace_tetra(Tetra_Geometry& geo);


	V3d get_standard_direction_table(int i, V3d p1, V3d p2, V3d p3, V3d p4);

	void redirect(double& e1, double& e2, V3d v1_old, V3d v2_old, V3d v1, V3d v2);

	double tetra_min_distance(V3d p0, long id, Tetra_Geometry& geo);


	class Tetra_domain:public RealDomain {
	public:
		Tetra_Geometry geo;
		Tetra_DataSpace dataspace;
		vector<V3d> tetra_center_list;
		vector<double> tetra_search_radius_list;

		bool check_if_ith_tetra_e_has_extern_surrounding_h(long eid) {
			for (int i = 0; i < 4; i++) {
				if (this->dataspace.Htable[this->dataspace.Etable[eid].ecoffs.h_ids[i]].tetra_id == -1) {
					return true;
				}
			}
			return false;
		}

		Parameter get_Parameter_from_Etable_or_Htable(long n, Field_type ft) {
			Parameter para;
			if (ft == E_field) {
				para.position = this->dataspace.Etable[n].eposition;
				para.direction = this->dataspace.Etable[n].edirection;
				para.field_type = E_field;
				/*
				if (this->check_if_ith_tetra_e_has_extern_surrounding_h(n)) {
					para.shader = BOUNDARY;
				}
				else {
					para.shader = NOT_BOUNDARY;
				}
				*/
				//set all the Es to be NOT_BOUNDARY
				para.shader = NOT_BOUNDARY;
				return para;
			}
			else {
				para.position = this->dataspace.Htable[n].hposition;
				para.direction = this->dataspace.Htable[n].hdirection;
				para.field_type = H_field;
				if (this->dataspace.Htable[n].tetra_id == -1) {
					para.shader = BOUNDARY;
				}
				else {
					/*
					if (this->Eparas[floor(n / 4.0)].shader == BOUNDARY) {
						para.shader = BOUNDARY;
					}
					else {
						para.shader = NOT_BOUNDARY;
					}
					*/
					para.shader = NOT_BOUNDARY;
				}
				return para;
			}
		}

		SparseVec get_Se_row_from_ith_tetra_element(long tetra_id, V3d position, V3d direction, Tetra_system::Tetra_Geometry& geo, Tetra_system::Tetra_DataSpace& ds) {
			SparseVec sv;
			V3d p1 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[0]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[0]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[0]].z };
			V3d p2 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[1]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[1]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[1]].z };
			V3d p3 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[2]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[2]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[2]].z };
			V3d p4 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[3]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[3]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[3]].z };
			for (long para_id = 0; para_id < BASENUM; para_id++) {
				sv.add(ds.tetra_e_table[tetra_id][para_id], NAN);
				//set evalues
				double evalues[BASENUM];
				for (int i = 0; i < BASENUM; i++) {
					evalues[i] = 0.0;
				}
				evalues[para_id] = 1.0;
				//get local e value
				for (int i = 0; i < EDGEBASENUM; i++) {
					evalues[i] = ds.Etable[ds.tetra_e_table[tetra_id][i]].edirection*Tetra_system::get_standard_direction_table(i + 1, p1, p2, p3, p4)*evalues[i];
				}
				for (int i = EDGEBASENUM; i < BASENUM; i = i + 2) {
					V3d v1_old = ds.Etable[ds.tetra_e_table[tetra_id][i]].edirection, v2_old = ds.Etable[ds.tetra_e_table[tetra_id][i + 1]].edirection;
					V3d v1 = Tetra_system::get_standard_direction_table(i + 1, p1, p2, p3, p4), v2 = Tetra_system::get_standard_direction_table(i + 2, p1, p2, p3, p4);
					Tetra_system::redirect(evalues[i], evalues[i + 1], v1_old, v2_old, v1, v2);
				}
				for (int i = 0; i < 4; i++) {
					long eid_1 = ds.tetra_e_table[tetra_id][EDGEBASENUM + 2 * i];
					long eid_2 = ds.tetra_e_table[tetra_id][EDGEBASENUM + 2 * i + 1];
					double a = Tetra_system::get_standard_direction_table(12 + 2 * i + 1 + 1, p1, p2, p3, p4)*Ni(12 + 2 * i + 1, ds.Etable[eid_2].eposition, p1, p2, p3, p4, ds.base_index_table);
					double b = Tetra_system::get_standard_direction_table(12 + 2 * i + 1, p1, p2, p3, p4)*Ni(12 + 2 * i + 1 + 1, ds.Etable[eid_1].eposition, p1, p2, p3, p4, ds.base_index_table);
					double e1 = evalues[EDGEBASENUM + 2 * i], e2 = evalues[EDGEBASENUM + 2 * i + 1];
					evalues[EDGEBASENUM + 2 * i] = (e1 - b*e2) / (1 - a*b);
					evalues[EDGEBASENUM + 2 * i + 1] = (e2 - a*e1) / (1 - a*b);
				}
				V3d curlE{ 0,0,0 };
				for (int i = 0; i < BASENUM; i++) {
					curlE = curlE + curlNi(i + 1, position, p1, p2, p3, p4, ds.base_index_table) * evalues[i];
				}
				sv.terms[para_id].value = curlE*direction;
			}
			return sv;
		}

		V3d get_tetra_center(long id) {
			V3d p1 = { geo.plist[geo.tetrahedral_list[id].vertices[0]].x,geo.plist[geo.tetrahedral_list[id].vertices[0]].y,geo.plist[geo.tetrahedral_list[id].vertices[0]].z };
			V3d p2 = { geo.plist[geo.tetrahedral_list[id].vertices[1]].x,geo.plist[geo.tetrahedral_list[id].vertices[1]].y,geo.plist[geo.tetrahedral_list[id].vertices[1]].z };
			V3d p3 = { geo.plist[geo.tetrahedral_list[id].vertices[2]].x,geo.plist[geo.tetrahedral_list[id].vertices[2]].y,geo.plist[geo.tetrahedral_list[id].vertices[2]].z };
			V3d p4 = { geo.plist[geo.tetrahedral_list[id].vertices[3]].x,geo.plist[geo.tetrahedral_list[id].vertices[3]].y,geo.plist[geo.tetrahedral_list[id].vertices[3]].z };
			return (p1 + p2 + p3 + p4) / 4.0;
		}

		double get_tetra_search_radius(long id, V3d base_point) {
			V3d p1 = { geo.plist[geo.tetrahedral_list[id].vertices[0]].x,geo.plist[geo.tetrahedral_list[id].vertices[0]].y,geo.plist[geo.tetrahedral_list[id].vertices[0]].z };
			V3d p2 = { geo.plist[geo.tetrahedral_list[id].vertices[1]].x,geo.plist[geo.tetrahedral_list[id].vertices[1]].y,geo.plist[geo.tetrahedral_list[id].vertices[1]].z };
			V3d p3 = { geo.plist[geo.tetrahedral_list[id].vertices[2]].x,geo.plist[geo.tetrahedral_list[id].vertices[2]].y,geo.plist[geo.tetrahedral_list[id].vertices[2]].z };
			V3d p4 = { geo.plist[geo.tetrahedral_list[id].vertices[3]].x,geo.plist[geo.tetrahedral_list[id].vertices[3]].y,geo.plist[geo.tetrahedral_list[id].vertices[3]].z };
			double d1 = norm(p1 - base_point);
			double d2 = norm(p2 - base_point);
			double d3 = norm(p3 - base_point);
			double d4 = norm(p4 - base_point);
			double value1 = d1 > d2 ? d1 : d2;
			double value2 = d3 > d4 ? d3 : d4;
			return value1 > value2 ? value1 : value2;
		}

		Tetra_domain() {}

		Tetra_domain(std::string nodefile, std::string tetrafile) {
			cout << "Build Tetra Domain..." << endl;
			this->geo=build_geo_tetra(nodefile, tetrafile);
			this->dataspace = geo_to_dataspace_tetra(this->geo);
			this->Ne = this->dataspace.Ne;
			this->Nh = this->dataspace.Nh;
			for (long i = 0; i < this->Ne; i++) {
				this->Eparas.push_back(this->get_Parameter_from_Etable_or_Htable(i, E_field));
			}
			for (long i = 0; i < this->Nh; i++) {
				this->Hparas.push_back(this->get_Parameter_from_Etable_or_Htable(i, H_field));
			}
			//set Sh
			this->Sh = vector<SparseVec>(this->Ne);
			for (long i = 0; i < this->Ne; i++) {
				//actually all the Es are NOT_BOUNDARY
				if (this->Eparas[i].shader == NOT_BOUNDARY) {
					SparseVec sv;
					for (long j = 0; j < 4; j++) {
						long surrounding_old_h_id = this->dataspace.Etable[i].ecoffs.h_ids[j];
						double h_coeff = this->dataspace.Etable[i].ecoffs.h_coff;
						sv.add(surrounding_old_h_id, h_coeff);
					}
					this->Sh[i] = sv;
				}
			}
			//set Se
			this->Se = vector<SparseVec>(this->Nh);
#pragma omp parallel for schedule(dynamic)
			for (long i = 0; i < this->Nh; i++) {
				if (this->Hparas[i].shader == NOT_BOUNDARY) {
					this->Se[i] = this->get_Se_row_from_ith_tetra_element(this->dataspace.Htable[i].tetra_id, this->Hparas[i].position, this->Hparas[i].direction, this->geo, this->dataspace);
				}
			}
			//set tetra_center_list & tetra_search_radius_list
			for (long i = 0; i < this->geo.tetrahedral_list.size(); i++) {
				V3d tetra_center = this->get_tetra_center(i);
				this->tetra_center_list.push_back(tetra_center);
				this->tetra_search_radius_list.push_back(this->get_tetra_search_radius(i, tetra_center));
			}
			cout << "Tetra Domain Finished!: Ne=" << this->Ne << " Nh=" << this->Nh << endl;
		}

		void output_mesh(string filename) {
			long N_edge = this->geo.edgelist.size();
			ofstream of;
			of.open(filename);
			for (long i = 0; i < N_edge; i++) {
				V3d p1 = V3d{ this->geo.plist[this->geo.edgelist[i].p1].x,this->geo.plist[this->geo.edgelist[i].p1].y,this->geo.plist[this->geo.edgelist[i].p1].z };
				V3d p2 = V3d{ this->geo.plist[this->geo.edgelist[i].p2].x,this->geo.plist[this->geo.edgelist[i].p2].y,this->geo.plist[this->geo.edgelist[i].p2].z };
				of << p1.x << " " << p1.y << " " << p1.z << " " << p2.x << " " << p2.y << " " << p2.z << endl;
			}
			of.close();
		}

		bool check_if_point_on_segment(V3d p0, V3d p1, V3d p2) {
			double bar = 1e-15;
			double ratio = (norm(p0 - p1) + norm(p0 - p2)) / norm(p1 - p2);
			if (abs(ratio - 1.0) < bar ) {
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
			double ratio = (area1 + area2 + area3) / area;
			double bar = 1e-15;
			if (abs(ratio - 1.0) < bar ) {
				return true;
			}
			else {
				return false;
			}
		}

		bool check_if_p_on_tetra_surface(V3d p0, V3d p1, V3d p2, V3d p3, V3d p4) {
			bool cond1 = this->check_if_point_on_segment(p0, p1, p2);
			bool cond2 = this->check_if_point_on_segment(p0, p1, p3);
			bool cond3 = this->check_if_point_on_segment(p0, p1, p4);
			bool cond4 = this->check_if_point_on_segment(p0, p2, p3);
			bool cond5 = this->check_if_point_on_segment(p0, p2, p4);
			bool cond6 = this->check_if_point_on_segment(p0, p3, p4);
			bool cond7 = this->check_if_point_on_triangle(p0, p2, p3, p4);
			bool cond8 = this->check_if_point_on_triangle(p0, p1, p3, p4);
			bool cond9 = this->check_if_point_on_triangle(p0, p1, p2, p4);
			bool cond10 = this->check_if_point_on_triangle(p0, p1, p2, p3);
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
			/*
			if (this->check_if_p_on_tetra_surface(p0, p1, p2, p3, p4)) {
				return true;
			}
			*/
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

		long get_tetra_fall_into(V3d p0, vector<long>& tetralist, Tetra_system::Tetra_Geometry& geo) {
			for (int i = 0; i < tetralist.size(); i++) {
				if (this->check_if_p_falls_into_tetra(p0, tetralist[i], geo)) {
					return tetralist[i];
				}
			}
			return -1;
			cout << "No tetra to fall into" << endl;
		}

		long search_which_tetra_fall_into(V3d position) {
			vector<long> tetralist;
			for (long i = 0; i < this->geo.tetrahedral_list.size(); i++) {
				if (norm(position - this->tetra_center_list[i]) < 1.5*this->tetra_search_radius_list[i]) {
					tetralist.push_back(i);
				}
			}
			return this->get_tetra_fall_into(position, tetralist, this->geo);
		}

		double bisection_for_partition(V3d start, V3d end, double tleft0, double tright0, long tetra_id_left, Tetra_system::Tetra_Geometry& geo) {
			double bar = 1e-15;
			double tleft = tleft0, tright = tright0;
			while ((tright - tleft) > bar) {
				double tmid = (tleft + tright) / 2.0;
				if (this->check_if_p_falls_into_tetra(start + (end - start)*tmid, tetra_id_left, geo)) {
					tleft = tmid;
				}
				else {
					tright = tmid;
				}
			}
			return tright;
		}

		bool check_if_p_on_and_along_edge(V3d p0, V3d dir, V3d p1, V3d p2) {
			bool on_edge = this->check_if_point_on_segment(p0, p1, p2);
			bool along_edge = abs(abs(((p1 - p2) / norm(p1 - p2)) * dir) - 1.0) < 1e-15;
			return on_edge && along_edge;
		}

		SparseVec refine_sv_on_and_along_edge(V3d p0, V3d dir, SparseVec sv, V3d p1, V3d p2) {
			SparseVec temp;
			for (long i = 0; i < sv.terms.size(); i++) {
				V3d e_position = this->Eparas[sv.terms[i].index].position;
				if (this->check_if_point_on_segment(e_position, p1, p2)) {
					temp.add(sv.terms[i].index, sv.terms[i].value);
				}
			}
			if (temp.terms.size() != 2) {
				cout << "refine_sv_on_edge: Error!" << endl;
			}
			V3d direction = this->Eparas[temp.terms[0].index].direction;
			V3d pos1 = this->Eparas[temp.terms[0].index].position;
			V3d pos2 = this->Eparas[temp.terms[1].index].position;
			double a = ((p0 - pos1) * direction) / ((pos2 - pos1) * direction);
			SparseVec new_sv;
			new_sv.add(temp.terms[0].index, 1 - a);
			new_sv.add(temp.terms[1].index, a);
			new_sv = (direction * dir) * new_sv;
			return new_sv;
		}

		SparseVec get_eSe_row_from_ith_tetra_element(long tetra_id, V3d position, V3d direction, Tetra_system::Tetra_Geometry& geo, Tetra_system::Tetra_DataSpace& ds) {
			SparseVec sv;
			V3d p1 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[0]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[0]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[0]].z };
			V3d p2 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[1]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[1]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[1]].z };
			V3d p3 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[2]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[2]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[2]].z };
			V3d p4 = { geo.plist[geo.tetrahedral_list[tetra_id].vertices[3]].x,geo.plist[geo.tetrahedral_list[tetra_id].vertices[3]].y,geo.plist[geo.tetrahedral_list[tetra_id].vertices[3]].z };
			for (long para_id = 0; para_id < BASENUM; para_id++) {
				sv.add(ds.tetra_e_table[tetra_id][para_id], NAN);
				//set evalues
				double evalues[BASENUM];
				for (int i = 0; i < BASENUM; i++) {
					evalues[i] = 0.0;
				}
				evalues[para_id] = 1.0;
				//get local e value
				for (int i = 0; i < EDGEBASENUM; i++) {
					evalues[i] = ds.Etable[ds.tetra_e_table[tetra_id][i]].edirection*Tetra_system::get_standard_direction_table(i + 1, p1, p2, p3, p4)*evalues[i];
				}
				for (int i = EDGEBASENUM; i < BASENUM; i = i + 2) {
					V3d v1_old = ds.Etable[ds.tetra_e_table[tetra_id][i]].edirection, v2_old = ds.Etable[ds.tetra_e_table[tetra_id][i + 1]].edirection;
					V3d v1 = Tetra_system::get_standard_direction_table(i + 1, p1, p2, p3, p4), v2 = Tetra_system::get_standard_direction_table(i + 2, p1, p2, p3, p4);
					Tetra_system::redirect(evalues[i], evalues[i + 1], v1_old, v2_old, v1, v2);
				}
				for (int i = 0; i < 4; i++) {
					long eid_1 = ds.tetra_e_table[tetra_id][EDGEBASENUM + 2 * i];
					long eid_2 = ds.tetra_e_table[tetra_id][EDGEBASENUM + 2 * i + 1];
					double a = Tetra_system::get_standard_direction_table(12 + 2 * i + 1 + 1, p1, p2, p3, p4)*Ni(12 + 2 * i + 1, ds.Etable[eid_2].eposition, p1, p2, p3, p4, ds.base_index_table);
					double b = Tetra_system::get_standard_direction_table(12 + 2 * i + 1, p1, p2, p3, p4)*Ni(12 + 2 * i + 1 + 1, ds.Etable[eid_1].eposition, p1, p2, p3, p4, ds.base_index_table);
					double e1 = evalues[EDGEBASENUM + 2 * i], e2 = evalues[EDGEBASENUM + 2 * i + 1];
					evalues[EDGEBASENUM + 2 * i] = (e1 - b*e2) / (1 - a*b);
					evalues[EDGEBASENUM + 2 * i + 1] = (e2 - a*e1) / (1 - a*b);
				}
				V3d E{ 0,0,0 };
				for (int i = 0; i < BASENUM; i++) {
					E = E + Ni(i + 1, position, p1, p2, p3, p4, ds.base_index_table) * evalues[i];
				}
				sv.terms[para_id].value = E*direction;
			}
			/*
			if (this->check_if_p_on_and_along_edge(position, direction, p1, p2)) {
				return this->refine_sv_on_and_along_edge(position, direction, sv, p1, p2);
			}
			if (this->check_if_p_on_and_along_edge(position, direction, p1, p3)) {
				return this->refine_sv_on_and_along_edge(position, direction, sv, p1, p3);
			}
			if (this->check_if_p_on_and_along_edge(position, direction, p1, p4)) {
				return this->refine_sv_on_and_along_edge(position, direction, sv, p1, p4);
			}
			if (this->check_if_p_on_and_along_edge(position, direction, p2, p3)) {
				return this->refine_sv_on_and_along_edge(position, direction, sv, p2, p3);
			}
			if (this->check_if_p_on_and_along_edge(position, direction, p2, p4)) {
				return this->refine_sv_on_and_along_edge(position, direction, sv, p2, p4);
			}
			if (this->check_if_p_on_and_along_edge(position, direction, p3, p4)) {
				return this->refine_sv_on_and_along_edge(position, direction, sv, p3, p4);
			}
			*/
			return sv;
		}

		SparseVec numerical_integration_Se_row(long tetra_id, V3d start, V3d end, V3d direction, Tetra_system::Tetra_Geometry& geo, Tetra_system::Tetra_DataSpace& ds) {
			SparseVec sv;
			V3d mid_point = (start + end)*0.5;
			sv = sv + this->get_eSe_row_from_ith_tetra_element(tetra_id, start, direction, geo, ds);
			sv = sv + 4.0*this->get_eSe_row_from_ith_tetra_element(tetra_id, mid_point, direction, geo, ds);
			sv = sv + this->get_eSe_row_from_ith_tetra_element(tetra_id, end, direction, geo, ds);
			sv = (1.0 / 6.0)*sv;
			return sv;
		}

		SparseVec make_weights_finer(V3d start, V3d end, SparseVec coarse_sv, double coarse_step, Tetra_system::Tetra_Geometry& geo) {
			if (coarse_sv.terms.size() == 1) {
				return coarse_sv;
			}
			SparseVec sv;
			for (long i = 0; i < (coarse_sv.terms.size() - 1); i++) {
				long tetra_id = coarse_sv.terms[i].index;
				double new_t = this->bisection_for_partition(start, end, coarse_sv.terms[i].value - coarse_step, coarse_sv.terms[i].value, tetra_id, geo);
				sv.add(coarse_sv.terms[i].index, new_t);
			}
			sv.add(coarse_sv.terms[coarse_sv.terms.size() - 1].index, 1.0);
			return sv;
		}

		SparseVec get_weights_for_interface_brick_E(V3d start, V3d end, Tetra_system::Tetra_Geometry& geo) {
			SparseVec sv;
			long N = 1e3;
			double step = 1.0 / N;
			long current_tetra = search_which_tetra_fall_into(start);
			long init = 1;
			while (current_tetra == -1) {
				current_tetra = search_which_tetra_fall_into(start + (end - start)*(init*step));
				init++;
				if (init > 5) {
					break;
				}
			}
			for (long i = init; i <= N; i++) {
				V3d p_now = start + (end - start)*(i*step);
				if (check_if_p_falls_into_tetra(p_now, current_tetra, geo) == false) {
					/*
					long next_tetra = search_which_tetra_fall_into(p_now);
					sv.add(current_tetra, i * step);
					current_tetra = next_tetra;
					*/
					
					long next_tetra= search_which_tetra_fall_into(p_now);
					long skipN = 0;
					while (next_tetra == -1 && i <= N) {
						i++;
						skipN++;
						next_tetra = search_which_tetra_fall_into(start + (end - start) * (i * step));
					}
					if (next_tetra != -1 && next_tetra != current_tetra) {
						sv.add(current_tetra, i * step);
						current_tetra = next_tetra;
					}
					if (skipN > 5) {
						cout << "Warning: big skip in getting eSe expansion! N="<<skipN << endl;
					}
					
				}
			}
			if (current_tetra != -1) {
				sv.add(current_tetra, 1.0);
			}
			sv = this->make_weights_finer(start, end, sv, step, geo);
			return sv;
		}

		SparseVec get_eSe_row_from_tetra_elements(V3d start, V3d end, V3d direction) {
			SparseVec weights = this->get_weights_for_interface_brick_E(start, end, this->geo);
			SparseVec eSe_row;
			for (long i = 0; i < weights.terms.size(); i++) {
				long tetra_id = weights.terms[i].index;
				V3d start_p, end_p, mid_point;
				double coeff;
				if (i == 0) {
					coeff = weights.terms[i].value;
					start_p = start;
					end_p = start + (end - start)*coeff;
					mid_point = start + (end - start)*coeff*0.5;
				}
				else {
					coeff = weights.terms[i].value - weights.terms[i - 1].value;
					start_p = start + (end - start)*weights.terms[i - 1].value;
					end_p = start + (end - start)*weights.terms[i].value;
					mid_point = start + (end - start)*(weights.terms[i].value + weights.terms[i - 1].value)*0.5;
				}
				eSe_row = eSe_row + coeff*this->numerical_integration_Se_row(tetra_id, start_p, end_p, direction, this->geo, this->dataspace);
			}
			return eSe_row;
		}
	};
	
	
	
}
