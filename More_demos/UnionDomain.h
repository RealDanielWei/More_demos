#pragma once
#include "Brick_system.h"
#include "Tetra_system.h"

namespace UnionDomain {

	enum Hitbox_tag
	{
		Outside_Hitbox, On_Hitbox, Inside_Hitbox
	};

	class Hitbox {
	public:
		double xmin, ymin, zmin;
		double xmax, ymax, zmax;
		double bar = 1e-15;

		Hitbox() {
			this->xmin = 0.0; this->ymin = 0.0; this->zmin = 0.0;
			this->xmax = 0.0; this->ymax = 0.0; this->zmax = 0.0;
		};

		double double_max(double a, double b) {
			return a > b ? a : b;
		}

		double double_min(double a, double b) {
			return a < b ? a : b;
		}

		Hitbox(Tetra_system::Tetra_domain& tdomain) {
			long Npoint = tdomain.geo.plist.size();
			double xplus = tdomain.geo.plist[0].x;
			double xminus = tdomain.geo.plist[0].x;
			double yplus = tdomain.geo.plist[0].y;
			double yminus = tdomain.geo.plist[0].y;
			double zplus = tdomain.geo.plist[0].z;
			double zminus = tdomain.geo.plist[0].z;
			for (long i = 0; i < Npoint; i++) {
				xplus = this->double_max(xplus, tdomain.geo.plist[i].x);
				xminus = this->double_min(xminus, tdomain.geo.plist[i].x);
				yplus = this->double_max(yplus, tdomain.geo.plist[i].y);
				yminus = this->double_min(yminus, tdomain.geo.plist[i].y);
				zplus = this->double_max(zplus, tdomain.geo.plist[i].z);
				zminus = this->double_min(zminus, tdomain.geo.plist[i].z);
			}
			this->xmax = xplus;
			this->xmin = xminus;
			this->ymax = yplus;
			this->ymin = yminus;
			this->zmax = zplus;
			this->zmin = zminus;
		}

		Hitbox(Brick_system::Brick_Domain& bdomain) {
			this->xmax = bdomain.setting.xmax;
			this->ymax = bdomain.setting.ymax;
			this->zmax = bdomain.setting.zmax;
			this->xmin = 0.0;
			this->ymin = 0.0;
			this->zmin = 0.0;
		}

		Hitbox_tag hitbox_test(V3d position) {
			bool not_outside_x = position.x > this->xmin - this->bar && position.x < this->xmax + this->bar;
			bool not_outside_y = position.y > this->ymin - this->bar && position.y < this->ymax + this->bar;
			bool not_outside_z = position.z > this->zmin - this->bar && position.z < this->zmax + this->bar;
			bool inside_x = position.x > this->xmin + this->bar && position.x < this->xmax - this->bar;
			bool inside_y = position.y > this->ymin + this->bar && position.y < this->ymax - this->bar;
			bool inside_z = position.z > this->zmin + this->bar && position.z < this->zmax - this->bar;
			if (not_outside_x && not_outside_y && not_outside_z) {
				if (inside_x && inside_y && inside_z) {
					return Inside_Hitbox;
				}
				else {
					return On_Hitbox;
				}
			}
			else {
				return Outside_Hitbox;
			}
		}	

	};

	class DomainInfo {
	public:
		Hitbox hitbox;
		V3d shift;
		vector<Shader> e_shader_list, h_shader_list;
		vector<long> e_ref_id_list, h_ref_id_list;
		vector<long> e_id_old_to_new, e_id_new_to_old;
		vector<long> h_id_old_to_new, h_id_new_to_old;

		DomainInfo(){}
		//get size
		long get_Ne_new() {
			return this->e_id_new_to_old.size();
		}
		long get_Ne_old() {
			return this->e_id_old_to_new.size();
		}
		long get_Nh_new() {
			return this->h_id_new_to_old.size();
		}
		long get_Nh_old() {
			return this->h_id_old_to_new.size();
		}
		//translation
		long e_old_to_new(long old_e_id) {
			return this->e_id_old_to_new[old_e_id];
		}
		long e_new_to_old(long new_e_id) {
			return this->e_id_new_to_old[new_e_id];
		}
		long h_old_to_new(long old_h_id) {
			return this->h_id_old_to_new[old_h_id];
		}
		long h_new_to_old(long new_h_id) {
			return this->h_id_new_to_old[new_h_id];
		}
	};

	class UnionDomain : public BasicDomain {
	public:
		vector<Brick_system::Brick_Domain> bg;
		vector<Tetra_system::Tetra_domain> tdomain_list;
		vector<DomainInfo> domain_info_list;

		UnionDomain(){}

		DomainInfo get_domain_info_for_tdomain(Tetra_system::Tetra_domain& tdomain, V3d shift) {
			DomainInfo info;
			info.hitbox = Hitbox(this->tdomain_list[this->tdomain_list.size() - 1]);
			info.shift = shift;
			for (long i = 0; i < tdomain.get_Ne(); i++) {
				Parameter para = tdomain.get_Eparas(i);
				if (para.shader == NOT_BOUNDARY) {
					info.e_shader_list.push_back(NOT_BOUNDARY);
					info.e_ref_id_list.push_back(-1);
					info.e_id_old_to_new.push_back(info.e_id_new_to_old.size());
					info.e_id_new_to_old.push_back(i);
				}
				else {
					cout << "Union::get_domain_info_for_tdomain: Invalid e shader!" << endl;
				}
			}
			for (long i = 0; i < tdomain.get_Nh(); i++) {
				Parameter para = tdomain.get_Hparas(i);
				if (para.shader == NOT_BOUNDARY) {
					info.h_shader_list.push_back(NOT_BOUNDARY);
					info.h_ref_id_list.push_back(-1);
					info.h_id_old_to_new.push_back(info.h_id_new_to_old.size());
					info.h_id_new_to_old.push_back(i);
				}
				else {
					if (para.shader == BOUNDARY) {
						info.h_shader_list.push_back(INTERFACE);
						info.h_ref_id_list.push_back(-1);
						info.h_id_old_to_new.push_back(-1);
					}
					else {
						cout << "Union::get_domain_info_for_tdomain: Invalid h shader!" << endl;
					}
				}
			}
			return info;
		}

		void add_tdomain(string node_file, string element_file, V3d shift) {
			this->tdomain_list.push_back(Tetra_system::Tetra_domain(node_file, element_file));
			this->domain_info_list.push_back(this->get_domain_info_for_tdomain(this->tdomain_list[this->tdomain_list.size() - 1], shift));
			
		}

		DomainInfo get_domain_info_for_bg(V3d shift) {
			DomainInfo info;
			info.hitbox = Hitbox(this->bg[0]);
			info.shift = shift;
			for (long i = 0; i < this->bg[0].get_Ne(); i++) {
				Parameter para = this->bg[0].Eparas[i];
				Hitbox_tag hit_result = Outside_Hitbox;
				long temp_ref_id = -1;
				//check whether hit or not
				for (long j = 0; j < this->domain_info_list.size(); j++) {
					V3d local_position = para.position + shift - this->domain_info_list[j].shift;
					Hitbox_tag local_hit_result = this->domain_info_list[j].hitbox.hitbox_test(local_position);
					if (local_hit_result != Outside_Hitbox) {
						hit_result = local_hit_result;
						temp_ref_id = j;
						break;
					}
				}
				//set shader & translations according to hit result
				if (hit_result == Outside_Hitbox) {
					info.e_shader_list.push_back(para.shader);
					info.e_id_old_to_new.push_back(info.e_id_new_to_old.size());
					info.e_id_new_to_old.push_back(i);
				}
				else {
					if (hit_result == On_Hitbox) {
						info.e_shader_list.push_back(INTERFACE);
					}
					else {
						info.e_shader_list.push_back(DISABLED);
					}
					info.e_id_old_to_new.push_back(-1);
				}
				//set ref
				info.e_ref_id_list.push_back(temp_ref_id);
			}
			for (long i = 0; i < this->bg[0].get_Nh(); i++) {
				Parameter para = this->bg[0].Hparas[i];
				Hitbox_tag hit_result = Outside_Hitbox;
				long temp_ref_id = -1;
				//check whether hit or not
				for (long j = 0; j < this->domain_info_list.size(); j++) {
					V3d local_position = para.position + shift - this->domain_info_list[j].shift;
					Hitbox_tag local_hit_result = this->domain_info_list[j].hitbox.hitbox_test(local_position);
					if (local_hit_result != Outside_Hitbox) {
						hit_result = local_hit_result;
						temp_ref_id = j;
						break;
					}
				}
				//set shader & translations according to hit result
				if (hit_result == Outside_Hitbox || hit_result == On_Hitbox) {
					info.h_shader_list.push_back(para.shader);
					info.h_id_old_to_new.push_back(info.h_id_new_to_old.size());
					info.h_id_new_to_old.push_back(i);
				}
				else {
					info.h_shader_list.push_back(DISABLED);
					info.h_id_old_to_new.push_back(-1);
				}
				//set ref
				info.h_ref_id_list.push_back(temp_ref_id);
			}
			return info;
		}

		void add_pml_bg(long Nx, long Ny, long Nz, long Npml) {
			Brick_system::brick_setting bsetting;
			bsetting.Nx = Nx + 2 * Npml;
			bsetting.Ny = Ny + 2 * Npml;
			bsetting.Nz = Nz + 2 * Npml;
			double deltax = (this->domain_info_list[0].hitbox.xmax - this->domain_info_list[0].hitbox.xmin) / Nx;
			double deltay = (this->domain_info_list[0].hitbox.ymax - this->domain_info_list[0].hitbox.ymin) / Ny;
			double deltaz = (this->domain_info_list[0].hitbox.zmax - this->domain_info_list[0].hitbox.zmin) / Nz;
			for (long i = 0; i < bsetting.Nx; i++) {
				bsetting.delta_x_list.push_back(deltax);
			}
			for (long i = 0; i < bsetting.Ny; i++) {
				bsetting.delta_y_list.push_back(deltay);
			}
			for (long i = 0; i < bsetting.Nz; i++) {
				bsetting.delta_z_list.push_back(deltaz);
			}
			bsetting.set();
			V3d shift = V3d{ this->domain_info_list[0].hitbox.xmin,this->domain_info_list[0].hitbox.ymin,this->domain_info_list[0].hitbox.zmin };
			shift = shift + this->domain_info_list[0].shift;
			shift = shift - V3d{ Npml * deltax,Npml * deltay,Npml * deltaz };
			this->bg.push_back(Brick_system::Brick_Domain(bsetting));
			this->domain_info_list.push_back(this->get_domain_info_for_bg(shift));
		}

		pair<long, long> eid_global_to_regional(long eid) {
			long region_id = 0;
			long Npre = 0;
			long Ntotal = this->domain_info_list[0].get_Ne_new();
			while (!(Ntotal > eid)) {
				region_id++;
				Npre = Ntotal;
				Ntotal += this->domain_info_list[region_id].get_Ne_new();
			}
			return make_pair(region_id, eid - Npre);
		}

		pair<long, long> hid_global_to_regional(long hid) {
			long region_id = 0;
			long Npre = 0;
			long Ntotal = this->domain_info_list[0].get_Nh_new();
			while (!(Ntotal > hid)) {
				region_id++;
				Npre = Ntotal;
				Ntotal += this->domain_info_list[region_id].get_Nh_new();
			}
			return make_pair(region_id, hid - Npre);
		}

		long eid_regional_to_global(pair<long, long> address_pair) {
			long eid = 0;
			for (long i = 0; i < address_pair.first; i++) {
				eid+= this->domain_info_list[i].get_Ne_new();
			}
			eid += address_pair.second;
			return eid;
		}

		long hid_regional_to_global(pair<long, long> address_pair) {
			long hid = 0;
			for (long i = 0; i < address_pair.first; i++) {
				hid += this->domain_info_list[i].get_Nh_new();
			}
			hid += address_pair.second;
			return hid;
		}

		long get_Ne() {
			long Ne = 0;
			for (long i = 0; i < this->domain_info_list.size(); i++) {
				Ne += this->domain_info_list[i].get_Ne_new();
			}
			return Ne;
		}

		long get_Nh(){
			long Nh = 0;
			for (long i = 0; i < this->domain_info_list.size(); i++) {
				Nh += this->domain_info_list[i].get_Nh_new();
			}
			return Nh;
		}

		Parameter get_Eparas(long i) {
			pair<long, long> regional_address = this->eid_global_to_regional(i);
			long region_id = regional_address.first;
			long regional_new_e_id = regional_address.second;
			long regional_old_e_id = this->domain_info_list[region_id].e_new_to_old(regional_new_e_id);
			Parameter para;
			if (region_id == this->domain_info_list.size() - 1) {
				para = this->bg[0].get_Eparas(regional_old_e_id);
			}
			else {
				para = this->tdomain_list[region_id].get_Eparas(regional_old_e_id);
			}
			para.shader = this->domain_info_list[region_id].e_shader_list[regional_old_e_id];
			para.position = para.position + this->domain_info_list[region_id].shift;
			return para;
		}

		Parameter get_Hparas(long i){
			pair<long, long> regional_address = this->hid_global_to_regional(i);
			long region_id = regional_address.first;
			long regional_new_h_id = regional_address.second;
			long regional_old_h_id = this->domain_info_list[region_id].h_new_to_old(regional_new_h_id);
			Parameter para;
			if (region_id == this->domain_info_list.size() - 1) {
				para = this->bg[0].get_Hparas(regional_old_h_id);
			}
			else {
				para = this->tdomain_list[region_id].get_Hparas(regional_old_h_id);
			}
			para.shader = this->domain_info_list[region_id].h_shader_list[regional_old_h_id];
			para.position = para.position + this->domain_info_list[region_id].shift;
			return para;
		}

		SparseVec e_regional_old_sv_to_global(long region_id, SparseVec raw_sv) {
			SparseVec sv;
			for (long i = 0; i < raw_sv.terms.size(); i++) {
				long old_regional_e_id = raw_sv.terms[i].index;
				long new_regional_e_id = this->domain_info_list[region_id].e_old_to_new(old_regional_e_id);
				long global_e_id = this->eid_regional_to_global(make_pair(region_id, new_regional_e_id));
				double value = raw_sv.terms[i].value;
				sv.add(global_e_id, value);
			}
			return sv;
		}

		SparseVec resolve_bg_interface_e_to_global_sv(long local_e_old_id) {
			long source_region = this->domain_info_list.size() - 1;
			Parameter para = this->bg[0].get_Eparas(local_e_old_id);
			V3d global_pos = para.position + this->domain_info_list[source_region].shift;
			V3d edge_direction = para.direction;
			double edge_length = this->bg[0].get_edge_length(local_e_old_id);
			long target_region = this->domain_info_list[source_region].e_ref_id_list[local_e_old_id];
			V3d regional_position = global_pos - this->domain_info_list[target_region].shift;
			V3d edge_start = regional_position - edge_direction * 0.5 * edge_length;
			V3d edge_end= regional_position + edge_direction * 0.5 * edge_length;
			SparseVec local_sv = this->tdomain_list[target_region].get_eSe_row_from_tetra_elements(edge_start, edge_end, edge_direction);
			return this->e_regional_old_sv_to_global(target_region, local_sv);
		}

		SparseVec h_bg_old_sv_to_global(SparseVec raw_sv) {
			SparseVec sv;
			long region_id = this->domain_info_list.size() - 1;
			for (long i = 0; i < raw_sv.terms.size(); i++) {
				long old_bg_h_id = raw_sv.terms[i].index;
				long new_bg_h_id = this->domain_info_list[region_id].h_old_to_new(old_bg_h_id);
				long global_h_id = this->hid_regional_to_global(make_pair(region_id, new_bg_h_id));
				double value = raw_sv.terms[i].value;
				sv.add(global_h_id, value);
			}
			return sv;
		}

		SparseVec resolve_interface_h_to_bg_global_sv(long source_region, long local_h_old_id) {
			Parameter para = this->tdomain_list[source_region].get_Hparas(local_h_old_id);
			V3d global_pos = para.position + this->domain_info_list[source_region].shift;
			V3d direction = para.direction;
			long target_region = this->domain_info_list.size() - 1;
			V3d regional_position = global_pos - this->domain_info_list[target_region].shift;
			SparseVec local_sv = this->bg[0].get_general_h_expansion(regional_position, direction, this->domain_info_list[target_region].h_shader_list);
			return this->h_bg_old_sv_to_global(local_sv);
		}

		SparseVec get_Sh_row(long eid){
			pair<long, long> regional_address = this->eid_global_to_regional(eid);
			long local_old_eid = this->domain_info_list[regional_address.first].e_new_to_old(regional_address.second);
			SparseVec local_old_sv;
			bool is_bg = false;
			if (regional_address.first == this->domain_info_list.size() - 1) {
				local_old_sv = this->bg[0].get_Sh_row(local_old_eid);
				is_bg = true;
			}
			else {
				local_old_sv = this->tdomain_list[regional_address.first].get_Sh_row(local_old_eid);
			}
			SparseVec sv;
			for (long i = 0; i < local_old_sv.terms.size(); i++) {
				long index = local_old_sv.terms[i].index;
				double value = local_old_sv.terms[i].value;
				Shader shader = this->domain_info_list[regional_address.first].h_shader_list[index];
				if (shader == NOT_BOUNDARY || shader == BOUNDARY) {
					long new_index = this->domain_info_list[regional_address.first].h_old_to_new(index);
					long global_index = this->hid_regional_to_global(make_pair(regional_address.first, new_index));
					sv.add(global_index, value);
				}
				else {
					sv = sv + value * this->resolve_interface_h_to_bg_global_sv(regional_address.first, index);
				}
			}
			return sv;
		}

		SparseVec get_Se_row(long hid){
			pair<long, long> regional_address = this->hid_global_to_regional(hid);
			long local_old_hid = this->domain_info_list[regional_address.first].h_new_to_old(regional_address.second);
			SparseVec local_old_sv;
			bool is_bg = false;
			if (regional_address.first == this->domain_info_list.size() - 1) {
				local_old_sv = this->bg[0].get_Se_row(local_old_hid);
				is_bg = true;
			}
			else {
				local_old_sv = this->tdomain_list[regional_address.first].get_Se_row(local_old_hid);
			}
			SparseVec sv;
			for (long i = 0; i < local_old_sv.terms.size(); i++) {
				long index = local_old_sv.terms[i].index;
				double value = local_old_sv.terms[i].value;
				Shader shader = this->domain_info_list[regional_address.first].e_shader_list[index];
				if (shader == NOT_BOUNDARY || shader == BOUNDARY) {
					long new_index = this->domain_info_list[regional_address.first].e_old_to_new(index);
					long global_index = this->eid_regional_to_global(make_pair(regional_address.first, new_index));
					sv.add(global_index, value);
				}
				else {
					sv = sv + value * this->resolve_bg_interface_e_to_global_sv(index);
				}
			}
			return sv;
		}

	};

	class RealizedUnionDomain :public RealDomain {
	public:
		UnionDomain unidomain;

		RealizedUnionDomain() {}

		void add_tdomain(string node_file, string element_file, V3d shift) {
			this->unidomain.add_tdomain(node_file, element_file, shift);
		}

		void add_pml_bg(long Nx, long Ny, long Nz, long Npml) {
			this->unidomain.add_pml_bg(Nx, Ny, Nz, Npml);
		}

		void realize() {
			this->Ne = this->unidomain.get_Ne();
			this->Nh = this->unidomain.get_Nh();
			this->Eparas = vector<Parameter>(this->Ne);
			this->Hparas = vector<Parameter>(this->Nh);
			this->Sh = vector<SparseVec>(this->Ne);
			this->Se = vector<SparseVec>(this->Nh);
			for (long i = 0; i < this->Ne; i++) {
				this->Eparas[i] = this->unidomain.get_Eparas(i);
				this->Sh[i] = this->unidomain.get_Sh_row(i);
			}
			for (long i = 0; i < this->Nh; i++) {
				this->Hparas[i] = this->unidomain.get_Hparas(i);
				this->Se[i] = this->unidomain.get_Se_row(i);
			}
		}
	};

}