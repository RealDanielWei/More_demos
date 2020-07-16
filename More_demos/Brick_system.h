#pragma once
#include <vector>
#include "LinearAlgebra.h"
#include "SparseVector.h"
#include "BasicDomain.h"

using namespace std;

namespace Brick_system {

	enum Direction
	{
		null, x, y, z, xy, yz, zx
	};

	struct spaceindex {
		Direction di;
		long nx, ny, nz;

		spaceindex() {
			this->nx = 0; this->ny = 0; this->nz = 0; this->di = null;
		}
		spaceindex(long ix, long iy, long iz, Direction Di) {
			this->nx = ix; this->ny = iy; this->nz = iz; this->di = Di;
		}
	};

	class brick_setting {
	public:
		long Nx, Ny, Nz;
		double xmax, ymax, zmax;
		vector<double> delta_x_list, delta_y_list, delta_z_list;
		vector<double> x_base_list, y_base_list, z_base_list;

		void set() {
			{
				//set x centers
				this->x_base_list = vector<double>(this->Nx + 1);
				this->x_base_list[0] = 0.0;
				for (long i = 1; i < this->Nx + 1; i++) {
					this->x_base_list[i] = this->x_base_list[i - 1] + this->delta_x_list[i - 1];
				}
			}
			{
				//set y centers
				this->y_base_list = vector<double>(this->Ny + 1);
				this->y_base_list[0] = 0.0;
				for (long i = 1; i < this->Ny + 1; i++) {
					this->y_base_list[i] = this->y_base_list[i - 1] + this->delta_y_list[i - 1];
				}
			}
			{
				//set z centers
				this->z_base_list = vector<double>(this->Nz + 1);
				this->z_base_list[0] = 0.0;
				for (long i = 1; i < this->Nz + 1; i++) {
					this->z_base_list[i] = this->z_base_list[i - 1] + this->delta_z_list[i - 1];
				}
			}
			this->xmax = this->x_base_list[this->Nx];
			this->ymax = this->y_base_list[this->Ny];
			this->zmax = this->z_base_list[this->Nz];
		}

		brick_setting() {

		}

		brick_setting(long Nx, long Ny, long Nz) {
			this->Nx = Nx;
			this->Ny = Ny;
			this->Nz = Nz;
			this->delta_x_list = vector<double>(this->Nx);
			this->delta_y_list = vector<double>(this->Ny);
			this->delta_z_list = vector<double>(this->Nz);
			this->x_base_list = vector<double>(this->Nx + 1);
			this->y_base_list = vector<double>(this->Ny + 1);
			this->z_base_list = vector<double>(this->Nz + 1);
		}

		spaceindex get_which_brick_it_falls_into(V3d position) {
			spaceindex spindex(-1, -1, -1, null);
			if (position.x<0.0 || position.x>this->xmax || position.y<0.0 || position.y>this->ymax || position.z<0.0 || position.z>this->zmax) {
				cout << "brick_setting::get_which_brick_it_falls_into: Invalid position!" << endl;
				return spindex;
			}
			//determine nx
			for (long i = 0; i < this->x_base_list.size() - 1; i++) {
				if (position.x >= this->x_base_list[i] && position.x <= this->x_base_list[i + 1]) {
					spindex.nx = i;
					break;
				}
			}
			//determine ny
			for (long i = 0; i < this->y_base_list.size() - 1; i++) {
				if (position.y >= this->y_base_list[i] && position.y <= this->y_base_list[i + 1]) {
					spindex.ny = i;
					break;
				}
			}
			//determine nz
			for (long i = 0; i < this->z_base_list.size() - 1; i++) {
				if (position.z >= this->z_base_list[i] && position.z <= this->z_base_list[i + 1]) {
					spindex.nz = i;
					break;
				}
			}
			return spindex;
		}

	};

	

	class Brick_Domain:public RealDomain {
	public:
		brick_setting setting;

		V3d brick_dir_to_V3d(Brick_system::Direction dir) {
			if (dir == Brick_system::x || dir == Brick_system::yz) {
				return V3d{ 1.0,0,0 };
			}
			if (dir == Brick_system::y || dir == Brick_system::zx) {
				return V3d{ 0,1.0,0 };
			}
			if (dir == Brick_system::z || dir == Brick_system::xy) {
				return V3d{ 0,0,1.0 };
			}
			cout << "brick_dir_to_V3d: Wrong Direction!" << endl;
		}

		long edgespaceindex_to_edgeindex(Direction di, spaceindex P, long Nx, long Ny, long Nz) {
			long Nxedge = Nx * (Ny + 1)*(Nz + 1);
			long Nyedge = (Nx + 1)*Ny*(Nz + 1);
			long Nzedge = (Nx + 1)*(Ny + 1)*Nz;
			switch (di)
			{
			case x:
				return P.nz*Nx*(Ny + 1) + P.ny*Nx + P.nx;
				break;
			case y:
				return Nxedge + P.nz*(Nx + 1)*Ny + P.ny*(Nx + 1) + P.nx;
				break;
			case z:
				return Nxedge + Nyedge + P.nz*(Nx + 1)*(Ny + 1) + P.ny*(Nx + 1) + P.nx;
				break;
			default:
				printf("invalid edge direction! ");
				return NAN;
				break;
			}

		}

		spaceindex edgeindex_to_edgespaceindex(long n, long Nx, long Ny, long Nz) {
			long Nxedge = Nx * (Ny + 1)*(Nz + 1);
			long Nyedge = (Nx + 1)*Ny*(Nz + 1);
			long Nzedge = (Nx + 1)*(Ny + 1)*Nz;
			Direction di;
			if (n < Nxedge) {
				di = x;
			}
			else {
				if (n < Nxedge + Nyedge) {
					di = y;
				}
				else {
					di = z;
				}
			}

			long ix, iy, iz;
			switch (di)
			{
			case x:
				iz = n / (Nx*(Ny + 1));
				n = n % (Nx*(Ny + 1));
				iy = n / Nx;
				ix = n % Nx;
				break;
			case y:
				n = n - Nxedge;
				iz = n / ((Nx + 1)*Ny);
				n = n % ((Nx + 1)*Ny);
				iy = n / (Nx + 1);
				ix = n % (Nx + 1);
				break;
			case z:
				n = n - Nxedge - Nyedge;
				iz = n / ((Nx + 1)*(Ny + 1));
				n = n % ((Nx + 1)*(Ny + 1));
				iy = n / (Nx + 1);
				ix = n % (Nx + 1);
				break;
			default:
				break;
			}
			return spaceindex(ix, iy, iz, di);
		}

		long pathchspaceindex_to_patchindex(Direction di, spaceindex P, long Nx, long Ny, long Nz) {
			long Nxypatch = Nx * Ny*(Nz + 1);
			long Nyzpatch = (Nx + 1)*Ny*Nz;
			long Nzxpatch = Nx * (Ny + 1)*Nz;

			switch (di)
			{
			case xy:
				return P.nz*Nx*Ny + P.ny*Nx + P.nx;
				break;
			case yz:
				return Nxypatch + P.nz*(Nx + 1)*Ny + P.ny*(Nx + 1) + P.nx;
				break;
			case zx:
				return Nxypatch + Nyzpatch + P.nz*Nx*(Ny + 1) + P.ny*Nx + P.nx;
				break;
			default:
				printf("invalid patch direction\n");
				return NAN;
				break;
			}
		}

		spaceindex patchindex_to_patchspaceidnex(long n, long Nx, long Ny, long Nz) {
			long Nxypatch = Nx * Ny*(Nz + 1);
			long Nyzpatch = (Nx + 1)*Ny*Nz;
			long Nzxpatch = Nx * (Ny + 1)*Nz;
			long ix, iy, iz;
			Direction di;
			if (n < Nxypatch) {
				di = xy;
			}
			else
			{
				if (n < Nxypatch + Nyzpatch) {
					di = yz;
				}
				else {
					di = zx;
				}
			}
			switch (di)
			{
			case xy:
				iz = n / ((Nx*Ny));
				n = n % ((Nx*Ny));
				iy = n / Nx;
				ix = n % Nx;
				break;
			case yz:
				n = n - Nxypatch;
				iz = n / ((Nx + 1)*Ny);
				n = n % ((Nx + 1)*Ny);
				iy = n / (Nx + 1);
				ix = n % (Nx + 1);
				break;
			case zx:
				n = n - Nxypatch - Nyzpatch;
				iz = n / (Nx*(Ny + 1));
				n = n % (Nx*(Ny + 1));
				iy = n / Nx;
				ix = n % Nx;
				break;
			default:
				break;
			}
			return spaceindex(ix, iy, iz, di);
		}

		V3d get_center(Direction di, spaceindex sp) {
			V3d position;
			position.x = this->setting.x_base_list[sp.nx];
			position.y = this->setting.y_base_list[sp.ny];
			position.z = this->setting.z_base_list[sp.nz];
			switch (di)
			{
			case Brick_system::null:
				cout << "get_center: invalid direciton!" << endl;
				break;
			case Brick_system::x:
				position.x += 0.5*this->setting.delta_x_list[sp.nx];
				break;
			case Brick_system::y:
				position.y += 0.5*this->setting.delta_y_list[sp.ny];
				break;
			case Brick_system::z:
				position.z += 0.5*this->setting.delta_z_list[sp.nz];
				break;
			case Brick_system::xy:
				position.x += 0.5*this->setting.delta_x_list[sp.nx];
				position.y += 0.5*this->setting.delta_y_list[sp.ny];
				break;
			case Brick_system::yz:
				position.y += 0.5*this->setting.delta_y_list[sp.ny];
				position.z += 0.5*this->setting.delta_z_list[sp.nz];
				break;
			case Brick_system::zx:
				position.z += 0.5*this->setting.delta_z_list[sp.nz];
				position.x += 0.5*this->setting.delta_x_list[sp.nx];
				break;
			default:
				cout << "get_center: invalid direciton!" << endl;
				break;
			}
			return position;
		}

		Shader get_shader(V3d position, Field_type ft) {
			double bar = 1e-12;
			if (position.x > bar && position.x<this->setting.xmax - bar && position.y>bar && position.y<this->setting.ymax - bar && position.z>bar && position.z < this->setting.zmax - bar) {
				return NOT_BOUNDARY;
			}
			else {
				if (ft == E_field) {
					return BOUNDARY;
				}
				else {
					return NOT_BOUNDARY;
				}
				
			}
		}

		Parameter get_Para(long n, Field_type ft) {
			Parameter para;
			spaceindex sp;
			if (ft == E_field) {
				sp = this->edgeindex_to_edgespaceindex(n, this->setting.Nx, this->setting.Ny, this->setting.Nz);
			}
			else {
				sp = this->patchindex_to_patchspaceidnex(n, this->setting.Nx, this->setting.Ny, this->setting.Nz);
			}
			para.position = this->get_center(sp.di, sp);
			para.direction = this->brick_dir_to_V3d(sp.di);
			para.field_type = ft;
			para.shader = this->get_shader(para.position, para.field_type);
			return para;
		}

		SparseVec get_loop(Direction di, long nx,long ny,long nz) {
			SparseVec sv;
			switch (di)
			{
			case Brick_system::null:
				cout << "get_loop: invalid direction!" << endl;
				break;
			case Brick_system::x:
				sv.add(this->pathchspaceindex_to_patchindex(xy, spaceindex(nx, ny, nz, xy), this->setting.Nx, this->setting.Ny, this->setting.Nz), 2.0 / (this->setting.delta_y_list[ny] + this->setting.delta_y_list[ny - 1]));
				sv.add(this->pathchspaceindex_to_patchindex(zx, spaceindex(nx, ny, nz, zx), this->setting.Nx, this->setting.Ny, this->setting.Nz), -2.0 / (this->setting.delta_z_list[nz] + this->setting.delta_z_list[nz - 1]));
				sv.add(this->pathchspaceindex_to_patchindex(xy, spaceindex(nx, ny - 1, nz, xy), this->setting.Nx, this->setting.Ny, this->setting.Nz), -2.0 / (this->setting.delta_y_list[ny] + this->setting.delta_y_list[ny - 1]));
				sv.add(this->pathchspaceindex_to_patchindex(zx, spaceindex(nx, ny, nz - 1, zx), this->setting.Nx, this->setting.Ny, this->setting.Nz), 2.0 / (this->setting.delta_z_list[nz] + this->setting.delta_z_list[nz - 1]));
				break;
			case Brick_system::y:
				sv.add(this->pathchspaceindex_to_patchindex(yz, spaceindex(nx, ny, nz, yz), this->setting.Nx, this->setting.Ny, this->setting.Nz), 2.0 / (this->setting.delta_z_list[nz] + this->setting.delta_z_list[nz - 1]));
				sv.add(this->pathchspaceindex_to_patchindex(xy, spaceindex(nx, ny, nz, xy), this->setting.Nx, this->setting.Ny, this->setting.Nz), -2.0 / (this->setting.delta_x_list[nx] + this->setting.delta_x_list[nx - 1]));
				sv.add(this->pathchspaceindex_to_patchindex(yz, spaceindex(nx, ny, nz - 1, yz), this->setting.Nx, this->setting.Ny, this->setting.Nz), -2.0 / (this->setting.delta_z_list[nz] + this->setting.delta_z_list[nz - 1]));
				sv.add(this->pathchspaceindex_to_patchindex(xy, spaceindex(nx - 1, ny, nz, xy), this->setting.Nx, this->setting.Ny, this->setting.Nz), 2.0 / (this->setting.delta_x_list[nx] + this->setting.delta_x_list[nx - 1]));
				break;
			case Brick_system::z:
				sv.add(this->pathchspaceindex_to_patchindex(zx, spaceindex(nx, ny, nz, zx), this->setting.Nx, this->setting.Ny, this->setting.Nz), 2.0 / (this->setting.delta_x_list[nx] + this->setting.delta_x_list[nx - 1]));
				sv.add(this->pathchspaceindex_to_patchindex(yz, spaceindex(nx, ny, nz, yz), this->setting.Nx, this->setting.Ny, this->setting.Nz), -2.0 / (this->setting.delta_y_list[ny] + this->setting.delta_y_list[ny - 1]));
				sv.add(this->pathchspaceindex_to_patchindex(zx, spaceindex(nx - 1, ny, nz, zx), this->setting.Nx, this->setting.Ny, this->setting.Nz), -2.0 / (this->setting.delta_x_list[nx] + this->setting.delta_x_list[nx - 1]));
				sv.add(this->pathchspaceindex_to_patchindex(yz, spaceindex(nx, ny - 1, nz, yz), this->setting.Nx, this->setting.Ny, this->setting.Nz), 2.0 / (this->setting.delta_y_list[ny] + this->setting.delta_y_list[ny - 1]));
				break;
			case Brick_system::xy:
				sv.add(this->edgespaceindex_to_edgeindex(y, spaceindex(nx + 1, ny, nz, y), this->setting.Nx, this->setting.Ny, this->setting.Nz), 1.0 / this->setting.delta_x_list[nx]);
				sv.add(this->edgespaceindex_to_edgeindex(x, spaceindex(nx, ny + 1, nz, x), this->setting.Nx, this->setting.Ny, this->setting.Nz), -1.0 / this->setting.delta_y_list[ny]);
				sv.add(this->edgespaceindex_to_edgeindex(y, spaceindex(nx, ny, nz, y), this->setting.Nx, this->setting.Ny, this->setting.Nz), -1.0 / this->setting.delta_x_list[nx]);
				sv.add(this->edgespaceindex_to_edgeindex(x, spaceindex(nx, ny, nz, x), this->setting.Nx, this->setting.Ny, this->setting.Nz), 1.0 / this->setting.delta_y_list[ny]);
				break;
			case Brick_system::yz:
				sv.add(this->edgespaceindex_to_edgeindex(z, spaceindex(nx, ny + 1, nz, z), this->setting.Nx, this->setting.Ny, this->setting.Nz), 1.0 / this->setting.delta_y_list[ny]);
				sv.add(this->edgespaceindex_to_edgeindex(y, spaceindex(nx, ny, nz + 1, y), this->setting.Nx, this->setting.Ny, this->setting.Nz), -1.0 / this->setting.delta_z_list[nz]);
				sv.add(this->edgespaceindex_to_edgeindex(z, spaceindex(nx, ny, nz, z), this->setting.Nx, this->setting.Ny, this->setting.Nz), -1.0 / this->setting.delta_y_list[ny]);
				sv.add(this->edgespaceindex_to_edgeindex(y, spaceindex(nx, ny, nz, y), this->setting.Nx, this->setting.Ny, this->setting.Nz), 1.0 / this->setting.delta_z_list[nz]);
				break;
			case Brick_system::zx:
				sv.add(this->edgespaceindex_to_edgeindex(x, spaceindex(nx, ny, nz + 1, x), this->setting.Nx, this->setting.Ny, this->setting.Nz), 1.0 / this->setting.delta_z_list[nz]);
				sv.add(this->edgespaceindex_to_edgeindex(z, spaceindex(nx + 1, ny, nz, z), this->setting.Nx, this->setting.Ny, this->setting.Nz), -1.0 / this->setting.delta_x_list[nx]);
				sv.add(this->edgespaceindex_to_edgeindex(x, spaceindex(nx, ny, nz, x), this->setting.Nx, this->setting.Ny, this->setting.Nz), -1.0 / this->setting.delta_z_list[nz]);
				sv.add(this->edgespaceindex_to_edgeindex(z, spaceindex(nx, ny, nz, z), this->setting.Nx, this->setting.Ny, this->setting.Nz), 1.0 / this->setting.delta_x_list[nx]);
				break;
			default:
				cout << "get_loop: invalid direction!" << endl;
				break;
			}
			return sv;
		}

		Brick_Domain(){}

		Brick_Domain(brick_setting setting) {
			cout << "Build Brick Domain..." << endl;
			this->setting = setting;
			this->Ne = setting.Nx * (setting.Ny + 1)*(setting.Nz + 1) + (setting.Nx + 1)*setting.Ny*(setting.Nz + 1) + (setting.Nx + 1)*(setting.Ny + 1)*setting.Nz;
			this->Nh = setting.Nx * setting.Ny*(setting.Nz + 1) + (setting.Nx + 1)*setting.Ny*setting.Nz + setting.Nx * (setting.Ny + 1)*setting.Nz;
			for (long i = 0; i < this->Ne; i++) {
				this->Eparas.push_back(this->get_Para(i, E_field));
			}
			for (long i = 0; i < this->Nh; i++) {
				this->Hparas.push_back(this->get_Para(i, H_field));
			}
			//set Sh
			this->Sh = vector<SparseVec>(this->Ne);
			for (long i = 0; i < this->Ne; i++) {
				if (this->Eparas[i].shader != BOUNDARY) {
					spaceindex sp = this->edgeindex_to_edgespaceindex(i, this->setting.Nx, this->setting.Ny, this->setting.Nz);
					this->Sh[i] = this->get_loop(sp.di, sp.nx, sp.ny, sp.nz);
				}
			}
			//set Se
			this->Se = vector<SparseVec>(this->Nh);
			for (long i = 0; i < this->Nh; i++) {
				if (this->Hparas[i].shader != BOUNDARY) {
					spaceindex sp = this->patchindex_to_patchspaceidnex(i, this->setting.Nx, this->setting.Ny, this->setting.Nz);
					this->Se[i] = this->get_loop(sp.di, sp.nx, sp.ny, sp.nz);
				}
			}
			cout << "Brick Domain Finished!: Ne=" << this->Ne << " Nh=" << this->Nh << endl;
		}

		double get_edge_length(long eid) {
			spaceindex sp = this->edgeindex_to_edgespaceindex(eid, this->setting.Nx, this->setting.Ny, this->setting.Nz);
			V3d edge_direction = this->brick_dir_to_V3d(sp.di);
			V3d edge_start = V3d{ this->setting.x_base_list[sp.nx],this->setting.y_base_list[sp.ny],this->setting.z_base_list[sp.nz] };
			double edge_length = 1.0;
			switch (sp.di)
			{
			case Brick_system::x:
				edge_length = this->setting.delta_x_list[sp.nx];
				break;
			case Brick_system::y:
				edge_length = this->setting.delta_y_list[sp.ny];
				break;
			case Brick_system::z:
				edge_length = this->setting.delta_z_list[sp.nz];
				break;
			default:
				cout << "Brick_Domain::get_edge_length: invalid direction!" << endl;
				break;
			}
			return edge_length;
		}

		long get_patch_id_by_spaceindex(Direction dir, long nx, long ny, long nz) {
			return this->pathchspaceindex_to_patchindex(dir, spaceindex(nx, ny, nz, dir), this->setting.Nx, this->setting.Ny, this->setting.Nz);
		}

		SparseVec get_simple_h_expansion(V3d position, V3d direction) {
			SparseVec sv;
			spaceindex spindex = this->setting.get_which_brick_it_falls_into(position);
			long idx = spindex.nx;
			long idy = spindex.ny;
			long idz = spindex.nz;
			double a = (position.x - this->setting.x_base_list[idx]) / this->setting.delta_x_list[idx];
			double b = (position.y - this->setting.y_base_list[idy]) / this->setting.delta_y_list[idy];
			double c = (position.z - this->setting.z_base_list[idz]) / this->setting.delta_z_list[idz];
			//set in x direction
			{
				SparseVec sv_x;
				sv_x.add(this->get_patch_id_by_spaceindex(Brick_system::yz, idx, idy, idz), 1 - a);
				sv_x.add(this->get_patch_id_by_spaceindex(Brick_system::yz, idx + 1, idy, idz), a);
				sv = sv + direction.x * sv_x;
			}
			//set in y direction
			{
				SparseVec sv_y;
				sv_y.add(this->get_patch_id_by_spaceindex(Brick_system::zx, idx, idy, idz), 1 - b);
				sv_y.add(this->get_patch_id_by_spaceindex(Brick_system::zx, idx, idy + 1, idz), b);
				sv = sv + direction.y * sv_y;
			}
			//set in z direction
			{
				SparseVec sv_z;
				sv_z.add(this->get_patch_id_by_spaceindex(Brick_system::xy, idx, idy, idz), 1 - c);
				sv_z.add(this->get_patch_id_by_spaceindex(Brick_system::xy, idx, idy, idz + 1), c);
				sv = sv + direction.z * sv_z;
			}
			return sv;
		}

		SparseVec get_h_expansion(spaceindex spindex, V3d position, V3d direction) {
			long idx = spindex.nx;
			long idy = spindex.ny;
			long idz = spindex.nz;
			/*
			double a = (2 * (position.x - this->setting.x_base_list[idx])) / (this->setting.delta_x_list[idx] + this->setting.delta_x_list[idx - 1]);
			double ax = (position.x - this->setting.x_base_list[idx]) / this->setting.delta_x_list[idx];
			double ac = this->setting.delta_x_list[idx] / (this->setting.delta_x_list[idx] + this->setting.delta_x_list[idx - 1]);
			double b = (2 * (position.y - this->setting.y_base_list[idy])) / (this->setting.delta_y_list[idy] + this->setting.delta_y_list[idy - 1]);
			double by = (position.y - this->setting.y_base_list[idy]) / this->setting.delta_y_list[idy];
			double bc = this->setting.delta_y_list[idy] / (this->setting.delta_y_list[idy] + this->setting.delta_y_list[idy - 1]);
			double c = (2 * (position.z - this->setting.z_base_list[idz])) / (this->setting.delta_z_list[idz] + this->setting.delta_z_list[idz - 1]);
			double cz = (position.z - this->setting.z_base_list[idz]) / this->setting.delta_z_list[idz];
			double cc = this->setting.delta_z_list[idz] / (this->setting.delta_z_list[idz] + this->setting.delta_z_list[idz - 1]);
			*/
			double a = (2 * (position.x - this->setting.x_base_list[idx]) + this->setting.delta_x_list[idx - 1]) / (this->setting.delta_x_list[idx] + this->setting.delta_x_list[idx - 1]);
			double a0 = (position.x - this->setting.x_base_list[idx]) / this->setting.delta_x_list[idx];
			double b = (2 * (position.y - this->setting.y_base_list[idy]) + this->setting.delta_y_list[idy - 1]) / (this->setting.delta_y_list[idy] + this->setting.delta_y_list[idy - 1]);
			double b0 = (position.y - this->setting.y_base_list[idy]) / this->setting.delta_y_list[idy];
			double c = (2 * (position.z - this->setting.z_base_list[idz]) + this->setting.delta_z_list[idz - 1]) / (this->setting.delta_z_list[idz] + this->setting.delta_z_list[idz - 1]);
			double c0 = (position.z - this->setting.z_base_list[idz]) / this->setting.delta_z_list[idz];
			SparseVec sv;
			{
				//interpolate in x
				SparseVec Tempsv_x;
				Tempsv_x.add(this->get_patch_id_by_spaceindex(Brick_system::yz, idx, idy - 1, idz - 1), (1 - a0)*(1 - b)*(1 - c));
				Tempsv_x.add(this->get_patch_id_by_spaceindex(Brick_system::yz, idx + 1, idy - 1, idz - 1), a0 * (1 - b) * (1 - c));
				Tempsv_x.add(this->get_patch_id_by_spaceindex(Brick_system::yz, idx, idy, idz - 1), (1 - a0) * b * (1 - c));
				Tempsv_x.add(this->get_patch_id_by_spaceindex(Brick_system::yz, idx + 1, idy, idz - 1), a0 * b * (1 - c));
				Tempsv_x.add(this->get_patch_id_by_spaceindex(Brick_system::yz, idx, idy, idz), (1 - a0) * b * c);
				Tempsv_x.add(this->get_patch_id_by_spaceindex(Brick_system::yz, idx + 1, idy, idz), a0 * b * c);
				Tempsv_x.add(this->get_patch_id_by_spaceindex(Brick_system::yz, idx, idy - 1, idz), (1 - a0)*(1 - b)*c);
				Tempsv_x.add(this->get_patch_id_by_spaceindex(Brick_system::yz, idx + 1, idy - 1, idz), a0*(1 - b)*c);
				sv = sv + direction.x*Tempsv_x;
			}
			{
				//interpolate in y
				SparseVec Tempsv_y;
				Tempsv_y.add(this->get_patch_id_by_spaceindex(Brick_system::zx, idx - 1, idy, idz - 1), (1 - a) * (1 - b0) * (1 - c));
				Tempsv_y.add(this->get_patch_id_by_spaceindex(Brick_system::zx, idx - 1, idy + 1, idz - 1), (1 - a) * b0 * (1 - c));
				Tempsv_y.add(this->get_patch_id_by_spaceindex(Brick_system::zx, idx - 1, idy, idz), (1 - a) * (1 - b0) * c);
				Tempsv_y.add(this->get_patch_id_by_spaceindex(Brick_system::zx, idx - 1, idy + 1, idz), (1 - a) * b0 * c);
				Tempsv_y.add(this->get_patch_id_by_spaceindex(Brick_system::zx, idx, idy, idz), a * (1 - b0) * c);
				Tempsv_y.add(this->get_patch_id_by_spaceindex(Brick_system::zx, idx, idy + 1, idz), a * b0 * c);
				Tempsv_y.add(this->get_patch_id_by_spaceindex(Brick_system::zx, idx, idy, idz - 1), a * (1 - b0) * (1 - c));
				Tempsv_y.add(this->get_patch_id_by_spaceindex(Brick_system::zx, idx, idy + 1, idz - 1), a * b0 * (1 - c));
				sv = sv + direction.y*Tempsv_y;
			}
			{
				//interpolate in z
				SparseVec Tempsv_z;
				Tempsv_z.add(this->get_patch_id_by_spaceindex(Brick_system::xy, idx, idy, idz), a * b * (1 - c0));
				Tempsv_z.add(this->get_patch_id_by_spaceindex(Brick_system::xy, idx, idy, idz + 1), a * b * c0);
				Tempsv_z.add(this->get_patch_id_by_spaceindex(Brick_system::xy, idx - 1, idy, idz), (1 - a) * b * (1 - c0));
				Tempsv_z.add(this->get_patch_id_by_spaceindex(Brick_system::xy, idx - 1, idy, idz + 1), (1 - a) * b * c0);
				Tempsv_z.add(this->get_patch_id_by_spaceindex(Brick_system::xy, idx - 1, idy - 1, idz), (1 - a) * (1 - b) * (1 - c0));
				Tempsv_z.add(this->get_patch_id_by_spaceindex(Brick_system::xy, idx - 1, idy - 1, idz + 1), (1 - a) * (1 - b) * c0);
				Tempsv_z.add(this->get_patch_id_by_spaceindex(Brick_system::xy, idx, idy - 1, idz), a * (1 - b) * (1 - c0));
				Tempsv_z.add(this->get_patch_id_by_spaceindex(Brick_system::xy, idx, idy - 1, idz + 1), a * (1 - b) * c0);
				sv = sv + direction.z*Tempsv_z;
			}
			
			return sv;
		}

		bool check_if_h_interpolation_in_brick_valid(long idx, long idy, long idz, vector<Shader>& h_shader_list) {
			//check in x
			{
				Shader shader_011 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::yz, idx, idy, idz)];
				Shader shader_001 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::yz, idx, idy - 1, idz)];
				Shader shader_010 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::yz, idx, idy, idz - 1)];
				Shader shader_000 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::yz, idx, idy - 1, idz - 1)];
				Shader shader_111 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::yz, idx + 1, idy, idz)];
				Shader shader_101 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::yz, idx + 1, idy - 1, idz)];
				Shader shader_110 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::yz, idx + 1, idy, idz - 1)];
				Shader shader_100 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::yz, idx + 1, idy - 1, idz - 1)];
				bool cond1 = shader_011 == DISABLED || shader_011 == INTERFACE;
				bool cond2 = shader_001 == DISABLED || shader_001 == INTERFACE;
				bool cond3 = shader_010 == DISABLED || shader_010 == INTERFACE;
				bool cond4 = shader_000 == DISABLED || shader_000 == INTERFACE;
				bool cond5 = shader_111 == DISABLED || shader_111 == INTERFACE;
				bool cond6 = shader_101 == DISABLED || shader_101 == INTERFACE;
				bool cond7 = shader_110 == DISABLED || shader_110 == INTERFACE;
				bool cond8 = shader_100 == DISABLED || shader_100 == INTERFACE;
				if (cond1 || cond2 || cond3 || cond4 || cond5 || cond6 || cond7 || cond8) {
					return false;
				}
			}
			//check in y
			{
				Shader shader_011 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::zx, idx, idy, idz)];
				Shader shader_001 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::zx, idx, idy, idz - 1)];
				Shader shader_010 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::zx, idx - 1, idy, idz)];
				Shader shader_000 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::zx, idx - 1, idy, idz - 1)];
				Shader shader_111 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::zx, idx, idy + 1, idz)];
				Shader shader_101 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::zx, idx, idy + 1, idz - 1)];
				Shader shader_110 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::zx, idx - 1, idy + 1, idz)];
				Shader shader_100 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::zx, idx - 1, idy + 1, idz - 1)];
				bool cond1 = shader_011 == DISABLED || shader_011 == INTERFACE;
				bool cond2 = shader_001 == DISABLED || shader_001 == INTERFACE;
				bool cond3 = shader_010 == DISABLED || shader_010 == INTERFACE;
				bool cond4 = shader_000 == DISABLED || shader_000 == INTERFACE;
				bool cond5 = shader_111 == DISABLED || shader_111 == INTERFACE;
				bool cond6 = shader_101 == DISABLED || shader_101 == INTERFACE;
				bool cond7 = shader_110 == DISABLED || shader_110 == INTERFACE;
				bool cond8 = shader_100 == DISABLED || shader_100 == INTERFACE;
				if (cond1 || cond2 || cond3 || cond4 || cond5 || cond6 || cond7 || cond8) {
					return false;
				}
			}
			//check in z
			{
				Shader shader_011 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::xy, idx, idy, idz)];
				Shader shader_001 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::xy, idx - 1, idy, idz)];
				Shader shader_010 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::xy, idx, idy - 1, idz)];
				Shader shader_000 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::xy, idx - 1, idy - 1, idz)];
				Shader shader_111 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::xy, idx, idy, idz + 1)];
				Shader shader_101 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::xy, idx - 1, idy, idz + 1)];
				Shader shader_110 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::xy, idx, idy - 1, idz + 1)];
				Shader shader_100 = h_shader_list[this->get_patch_id_by_spaceindex(Brick_system::xy, idx - 1, idy - 1, idz + 1)];
				bool cond1 = shader_011 == DISABLED || shader_011 == INTERFACE;
				bool cond2 = shader_001 == DISABLED || shader_001 == INTERFACE;
				bool cond3 = shader_010 == DISABLED || shader_010 == INTERFACE;
				bool cond4 = shader_000 == DISABLED || shader_000 == INTERFACE;
				bool cond5 = shader_111 == DISABLED || shader_111 == INTERFACE;
				bool cond6 = shader_101 == DISABLED || shader_101 == INTERFACE;
				bool cond7 = shader_110 == DISABLED || shader_110 == INTERFACE;
				bool cond8 = shader_100 == DISABLED || shader_100 == INTERFACE;
				if (cond1 || cond2 || cond3 || cond4 || cond5 || cond6 || cond7 || cond8) {
					return false;
				}
			}
			return true;
		}

		bool move_h_interpolation_brick_if_needed(long& idx, long& idy, long& idz, vector<Shader>& h_shader_list) {
			if (this->check_if_h_interpolation_in_brick_valid(idx, idy, idz, h_shader_list)) {
				return true;
			}
			if (this->check_if_h_interpolation_in_brick_valid(idx, idy + 1, idz, h_shader_list)) {
				idy = idy + 1;
				return true;
			}
			if (this->check_if_h_interpolation_in_brick_valid(idx, idy, idz + 1, h_shader_list)) {
				idz = idz + 1;
				return true;
			}
			if (this->check_if_h_interpolation_in_brick_valid(idx, idy + 1, idz + 1, h_shader_list)) {
				idy = idy + 1;
				idz = idz + 1;
				return true;
			}
			if (this->check_if_h_interpolation_in_brick_valid(idx + 1, idy, idz, h_shader_list)) {
				idx = idx + 1;
				return true;
			}
			if (this->check_if_h_interpolation_in_brick_valid(idx + 1, idy + 1, idz, h_shader_list)) {
				idx = idx + 1;
				idy = idy + 1;
				return true;
			}
			if (this->check_if_h_interpolation_in_brick_valid(idx + 1, idy, idz + 1, h_shader_list)) {
				idx = idx + 1;
				idz = idz + 1;
				return true;
			}
			if (this->check_if_h_interpolation_in_brick_valid(idx + 1, idy + 1, idz + 1, h_shader_list)) {
				idx = idx + 1;
				idy = idy + 1;
				idz = idz + 1;
				return true;
			}
			return false;
		}

		SparseVec get_general_h_expansion(V3d position, V3d direction, vector<Shader>& h_shader_list) {
			Brick_system::spaceindex sp = this->setting.get_which_brick_it_falls_into(position);
			long idx = sp.nx;
			long idy = sp.ny;
			long idz = sp.nz;
			if (this->move_h_interpolation_brick_if_needed(idx, idy, idz, h_shader_list)) {
				return this->get_h_expansion(Brick_system::spaceindex(idx, idy, idz, Brick_system::null), position, direction);
			}
			else {
				cout << "Brick_Domain::get_general_h_expansion: No useful brick" << endl;
			}
		}
	};
}