#pragma once
#include <map>
#include <functional>
#include <vector>
#include "BasicDomain.h"
#include "LinearAlgebra.h"
#include "SparseVector.h"
using namespace std;
namespace Computationspace {

	typedef function<double(double, V3d, V3d, Field_type)> AnaFunc;

	struct Tag
	{
		Shader shader;
		long additional_index;

		Tag(Shader shader, long index) {
			this->shader = shader;
			this->additional_index = index;
		}
	};

	class BoundaryCondition {
	public: 
		bool IsInterface;
		BoundaryCondition() {
			this->IsInterface = false;
		}
		BoundaryCondition(bool value) {
			this->IsInterface = value;
		}
	};

	class AnalytiacalBC :public BoundaryCondition {
	public:
		AnaFunc func;
		AnalytiacalBC(AnaFunc func) {
			this->func = func;
		}

	};

	class ComputationSpace {
	public:
		BasicDomain& domain;
		AnalytiacalBC anaBC;
		vector<BoundaryCondition> bcs;
		vector<double> Evector, Hvector;
		vector<double> curle_pre;
		vector<Tag> Etags, Htags;
		vector<double> Ecompare, Hcompare;

		ComputationSpace(BasicDomain& anydomain, AnaFunc anafunction):domain(anydomain),anaBC(AnalytiacalBC(anafunction)) {
			this->Evector = vector<double>(this->domain.get_Ne());
			this->Hvector = vector<double>(this->domain.get_Nh());
			this->curle_pre = vector<double>(this->domain.get_Nh());
			this->Ecompare = vector<double>(this->domain.get_Ne());
			this->Hcompare = vector<double>(this->domain.get_Nh());
			for (long i = 0; i < this->domain.get_Ne(); i++) {
				Parameter para = this->domain.get_Eparas(i);
				this->Etags.push_back(Tag(para.shader, -1));
			}
			for (long i = 0; i < this->domain.get_Nh(); i++) {
				Parameter para = this->domain.get_Hparas(i);
				this->Htags.push_back(Tag(para.shader, -1));
			}
		}

		double get_curlh(long n) {
			double value = 0.0;
			SparseVec sv = this->domain.get_Sh_row(n);
			for (long i = 0; i < sv.terms.size(); i++) {
				long index = sv.terms[i].index;
				double coeff = sv.terms[i].value;
				value += this->Hvector[index] * coeff;
			}
			return value;
		}

		double get_curle(long n) {
			double value = 0.0;
			SparseVec sv = this->domain.get_Se_row(n);
			for (long i = 0; i < sv.terms.size(); i++) {
				long index = sv.terms[i].index;
				double coeff = sv.terms[i].value;
				value += this->Evector[index] * coeff;
			}
			return value;
		}

		void set_E(double time) {
			for (long i = 0; i < this->Evector.size(); i++) {
				Parameter para = this->domain.get_Eparas(i);
				this->Evector[i] = this->anaBC.func(time, para.position, para.direction, E_field);
				this->Ecompare[i] = this->Evector[i];
			}
		}

		void set_H(double time) {
			//must be called after set_E!
			for (long i = 0; i < this->Hvector.size(); i++) {
				Parameter para = this->domain.get_Hparas(i);
				this->Hvector[i] = this->anaBC.func(time, para.position, para.direction, H_field);
				this->curle_pre[i] = this->get_curle(i);
				this->Hcompare[i] = this->Hvector[i];
			}
		}

		void update_E(double deltat, double time) {
			double epsi0 = 8.854187817e-12;
			for (long i = 0; i < this->Evector.size(); i++) {
				Parameter para = this->domain.get_Eparas(i);
				if (this->Etags[i].shader == BOUNDARY) {
					this->Evector[i] = this->anaBC.func(time, para.position, para.direction, E_field);
				}
				else {
					double curlh = this->get_curlh(i);
					this->Evector[i] += (deltat / epsi0)*curlh;
				}
				this->Ecompare[i] = this->anaBC.func(time, para.position, para.direction, E_field);
			}
		}

		void update_H(double deltat, double time) {
			double miu = 4 * 3.14159265358979*1e-7;
			for (long i = 0; i < this->Hvector.size(); i++) {
				Parameter para = this->domain.get_Hparas(i);
				if (this->Htags[i].shader == BOUNDARY) {
					this->Hvector[i] = this->anaBC.func(time, para.position, para.direction, H_field);
				}
				else {
					double curle = this->get_curle(i);
					this->Hvector[i] -= (deltat / miu)* (2 * curle - curle_pre[i]);
					curle_pre[i] = curle;
				}
				this->Hcompare[i] = this->anaBC.func(time, para.position, para.direction, H_field);
			}
		}

		double get_double_vector_err(vector<double>& v, vector<double>& vcompare) {
			double sum_err = 0.0, sum = 0.0;
			for (long i = 0; i < v.size(); i++) {
				sum_err += (v[i] - vcompare[i])*(v[i] - vcompare[i]);
				sum += vcompare[i] * vcompare[i];
			}
			return sqrt(sum_err) / sqrt(sum);
		}

		void time_march(double deltat, long n_step) {
			ofstream err_out;
			ofstream sampleline;
			ofstream tetra_sample;
			err_out.open("err_out.txt");
			sampleline.open("sampleline.txt");
			tetra_sample.open("tetra_sample.txt");
			long watch_point = 132;//this->domain.edgespaceindex_to_edgeindex(Brick_system::x, Brick_system::spaceindex(5, 5, 5, Brick_system::x), this->domain.setting.Nx, this->domain.setting.Ny, this->domain.setting.Nz);
			///////////////   real time_marching  /////////////////////
			this->set_E(0 * deltat);
			this->set_H(0.5*deltat);
			cout << 0 << " E_err=" << this->get_double_vector_err(this->Evector, this->Ecompare) << " H_err=" << this->get_double_vector_err(this->Hvector, this->Hcompare) << endl;
			err_out << 0 << " " << this->get_double_vector_err(this->Evector, this->Ecompare) << " " << this->get_double_vector_err(this->Hvector, this->Hcompare) << endl;
			for (long timestep = 1; timestep < (n_step)+1; timestep++) {
				this->update_E(deltat, timestep*deltat);
				this->update_H(deltat, (timestep + 0.5)*deltat);
				cout << timestep << " E_err=" << this->get_double_vector_err(this->Evector, this->Ecompare) << " H_err=" << this->get_double_vector_err(this->Hvector, this->Hcompare) << "\t";
				cout << "Ewatch=" << watch_point << " value= " << this->Evector[watch_point] << " ";
				this->domain.printPara(this->domain.get_Eparas(watch_point));
				tetra_sample << this->Evector[watch_point] << endl;
				err_out << timestep << " " << this->get_double_vector_err(this->Evector, this->Ecompare) << " " << this->get_double_vector_err(this->Hvector, this->Hcompare) << endl;
			}
			err_out.close();
			sampleline.close();
			tetra_sample.close();
		}
	};

}