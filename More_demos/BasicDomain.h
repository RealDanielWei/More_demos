#pragma once
#include <fstream>
#include <iomanip>
#include <limits>
#include <complex>

typedef std::complex<double> Phaser;
typedef std::numeric_limits< double > dbl;

enum Shader {
	BOUNDARY, NOT_BOUNDARY, INTERFACE
};

enum Field_type {
	E_field, H_field
};

struct Parameter
{
	V3d position;
	V3d direction;
	Field_type field_type;
	Shader shader;
};

class BasicDomain {
public:
	virtual long get_Ne() = 0;
	virtual long get_Nh() = 0;
	virtual Parameter get_Eparas(long i) = 0;
	virtual Parameter get_Hparas(long i) = 0;
	virtual SparseVec get_Sh_row(long i) = 0;
	virtual SparseVec get_Se_row(long i) = 0;
	
	void printPara(Parameter para) {
		if (para.field_type == E_field) {
			cout << "E ";
		}
		else {
			cout << "H ";
		}
		switch (para.shader)
		{
		case NOT_BOUNDARY:
			cout << "not_boundary ";
			break;
		case BOUNDARY:
			cout << "boundary ";
			break;
		default:
			cout << "unkown";
			break;
		}
		cout << "pos=(" << para.position.x << "," << para.position.y << "," << para.position.z << ") ";
		cout << " dir=(" << para.direction.x << "," << para.direction.y << "," << para.direction.z << ") ";
		cout << endl;
	}

	//////////////////////    output functions    /////////////////////////

	void output_Se() {
		ofstream Se_in_in_out, Se_b_in_out;
		Se_in_in_out.open("Se_in_in.txt");
		Se_b_in_out.open("Se_b_in.txt");
		for (long row_id = 0; row_id < this->get_Nh(); row_id++) {
			SparseVec sv = this->get_Se_row(row_id);
			for (long j = 0; j < sv.terms.size(); j++) {
				long column_id = sv.terms[j].index;
				double term_value = sv.terms[j].value;
				Parameter termpara = this->get_Eparas(column_id);
				if (termpara.shader == BOUNDARY) {
					Se_b_in_out << row_id << " " << column_id << " " << setprecision(dbl::max_digits10) << term_value << endl;
				}
				else {
					if (termpara.shader == NOT_BOUNDARY) {
						Se_in_in_out << row_id << " " << column_id << " " << setprecision(dbl::max_digits10) << term_value << endl;
					}
					else {
						cout << "Wrong in Se output!" << endl;
					}
				}
			}
		}

		Se_in_in_out.close();
		Se_b_in_out.close();
	}

	void output_Sh() {
		ofstream Sh_in_in_out, Sh_b_in_out;
		Sh_in_in_out.open("Sh_in_in.txt");
		Sh_b_in_out.open("Sh_b_in.txt");
		for (long row_id = 0; row_id < this->get_Ne(); row_id++) {
			SparseVec sv = this->get_Sh_row(row_id);
			for (long j = 0; j < sv.terms.size(); j++) {
				long column_id = sv.terms[j].index;
				double term_value = sv.terms[j].value;
				Parameter termpara = this->get_Hparas(column_id);
				if (termpara.shader == BOUNDARY) {
					Sh_b_in_out << row_id << " " << column_id << " " << setprecision(dbl::max_digits10) << term_value << endl;
				}
				else {
					if (termpara.shader == NOT_BOUNDARY) {
						Sh_in_in_out << row_id << " " << column_id << " " << setprecision(dbl::max_digits10) << term_value << endl;
					}
					else {
						cout << "Wrong in Sh output!" << endl;
					}

				}
			}
		}
		Sh_in_in_out.close();
		Sh_b_in_out.close();
	}
};

class RealDomain :public BasicDomain {
public:
	long Ne, Nh;
	vector<Parameter> Eparas, Hparas;
	vector<SparseVec> Se, Sh;

	long get_Ne() {
		return this->Ne;
	}

	long get_Nh() {
		return this->Nh;
	}

	Parameter get_Eparas(long i) {
		return this->Eparas[i];
	}

	Parameter get_Hparas(long i) {
		return this->Hparas[i];
	}

	SparseVec get_Sh_row(long i) {
		return this->Sh[i];
	}

	SparseVec get_Se_row(long i) {
		return this->Se[i];
	}

};

class Test_data_set {
public:
	vector<Phaser> e_in;
	vector<Phaser> e_b;
	vector<Phaser> h_in;
	vector<Phaser> h_b;
	vector<Phaser> curle_in;
	vector<Phaser> curle_b;
	vector<Phaser> curlh_in;

	Test_data_set(BasicDomain& bd) {
		long Ne = bd.get_Ne();
		long Nh = bd.get_Nh();

		Phaser zero = 0.0 + 0.0i;
		this->e_in = vector<Phaser>(Ne, zero);
		this->e_b = vector<Phaser>(Ne, zero);
		this->curlh_in = vector<Phaser>(Ne, zero);
		this->h_in = vector<Phaser>(Nh, zero);
		this->h_b = vector<Phaser>(Nh, zero);
		this->curle_in = vector<Phaser>(Nh, zero);
		this->curle_b = vector<Phaser>(Nh, zero);

		//parameters
		double phase = 0.0, k = 1.0, c = 3.0e8, epsi = 8.854187817e-12, miu = 4 * 3.14159265358979*1e-7, ita = sqrt(epsi / miu);
		double Ax = 100;
		double Ay = 100;
		double Az = 100;
		V3d vx{ 1.0,0.0,0.0 };
		V3d vy{ 0.0,1.0,0.0 };
		V3d vz{ 0.0,0.0,1.0 };
		//set e_in & e_b & curlh
		for (long i = 0; i < Ne; i++) {
			Parameter para = bd.get_Eparas(i);
			Phaser Ex = Ax * exp(-k * para.position.z * 1.0i + phase * 1.0i);
			Phaser Ey = Ay * exp(-k * para.position.x * 1.0i + phase * 1.0i);
			Phaser Ez = Az * exp(-k * para.position.y * 1.0i + phase * 1.0i);
			Phaser curlh_x = k * 1.0i * Ax * ita * exp(-k * para.position.z * 1.0i + phase * 1.0i);
			Phaser curlh_y = k * 1.0i * Ay * ita * exp(-k * para.position.x * 1.0i + phase * 1.0i);
			Phaser curlh_z = k * 1.0i * Az * ita * exp(-k * para.position.y * 1.0i + phase * 1.0i);
			if (para.shader == BOUNDARY) {
				e_b[i] = para.direction * vz * Ez + para.direction * vx * Ex + para.direction * vy * Ey;
			}
			else {
				e_in[i] = para.direction * vz * Ez + para.direction * vx * Ex + para.direction * vy * Ey;
				curlh_in[i] = para.direction * vz * curlh_z + para.direction * vx * curlh_x + para.direction * vy * curlh_y;
			}
		}
		//set h_in & curle
		for (long i = 0; i < Nh; i++) {
			Parameter para = bd.get_Hparas(i);
			Phaser Hx = Az * ita * exp(-k * para.position.y * 1.0i + phase * 1.0i);
			Phaser Hy = Ax * ita * exp(-k * para.position.z * 1.0i + phase * 1.0i);
			Phaser Hz = Ay * ita * exp(-k * para.position.x * 1.0i + phase * 1.0i);
			Phaser curle_x = -k * 1.0i * Az * exp(-k * para.position.y * 1.0i + phase * 1.0i);
			Phaser curle_y = -k * 1.0i * Ax * exp(-k * para.position.z * 1.0i + phase * 1.0i);
			Phaser curle_z = -k * 1.0i * Ay * exp(-k * para.position.x * 1.0i + phase * 1.0i);
			if (para.shader == BOUNDARY) {
				h_b[i] = para.direction * vx * Hx + para.direction * vy * Hy + para.direction * vz * Hz;
				curle_b[i] = para.direction * vx * curle_x + para.direction * vy * curle_y + para.direction * vz * curle_z;
			}
			else {
				h_in[i] = para.direction * vx * Hx + para.direction * vy * Hy + para.direction * vz * Hz;
				curle_in[i] = para.direction * vx * curle_x + para.direction * vy * curle_y + para.direction * vz * curle_z;
			}
		}
	}

	void output() {
		ofstream e_in_out, e_b_out, h_in_out, h_b_out, curle_in_out, curle_b_out, curlh_in_out;
		e_in_out.open("e_in.txt");
		e_b_out.open("e_b.txt");
		h_in_out.open("h_in.txt");
		h_b_out.open("h_b.txt");
		curle_in_out.open("curle_in.txt");
		curle_b_out.open("curle_b.txt");
		curlh_in_out.open("curlh_in.txt");

		for (long i = 0; i < e_in.size(); i++) {
			e_in_out << e_in[i].real() << " " << e_in[i].imag() << endl;
			e_b_out << e_b[i].real() << " " << e_b[i].imag() << endl;
			curlh_in_out << curlh_in[i].real() << " " << curlh_in[i].imag() << endl;
		}
		for (long i = 0; i < h_in.size(); i++) {
			h_in_out << h_in[i].real() << " " << h_in[i].imag() << endl;
			h_b_out << h_b[i].real() << " " << h_b[i].imag() << endl;
			curle_in_out << curle_in[i].real() << " " << curle_in[i].imag() << endl;
			curle_b_out << curle_b[i].real() << " " << curle_b[i].imag() << endl;
		}

		e_in_out.close();
		e_b_out.close();
		curlh_in_out.close();

		h_in_out.close();
		h_b_out.close();
		curle_in_out.close();
		curle_b_out.close();
	}
};