#include "UnionDomain.h"
#include "ComputationSpace.h"

double fshape(double t) {
	double tao = 2*1.2e-8, t0 = 4 * tao;
	return (t - t0) * exp(-((t - t0) * (t - t0)) / (tao * tao));
}

double ana_field(double time, V3d position, V3d direction, Field_type field_marker) {
	double epsi = 8.854187817e-12, miu = 4 * 3.14159265358979 * 1e-7;
	double ita = sqrt(miu / epsi);
	ita = 1.0 / ita;
	double c_light = 3.0e8;
	double zmin = -0.75;
	double z0 = 3.5 + 1.0;
	double f_minusz = fshape(time + (position.z - z0) / c_light);
	if (field_marker == E_field) {
		V3d E{ 100 * 2 * f_minusz,0.0,0.0 };
		return direction * E;
	}
	else {
		V3d H{ 0.0,-1 * ita * 100 * 2 * f_minusz,0.0 };
		return direction * H;
	}
}

int main() {
	
	UnionDomain::RealizedUnionDomain unidomain= UnionDomain::RealizedUnionDomain();
	unidomain.add_tdomain("nodes_scatter.txt", "elements_scatter.txt", V3d{ 0.0,0.0,0.0 });
	unidomain.add_pml_bg(5, 5, 5, 5);
	unidomain.realize();
	
	/*
	cout << "output Sh...Ne=" << unidomain.get_Ne() << endl;
	unidomain.output_Sh();
	cout << "output Se...Nh=" << unidomain.get_Nh() << endl;
	unidomain.output_Se();
	Test_data_set testdata(unidomain);
	testdata.output();
	*/
	Computationspace::ComputationSpace cs(unidomain, ana_field);
	cs.time_march(2.5e-11, 1000000);
	/*
	Tetra_system::Tetra_domain tdomain("nodes_scatter.txt", "elements_scatter.txt");
	
	cout << "output Sh...Ne=" << tdomain.get_Ne() << endl;
	tdomain.output_Sh();
	cout << "output Se...Nh=" << tdomain.get_Nh() << endl;
	tdomain.output_Se();
	Test_data_set testdata(tdomain);
	testdata.output();
	*/
	/*
	Computationspace::ComputationSpace cs(tdomain, ana_field);
	cs.time_march(4e-11, 100000);
	*/
	system("pause");
}