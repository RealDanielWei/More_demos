#pragma once
#include <vector>

using namespace std;

struct Term
{
	long index;
	double value;
};

class SparseVec {
public:
	vector<Term> terms;
	void add(long index,double value) {
		Term t{ index,value };
		this->terms.push_back(t);
	}
	void show() {
		for (long i = 0; i < this->terms.size(); i++) {
			cout << "(" << this->terms[i].index << "," << this->terms[i].value << ") ";
		}
		cout << endl;
	}
	friend SparseVec operator +(const SparseVec& A, const SparseVec& B) {
		SparseVec sv;
		//add A
		for (long i = 0; i < A.terms.size(); i++) {
			sv.terms.push_back(Term{ A.terms[i].index,A.terms[i].value });
		}
		//add B
		for (long i = 0; i < B.terms.size(); i++) {
			bool flag = false; long id = -1;
			for (long j = 0; j < A.terms.size(); j++) {
				if (B.terms[i].index == A.terms[j].index) {
					flag = true;
					id = j;
				}
			}
			if (flag == true) {
				sv.terms[id].value = sv.terms[id].value + B.terms[i].value;
			}
			else {
				sv.terms.push_back(Term{ B.terms[i].index,B.terms[i].value });
			}
		}
		return sv;
	}
	friend SparseVec operator *(double k, const SparseVec& A) {
		SparseVec sv;
		for (long i = 0; i < A.terms.size(); i++) {
			sv.terms.push_back(Term{ A.terms[i].index,k*A.terms[i].value });
		}
		return sv;
	}
};