class friction;

#ifndef FRICTION_H
#define FRICTION_H

#include <gsl/gsl_vector.h>
#include <iostream>
#include <math.h>
#include "dl_exception.h"
#include "ltraj_param.h"
using namespace std;

class friction:public ltraj_param{
public:
	friction();
	~friction();
	friction(const friction& nl);
	friction(double gamma);
	friction(double gamma, double nu1, double nu2, int N);

	double g(double x);
	double g(uint n);
	double x(uint n);
	double get_dx();
	using ltraj_param::print;
	void print(std::ostream& fout);
	void update_from_params();

	double LLR(ltraj* wrt, uint i1, uint i2);

	friction* clone();

private:
	void init();
	void free_gdata();
	void setup(int size);

	gsl_vector* gdata;
};

#endif
