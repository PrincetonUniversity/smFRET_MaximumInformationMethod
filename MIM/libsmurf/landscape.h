class landscape;

#ifndef LANDSCAPE_H
#define LANDSCAPE_H

#include <gsl/gsl_vector.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "dl_exception.h"
#include "ltraj_param.h"
using namespace std;

class landscape:public ltraj_param{
public:
	landscape();
	~landscape();
	landscape(const landscape& nl);
	landscape(double x1, double x2, uint num=50);

	double pmf(double x);
	double pmf(uint n);
	double pdf(double x);
	double pdf(uint n);
	double pmfp(double x);
	double pmfp(uint n);
	double x(uint n);
	double get_dx();
	void set_pdf(gsl_vector* pdfn, double x1n, double dxn, double beta);
	void set_pdf_constant(double x1, double x2, int N);
	void set_pmf(gsl_vector* pmfn, double x1n, double dxn);
	void set_pmfp(gsl_vector* pmfpn, double x1n, double dxn);
	void load_pmf_from_file(const char* infile);
	using ltraj_param::print;
	void print(std::ostream& fout);
	void print_pdf(char* file);
	void print_pmf(std::ostream& pout);
	void print_pmf(char* file);
	void print_pmfp(std::ostream& pout);
	void update_from_params();

	double LLR(ltraj* wrt, uint i1, uint i2);
	
	landscape* clone();

private:
	void init();
	void free_pdata();
	void free_params();
	void setup(int size);
	void setup_pdata(int size);
	void cache_xint(ltraj**, uint);
	void free_xint(int);

	gsl_vector* pmfpdata;
};

#endif
