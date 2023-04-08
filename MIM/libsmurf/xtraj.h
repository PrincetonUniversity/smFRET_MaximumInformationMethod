#ifndef XTRAJ_H
#define XTRAJ_H

#include <gsl/gsl_vector.h>
#include <fstream>
#include "fretdata.h"
#include "dl_exception.h"
#include "ltraj_param.h"
class ltraj;

class xtraj:public ltraj_param{
public:
	xtraj();
	xtraj(double tmax, double d);
	xtraj(gsl_vector* nx, double d);
	xtraj(const xtraj &nx);
	~xtraj();

	double x(double t) const;
	double x(uint n) const;
	double v(double t) const;
	double v(uint n) const;
	double l(double t) const;
	double l(uint n) const;
	double t(double t) const;
	double t(uint n) const;
	double get_dt() const;
	double get_T() const;
	gsl_vector* get_xdata() const;
	gsl_vector* get_vdata() const;
	void set_x(uint n, double xn);
	void set_v(uint n, double vn);
	void set_l(uint n, double ln);
	void set_xdata(gsl_vector* xn, double d);
	void set_xdata(double* xn);
	void set_xdata_zero(double T, double d);
	void set_xdata_zero();
	void set_xdata_constant(double T, double d, double c);
	void set_xdata_constant(double c);
	void set_xdata_from_data_const(s_fretdata* photons, double alpha, double dt);
	void set_xdata_from_data_interp(s_fretdata* photons, double alpha, double dt);
	void set_xdata_from_xdata(xtraj* nxt);
	void set_xdata_from_file(char* file, int N);
	void set_ltraj(ltraj* assocn){assoc=assocn;}
	void add_dx(gsl_vector* dx);
	void subtract_dx(gsl_vector* dx);
	void augment(double factor);
	using ltraj_param::print;
	void print(ostream &xout);
	xtraj restrict(int N=1) const;
	void normalize();
	void calc_xdata_from_vdata(uint beg, uint end);
	void calc_xdata_from_vdata(uint i);
	void calc_vdata_from_xdata(uint beg, uint end);
	void calc_vdata_from_xdata(uint i=0);
	void calc_from_ldata(uint beg, uint end);
	void calc_from_ldata(uint i);
	void update_from_params(uint i);
	void update_from_params(uint i1, uint i2);
	void get_x_minmax(double* min, double* max);
	void get_v_minmax(double* min, double* max);

	double LLR(ltraj* wrt, uint i1, uint i2);

	enum L_type{POSITION, VELOCITY};
	void L_select(L_type kind);
	L_type L_get();

	xtraj& operator=(const xtraj& xa);
	xtraj* clone();

private:
	void init();
	void free_xdata();
	void alloc_xdata(int size);
	void setup(int N=0);

	gsl_vector* xdata;
	gsl_vector* vdata;
	ltraj* assoc;

	double dt1;
	L_type L_chosen;
};

#endif
