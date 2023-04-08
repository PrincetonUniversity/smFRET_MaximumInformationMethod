/*

declaration of trajectory class
Describes the distance as a function of time in a single molecule trajectory.
NOTE: Assumes the time vector has constant dt.

*/
class trajectory;

#ifndef TRAJ_H
#define TRAJ_H
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ostream>
#include <time.h>
#include "fretdata.h"

class trajectory{
public:
	friend class s_fretdata;
	trajectory(uint32_t);
	double calc_ll();
	double get_x(uint32_t n){return gsl_vector_get(traj,n);}
//	void set_x(uint32_t n, double x){gsl_vector_set(traj,n,x);}
	void fill(gsl_vector* t,gsl_vector* x, gsl_vector* s);
	void print(std::ostream& dout);
	void set_ll(double lln){ll=lln;}
	double distance(double t);
	void set_x(gsl_vector* nx){gsl_vector_memcpy(traj,nx);}

private:
	uint32_t findtime(double t);

	gsl_rng* rnd;
	double T;
	double dt;
	uint32_t Nt;
	uint32_t Nx;
	double ll;
	gsl_vector* time;
	gsl_vector* traj;
};

class tparams{
	public:
	tparams(double g,double b,double c){
		gamma=g;
		beta=b;
		c1=c;
	}

	double gamma;
	double beta;
	double c1;
};
#endif
