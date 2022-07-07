#include <fstream>
#include <string>
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
using namespace std;

enum dyn_enum {FAST_FRICTION,SLOW_FRICTION};

class langsim{
public:
	langsim();
	~langsim();
	void Evolve();
	void Reset();
	void Register_Vp(double (*NVp)(double)){Vp=NVp;}
	void Register_gam(double (*Ng)(double,double)){gam=Ng;}
	void Print(ostream& f);

	void set_x0(double xn){x0=xn;}
	void set_v0(double vn){v0=vn;}
	void set_beta(double b){beta=b;}
	void set_gamma(double g){dgam=g;}
	void set_mass(double massn){mass=mass;}
	void set_dt(double dtn){dt=dtn;}
	void set_dyn(dyn_enum nd){dyn=nd;}
	void set_markovian(bool m){markovian=m;}

	double get_x0(){return x0;}
	double get_v0(){return v0;}
	double get_x(){return x;}
	double get_v(){return v;}
	double get_beta(){return beta;}
	double get_gamma(){return dgam;}
	double get_mass(){return mass;}
	double get_dt(){return dt;}
	double get_dyn(){return dyn;}
	double get_markovian(){return markovian;}

private:
	void Evolve_FF();
	void Evolve_SF();
	static double gam_default(double x,double b);
	static double Vp_default(double x);
	static double K(double t);

	double (*Vp)(double);
	double (*gam)(double,double);
	gsl_rng *rnd;
	dyn_enum dyn;
	double dt;
	double x,x0;
	double v,v0;
	double beta;
	double mass;
	double t;
	double dgam;
	gsl_vector* mem;
	int memlength;
	bool markovian;
};
