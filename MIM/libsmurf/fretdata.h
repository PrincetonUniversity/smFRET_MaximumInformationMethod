/*

declaration of fretdata class
Manages .fr files and their data.

*/

//class s_fretdata;

#ifndef FRETDATA_H
#define FRETDATA_H

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include "libsmurf.h"
#include "traj.h"
using namespace std;

class s_fretdata{
public:
	s_fretdata(string datafile,string procfile);
	s_fretdata(string datafile);
	~s_fretdata();
	string getrawfile() {return rawfilename;}
	s_fr get_photon(uint32_t);
	s_fr get_photon(double t, bool before=1);
	uint32_t findtime(double t, bool before=1);
	s_x mistep(double t, double a);
	s_x mistepe(double t, double a);
	s_x mistepm(double tm, double a);
	s_x mipstep(double t, double a, int n=1);
	s_x miestep(double t, double a, int n=1);
	s_x box(double t1, double t2);
	s_x boxe(double t1, double t2);
	bool KS(double t1, double t2, int chan,int mul);
        double KSVal(double t1, double t2, int chan, int mul);
        double KSExp(double t1, double t2, int chan);
        double ChiSq(double t1, double t2, int chan, string testnum);
        double coxOakes(double t1, double t2, int chan, string testnum);
        double coxOakes2(double t1, double t2, int chan, string testnum);
        void buildDist(double t1, double t2, int chan, string testnum);
        double Variance(double t1, double t2, int chan);
	bool coxOakesBin(double t1, double t2,int chan, double confid, double mean, double mean2 ,uint32_t photons, string testnum);
	double LL(trajectory traj);
	bool isvalid();
	void mim(double**,uint32_t&,double a);
	void print_log(string);
	void print_stats();
	s_cal b;
        double coxOakesI,coxOakesII,coxOakesIII,coxOakesIV,CONFIDENCE1, CONFIDENCE2;
	double GaussianProb(double x, double variance);
        s_x xpgen(double t1, double t2, int n);
	void cross(double);
	void gen_cache(gsl_vector* t,gsl_vector* x,gsl_vector* s);
	double Ia(double);
	double Id(double);
        double chiSqInv(double prob, double d);
	string procfilename;
	fstream processed;
private:
	void init();
        gsl_rng *rnd;
	string rawfilename;
	fstream rawdata;

};

double MLEIntegrand(double t, void * params);
double JIntegrand(double t, void * params);
double gas_f(double x, void * params);
int drdt_f (const gsl_vector * dd, void *params, gsl_vector * f);

//Used to contain the data whose likelihood is being maximized
class xpparams{
	public:
	xpparams(int nd,int np,double t,s_cal nb){
		Np=np;
		Ndim=nd;
		ts=new s_fr[Np];
		cs=new double[Ndim];
		b=nb;
		T=t;
	}

	~xpparams(){
		delete[] ts;
		delete[] cs;
	}

	s_fr* ts;
	s_cal b;
	uint32_t Np;
	uint32_t Ndim;
	uint32_t currdim;
	double* cs;
	double T;
};

class gasParams{
    public:
    gasParams(double a){
        A=a;
    }
    
    double A;
};
#endif













