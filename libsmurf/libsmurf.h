#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <stdint.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
using namespace std;

class s_fretdata;
class s_fret2;
class s_tt;
class s_x;
class s_cal;
class s_p;
class s_fr;

#ifndef LIBSMURF_H
#define LIBSMURF_H

const uint32_t MAGIC=0xd3edf5f2;
const uint16_t VERSION=1;
const uint16_t ft_image=1;
const uint16_t ft_time=2;
const uint16_t ft_corr=3;
const uint16_t ft_intensity=4;
const uint16_t ft_fret=6;

const uint32_t NOT_SMURF=0x01;
const uint32_t NOT_FRET=0x02;
const uint32_t NOT_TT=0x04;

const double PHI=0.61803399;

const uint16_t DONOR=0x01;
const uint16_t ACCEPTOR=0x02;

enum xt_type {UNSET,NONE,LPW,JANDM,XUN,CY35};
extern xt_type SMURF_XT;

inline void read_Small_Endian(fstream& file, char* data, int size);
inline void write_Small_Endian(fstream& file, char* data, int size);

uint32_t smurf_loadwhich(const string fname);

uint32_t smurf_loadtt(const string fname, s_tt* data);
uint32_t smurf_savett(const string fname, s_tt* data);

uint32_t smurf_loadfret2(const string fname, s_fret2* data);
uint32_t smurf_savefret2(const string fname, s_fret2* data);

s_cal smurf_process_fret(string fin, fstream &tout);
s_cal smurf_process_freta(string infile, fstream &tout);
s_cal smurf_process_fretd(string infile, fstream &tout);
s_cal smurf_fret_calibrate(fstream& tin,uint32_t lt,string infile);
s_cal smurf_fret_cross(fstream& tin,uint32_t lt,double Ebar);
xt_type smurf_fret_xt_parse(char* xtv);
s_fr read_num(fstream& f, long n);
double LLT(fstream &tin, uint32_t ns, uint32_t n, uint32_t nt,int chan);
double tfind(fstream &data,uint32_t ns,uint32_t n,int chan);
uint32_t findtime(fstream& f, uint32_t N, double t,bool before=1);
int smurf_isfret(const string fname);
uint32_t smurf_tfind(uint32_t[],uint32_t);
double smurf_LLT(uint32_t[],uint32_t, uint32_t);
uint32_t smurf_tfind(uint32_t* data, uint32_t n,uint32_t* sampdata, uint32_t nsamp);
void spinner(uint32_t, uint32_t, int nskip);

double xcorr_calc(double tv[], double xcv[], s_fretdata* farr[], int nin, int N, double T, double alpha=0.1, int dim=1, double x1=-1, double x2=-1);
double xcorr_e(double tc[], double xct[], int N);
double xcorr_int(double tc[], double xct[], int N);

inline void SHFT2(double &a,double &b,double c){
	a=b;
	b=c;
}

inline void SHFT3(uint32_t &a,uint32_t &b,uint32_t &c,uint32_t d){
	a=b;
	b=c;
	c=d;
}

double smurf_MLE2B(double nd, double na, double bd, double ba, double IdbT, double IabT);
double smurf_MLEEB(double nd, double na, double bd, double ba, double IdbT, double IabT);
double smurf_STD2B(double x, double bd, double ba, double idbT, double iabT);
double smurf_STDEB(double E, double bd, double ba, double idbT, double iabT);

double smurf_STD2B(double x, double T, s_cal b);
double smurf_STDEB(double E, double T, s_cal b);
double smurf_MLE2B(double nd, double na, s_cal b);
double smurf_MLEEB(double nd, double na, s_cal b);

void PrintMatrix(gsl_matrix* toPrint);
void PrintMatrix(gsl_matrix* toPrint,string filename);

extern "C" void eig(gsl_matrix** W, gsl_matrix** d,gsl_matrix* g,gsl_matrix* M);

class s_fret2{
public:
	s_fret2(){
		power=clock=wc1=wc2=ps1=ps2=0;
		l1=l2=0;
	}

	double power,clock,wc1,wc2,ps1,ps2;
	uint32_t *data1,*data2;
	uint32_t l1,l2;
};

class s_tt{
public:
	double power,clock,wc,ps;
	uint32_t *data;
	uint32_t l;
};

class s_x{
public:
	s_x(int nset=1){
		tmin=0;
		tmax=0;
		x=0;
		std=0;
		xp=0;
		stdp=0;
		r=0;
		rs=0;
		rp=0;
		rps=0;
		n=nset;
		c=gsl_vector_alloc(n);
		cS=gsl_matrix_alloc(n,n);
		allocated=true;
		error=0;
	}

	~s_x(){
		gsl_vector_free(c);
		gsl_matrix_free(cS);
	}

	ostream& operator<<(ostream& f){
		print(f);
		return f;
	}

	void print(ostream& f){
		f.setf(ios_base::left|ios_base::fixed);
		f.precision(6);
		f.width(12);
		f << tmin;
		f.width(12);
		f  << tmax;
		f.width(12);
		f << x;
		f.width(12);
		f << std << endl;
	}

	void printp(ostream& f){
//		f.setf(ios_base::left|ios_base::fixed|ios_base::scientific);
		int w=12;
		f.precision(4);
		f.width(w);
		f << scientific << left << tmin;
		f.width(w);
		f  << tmax;
		for(int i=0;i<n;i++){
			f.width(w);
			f << gsl_vector_get(c,i);
			f.width(w);
			f << sqrt(gsl_matrix_get(cS,i,i));
		}
		f << endl;

		return;
		f.width(11);
		f << scientific << x;
		f.width(11);
		f << std;
		f.width(11);
		f << xp;
		f.width(11);
		f << stdp << endl;
	}

	s_x(const s_x& a){
		tmin=a.tmin;
		tmax=a.tmax;
		x=a.x;
		std=a.std;
		r=a.r;
		rs=a.rs;
		xp=a.xp;
		stdp=a.stdp;
		n=a.n;
		c=gsl_vector_alloc(a.n);
		cS=gsl_matrix_alloc(a.n,a.n);
		error=a.error;
		gsl_vector_memcpy(c,a.c);
		gsl_matrix_memcpy(cS,a.cS);
	}

	s_x operator=(s_x a){
		tmin=a.tmin;
		tmax=a.tmax;
		x=a.x;
		std=a.std;
		r=a.r;
		rs=a.rs;
		xp=a.xp;
		stdp=a.stdp;
		n=a.n;
		if(allocated){
			gsl_vector_free(c);
			gsl_matrix_free(cS);
		}
		c=gsl_vector_alloc(a.n);
		cS=gsl_matrix_alloc(a.n,a.n);
		gsl_vector_memcpy(c,a.c);
		gsl_matrix_memcpy(cS,a.cS);
		error=a.error;
		return *this;
	}

	bool operator!=(s_x a){
		if(fabs(x-a.x)>1.96*(std+a.std))
			return true;
		else
			return false;
	}

	bool isinvalid(){
		if((isnan(x))||(isnan(std))||(isnan(xp))||(isnan(stdp))||(tmin>tmax))
			return true;
		else
			return false;
	}

	void adjust_mean(double dx,int dim=1){
		if(dim==1) x+=dx;
		if(dim==2) xp+=dx;
		gsl_vector_set(c,dim-1,gsl_vector_get(c,dim-1)+dx);
	}

	double tmin,tmax,x,std,xp,stdp,r,rs,rp,rps;
	int error;
	int n;
	gsl_vector* c;
	gsl_matrix* cS;
	bool allocated;
};

class s_p{
public:
	s_p(){time=type=0;}

	double time;
	uint64_t itime;
	int type;
//	float a;
//	float d;
};

class s_cal{
public:
	s_cal(){
		ta_bleach=0;	//Acceptor bleaching time
		td_bleach=0; 	//Donor bleaching time
		ba=0;			//Acceptor signal to background ratio
		vba=0;
		bd=0;			//Donor signal to background ratio
		vbd=0;
		IdB=0;			//Donor intensity w/ background at x=\inf
		vIdB=0;
		IaB=0;			//Acceptor intensity w/ background at x=0
		vIaB=0;
		nt=0;			//Last photon before bleaching
		ntl=0;			//Last photon of the trajectory
		P=0;
                Vers=2;
	}

	s_cal operator=(s_cal b){
		ta_bleach=b.ta_bleach;
		td_bleach=b.td_bleach;
		ba=b.ba;
		vba=b.vba;
		bd=b.bd;
		vbd=b.vbd;
		IdB=b.IdB;
		vIdB=b.vIdB;
		IaB=b.IaB;
		vIaB=b.vIaB;
		nt=b.nt;
		ntl=b.ntl;
		P=b.P;
                Vers = b.Vers;
		return *this;
	}

	double ta_bleach;
	double td_bleach;
	double ba;
	double vba;
	double bd;
	double vbd;
	double IdB;
	double vIdB;
	double IaB;
	double vIaB;
	double P;
	uint32_t nt;
	uint32_t ntl;
        int Vers;
};

class s_fr{
public:
	s_fr(){
		cpca=0;
		cpcd=0;
		type=0;
		time=0;
	}

	s_fr operator=(s_fr b){
		type=b.type;
		time=b.time;
		cpca=b.cpca;
		cpcd=b.cpcd;
		return *this;
	}

	uint32_t cpc(int chan=0){
		switch(chan){
		case 0:
			return cpca+cpcd;
			break;
		case DONOR:
			return cpcd;
			break;
		case ACCEPTOR:
			return cpca;
			break;
		}
		return cpca+cpcd;
	}

	int type;
	double time;
	uint32_t cpca;
	uint32_t cpcd;
};

#endif //LIBSMURF_H
