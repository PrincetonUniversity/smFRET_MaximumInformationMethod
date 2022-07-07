#include <fstream>
#include <string>
#include <iostream>
#include <unistd.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "libsmurf.h"
#include "fretdata.h"
using namespace std;

class xcorr{
public:
	xcorr(s_fretdata** f, int n, int d=1);

	double xcorr_calc(double[],double[],double,double);
	double xcorr_calc(double[],double[]);
	double xcorr_calc_ct(double[],double[],double,double);
	double xcorr_calc_ct(double[],double[],double[]);
	double xcorr_calc_ct_xx(double[],double *,double *, double, double); // (HY-20060827), compute <x1(0)x2(t)>
	double xcorr_calc_ct_pos(double[],double[],double,double);
	void set_data(s_fretdata** d, int nin);
	void set_times(double T, double t, int N); //(HY-20060712)
	void set_times2(double T, double t, int Nt, int nx); // (HY-20060827)
	void set_alpha(double a){alpha=a;}
	double xcorr_e(double tc[], double xct[], int N);
	double xcorr_int(double tc[], double xct[], int N);
	double xcorr_int2(double tc[], double xct[], int N);

//	double xcorr::xcorr_e(double tc[], double xct[], int N);
//	double xcorr::xcorr_int(double tc[], double xct[], int N);
//	double xcorr::xcorr_int2(double tc[], double xct[], int N);

private:
	int ct_fill(s_fretdata*, s_x*, double, double);
	void ctc(double xc[], double xcn[], double xce[], int Nt,const s_x xd1[],const s_x xd2[],int N,int dim); // (HY-20060717)
	void ctc_xx(double xc[], double xcn[], double xce[], int Nt,const s_x xd1[],const s_x xd2[],int N,int dim, double x1, double x2); // (HY-20060827)
	void ctc_pos(double tau[],double xc[],int Nt,const s_x xd1[],const s_x xd2[],int N,int dim,double,double);
	double xcc(double tau[],double xc[],int Nt,const s_x xd1[],const s_x xd2[],int N,int dim);
	void xci(double xc[],double tau[],int Nt, s_x x1, s_x x2, int dim);
	double Gaussian(double,double);

	s_fretdata** farr;
	int Ntraj;
	int Ncorr;
	int Nxcorr; // number of positions to be correlated
	double Tcorr;
	double tcorr;
	int dim;
	double alpha;
};

