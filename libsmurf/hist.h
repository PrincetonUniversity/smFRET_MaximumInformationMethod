#ifndef HIST_H
#define HIST_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_math.h>
#include <unistd.h>
#include "libsmurf.h"

#define MAX_MAXENT_ITER 100000

const double HIST_C=1.357; //For Gaussian Kernel, Eq. 2.35 in Chap. 7 of Eggermont and LaRiccia
const double HIST_R=0.776388; //For Gaussian Kernel, Eq. 2.18 in Chap. 7 of Eggermont and LaRiccia
const int HIST_SCALE_LINEAR=0;
const int HIST_SCALE_LOG=1;

class hist{
public:
	hist(int ndim,int nh,int ns);
	~hist();
	void MakeHist(double dev=-1);
	void MakeHist(int);
	void Bootstrap(int nb);
	void DoMEM(int flag);
	void set_lambda(double p) {lambdaM=p;}
	void Print(const char*, bool bmem=false);
	void BootPrint(vector<double>* hxv, vector< vector<double> >* hyv, vector< vector<double> >* hdcv, vector< vector<double> >* hytv);
	void Print_MEM(const char*);
	void Print_MEM();
	void Print_data(ostream&);
	void gausskern(double dev,int nk);
	void gausshist(bool boostrapping=false);
	void normalize(gsl_vector*);
	void toProb(gsl_vector*);
	void adddatum(double w,gsl_vector* val,gsl_vector* dev);
	void adddatum(double w,gsl_vector* val,gsl_vector* dev,double c);
	void adddatum(s_x);
	double std();
	double mean();
	uint32_t getlength() {return length;}
	double get_ave_dev(int dim);
	void clear();
	void set_stepsize(double p) {step_size=p;}
	void set_tolerance(double p) {tolerance=p;}
	void setNumHist(int num);
	void setHistX(gsl_vector* x,int dim=0);
	void setHistXAveDev(double x,int dim=0);
	void setHistRaw(gsl_vector* y);
	void getHistRaw(gsl_vector* y);
	void getHistDC(gsl_vector* y);
	void getHistX(gsl_vector* x,int dim=0);
	double getHistX(int index,int dim=0);
	int NumHist(int n=0){return numhist;}

private:
	void MatrixMul(gsl_matrix* A,CBLAS_TRANSPOSE trA,gsl_matrix* B,CBLAS_TRANSPOSE trB,gsl_matrix** S);
	void rawhist();
	void augmentlavail();
	double HKernel(gsl_vector*,gsl_vector*,gsl_vector*,double);
	//This method is added for debugging...It makes calls to functions to make sure they work properly...
	void test();
	void fill_hx();
	void hist_alloc_space();
	void hist_free_space();
	void hist_alloc_space_aug();
	void hist_free_space_aug();
	inline double N(double x, double h2);
	void initHKernel2d();

	//MEM functions
	//double X2();
	//double H();
	void init_hdc();
	bool is_hdc_valid();
	double  MaxEnt_iterate(int flag=0);
	void buildResponse();
	gsl_vector* k;
	double lambdaM;
	double step_size;
	double tolerance;
	double testparam;

	//Histogram data
	int scale;
	int dim;			//dimensionality of histogram
	int n;				//number of histogram points counted
	int numhist;
	gsl_vector** hx;	//x values of histogram
	gsl_vector* hy;		//probability densities
	gsl_vector* sigma;  //This contains the histogram errors for each bin as determined by bootstrapping... 
	gsl_vector* kernel;	//kernel
	gsl_vector* xk;
	gsl_vector* meanAlpha; //This is the mean alpha value for every bin which should be the chosen alpha in boots
	//but if we want to use say constant-time binning then we will be determining the mean for every bin...
	gsl_vector *minx,*maxx,*binsize;
	gsl_vector* kvp;		//2d: variable variance data for x'
    gsl_vector* thy;  //This is histogram produced by function gausshist for bootstrapping or not...
	gsl_vector* hdc; //deconvoluted histogram which is changing during every itertion of MaxEnt_iterate...
	gsl_vector* hyt; //trial histogram which is hdc reconvoluted using the response matrix...
	bool bootInit;

	gsl_rng* rnd;

	uint32_t length;	//number of data points
	uint32_t lavail;	//size of vectors available for data points
	gsl_vector* weight;		//weight factors for data points
	gsl_vector** value;		//values of data points
	gsl_vector** dev;		//standard deviation of data points
	gsl_vector* cov;		//2d: covariance between x and x'
	double h,lambda;
	double* cache;
	
	//Some Skilling stuff....
	gsl_matrix* Response;

};

#endif
