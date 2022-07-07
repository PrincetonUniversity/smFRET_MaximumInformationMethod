#ifndef LTRAJ_PARAM_H
#define LTRAJ_PARAM_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "dl_exception.h"
using namespace std;
class ltraj;

class ltraj_param{
 public:
	ltraj_param();
	virtual ~ltraj_param();

	double get_param(double u) const;
	double get_param(uint i) const;
	void set_param(uint i, double val);
	void set_param(gsl_vector* nparam, double du=0, double ustart=0);

	void set_constant(double c, double u1n, double u2n, int nn);
	//virtual void update_from_params()=0;
	virtual void update_from_params(uint i){};
	virtual void update_from_params(uint i1, uint i2){};

	uint n(double ui) const;
	double u(uint ni) const;
	uint length() const;

	virtual void print(const char* fname);
	virtual void print(std::ostream& fout)=0;

  	void mcstep(vector<ltraj*> &wrt);
	void mcstep(vector<ltraj*> &wrt, uint i);
 	void mcstep(vector<ltraj*> &wrt, uint i1, uint i2, uint i3);
	bool accept_move(double llr);
	virtual double LLR(vector<ltraj*> &wrt, uint i1, uint i2);
	virtual double LLR(ltraj* wrt, uint i1, uint u2)=0;

	void set_mcstep_d_lin(double val) {mcstep_d_lin=val;}
	void set_mcstep_d_ind(double val) {mcstep_d_ind=val;}
	double get_mcstep_d_lin() {return mcstep_d_lin;}
	double get_mcstep_d_ind() {return mcstep_d_ind;}
	void set_nindmc(int val) {nindmc=val;}
	double get_nindmc() {return nindmc;}
	void set_nlinmc(int val) {nlinmc=val;}
	double get_nlinmc() {return nlinmc;}
	void set_accept_ratio_target_ind(double val) {accept_ratio_target_ind=val;}
	double get_accept_ratio_target_ind(){return accept_ratio_target_ind;}
	void set_accept_ratio_target_lin(double val) {accept_ratio_target_lin=val;}
	double get_accept_ratio_target_lin(){return accept_ratio_target_lin;}
	double accept_ratio_ind() const;
	void reset_accept_ratio_ind();
	double accept_ratio_lin() const;
	void reset_accept_ratio_lin();
	uint get_xint(uint ind){return xint[current_xint][ind];}
	bool xint_cached(){if(xint) return true;else return false;}

	void print_params_history_mean_std(const char* fname);
	void print_params_history_mean_std(std::ostream& fout);
	void print_params_history(const char* fname);
	void print_params_history(std::ostream& fout);
	void get_params_history_mean(gsl_vector* pmean);
	void get_params_history_std(gsl_vector* pstd);
	void get_params_history_minmax(double* min, double* max);
	void get_params_minmax(double* min, double* max);
	void get_params_history_mean_std(gsl_vector* pm, gsl_vector* ps);
	void gen_params_history_covar();
	void gen_params_history_cdf(uint N=500);
	void gen_params_history_pdf(uint N=500, uint mini=0, uint maxi=0);
	void print_params_history_covar(const char* fname);
	void print_params_history_covar(std::ostream& fout);
	void print_params_history_cdf(const char* fname, uint N=100);
	void print_params_history_cdf(std::ostream& fout, uint N=100);
	void print_params_history_pdf(const char* fname, uint N=100, uint mini=0, uint maxi=0);
	void print_params_history_pdf(std::ostream& fout, uint N=100, uint mini=0, uint maxi=0);

	void start_collecting(){collecting_history=true;}
	void stop_collecting() {collecting_history=false;}
	void reset_params_history();

	enum P_type{STRAIGHT,TAYLOR};
	void P_select(P_type);
	P_type P_get(){return P_chosen;}
	std::string getname(){return name;}
	double temp;

 protected:
	virtual void setup(int nl)=0;
	virtual ltraj_param* clone()=0;
	void copy_assist(const ltraj_param* nl);
	void update_from_taylor_params(gsl_vector* dest);

	void alloc_params(int n);
	void free_params();
	gsl_vector* params;
	bool lp_managing_params;

	void addto_params_history();
	void inc_params_history(int n);
	void free_params_history();
	void free_params_history_covar();
	void free_params_history_cdf();
	void free_params_history_pdf();
	gsl_matrix* params_history;
	gsl_matrix* params_history_covar;
	gsl_matrix* params_history_cdf;
	gsl_vector* params_history_cdfu;
	double min_cdfu, max_cdfu, dcdfu;
	gsl_matrix* params_history_pdf;
	gsl_vector* params_history_pdfu;
	double min_pdfu, max_pdfu, dpdfu;
	uint n_history_avail;
	uint n_history_used;
	bool collecting_history;

	void cache_xint(vector<ltraj*>& wrt);
	void free_xint();
	void set_caching_xint(bool val){caching_xint=val;}
	bool get_caching_xint(){return caching_xint;}
	uint** xint;
	uint* xint_sizes;
	uint xint_size;
	int current_xint; //TODO:HACK!!!
	bool caching_xint;

	double du;
	double du1;
	double u1;
	double mcstep_d_lin;
	double mcstep_d_ind;
	double accept_ratio_target_lin;
	double accept_ratio_target_ind;
	int nindmc;
	int nlinmc;
	bool mcstep_ind;
	bool mcstep_lin;
	gsl_rng *r;
	ltraj_param* ltt;
	P_type P_chosen;

	//acceptance statistics
	int num_ind_tries;
	int num_ind_accept;
	int num_lin_tries;
	int num_lin_accept;

	std::string name;
};

#endif
