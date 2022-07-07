class ltraj;

#ifndef LTRAJ_H
#define LTRAJ_H

#include <gsl/gsl_vector.h>
#include <string.h>
#include <stdio.h>
#include "libsmurf.h"
#include "fretdata.h"
#include "xtraj.h"
#include "landscape.h"
#include "friction.h"

class ltraj{
public:
	friend double landscape::LLR(ltraj*,uint,uint);
	ltraj(string datafile);
	ltraj();
	~ltraj();

	void set_datafile(string datafile);
	void load_burst(string datafile, string burstfile, int burstnumber);

	double L(xtraj *xt);
	double LL(xtraj *xt);
	double LLR(xtraj *xt1, xtraj *xt2, uint i1, uint i2);
	double LLR(xtraj *xt, double xn, uint index);
	double LL_Intensity(xtraj *xt);
	double LLR_Intensity(xtraj *xt1, xtraj *xt2, uint i1, uint i2);
	double LLR_Intensity(xtraj *xt, double xn, uint index);
 	double LL_Langevin(xtraj *xt);
 	double LLR_Langevin(xtraj *xt1, xtraj *xt2, uint i1, uint i2);
	double LLR_Langevin(friction *g1, friction *g2, uint i1, uint i2);
	double LLR_Langevin(xtraj *xt, double dn, uint index);

	void dLLda(double* grad);
	void dLLda(xtraj& xt, double* grad);
	void mean_dLLda(double* grad, double* gradstd, double error);
	void mean_dLLda(xtraj& xt, double* grad, double* gradstd, double error);

	void gibbs();

	double D(int order, uint i=0, double t=0);
	void D_p(int order, uint i, double t, double* grad);
	void D_adjust_params(double* grad, double ss, double sl);
	enum D_type{SIMPLE, SIMPLEV, PMF, PMFV, MEMORY};
	void D_select(D_type kind);
	void D_select(const char* kind);
	D_type get_D();

	double get_time();

	void set_photons(string datafile);
	s_fretdata* get_photons();

	void set_xtraj(xtraj* nx);
	xtraj* get_xtraj();

	void set_landscape(landscape* nl);
	landscape* get_landscape();

	void set_friction(friction* ng);
	friction* get_friction();

	double get_beta();
	void set_beta(double nb);

	double get_mass();
	void set_mass(double nm);

	void set_ncalc(int ncn);
	void set_neq(int neq);

	int param_length();

	void print(char* file);
	void print_xdata(char* file);
	void print_pmf(char* file);
	void print_pdf(char* file);
	void print_pmfp(char* file);

	void optimize_gamma(double g1, double g2, double tol=1e-2);
	int optimize_xdata();
	int optimize_xdata_dir(xtraj dir);

private:
	void init();
	double D_simp(int order, uint i=0, double t=0);
	void D_simp_p(int order, uint i, double t, double* grad);
	void D_simp_adjust_params(double* grad, double ss, double sl);
	double D_simpv(int order, uint i=0, double t=0);
	void D_simpv_p(int order, uint i, double t, double* grad);
	void D_simpv_adjust_params(double* grad, double ss, double sl);
	double D_pmf(int order, uint i=0, double t=0);
	void D_pmf_p(int order, uint i, double t, double* grad);
	void D_pmf_adjust_params(double* grad, double ss, double sl);
	double D_pmfv(int order, uint i=0, double t=0);
	void D_pmfv_p(int order, uint i, double t, double* grad);
	void D_pmfv_adjust_params(double* grad, double ss, double sl);
	double D_memory(int order, uint i=0, double t=0);
	void D_memory_p(int order, uint i, double t, double* grad);
	void D_memory_adjust_params(double* grad, double ss, double sl);
	bool accept_move(double llr);
	inline double zeta_d(double x);
	double zeta_d_p(double x);
	inline double zeta_a(double x);
	double zeta_a_p(double x);
	inline double zeta_xyz(double x, double y, double z);
	void cache_zetas(double dxn);
	void free_zeta_caches();
	void cache_photons();
	void free_photons();
	void cache_pn();
	void free_pn();

	xtraj* xdata;
	landscape *pdata;
	s_fretdata *photons;
	friction *gamma;
	double beta;
	double mass;
	double bdi;
	double bai;
	double wx2;
	double wz2;
	s_fr* pcache;
	uint* pncache;
	uint* npcache;
	uint* nacache;
	uint* ndcache;
	double* zeta_a_cache;
	double* zeta_d_cache;
	double zeta_cache_dx;
	uint32_t np;
	s_cal cal;
	gsl_rng *r;
	bool D2CONST;
	bool D1ZERO;
	bool use_cached_x;
	D_type D_chosen;
	int N_PARAM;
	int N_CALC;
	int N_EQ;
	uint param_gamma;
	uint param_mass;
	uint param_pmf_start;
	uint param_pmf_end;
	uint param_memory_start;
	uint param_memory_end;
};

#endif
