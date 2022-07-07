#include "ltraj_param.h"
#include "ltraj.h"
#include <time.h>
#include <gsl/gsl_sf.h>

ltraj_param::ltraj_param(){
 	r=gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r,time(NULL));
	params=NULL;
	params_history=NULL;
	params_history_covar=NULL;
	params_history_cdf=NULL;
	params_history_cdfu=NULL;
	params_history_pdf=NULL;
	params_history_pdfu=NULL;
	n_history_avail=0;
	n_history_used=0;
	collecting_history=true;
	lp_managing_params=false;
	P_chosen=STRAIGHT;
	u1=0;
	du=0;
	du1=0;
	mcstep_d_lin=0.1;
	mcstep_d_ind=0.1;
	accept_ratio_target_ind=0.6;
	accept_ratio_target_lin=0.6;

	nindmc=2;
	nlinmc=7;
	temp=1;

	num_ind_tries=0;
	num_ind_accept=0;
	num_lin_tries=0;
	num_lin_accept=0;

	xint=NULL;
	caching_xint=false;
}

ltraj_param::~ltraj_param(){
	gsl_rng_free(r);
	free_params_history();
	free_params_history_cdf();
	free_params_history_pdf();
	if(lp_managing_params)
		free_params();
}

double ltraj_param::get_param(uint i) const {
	return gsl_vector_get(params,i);
}

void ltraj_param::set_param(uint i, double val) {
	gsl_vector_set(params,i,val);
}

uint ltraj_param::n(double ui) const {
	if(ui<u1) return 0;
	uint r=uint((ui-u1)/du);
	if(r>length()-1) r=length()-1;
	return r;
}

double ltraj_param::u(uint ni) const {
	return u1+ni*du;
}

uint ltraj_param::length() const {
	return params->size;
}

double ltraj_param::accept_ratio_ind() const {
	return double(num_ind_accept)/double(num_ind_tries);
}

void ltraj_param::reset_accept_ratio_ind(){
	num_ind_accept=0;
	num_ind_tries=0;
}

double ltraj_param::accept_ratio_lin() const {
	return double(num_lin_accept)/double(num_lin_tries);
}

void ltraj_param::reset_accept_ratio_lin(){
	num_lin_accept=0;
	num_lin_tries=0;
}

void ltraj_param::set_constant(double c, double u1n, double u2n, int nn){
	du=(u2n-u1n)/(nn);
	du1=1/du;
	u1=u1n;
	setup(nn);
	gsl_vector_set_all(params,c);
}

void ltraj_param::print(const char* file){
	ofstream fout(file);
	print(fout);
}

void ltraj_param::mcstep(vector<ltraj*> &wrt){
	int stop=length();
	int size=stop;
	int fac=4;
	int nit=nlinmc;
	int mnit=(int)floor(log(stop/4.0)/log(double(fac)));
	if(mnit<nit) nit=mnit;
	int nmc=length()*nindmc;

	cache_xint(wrt);

	reset_accept_ratio_lin();
	for(int i=0;i<nit;i++){
		size=size/fac;
		mcstep(wrt,0,0,size);
		for(int j=size;j<stop-size;j+=size){
			mcstep(wrt,j-size,j,j+size);
		}
		mcstep(wrt,stop-size,stop,stop);
	}
	if(accept_ratio_lin()>accept_ratio_target_lin) mcstep_d_lin*=1.3;
	else if(accept_ratio_lin()<accept_ratio_target_lin) mcstep_d_lin*=0.8;

	reset_accept_ratio_ind();
 	for(int j=0;j<nmc;j++)
		mcstep(wrt,(int)(gsl_ran_flat(r,0,stop)));

	if(accept_ratio_ind()>accept_ratio_target_ind) mcstep_d_ind*=1.3;
	else if(accept_ratio_ind()<accept_ratio_target_ind) mcstep_d_ind*=0.8;

	if(collecting_history)
		addto_params_history();

	free_xint();
}

void ltraj_param::mcstep(vector<ltraj*> &wrt, uint i){
	double uo=get_param(i);
	double dum=gsl_ran_flat(r,-mcstep_d_ind,mcstep_d_ind);
	num_ind_tries++;

 	ltt=clone();
 	ltt->set_param(i,uo+dum);
 	ltt->update_from_params(i);
	bool am=false;
 	if(i==0) am=accept_move(temp*LLR(wrt,0,i+1));
 	else if(i==length()-1) am=accept_move(temp*LLR(wrt,i-1,i));
 	//else am=am=accept_move(temp*LLR(wrt,i-1,i+1));
        else am=accept_move(temp*LLR(wrt,i-1,i+1));

 	if(am){
		num_ind_accept++;
 		set_param(i,uo+dum);
		update_from_params(i);
	}

 	delete ltt;
	ltt=NULL;
}

void ltraj_param::mcstep(vector<ltraj*> &wrt, uint i1, uint i2, uint i3){
	double dum=gsl_ran_flat(r,-mcstep_d_lin,mcstep_d_lin);
	num_lin_tries++;

	ltt=clone();

	if(i2-i1) for(uint i=i1;i<i2;i++){
		ltt->set_param(i,get_param(i)+((double)(i-i1))/((double)(i2-i1))*dum);
	}
	if(i3-i2) for(uint i=i2;i<i3;i++){
		ltt->set_param(i,get_param(i)+((double)(i3-i))/((double)(i3-i2))*dum);
	}
	ltt->update_from_params(i1,i3);
	//ltt->print("ltt.out");

	if(accept_move(temp*LLR(wrt,i1,i3))){
 		num_lin_accept++;
		gsl_vector_memcpy(params,ltt->params);
		update_from_params(i1,i3);
	}

	delete ltt;
	ltt=NULL;
}

bool ltraj_param::accept_move(double llr){
	if((llr<=0)||(gsl_ran_flat(r,0,1)<exp(-llr)))
		return true;
	else
		return false;
}

double ltraj_param::LLR(vector<ltraj*> &wrt, uint i1, uint i2){
	double llr=0;
	for(uint i=0;i<wrt.size();i++){
		current_xint=i;
		llr+=LLR(wrt[i],i1,i2);
	}
	return llr;
}

void ltraj_param::print_params_history_mean_std(const char* fname){
	ofstream fout(fname);
	print_params_history_mean_std(fout);
}

void ltraj_param::print_params_history_mean_std(std::ostream& fout){
	gsl_vector* pm=gsl_vector_calloc(params->size);
	gsl_vector* ps=gsl_vector_calloc(params->size);

	get_params_history_mean_std(pm,ps);

	double cm=0;
	double cs=0;
	for(uint i=0;i<pm->size;i++){
		cm=gsl_vector_get(pm,i);
		cs=gsl_vector_get(ps,i);
		fout << scientific << setiosflags(ios::right) << setprecision(4);
		fout << setw(9) << u(i);
		if(fabs(cm)<1e10){
			fout << setw(15) << cm;
			fout << setw(15) << cs;
		} else{
			fout << setw(15) << 0;
			fout << setw(15) << 0;
		}
		fout << endl;
	}

	gsl_vector_free(pm);
	gsl_vector_free(ps);
}

void ltraj_param::print_params_history(const char* fname){
	ofstream fout(fname);
	print_params_history(fout);
}

void ltraj_param::print_params_history(std::ostream& fout){
	if(!params_history) return;

	double cv=0;
	fout << scientific << setiosflags(ios::right) << setprecision(4);
	for(uint i=0;i<params_history->size1;i++){
		for(uint j=0;j<n_history_used;j++){
			cv=gsl_matrix_get(params_history,i,j);
			if(fabs(cv)>1e10) cv=0;
			fout << setw(15) << cv;
		}
		fout << endl;
	}
}

void ltraj_param::get_params_history_mean(gsl_vector* pm){
	if(!params_history){
		gsl_vector_set_all(pm,0);
		return;
	}

	if(pm->size!=params_history->size1)
		throw dle::LTRAJ_PARAM_LENGTH_MISMATCH;

	for(uint i=0;i<pm->size;i++){
		double tval=0;
		for(uint j=0;j<n_history_used;j++){
			tval+=gsl_matrix_get(params_history,i,j);
		}
		tval=tval/n_history_used;
		gsl_vector_set(pm,i,tval);
	}
}

void ltraj_param::get_params_history_std(gsl_vector* ps){
	if(!params_history){
		gsl_vector_set_all(ps,0);
		return;
	}

	if(ps->size!=params_history->size1)
		throw dle::LTRAJ_PARAM_LENGTH_MISMATCH;

	double tval,tval2;
	for(uint i=0;i<ps->size;i++){
		tval=0;
		tval2=0;
		for(uint j=0;j<n_history_used;j++){
			tval+=gsl_matrix_get(params_history,i,j);
			tval2+=gsl_pow_2(gsl_matrix_get(params_history,i,j));
		}
		tval2=sqrt(tval2/n_history_used-gsl_pow_2(tval/n_history_used));
		gsl_vector_set(ps,i,tval2);
	}
}

void ltraj_param::get_params_history_mean_std(gsl_vector* pm, gsl_vector* ps){
	if(!params_history){
		gsl_vector_set_all(pm,0);
		gsl_vector_set_all(ps,0);
		return;
	}

	if(pm->size!=params_history->size1)
		throw dle::LTRAJ_PARAM_LENGTH_MISMATCH;

	if(ps->size!=params_history->size1)
		throw dle::LTRAJ_PARAM_LENGTH_MISMATCH;

	for(uint i=0;i<pm->size;i++){
		double tval=0;
		double tval2=0;
		for(uint j=0;j<n_history_used;j++){
			tval+=gsl_matrix_get(params_history,i,j);
			tval2+=gsl_pow_2(gsl_matrix_get(params_history,i,j));
		}
		tval2=sqrt(tval2/n_history_used-gsl_pow_2(tval/n_history_used));
		tval=tval/n_history_used;
		gsl_vector_set(pm,i,tval);
		gsl_vector_set(ps,i,tval2);
	}
}

void ltraj_param::get_params_history_minmax(double* min, double* max){
	if(!params_history){
		*min=0;
		*max=0;
	}

	gsl_matrix_minmax(params_history,min,max);
}

void ltraj_param::get_params_minmax(double* min, double* max){
	gsl_vector_minmax(params,min,max);
}

void ltraj_param::gen_params_history_covar(){
	if(!params_history) return;
	free_params_history_covar();

	uint N=params_history->size1;
	params_history_covar=gsl_matrix_calloc(N,N);

	double sx1, sx2, sx1x2;
	double cvt;
	for(uint i=0;i<N;i++){
		for(uint j=0;j<N;j++){
			sx1=0;
			sx2=0;
			sx1x2=0;
			for(uint k=0;k<n_history_used;k++){
				sx1+=gsl_matrix_get(params_history,i,k);
				sx2+=gsl_matrix_get(params_history,j,k);
				sx1x2+=gsl_matrix_get(params_history,i,k)*gsl_matrix_get(params_history,j,k);
			}
			cvt=sx1x2/n_history_used-sx1*sx2/gsl_pow_2(n_history_used);
			gsl_matrix_set(params_history_covar,i,j,cvt);
		}
	}
}

void ltraj_param::gen_params_history_cdf(uint N){
	if(!params_history) return;
	free_params_history_cdf();


	params_history_cdfu=gsl_vector_calloc(N);
	params_history_cdf=gsl_matrix_calloc(params_history->size1,N);

	gsl_matrix_minmax(params_history,&min_cdfu,&max_cdfu);
	dcdfu=(max_cdfu-min_cdfu)/N;
	double dcdf=1.0/params_history->size2;

	for(uint i=0;i<N;i++)
		gsl_vector_set(params_history_cdfu,i,min_cdfu+i*dcdfu);

	for(uint i=0;i<params_history->size1;i++){
		double tmpval;
		uint tmpind;
		double newcdf;
		for(uint j=0;j<n_history_used;j++){
			tmpval=gsl_matrix_get(params_history,i,j);
			tmpind=(uint)floor((tmpval-min_cdfu)/N);
			if(tmpind>=N) tmpind=N-1;
			newcdf=gsl_matrix_get(params_history_cdf,i,tmpind)+dcdf;
			gsl_matrix_set(params_history_cdf,i,tmpind,newcdf);
		}

		newcdf=0;
		for(uint j=0;j<params_history_cdf->size2;j++){
			newcdf+=gsl_matrix_get(params_history_cdf,i,j);
			gsl_matrix_set(params_history_cdf,i,j,newcdf);
		}
	}
}

void ltraj_param::gen_params_history_pdf(uint N, uint mini, uint maxi){
	if(!params_history) return;
	free_params_history_pdf();

	if(maxi==0) maxi=length();

	params_history_pdfu=gsl_vector_calloc(N);
	params_history_pdf=gsl_matrix_calloc(maxi-mini,N);

	if(getname()=="xtraj"){
		max_pdfu=1.6;
		min_pdfu=0.4;
	} else{
		double phm=0;
		double phs=0;
		long unsigned int tnum=0;
		uint i=params_history->size1/2;	
		for(uint j=0;j<params_history->size2;j++){
			tnum++;
			phm+=gsl_matrix_get(params_history,i,j);
			phs+=gsl_pow_2(gsl_matrix_get(params_history,i,j));
		}
		phs=sqrt(phs/tnum-gsl_pow_2(phm/tnum));
		phm=phm/tnum;
		max_pdfu=phm+10*phs;
		min_pdfu=phm-10*phs;
	}

	dpdfu=(max_pdfu-min_pdfu)/N;
	double dpdf=1.0/n_history_used/dpdfu;

	for(uint i=0;i<N;i++)
		gsl_vector_set(params_history_pdfu,i,min_pdfu+i*dpdfu);

	for(uint i=mini;i<maxi;i++){//params_history->size1;i++){
		double tmpval;
		uint tmpind;
		double newpdf;
		for(uint j=0;j<n_history_used;j++){
			tmpval=gsl_matrix_get(params_history,i,j);
			tmpind=(uint)floor((tmpval-min_pdfu)/dpdfu);
			if(tmpind>=N) tmpind=N-1;
			newpdf=gsl_matrix_get(params_history_pdf,i-mini,tmpind)+dpdf;
			gsl_matrix_set(params_history_pdf,i-mini,tmpind,newpdf);
		}
	}
}

void ltraj_param::print_params_history_cdf(const char* file, uint N){
	ofstream fout(file);
	print_params_history_cdf(fout, N);
}

void ltraj_param::print_params_history_cdf(std::ostream& fout, uint N){
	gen_params_history_cdf(N);
	if(!params_history_cdf) return;

	fout << scientific << setiosflags(ios::right) << setprecision(4);
	for(uint i=0;i<params_history_cdf->size1;i++){
		for(uint j=0;j<params_history_cdf->size2;j++){
			fout << setw(15) << gsl_matrix_get(params_history_cdf,i,j);
		}
		fout << endl;
	}
}

void ltraj_param::print_params_history_pdf(const char* file, uint N, uint mini, uint maxi){
	ofstream fout(file);
	fout << "% " << getname() << " mc-history pdf" << endl
		 << "% line 1:   values for dV/dx" << endl
		 << "% column 1: values for x" << endl
		 << "% the rest: probability densities" << endl;
	print_params_history_pdf(fout, N, mini, maxi);
}

void ltraj_param::print_params_history_pdf(std::ostream& fout, uint N, uint mini, uint maxi){
	gen_params_history_pdf(N, mini, maxi);
	if(!params_history_pdf) return;

	fout << scientific << setiosflags(ios::right) << setprecision(4);

	fout << setw(10) << 0;
	for(uint i=0;i<params_history_pdf->size2;i++){
		fout << setw(15) << gsl_vector_get(params_history_pdfu,i);
	}
	fout << endl;

	double cv=0;
	for(uint i=0;i<params_history_pdf->size1;i++){
		fout << setw(10) << u(i);
		for(uint j=0;j<params_history_pdf->size2;j++){
			cv=gsl_matrix_get(params_history_pdf,i,j);
			if(fabs(cv)>1e10) fout << setw(15) << 0;
			else fout << setw(15) << cv;
		}
		fout << endl;
	}
}

void ltraj_param::reset_params_history(){
	n_history_used=0;
}

void ltraj_param::copy_assist(const ltraj_param* lpn){
	num_ind_tries=lpn->num_ind_tries;
	num_ind_accept=lpn->num_ind_accept;
	num_lin_tries=lpn->num_lin_tries;
	num_lin_accept=lpn->num_lin_accept;

	mcstep_d_lin=lpn->mcstep_d_lin;
	mcstep_d_ind=lpn->mcstep_d_ind;
	du=lpn->du;
	du1=1/du;
	u1=lpn->u1;
	P_select(lpn->P_chosen);
	if(lp_managing_params) gsl_vector_memcpy(params,lpn->params);
}

void ltraj_param::update_from_taylor_params(gsl_vector* dest){
	if(!dest) return;
	gsl_vector_set_all(dest,0);
	double val;
	double tay=1;
	double xi=0;
	double xc=get_param(uint(0));
	for(uint i=0;i<dest->size;i++){
		tay=1;
		val=tay*get_param(uint(1));
		xi=u(i)-xc;
		for(uint j=2;j<params->size;j++){
			tay*=xi/(j-1);
			val+=tay*get_param(j);
		}
		gsl_vector_set(dest,i,val);
	}
}

void ltraj_param::P_select(P_type P_new){
	P_chosen=P_new;
	if(P_chosen==STRAIGHT){
		if(lp_managing_params){
			free_params();
			lp_managing_params=false;
		}
	} else if(P_chosen==TAYLOR){
		if(lp_managing_params){
		} else{
			lp_managing_params=true;
			alloc_params(15);
		}
	}
	update_from_params(0);
}

void ltraj_param::alloc_params(int n){
	if(!lp_managing_params) return;

	params=gsl_vector_calloc(n);
}

void ltraj_param::addto_params_history(){
	if(n_history_used>=n_history_avail)
		inc_params_history(100);

	gsl_vector_view hview=gsl_matrix_column(params_history,n_history_used);
	gsl_vector_memcpy(&(hview.vector),params);
	n_history_used++;
}

void ltraj_param::inc_params_history(int n){
	if(!params_history){
		params_history=gsl_matrix_calloc(params->size,n);
		n_history_avail=n;
		return;
	}

	int s1=params_history->size1;
	int s2=params_history->size2;
	gsl_matrix* tmpmat=gsl_matrix_calloc(s1,s2+n);

	gsl_matrix_view cview=gsl_matrix_submatrix(tmpmat,0,0,s1,s2);
	gsl_matrix_memcpy(&(cview.matrix),params_history);

	gsl_matrix_free(params_history);
	params_history=tmpmat;
	n_history_avail=s2+n;
}

void ltraj_param::free_params_history(){
	if(params_history)
		gsl_matrix_free(params_history);
	params_history=NULL;
}

void ltraj_param::free_params_history_covar(){
	if(params_history_covar)
		gsl_matrix_free(params_history_covar);
	params_history_covar=NULL;
}

void ltraj_param::free_params_history_cdf(){
	if(params_history_cdf)
		gsl_matrix_free(params_history_cdf);
	params_history_cdf=NULL;

	if(params_history_cdfu)
		gsl_vector_free(params_history_cdfu);
	params_history_cdfu=NULL;
}

void ltraj_param::free_params_history_pdf(){
	if(params_history_pdf)
		gsl_matrix_free(params_history_pdf);
	params_history_pdf=NULL;

	if(params_history_pdfu)
		gsl_vector_free(params_history_pdfu);
	params_history_pdfu=NULL;
}

void ltraj_param::free_params(){
	if(!lp_managing_params)
		return;

	if(params){
		gsl_vector_free(params);
	}

	params=NULL;
}

void ltraj_param::cache_xint(vector<ltraj*> &wrt){
	if(!caching_xint) return;
	xint=new uint*[wrt.size()];
	xint_sizes=new uint[wrt.size()];
	xint_size=wrt.size();

	for(uint i=0;i<wrt.size();i++){
		xint_sizes[i]=wrt[i]->get_xtraj()->length();
		xint[i]=new uint[xint_sizes[i]];

		for(uint j=0;j<xint_sizes[i];j++)
			xint[i][j]=n(wrt[i]->get_xtraj()->x(j));
	}
}

void ltraj_param::free_xint(){
	if(!xint) return;
	for(uint i=0;i<xint_size;i++){
		delete[] xint[i];
		xint[i]=NULL;
	}
	delete[] xint;
	xint=NULL;
}

