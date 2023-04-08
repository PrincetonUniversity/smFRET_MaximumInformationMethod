#include "landscape.h"
#include "ltraj.h"

landscape::landscape(){
	init();
}

landscape::~landscape(){
	free_pdata();
}

landscape::landscape(const landscape& nl){
	init();
	copy_assist(&nl);
	setup(nl.pmfpdata->size);
	if(nl.pmfpdata)
		gsl_vector_memcpy(pmfpdata,nl.pmfpdata);
}

landscape::landscape(double x1, double x2, uint n){
	init();
	set_pdf_constant(x1,x2,n);
}

double landscape::pmfp(double xi){
	return pmfp(n(xi));
}

double landscape::pmfp(uint ni){
	if(ni>pmfpdata->size-1) ni=pmfpdata->size-1;
	if(ni<0) ni=0;
	return gsl_vector_get(pmfpdata,ni);
}

double landscape::x(uint ni){
	return u(ni);
}

double landscape::get_dx(){
	return du;
}

void landscape::set_pmf(gsl_vector* pmfn, double x1n, double dxn){
	gsl_vector* pmfpn=gsl_vector_calloc(pmfn->size);	
	double pmfpt;

	for(uint i=0;i<pmfn->size-1;i++){
		pmfpt=(gsl_vector_get(pmfn,i+1)-gsl_vector_get(pmfn,i))/dxn;
		gsl_vector_set(pmfpn,i,pmfpt);
	}

	set_pmfp(pmfpn,x1n,dxn);
	gsl_vector_free(pmfpn);
}

void landscape::set_pmfp(gsl_vector* pmfpn, double x1n, double dxn){
	setup(pmfpn->size);
	gsl_vector_memcpy(pmfpdata,pmfpn);
	u1=x1n;
	du=dxn;
}

void landscape::set_pdf_constant(double x1, double x2, int N){
	//TODO: resolution? get from env variables like hist
	gsl_vector* paramn=gsl_vector_calloc(N);
	gsl_vector_set_all(paramn,0);
	set_pmfp(paramn, x1, (x2-x1)/N);
	update_from_params();
	gsl_vector_free(paramn);
}

void landscape::load_pmf_from_file(const char* file){
 	//TODO: Make this prettier
 	ifstream pintest(file);

 	double a[2];
	int num=-1;
	while(!pintest.eof()){
		pintest >> a[0] >> a[1];
		num++;
	}
	pintest.close();

	int ind=1;

	gsl_vector* pmfx=gsl_vector_alloc(num);	
	gsl_vector* pmfy=gsl_vector_alloc(num);	
	ifstream pin(file);
	for(int i=0;i<num;i++){
		pin >> a[0] >> a[1];
		gsl_vector_set(pmfx,i,a[0]);
		gsl_vector_set(pmfy,i,a[ind]);
	}
	pin.close();

	double x1=gsl_vector_get(pmfx,0);
	double dx=gsl_vector_get(pmfx,1)-x1;

	set_pmf(pmfy,x1,dx);
}

void landscape::print(std::ostream& fout){
	print_pmfp(fout);
}

void landscape::print_pmf(std::ostream& pout){
	double pi=0;
	pout << x(0) << '\t' << pi << endl;
 	for(uint i=1;i<pmfpdata->size;i++){
		pi+=du*pmfp(i-1);
 		pout << x(i) << '\t' << pi << endl;
	}
}

void landscape::print_pmf(char* file){
 	ofstream pout(file);
	print_pmf(pout);
}

void landscape::print_pmfp(std::ostream& pout){
	for(uint i=0;i<pmfpdata->size;i++)
		pout << x(i) << '\t' << pmfp(i) << endl;
}

void landscape::init(){
	name="landscape";
	pmfpdata=NULL;
	nindmc=3;
	nlinmc=0;
	mcstep_d_lin=1.0;
	mcstep_d_ind=1.0;
	caching_xint=false;
}

void landscape::free_pdata(){
	if(pmfpdata) gsl_vector_free(pmfpdata);
	pmfpdata=NULL;
	params=NULL;
}

void landscape::setup(int size){
	if(!pmfpdata){
		pmfpdata=gsl_vector_calloc(size);
	} else if(pmfpdata->size!=uint(size)){
		free_pdata();
		pmfpdata=gsl_vector_calloc(size);
	}

	if(P_chosen==STRAIGHT)
		params=pmfpdata;
}		

void landscape::update_from_params(){
	if(P_chosen==STRAIGHT){
		return;
	} else if(P_chosen==TAYLOR){
		update_from_taylor_params(pmfpdata);
	}

	if(!lp_managing_params)
		params=pmfpdata;
}

double landscape::LLR(ltraj* wrt, uint i1, uint i2){
	if(this->length()!=ltt->length()) throw dle(dle::LTRAJ_LLP_SIZE);

	double lle1=0;
	double lle2=0;
	double lll1=0;
	double lll2=0;
	double d1_1=0;
	double d1_2=0;
	double d2_1=1;
	double d2_2=1;
	double dt=wrt->xdata->get_dt();
	landscape* pdatahold=wrt->pdata;
	uint stop=wrt->xdata->length();

	for(uint i=0;i<stop-1;i++){
 		if((caching_xint)&&((xint[current_xint][i]<i1)||(xint[current_xint][i]>i2))){
 			continue;
 		}
		if(!wrt->D1ZERO){		
			wrt->pdata=this;
			d1_1=wrt->D(1,i);
			wrt->pdata=static_cast<landscape*>(ltt);
			d1_2=wrt->D(1,i);
		}
		if(!wrt->D2CONST){
			wrt->pdata=this;
			d2_1=wrt->D(2,i);
			wrt->pdata=static_cast<landscape*>(ltt);;
			d2_2=wrt->D(2,i);
			lll1-=log(d2_1);
			lll2-=log(d2_2);
		}
		lle1-=gsl_pow_2(wrt->xdata->l(i+1)-wrt->xdata->l(i)-d1_1*dt)/d2_1;
		lle2-=gsl_pow_2(wrt->xdata->l(i+1)-wrt->xdata->l(i)-d1_2*dt)/d2_2;
 	}

 	if(wrt->D2CONST){
		wrt->pdata=this;
 		lle1/=4*wrt->D(2)*dt;
		lll1=0;//-((stop-1)/2)*log(4*M_PI*D(2,0,0)*dt);
		wrt->pdata=static_cast<landscape*>(ltt);
 		lle2/=4*wrt->D(2)*dt;
		lll2=0;//-((stop-1)/2)*log(4*M_PI*D(2,0,0)*dt); These cancel
 	} else{
		lle1/=4*dt;
		lle2/=4*dt;
		lll1=lll1/2;//-((stop-1)/2)*log(4*M_PI*dt);
		lll2=lll2/2;//-((stop-1)/2)*log(4*M_PI*dt); These cancel!
 	}

	wrt->pdata=pdatahold;
	return lle1-lle2+lll1-lll2;
}

landscape* landscape::clone(){
	return new landscape(*this);
}

// void landscape::set_pdf(gsl_vector* pdfn, double x1n, double dxn, double beta){
//  	setup_pdata(pdfn->size);

//  	gsl_vector_memcpy(pdfdata,pdfn);
//  	x1=x1n;
//  	dx=dxn;
// }


// void landscape::calc_pmf_from_pdf(double beta){
// 	//TODO: units?!
// 	for(int i=0;i<length();i++)
// 		gsl_vector_set(pmfdata,i,-log(gsl_vector_get(pdfdata,i))/beta);
// }

// void landscape::calc_pdf_from_pmf(double beta){
// 	//TODO: units?!
// 	for(int i=0;i<length();i++)
// 		gsl_vector_set(pdfdata,i,exp(-beta*gsl_vector_get(pmfdata,i)));
// }

// void landscape::calc_pmfp(){
// 	double diff;
// 	double dx2=dx*2;

// 	diff=gsl_vector_get(pmfdata,1)-gsl_vector_get(pmfdata,0);
// 	gsl_vector_set(pmfpdata,0,diff/dx);

// 	for(int i=1;i<length()-1;i++){
// 		diff=gsl_vector_get(pmfdata,i+1)-gsl_vector_get(pmfdata,i-1);
// 		gsl_vector_set(pmfpdata,i,diff/dx2);
// 	}

// 	diff=gsl_vector_get(pmfdata,length()-1)-gsl_vector_get(pmfdata,length()-2);
// 	gsl_vector_set(pmfpdata,length()-1,diff/dx);
// }

// void landscape::set_param(gsl_vector* paramn, double x1n, double dxn){
// 	if((x1n>=0)&&(dxn>=0)){
// 		x1=x1n;
// 		dx=dxn;
// 		setup_params(paramn->size);
// 	} else if(paramn->size != params->size){
// 		throw dle::LANDSCAPE_LENGTH_MISMATCH;
// 	}
// 	gsl_vector_memcpy(params,paramn);
// 	update_from_params();
// }

// void landscape::set_param(int index, double param){
// 	gsl_vector_set(params, index, param);
// 	update_from_params();
// }

// void landscape::print_pdf(char* file){
// 	ofstream pout(file);
// 	for(int i=0;i<length();i++)
// 		pout << x(i) << '\t' << pdf(i) << endl;
// }

// double landscape::pmf(double xi){
// 	return pmf(n(xi));
// }

// double landscape::pmf(int ni){
//  if(length()==0) throw dle(dle::LANDSCAPE_NO_DATA);
// 	if((ni>=length())||(ni<0)) throw dle(dle::LANDSCAPE_RANGE);
// 	return gsl_vector_get(pmfdata,ni);
// }

// double landscape::pdf(double xi){
// 	return pdf(n(xi));
// }

// double landscape::pdf(int ni){
// 	if(length()==0) throw dle(dle::LANDSCAPE_NO_DATA);
// 	if((ni>=length())||(ni<0)) throw dle(dle::LANDSCAPE_RANGE);
// 	return gsl_vector_get(pdfdata,ni);
// // // }

// void ltraj::set_pdata_params_from_file(char* file){
// 	//TODO: Make this prettier
// 	ifstream pintest(file);

// 	double a;
// 	int num=-1;
// 	while(!pintest.eof()){
// 		pintest >> a;
// 		num++;
// 	}
// 	pintest.close();

// 	gsl_vector* pdfy=gsl_vector_alloc(num);	
// 	ifstream pin(file);
// 	for(int i=0;i<num;i++){
// 		pin >> a;
// 		gsl_vector_set(pdfy,i,a);
// 	}
// 	pin.close();

// 	pdata->set_param(pdfy);
// }

