#include "ltraj.h"
#include <time.h>
#include <fstream>
#include <typeinfo>
#include <gsl/gsl_sf_gamma.h>

ltraj::ltraj(){
	init();
}

ltraj::ltraj(string datafile){
	init();
	set_datafile(datafile);
}

ltraj::~ltraj(){
	gsl_rng_free(r);
	if(photons){
		delete photons;
		photons=NULL;
	}
	free_photons();
	free_zeta_caches();
}

void ltraj::init(){
	beta=1;
	mass=1;
	np=0;
	photons=NULL;
	pcache=NULL;
	pncache=NULL;
	npcache=NULL;
	ndcache=NULL;
	nacache=NULL;
	xdata=NULL;
	pdata=NULL;
	gamma=NULL;
	zeta_a_cache=NULL;
	zeta_d_cache=NULL;
	r=gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r,time(NULL));
	D_select(SIMPLE);
	use_cached_x=false;
}

void ltraj::set_datafile(string datafile){
	set_photons(datafile);
	cal=photons->b;
	cache_photons();
	D_select(D_chosen);
}

void ltraj::load_burst(string datafile, string burstfile, int burstnum){
	//find limits of selected burst
	ifstream data(datafile.c_str());
	ifstream bursts(burstfile.c_str());

	double bstart=0, bstop=0;

	for(int i=0;i<burstnum;i++){
		bursts >> bstart >> bstart >> bstart >>  bstop;
	}

	np=(long)(bstop-bstart);
	double* bphotons=new double[np];

	for(int i=0;i<bstart;i++){
		data >> bphotons[0];
		data.ignore(1000,'\n');
	}

	for(uint i=0;i<np;i++){
		data >> bphotons[i];
		data.ignore(1000,'\n');
	}

	for(uint i=1;i<np;i++){
		bphotons[i]=bphotons[i]-bphotons[0];
	}
	bphotons[0]=0;

	//setup the photon cache
	free_photons();
	pcache=new s_fr[np];
	s_fr cp;

	for(uint i=0;i<np;i++){
		cp.cpca=0;
		cp.cpcd=i;
		cp.time=bphotons[i];
		cp.type=DONOR;
		pcache[i]=cp;
	}

	cal.ba=0;
	cal.bd=21;
	cal.IaB=0;
	cal.IdB=4.2e4;
	cal.ta_bleach=bphotons[np-1];
	wx2=gsl_pow_2(0.51);
	wz2=gsl_pow_2(2.04);

	delete[] bphotons;
	//	xdata.start_moving();
}

double ltraj::L(xtraj *xt){
	return exp(LL(xt));
}

double ltraj::LL(xtraj *xt){
	return LL_Langevin(xt) + LL_Intensity(xt);
}

// double ltraj::LLR(const ltraj_param &lp1, const ltraj_param &lp2, int i1, int i2){
// 	cout << typeid(lp1).name() << endl;
// 	cout << typeid(&lp1).name() << endl;
// 	exit(1);
// }

double ltraj::LLR(xtraj *xt1, xtraj *xt2, uint i1, uint i2){
	return LLR_Langevin(xt1,xt2,i1,i2) + LLR_Intensity(xt1,xt2,i1,i2);
}

inline double ltraj::LLR(xtraj *xt, double xn, uint i){
	return LLR_Langevin(xt,xn,i) + LLR_Intensity(xt,xn,i);
}

double ltraj::LL_Intensity(xtraj *xt){
	double ll=0;
	double lni=0;
	double xtx=0;

	for(uint i=0;i<np;i++){
		xtx=xt->x(pncache[i]);
		if(pcache[i].type==DONOR){
			lni=log(cal.IdB*zeta_d(xtx));
		}
		else if(pcache[i].type==ACCEPTOR){
			lni=log(cal.IaB*zeta_a(xtx));
		}
		else{
			cerr << "Invalid photon type!" << endl;				
		}
// 		if(xt.is_moving()){
// 			lni*=zeta_xyz(xt.px(pncache[i]),xt.py(pncache[i]),xt.pz(pncache[i]));
// 		}
		ll+=lni;
	}

	double intd=0, inta=0, intdi=0, intai=0;
	uint stop=xt->length();
	for(uint i=0;i<stop;i++){
		intdi+=zeta_d(xt->x(i));
		intai+=zeta_a(xt->x(i));
// 		if(xt.is_moving()){
// 			intdi*=zeta_xyz(xt.px(pncache[i]),xt.py(pncache[i]),xt.pz(pncache[i]));
// 			intai*=zeta_xyz(xt.px(pncache[i]),xt.py(pncache[i]),xt.pz(pncache[i]));
// 		}
		intd+=intdi;
		inta+=intai;
	}
	ll-=cal.IdB*xt->get_dt()*log(intd);
	ll-=cal.IaB*xt->get_dt()*log(inta);

	return ll;
}

double ltraj::LLR_Intensity(xtraj *xt1, xtraj *xt2, uint i1, uint i2){
	if(xt1->length()!=xt2->length()) throw dle(dle::LTRAJ_LLP_SIZE);

	double xo=0;
	double xn=0;
	double zdo=0;
	double zdn=0;
	double zao=0;
	double zan=0;
	double lnzdxon=0;
	double lnzaxon=0;
	double llp=0;
	double llid=0;
	double llia=0;
	double zxyzo=1;
	double zxyzn=1;

	for(uint i=i1;i<i2;i++){
		xo=xt1->x(i);
		xn=xt2->x(i);
// 		if(xt1.is_moving()){
// 			zxyzo=zeta_xyz(xt1.px(i),xt1.py(i),xt1.pz(i));
// 			zxyzn=zeta_xyz(xt2.px(i),xt2.py(i),xt2.pz(i));
// 		}
		zdo=zeta_d(xo)*zxyzo;
		zdn=zeta_d(xn)*zxyzn;
		zao=zeta_a(xo)*zxyzo;
		zan=zeta_a(xn)*zxyzn;

		lnzdxon=log(zdo/zdn);
		lnzaxon=log(zao/zan);

		llp+=ndcache[i]*log(zdo/zdn)+nacache[i]*log(zao/zan);
		llia+=zao-zan;
		llid+=zdo-zdn;
	}

	return llp-xt1->get_dt()*(cal.IaB*llia+cal.IdB*llid);
}

inline double ltraj::LLR_Intensity(xtraj *xt, double dx, uint index){
	cerr << "b0rken!" << endl;
	exit(1);
	double xo=xt->x(index);
	double xn=xo+dx;
	double zdo=zeta_d(xo);
	double zdn=zeta_d(xn);
	double zao=zeta_a(xo);
	double zan=zeta_a(xn);

	return ndcache[index]*log(zdo/zdn)+nacache[index]*log(zao/zan)
		-(cal.IdB*(zdo-zdn)+cal.IaB*(zao-zan))*xt->get_dt();
}

double ltraj::LL_Langevin(xtraj* xt){
	double lle=0;
	double lll=0;
	double d1=0;
	double d2=1;
	double xl=xt->length()-1;
	xtraj* xdatahold=xdata;
	xdata=xt;

	for(uint i=1;i<xl-1;i++){
		if(!D1ZERO){
			d1=D(1,i,xt->t(i));
		}
		if(!D2CONST){
			d2=D(2,i,xt->t(i));
			lll-=log(d2);
		}
		lle-=gsl_pow_2(xt->x(i)-xt->x(i-1)-d1)/d2;
	}

	if(D2CONST){
		lle/=4*D(2,0,0)*xt->get_dt();
		lll=-((xt->length()-1)/2)*log(4*M_PI*D(2,0,0)*xt->get_dt());
	} else{
		lle/=(4*xt->get_dt());
		lll=lll/2-((xt->length()-1)/2)*log(4*M_PI*xt->get_dt());
	}

	xdata=xdatahold;
	return lle+lll;
}

double ltraj::LLR_Langevin(xtraj *xt1, xtraj *xt2, uint i1, uint i2){
	if(xt1->length()!=xt2->length()) throw dle(dle::LTRAJ_LLP_SIZE);

	double lle=0;
	double lll=0;
	double d1_1=0;
	double d1_2=0;
	double d2_1=1;
	double d2_2=1;
	double dt=xt1->get_dt();
	if(i2==xt1->length()) i2--;
	xtraj* xtrajhold=xdata;

	for(uint i=i1;i<i2;i++){
		if(!D1ZERO){
			xdata=xt1;
			d1_1=D(1,i,xdata->t(i));
			xdata=xt2;
			d1_2=D(1,i,xdata->t(i));
		}
		if(!D2CONST){
			xdata=xt1;
			d2_1=D(2,i,xdata->t(i));
			xdata=xt2;
			d2_2=D(2,i,xdata->t(i));
			lll-=log(d2_1)-log(d2_2);
		}
		lle-=gsl_pow_2(xt1->l(i+1)-xt1->l(i)-d1_1*dt)/d2_1
			-gsl_pow_2(xt2->l(i+1)-xt2->l(i)-d1_2*dt)/d2_2;
 	}

	if(D2CONST){
		lle/=4*D(2,0,0)*dt;
		lll=0;//((i2-i1-1)/2)*log(4*M_PI*D(2,0,0)*xt1.get_dt(N));
	} else{
		lle/=(4*dt);
		lll=lll/2-((i2-i1-1)/2)*log(4*M_PI*dt);
	}

	xdata=xtrajhold;
	return lle+lll;
}

double ltraj::LLR_Langevin(friction *g1, friction *g2, uint i1, uint i2){
	if(g1->length()!=g2->length()) throw dle(dle::LTRAJ_LLP_SIZE);
	double lle1=0;
	double lle2=0;
	double lll1=0;
	double lll2=0;
	double d1_1=0;
	double d1_2=0;
	double d2_1=1;
	double d2_2=1;
	double dt=xdata->get_dt();
	friction* gammahold=gamma;
	uint stop=xdata->length();

	for(uint i=0;i<stop-1;i++){
		if(!D1ZERO){
			gamma=g1;
			d1_1=D(1,i,xdata->t(i));
			gamma=g2;
			d1_2=D(1,i,xdata->t(i));
		}
		if(!D2CONST){
			gamma=g1;
			d2_1=D(2,i,xdata->t(i));
			gamma=g2;
			d2_2=D(2,i,xdata->t(i));
			lll1-=log(d2_1);
			lll2-=log(d2_2);
		}
		lle1-=gsl_pow_2(xdata->l(i+1)-xdata->l(i)-d1_1*dt)/d2_1;
		lle2-=gsl_pow_2(xdata->l(i+1)-xdata->l(i)-d1_2*dt)/d2_2;
 	}

 	if(D2CONST){
		gamma=g1;
 		lle1/=4*D(2,0,0)*dt;
		lll1=-((stop-1)/2.0)*log(4*M_PI*D(2,0,0)*dt);
		gamma=g2;
 		lle2/=4*D(2,0,0)*dt;
		lll2=-((stop-1)/2.0)*log(4*M_PI*D(2,0,0)*dt);
 	} else{
		lle1/=4.0*dt;
		lle2/=4.0*dt;
		lll1=lll1/2.0;//-((stop-1)/2)*log(4*M_PI*dt);
		lll2=lll2/2.0;//-((stop-1)/2)*log(4*M_PI*dt); These cancel!
 	}

	gamma=gammahold;
	return lle1-lle2+lll1-lll2;
}

// inline double ltraj::LLR_Langevin(const xtraj &xt, double dx, int index){
// 	cout << "broken!!" << endl;
// 	exit(1);
// 	//TODO: PROBLEM WITH VELOCITY!!!
// 	double lle=0;
// 	double lll=0;
// 	double d1_1=0;
// 	double d1_2=0;
// 	double d2_1=1;
// 	double d2_2=1;
// 	//	xtraj xt=xtp;
// 	double dt=xt.get_dt();
// 	double xtt=xt.t(index);
// 	double xo;
// 	double xn;

// 	xo=xt.l(index);
// 	xn=xo+dx;

// 	if(index>0){
// 		double xm;
// 		xm=xt.l(index-1);

// 		if(!D1ZERO){
// 			d1_1=D(1,xt,pdata,index-1,xtt-dt);
// 			d1_2=d1_1;
// 		}
// 		if(!D2CONST){
// 			d2_1=D(2,xt,pdata,index-1,xtt-dt);
// 			d2_2=d2_1;
// 		}
// 		lle-=gsl_pow_2(xo-xm-d1_1*dt)/d2_1
// 			-gsl_pow_2(xn-xm-d1_2*dt)/d2_2;
// 	}

// 	if(index<xt.length()-1){
// 		double xp;
// 		xp=xt.l(index+1);

// 		if(!D1ZERO){
// 			d1_2=D(1,xt,pdata,index,xtt);
// 			d1_1=d1_2;//D(1,xt,index,xtt);
// 		}
// 		if(!D2CONST){
// 			d2_2=D(2,xt,pdata,index,xtt);
// 			d2_1=d2_2;//D(2,xt,index,xtt);
// 			lll-=log(d2_1/d2_2);
// 		}
// 		lle-=gsl_pow_2(xp-xo-d1_1*dt)/d2_1
// 			-gsl_pow_2(xp-xn-d1_2*dt)/d2_2;
// 	}

// 	if(D2CONST){
// 		lle=lle/(4*D(2,xt,pdata,0,0)*dt);
// 		lll=0;//(2/2)*log(4*M_PI*D(2,0,0)*xt.get_dt(N));
// 	} else{
// 		lle/=(4*dt);
// 		lll=lll/2;//-(2/2)*log(4*M_PI*xt.get_dt(N));
// 	}

// 	return lle+lll;	
// }

void ltraj::dLLda(double* grad){
	dLLda(*xdata, grad);
	return;
}

void ltraj::dLLda(xtraj& xt, double* grad){
	uint xl=xt.length();
	double xtdt=xt.get_dt();
	int di=1;
	double d1=0;
	double dlldd1=0;
	double* d1s=new double[N_PARAM];
	double* d1da=new double[N_PARAM];
	double d2=1;
	double dlldd2_1=0;
	double dlldd2_2=0;
	double* d2s1=new double[N_PARAM];
	double* d2s2=new double[N_PARAM];
	double* d2da=new double[N_PARAM];

	for(int j=0;j<N_PARAM;j++){
		d1s[j]=0;
		d1da[j]=0;
		d2s1[j]=0;
		d2s2[j]=0;
		d2da[j]=0;
	}

	for(uint i=0;i<xl-di;i++){
		d1=D(1,i,xt.t(i));
		d2=D(2,i,xt.t(i));
		if(!D1ZERO){
			dlldd1=(xt.l(i+di)-xt.l(i)-d1*xt.get_dt()*di)/(2*d2);
			D_p(1,i,xt.t(i),d1da);
			for(int j=0;j<N_PARAM;j++){
				d1s[j]+=dlldd1*d1da[j];
			}
		}

		dlldd2_1=-1/(2*d2);
		dlldd2_2=gsl_pow_2((xt.l(i+di)-xt.l(i)-d1*xtdt*di)/d2)/(4*xt.get_dt()*di);
		D_p(2,i,xt.t(i),d2da);
		for(int j=0;j<N_PARAM;j++){
			d2s1[j]+=dlldd2_1*d2da[j]+dlldd2_2*d2da[j];
		}
	}

	//TODO: Makes all trajectories same weight... good? bad?
	for(int j=0;j<N_PARAM;j++){
		grad[j]=d1s[j]+d2s1[j];
	}

	delete[] d1s;
	delete[] d1da;
	delete[] d2s1;
	delete[] d2s2;
	delete[] d2da;
}

void ltraj::mean_dLLda(double* grad, double* gradstd, double error){
	mean_dLLda(*xdata, grad, gradstd, error);
	return;
}

void ltraj::mean_dLLda(xtraj& xt, double* grad, double* gradstd, double error){
// 	double* gs=new double[N_PARAM];
// 	double* gs2=new double[N_PARAM];
// 	double* mg=new double[N_PARAM];
// 	double* rsg=new double[N_PARAM];

// 	for(int j=0;j<N_PARAM;j++){
// 		gs[j]=0;
// 		gs2[j]=0;
// 		mg[j]=0;
// 		rsg[j]=0;
// 	}

// 	ofstream dlout("dll.txt");
// 	for(int i=0;i<N_EQ+N_CALC;i++){
// 		spinner(i,N_EQ+N_CALC,1);
// 		xt.mcstep(&this);
// 		xt.calc_from_ldata(0,xt.length());
// 		print("dl.out");
// 		dLLda(grad);
// 		if(i%1==0){
//  			for(int j=0;j<N_PARAM;j++)
//  				dlout << grad[j] << '\t';
//  			dlout << endl << flush;
// 		}
// 		if(i<N_EQ) continue;
// 		for(int j=0;j<N_PARAM;j++){
// 			gs[j]+=grad[j];
// 			gs2[j]+=gsl_pow_2(grad[j]);
// 		}
// 	}

// 	for(int j=0;j<N_PARAM;j++){
// 		grad[j]=gs[j]/N_CALC;
// 		gradstd[j]=sqrt((gs2[j]/N_CALC-gsl_pow_2(mg[j]))/N_CALC);
// 	}

// 	delete[] gs;
// 	delete[] gs2;
// 	delete[] mg;
// 	delete[] rsg;

// 	return;
}

double ltraj::D(int order, uint i, double t){
	switch(D_chosen){
	case SIMPLE:
		return D_simp(order, i, t);
		break;
	case SIMPLEV:
		return D_simpv(order, i, t);
		break;
	case PMF:
		return D_pmf(order, i, t);
		break;
	case PMFV:
		return D_pmfv(order, i, t);
		break;
	case MEMORY:
		return D_memory(order, i, t);
		break;
	default:
		throw dle(dle::LTRAJ_NO_SUCH_MODEL);
		break;
	}
}

void ltraj::D_p(int order, uint i, double t, double* grad){
	switch(D_chosen){
	case SIMPLE:
		D_simp_p(order, i, t, grad);
		break;
	case SIMPLEV:
		D_simpv_p(order, i, t, grad);
		break;
	case PMF:
		D_pmf_p(order, i, t, grad);
		break;
	case PMFV:
		D_pmfv_p(order, i, t, grad);
		break;
	case MEMORY:
		return D_memory_p(order, i, t, grad);
		break;
	default:
		throw dle(dle::LTRAJ_NO_SUCH_MODEL);
		break;
	}
	return;
}

void ltraj::D_adjust_params(double* grad, double ss, double sl){
	switch(D_chosen){
	case SIMPLE:
		D_simp_adjust_params(grad, ss, sl);
		break;
	case SIMPLEV:
		D_simpv_adjust_params(grad, ss, sl);
		break;
	case PMF:
		D_pmf_adjust_params(grad, ss, sl);
		break;
	case PMFV:
		D_pmfv_adjust_params(grad, ss, sl);
		break;
	case MEMORY:
		return D_memory_adjust_params(grad, ss, sl);
		break;
	default:
		throw dle(dle::LTRAJ_NO_SUCH_MODEL);
		break;
	}
	return;
}

void ltraj::D_select(D_type kind){
	D_chosen=kind;
	if(D_chosen==SIMPLE){
		D1ZERO=true;
		D2CONST=false;
		param_gamma=0;
		N_PARAM=1;
		if(xdata) xdata->L_select(xtraj::POSITION);
	} else if(D_chosen==SIMPLEV){
		D1ZERO=false;
		D2CONST=true;
		param_gamma=0;
		N_PARAM=1;
		if(xdata) xdata->L_select(xtraj::VELOCITY);
	}else if(D_chosen==PMF){
		D1ZERO=false;
		D2CONST=false;
		param_gamma=0;
		param_pmf_start=1;
		param_pmf_end=1+pdata->length();
		N_PARAM=param_pmf_end;
		if(xdata) xdata->L_select(xtraj::POSITION);
	}else if(D_chosen==PMFV){
		D1ZERO=false;
		D2CONST=false;
		param_gamma=0;
		param_pmf_start=1;
		param_pmf_end=param_pmf_start+pdata->length();
		N_PARAM=param_pmf_end;
		if(xdata) xdata->L_select(xtraj::VELOCITY);
	} else if(D_chosen==MEMORY){
// 		D1ZERO=false;
// 		D2CONST=true;
// 		param_gamma=0;
// 		param_memory_start=1;
// 		param_memory_end=param_memory_start+mkernel.length();
// 		N_PARAM=param_memory_end;
	}
}

void ltraj::D_select(const char* kind){
	if(strcmp(kind,"simple")==0){
		cout << "Langevin model: SIMPLE" << endl;
		D_select(SIMPLE);
	} else if(strcmp(kind,"simplev")==0){
		cout << "Langevin model: SIMPLEV" << endl;
		D_select(SIMPLEV);
	} else if(strcmp(kind,"pmf")==0){
		cout << "Langevin model: PMF" << endl;
		D_select(PMF);
	} else if(strcmp(kind,"pmfv")==0){
		cout << "Langevin model: PMFV" << endl;
		D_select(PMFV);
	} else if(strcmp(kind,"memory")==0){
		cout << "Langevin model: MEMORY" << endl;
		D_select(MEMORY);
	} else{
		cout << "No such model: " << kind << endl;
		throw dle::LTRAJ_NO_SUCH_MODEL;
	}
}

ltraj::D_type ltraj::get_D(){
	return D_chosen;
}



// void ltraj::D_select_from_env(){
// 	char* kind=getenv("SMURF_LMODEL");
// 	if(kind){

// }

double ltraj::get_time(){
	return pcache[np-1].time;
}

void ltraj::set_photons(string datafile){
	photons=new s_fretdata(datafile);
}

s_fretdata* ltraj::get_photons(){
	return photons;
}

// void ltraj::set_xtraj_from_data(double alpha, double dt){
// 	s_x tx, txold;
// 	double cx=1;
// 	double ct=0;
// 	int i;

// 	xtt.set_xdata_zero(cal.ta_bleach, dt);
// 	double oldx=1;
// 	double oldt=0;
// 	double a;
// 	double b;
// 	double tmin=0,tmax=0;

// 	for(tx=photons->mipstep(0,alpha,1);
// 		tx.tmax<cal.ta_bleach;
// 		tx=photons->mipstep(tx.tmax,alpha,1))
// 	{
// 		if(!tx.error){
// 			cx=gsl_vector_get(tx.c,0);
// 			tmax=tx.tmax/2+tx.tmin/2;
// 		} else{
// 			continue;
// 		}

// 		for(i=xtt.n(tmin);i<xtt.n(tmax);i++){
// 			a=(xtt.t(i)-tmin)/(tmax-tmin);
// 			b=(tmax-xtt.t(i))/(tmax-tmin);
// 			xtt.set_x(i, a*cx+b*oldx);
// 		}
// 		oldx=cx;
// 		oldt=tx.tmax;
// //		txold=tx;
// 		tmin=tmax;
// 	}

// 	for(;i<xtt.length();i++)
// 		xtt.set_x(i,cx);

// 	xdata=xtt;
// 	cache_pn();
// 	return;
// }

void ltraj::set_xtraj(xtraj* nx){
	xdata=nx;
	xdata->set_ltraj(this);
	cache_pn();
}

xtraj* ltraj::get_xtraj(){
	return xdata;
}

void ltraj::set_landscape(landscape* nl){
	pdata=nl;
}

landscape* ltraj::get_landscape(){
	return pdata;
}

void ltraj::set_friction(friction* ng){
	gamma=ng;
}

friction* ltraj::get_friction(){
	return gamma;
}

void ltraj::set_ncalc(int ncn){
	N_CALC=ncn;
}

void ltraj::set_neq(int neq){
	N_EQ=neq;
}

double ltraj::get_beta(){
	return beta;
}

void ltraj::set_beta(double nb){
	beta=nb;
}

double ltraj::get_mass(){
	return mass;
}

void ltraj::set_mass(double nm){
	mass=nm;
}

int ltraj::param_length(){
	return N_PARAM;
}

void ltraj::print(char* file){
	ofstream lout(file);
	lout << scientific;
	lout << gamma->length() << '\t' << xdata->length() << '\t' << pdata->length() << endl;
	gamma->print(lout);
	xdata->print(lout);
	pdata->print_pmf(lout);
}

void ltraj::optimize_gamma(double g1, double g4, double tol){
	cerr << "Broken" << endl;
	exit(1);
// 	//Optimize gamma given the current x trajectory, based on LL_Langevin
// 	double g2=0.3*g4+0.7*g1;
// 	double g3=0.7*g4+0.3*g1;
// 	//xtt=xdata->restrict(1000);

// 	set_gamma(g1);
// 	double ll1=LL_Langevin(&xtt);
// 	set_gamma(g2);
// 	double ll2=LL_Langevin(&xtt);
// 	set_gamma(g3);
// 	double ll3=LL_Langevin(&xtt);
// 	set_gamma(g4);
// 	double ll4=LL_Langevin(&xtt);

// 	while(fabs((g3-g2)/g3)>tol){
// 		if(ll2>ll3){
// 			g4=g3;
// 			ll4=ll3;
// 			g3=g2;
// 			ll3=ll2;
// 			g2=0.4*g1+0.6*g2;
// 			set_gamma(g2);
// 			ll2=LL_Langevin(&xtt);
// 		} else{
// 			g1=g2;
// 			ll1=ll2;
// 			g2=g3;
// 			ll2=ll3;
// 			g3=0.4*g4+0.6*g3;
// 			set_gamma(g3);
// 			ll3=LL_Langevin(&xtt);
// 		}
// 	}

// 	if(ll2>ll3) set_gamma(g2);
// 	else set_gamma(g3);
}

int ltraj::optimize_xdata(){
// 	double step_size;

// 	step_size=1e-3;
// 	double beta_old=beta;
// 	int N_OPT=4;
// 	int N_MC=50;
// 	double beta_start=beta;

// 	for(int i=0;i<N_OPT;i++){
// 		set_beta(exp((i*log(beta_old)+(N_OPT-i-1)*log(beta_start))/(N_OPT-1)));
//   		for(int j=0;j<N_MC;j++){
// 			xdata->mcstep(this);
//   		}
// 		cerr << i+1 << '/' << N_OPT << '\t' << beta << '\t' << LL(xdata) << endl;
// 	}

// 	set_beta(beta_old);
// 	return 0;
	return 0;
}

double ltraj::D_simp(int order, uint i, double t){
	if(order==1){
		return 0;
	} else if(order==2){
		uint gint=0;
		if(gamma->xint_cached())
			gint=gamma->get_xint(i);
		else
			gint=gamma->n(xdata->x(i));

		return 1/gamma->g(gint)/beta;
	} else{
		throw dle(dle::INVALID_ORDER);
	}
}

void ltraj::D_simp_p(int order, uint i, double t, double* grad){
	if(order==1){
		grad[param_gamma]=0;
	} else if(order==2){
		uint gint=0;
		if(gamma->xint_cached())
			gint=gamma->get_xint(i);
		else
			gint=gamma->n(xdata->x(i));

		grad[param_gamma]=-1/gsl_pow_2(gamma->g(gint))/beta;
	} else{
		throw dle(dle::INVALID_ORDER);
	}
	return;
}

void ltraj::D_simp_adjust_params(double* grad, double ss, double sl){
	if(fabs(grad[param_gamma]*ss)>fabs(sl*gamma->g(uint(0)))){
		ss=sl*gamma->g(uint(0))/grad[param_gamma]*grad[param_gamma]/fabs(grad[param_gamma]);
	}

	gamma->set_param(uint(0),gamma->get_param(uint(0))+grad[param_gamma]*ss);
}

double ltraj::D_simpv(int order, uint i, double t){
	if(order==1){
		uint gint=0;
		if(gamma->xint_cached())
			gint=gamma->get_xint(i);
		else
			gint=gamma->n(xdata->x(i));

		return -gamma->g(gint)*xdata->v(i);
	} else if(order==2){
		uint gint=0;
		if(gamma->xint_cached())
			gint=gamma->get_xint(i);
		else
			gint=gamma->n(xdata->x(i));

		return gamma->g(gint)/beta/mass;
	} else{
		throw dle(dle::INVALID_ORDER);
	}
}

void ltraj::D_simpv_p(int order, uint i, double t, double* grad){
	exit(1);
	if(order==1){
		grad[param_gamma]=-xdata->v(i);
		grad[param_mass]=0;
	} else if(order==2){
		grad[param_gamma]=1/beta/mass;
		grad[param_mass]=-gamma->g(uint(0))/beta/gsl_pow_2(mass);
	} else{
		throw dle(dle::INVALID_ORDER);
	}
}

void ltraj::D_simpv_adjust_params(double* grad, double ss, double sl){
	exit(1);
	if(fabs(grad[param_gamma]*ss)>fabs(sl*gamma->g(uint(0)))){
		ss=sl*gamma->g(uint(0))/grad[param_gamma]*grad[param_gamma]/fabs(grad[param_gamma]);
	}

	gamma->set_param(uint(0),gamma->g(uint(0))+grad[param_gamma]*ss);
	mass+=grad[param_mass]*ss;
}

double ltraj::D_pmf(int order, uint i, double t){
	if(order==1){
		uint pint=0;
		if(pdata->xint_cached())
			pint=pdata->get_xint(i);
		else
			pint=pdata->n(xdata->x(i));

		uint gint=0;
		if(gamma->xint_cached())
			gint=gamma->get_xint(i);
		else
			gint=gamma->n(xdata->x(i));
		return -pdata->pmfp(pint)/gamma->g(gint);
	} else if(order==2){
		uint gint=0;
		if(gamma->xint_cached())
			gint=gamma->get_xint(i);
		else
			gint=gamma->n(xdata->x(i));
		return 1/gamma->g(gint)/beta;
	} else{
		throw dle(dle::INVALID_ORDER);
	}
}

void ltraj::D_pmf_p(int order, uint i, double t, double* grad){
	exit(1);
	if(order==1){
		uint pdn=pdata->n(xdata->x(i));
		grad[param_gamma]=pdata->pmfp(pdn)/gsl_pow_2(gamma->g(uint(0)));
		for(uint i=param_pmf_start;i<param_pmf_end;i++){
			if(i-param_pmf_start == pdn)
				grad[i]=-1/gamma->g(uint(0));
			else
				grad[i]=0;
		}
	} else if(order==2){
		grad[param_gamma]=-1/gsl_pow_2(gamma->g(uint(0)))/beta;
		for(uint i=param_pmf_start;i<param_pmf_end;i++){
			grad[i]=0;
		}
	}
	return;
}

void ltraj::D_pmf_adjust_params(double* grad, double ss, double sl){
	exit(1);
	if(fabs(grad[param_gamma]*ss)>fabs(sl*gamma->g(uint(0)))){
		ss=sl*gamma->g(uint(0))/grad[param_gamma]*grad[param_gamma]/fabs(grad[param_gamma]);
	}

	gamma->set_param(uint(0),gamma->get_param(uint(0))+grad[param_gamma]*ss);

	int pmfi;
	double pdx=pdata->get_dx();
	for(uint i=0;i<param_pmf_end-param_pmf_start;i++){
		pmfi=param_pmf_start+i;
		pdata->set_param(i,pdata->get_param(i)+grad[pmfi]/pdx*ss);
	}
}

double ltraj::D_pmfv(int order, uint i, double t){
	if(order==1){
		return -gamma->g(uint(0))*xdata->v(i)-pdata->pmfp(xdata->x(i))/mass;
	} else if(order==2){
		return gamma->g(uint(0))/beta/mass;
	} else{
		throw dle(dle::INVALID_ORDER);
	}
}

void ltraj::D_pmfv_p(int order, uint i, double t, double* grad){
	if(order==1){
		uint pdn=pdata->n(xdata->x(i));
		grad[param_gamma]=-xdata->v(i);
		for(uint i=param_pmf_start;i<param_pmf_end;i++){
			if(i-param_pmf_start == pdn)
				grad[i]=-1/mass;
			else
				grad[i]=0;
		}
	} else if(order==2){
		grad[param_gamma]=1/beta/mass;
		for(uint i=param_pmf_start;i<param_pmf_end;i++){
			grad[i]=0;
		}
	}
	return;
}

void ltraj::D_pmfv_adjust_params(double* grad, double ss, double sl){
	if(fabs(grad[param_gamma]*ss)>fabs(sl*gamma->g(uint(0)))){
		ss=sl*gamma->g(uint(0))/grad[param_gamma]*grad[param_gamma]/fabs(grad[param_gamma]);
	}

	gamma->set_param(uint(0),gamma->get_param(uint(0))+grad[param_gamma]*ss);

	int pmfi;
	double pdx=pdata->get_dx();
	for(uint i=0;i<param_pmf_end-param_pmf_start;i++){
		pmfi=param_pmf_start+i;
		pdata->set_param(i,pdata->get_param(i)+grad[pmfi]/pdx*ss);
	}
}

//int ltraj::D_pmf_stop(double* grad, double* gradstd){
	//Check to see if we can stop with this gradient calculation and move on to the next.
//}

double ltraj::D_memory(int order, uint i, double t){
	return 0;
// 	if(order==1){
// 		double fric=0;
// 		for(int j=0;(j<mkernel.length())&&(i-j>0);j--){
// 			fric+=xdata->v(i-j)*mkernel.x(j);
// 		}
// 		return -gamma->g(uint(0))*fric;
// 	} else if(order==2){
// 		return gamma->g(uint(0))/beta/mass; //TODO: Is this right with the memory kernel?!
// 	} else{
// 		throw dle(dle::INVALID_ORDER);
// 	}
}

void ltraj::D_memory_p(int order, uint i, double t, double* grad){
// 	if(order==1){
// 		double fric=0;
// 		for(int k=param_memory_start;k<param_memory_end;k++){
// 			grad[k]=-gamma->g(uint(0))*xdata->v(i-(k-param_memory_start));
// 			fric+=mkernel.x(k-param_memory_start)*xdata->v(i-(k-param_memory_start));
// 		}
// 		grad[param_gamma]=fric;
// 		grad[param_mass]=0;
// 	} else if(order==2){
// 		grad[param_gamma]=1/beta/mass;
// 		grad[param_mass]=-gamma->g(uint(0))/beta/gsl_pow_2(mass);
// 		for(int k=param_memory_start;k<param_memory_end;k++){
// 			grad[k]=0;
// 		}
// 	}
}

void ltraj::D_memory_adjust_params(double* grad, double ss, double sl){
	if(fabs(grad[param_gamma]*ss)>fabs(sl*gamma->g(uint(0)))){
		ss=sl*gamma->g(uint(0))/grad[param_gamma]*grad[param_gamma]/fabs(grad[param_gamma]);
	}

	gamma->set_param(uint(0),gamma->get_param(uint(0))+grad[param_gamma]*ss);

	int pmfi;
	double pdx=pdata->get_dx();
	for(uint i=0;i<param_pmf_end-param_pmf_start;i++){
		pmfi=param_pmf_start+i;
		pdata->set_param(i,pdata->get_param(i)+grad[pmfi]/pdx*ss);
	}
}

double ltraj::zeta_d(double xi){
	double x6=gsl_pow_6(xi);
	return (1-bdi)*x6/(1+x6)+bdi;
}

double ltraj::zeta_a(double xi){
	double x6=gsl_pow_6(xi);
	return (1-bai)/(1+x6)+bai;
}

double ltraj::zeta_xyz(double xi, double yi, double zi){
	return exp((-gsl_pow_2(xi)+gsl_pow_2(yi))/wx2-gsl_pow_2(zi)/wz2);
}

// inline double ltraj::zeta_d(double xi){
// 	return zeta_d_cache[(int)(xi/zeta_cache_dx)];
// }

// inline double ltraj::zeta_a(double xi){
// 	return zeta_a_cache[(int)(xi/zeta_cache_dx)];
// }

double ltraj::zeta_d_p(double xi){
	double x5=gsl_pow_5(xi);
	return (1-bdi)*6*x5/gsl_pow_2(1+xi*x5);
}

void ltraj::cache_zetas(double dxn){
	zeta_cache_dx=dxn;
	free_zeta_caches();
	zeta_d_cache=new double[(int)(3.0/zeta_cache_dx)];
	zeta_a_cache=new double[(int)(3.0/zeta_cache_dx)];

	double x6;
	for(int i=0;i<(int)(3.0/zeta_cache_dx);i++){
		x6=gsl_pow_6(i*zeta_cache_dx);
		zeta_a_cache[i]=(1-bai)/(1+x6)+bai;
		zeta_d_cache[i]=(1-bdi)*x6/(1+x6)+bdi;
	}
}

void ltraj::free_zeta_caches(){
	if(zeta_d_cache) delete[] zeta_d_cache;
	zeta_d_cache=NULL;

	if(zeta_a_cache) delete[] zeta_a_cache;
	zeta_a_cache=NULL;
}

double ltraj::zeta_a_p(double xi){
	double x5=gsl_pow_5(xi);
	return -(1-bai)*6*x5/gsl_pow_2(1+xi*x5);
}

void ltraj::cache_photons(){
	free_photons();
	np=cal.nt-1;
	pcache=new s_fr[np];
	bdi=1/cal.bd;
	bai=1/cal.ba;

	for(uint i=0;i<np;i++){
		pcache[i]=photons->get_photon((uint32_t)(i));
	}
}

void ltraj::free_photons(){
	if(pcache)
		delete[] pcache;
	free_pn();
}

void ltraj::cache_pn(){
	free_pn();

	pncache=new uint[np];
	npcache=new uint[xdata->length()+1];
	nacache=new uint[xdata->length()+1];
	ndcache=new uint[xdata->length()+1];

	for(uint i=0;i<xdata->length();i++){
		nacache[i]=0;
		ndcache[i]=0;
	}

	for(uint i=0;i<np;i++){
		pncache[i]=xdata->n(pcache[i].time);
		npcache[pncache[i]]=i;

		if(pcache[i].type==DONOR){
			ndcache[pncache[i]]++;
		}
		else if(pcache[i].type==ACCEPTOR){
			nacache[pncache[i]]++;
		}
		else{
			cerr << "Invalid photon type!" << endl;				
		}
	}

	int nx=xdata->length();
	for(int i=nx;i>0;i--){
		npcache[i]=npcache[i-1];
	}
	npcache[0]=0;
}	

void ltraj::free_pn(){
	if(pncache){
		delete[] pncache;
		pncache=NULL;
	}

	if(npcache){
		delete[] npcache;
		npcache=NULL;
	}

	if(ndcache){
		delete[] ndcache;
		ndcache=NULL;
	}

	if(nacache){
		delete[] nacache;
		nacache=NULL;
	}
}
