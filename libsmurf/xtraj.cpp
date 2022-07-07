#include <iostream>
#include <math.h>
#include "xtraj.h"
#include "ltraj.h"

xtraj::xtraj(){
	init();
}

xtraj::xtraj(double T, double d){
	init();
	set_xdata_constant(T,d,1.0);
}

xtraj::xtraj(gsl_vector* nx, double d){
	init();
	set_xdata(nx, d);
}

xtraj::xtraj(const xtraj &xa){
	init();
	du=xa.du;
	dt1=1/du;
	u1=xa.u1;
	L_select(xa.L_chosen);
	setup(xa.length());
	if(xa.xdata)
		gsl_vector_memcpy(xdata,xa.xdata);
	if(xa.vdata)
		gsl_vector_memcpy(vdata,xa.vdata);
}

xtraj::~xtraj(){
	free_xdata();
}

double xtraj::x(double t) const{
	return x(n(t));
}

double xtraj::x(uint n) const{
	return gsl_vector_get(xdata,n);
}

double xtraj::v(double t) const{
	return v(n(t));
}

double xtraj::v(uint n) const{
	return gsl_vector_get(vdata,n);
}

double xtraj::l(double t) const{
	return l(n(t));
}

double xtraj::l(uint n) const{
	return gsl_vector_get(params,n);
}

double xtraj::t(uint n) const{
	return u(n);
}

double xtraj::get_dt() const{
	return du;
}

double xtraj::get_T() const{
	return du*length();
}

gsl_vector* xtraj::get_xdata() const{
	return xdata;
}

void xtraj::set_x(uint n, double xn){
	gsl_vector_set(xdata,n,xn);
}

void xtraj::set_v(uint n, double vn){
	gsl_vector_set(vdata,n,vn);
}

void xtraj::set_l(uint n, double ln){
	gsl_vector_set(params,n,ln);
}

void xtraj::set_xdata(gsl_vector* nx, double d){
	setup(nx->size);
	gsl_vector_memcpy(xdata,nx);
	du=d;
	dt1=1/d;
	calc_vdata_from_xdata();
}

void xtraj::set_xdata_zero(double T, double d){
	u1=0;
	du=d;
	dt1=1/d;
	setup((int)(T/du)+1);
}

void xtraj::set_xdata_zero(){
	gsl_vector_set_all(xdata,0);
	calc_vdata_from_xdata();
}

void xtraj::set_xdata_constant(double T, double d, double c){
	set_xdata_zero(T,d);
	gsl_vector_set_all(xdata,c);
	calc_vdata_from_xdata();
}

void xtraj::set_xdata_constant(double c){
	gsl_vector_set_all(xdata,c);
	calc_vdata_from_xdata();
}

void xtraj::set_xdata_from_data_const(s_fretdata* photons, double alpha, double dx){
    s_x tx, txold;
    double cx=1;
    uint i;

    set_xdata_zero(photons->b.ta_bleach, dx);
    double tmin=0,tmax=0;

    for(tx=photons->mipstep(0,alpha,1);
        tx.tmax<photons->b.ta_bleach;
        tx=photons->mipstep(tx.tmax,alpha,1))
    {
        if(!tx.error){
            cx=gsl_vector_get(tx.c,0);
            tmax=tx.tmax;
        } else{
            continue;
        }

        for(i=n(tmin);i<n(tmax);i++){
            set_x(i,cx);
        }
        tmin=tmax;
    }

    for(;i<length();i++)
        set_x(i,cx);
	calc_vdata_from_xdata();
}

void xtraj::set_xdata_from_data_interp(s_fretdata* photons, double alpha, double dx){
    s_x tx, txold;
    double cx=1;
    uint i;

    set_xdata_zero(photons->b.ta_bleach, dx);
    double oldx=1;
    double oldt=0;
    double a;
    double b;
    double tmin=0,tmax=0;

    for(tx=photons->mipstep(0,alpha,1);
        tx.tmax<photons->b.ta_bleach;
        tx=photons->mipstep(tx.tmax,alpha,1))
    {
        if(!tx.error){
            cx=gsl_vector_get(tx.c,0);
            tmax=tx.tmax/2+tx.tmin/2;
        } else{
            continue;
        }

        for(i=n(tmin);i<n(tmax);i++){
            a=(t(i)-tmin)/(tmax-tmin);
            b=(tmax-t(i))/(tmax-tmin);
            set_x(i, a*cx+b*oldx);
        }
        oldx=cx;
        oldt=tx.tmax;
        tmin=tmax;
    }

    for(;i<length();i++)
        set_x(i,cx);
	calc_vdata_from_xdata();
}

void xtraj::set_xdata_from_file(char* file, int N){
	//TODO: Make this prettier
	ifstream xintest(file);

	if(!xintest.is_open()) return;

	double a[2];
	int num=-1;
	while(!xintest.eof()){
		xintest >> a[0] >> a[1];
		xintest.ignore(1000,'\n');
		num++;
	}
	xintest.close();
	if(num<=0) return;

	gsl_vector* nxd=gsl_vector_alloc(num);	
	gsl_vector* ntd=gsl_vector_alloc(num);	
	ifstream xin(file);
	for(int i=0;i<num;i++){
		xin >> a[0] >> a[1];
		xin.ignore(1000,'\n');
		gsl_vector_set(ntd,i,a[0]);
		gsl_vector_set(nxd,i,a[1]);
	}
	xin.close();

	double dx=gsl_vector_get(ntd,1)-gsl_vector_get(ntd,0);

	xtraj xtt;
	xtt.set_xdata(nxd,dx);
	xtt=xtt.restrict(N);
	*this=xtt;
	gsl_vector_free(nxd);
	gsl_vector_free(ntd);
	return;
}

// void xtraj::set_xdata_from_xdata(xtraj* nxt){
// //Set trajectory from another xtraj of different time resolution,
// //interpolating or averaging as required.
// 	if(tmax!=nxt->tmax) throw dle(dle::XTRAJ_RANGE);

// 	if(dt<nxt->dt){	//Need to interpolate
// 		int ni=0;
// 		double ipx=0;
// 		double tb=0;
// 		double ta=0;
// 		double ct=0;
// 		for(int i=0;i<length();i++){
// //TODO: Clean this up
// 			ct=t(i);
// 			try{
// 				ni=nxt->n(ct);
// 				ta=nxt->t(ni);
// 				tb=nxt->t(ni+1);
// 				ipx=((ct-ta)*nxt->x(ni+1)+(tb-ct)*nxt->x(ni))/(tb-ta);
// 			} catch(...){
// 				gsl_vector_set(xdata,i,nxt->x(nxt->length()-1));
// 				continue;
// 			}
// 			gsl_vector_set(xdata,i,ipx);
// 		}
// 	} else if(dt>nxt->dt){ //TODO: Need to average
// 		throw dle(dle::NOT_READY);
// 	}
// }

// void xtraj::add_dx(gsl_vector* dx){
// 	if(dx->size!=length()) throw dle(dle::XTRAJ_RANGE);
// 	gsl_vector_add(xdata,dx);
// }

// void xtraj::subtract_dx(gsl_vector* dx){
// 	if(dx->size!=length()) throw dle(dle::XTRAJ_RANGE);
// 	gsl_vector_sub(xdata,dx);
// }

// void xtraj::augment(double factor){
// 	cout << "NO!" << endl;
// 	exit(1);
// //	xtraj* old=new xtraj(this);
// //	dt*=factor;
// //	setup((int)(tmax/dt)+1);
// //	set_traj_from_xtraj(old);
// //	delete old;
// }

void xtraj::print(ostream &xout){
	for(uint i=0;i<length();i++){
		xout << i*du << '\t' << x(i) << '\t' << v(i);
		xout << endl;
	}
}

xtraj xtraj::restrict(int N) const{
	if(N<=0) throw dle(dle::XTRAJ_RANGE);
	if(N==1) return *this;

	xtraj temp=*this;
	temp.du=N*get_dt();
	temp.dt1=1/temp.du;
	temp.setup(length()/N);

	for(uint i=0;i<temp.length();i++){
		temp.set_x(i,x(i*N));
	}
	return temp;
}

void xtraj::calc_xdata_from_vdata(uint beg, uint end){
	if(L_chosen==POSITION) return;
	//if(uint(end)>xdata->size) end=xdata->size;
	end=xdata->size;
	for(uint i=beg+1;i<end;i++){
		set_x(i,x(i-1)+v(i-1)*du);
	}
}

void xtraj::calc_xdata_from_vdata(uint i){
	if(L_chosen==POSITION) return;
	if(i>=xdata->size-1) return;
	double dx=v(i)*du-(x(i+1)-x(i));
	for(uint j=i;uint(j)<xdata->size;j++){
		set_x(j,x(j)+dx);
	}
}

void xtraj::calc_vdata_from_xdata(uint beg, uint end){
	if(L_chosen==POSITION) return;
	if(uint(end)==xdata->size) end--;
	for(uint i=beg;i<end;i++)
		set_v(i,(x(i+1)-x(i))/du);
	if(uint(end)==xdata->size-1)
		set_v(xdata->size-1,v(uint(xdata->size-2)));
}

void xtraj::calc_vdata_from_xdata(uint i){
	if(L_chosen==POSITION) return;
	set_v(i,(x(i+1)-x(i))/du);
}

void xtraj::calc_from_ldata(uint beg, uint end){
	if(L_chosen==POSITION)
		calc_vdata_from_xdata(beg,end);
	else if(L_chosen==VELOCITY)
		calc_xdata_from_vdata(beg,end);
}

void xtraj::calc_from_ldata(uint i){
	if(L_chosen==POSITION)
		calc_vdata_from_xdata(i);
	else if(L_chosen==VELOCITY)
		calc_xdata_from_vdata(i);
}

void xtraj::update_from_params(uint i){
	calc_from_ldata(i);
}

void xtraj::update_from_params(uint i1, uint i2){
	calc_from_ldata(i1,i2);
}

void xtraj::get_x_minmax(double* min, double* max){
	gsl_vector_minmax(xdata,min,max);
}

void xtraj::get_v_minmax(double* min, double* max){
	gsl_vector_minmax(vdata,min,max);
}


// void xtraj::mcstep(ltraj* wrt){
// 	int stop=length();
// 	int rn=0;
// 	int size=stop;
// 	int fac=4;
// 	int nit=7;
// 	int mnit=(int)floor(log(stop/4)/log(fac));
// 	if(mnit<nit) nit=mnit;

//   	for(int i=0;i<nit;i++){
//   		size=size/fac;
// 		mcstep_lin(wrt,0,0,size);
// 		for(int j=size;j<stop-size;j+=size){
// 			mcstep_lin(wrt,j-size,j,j+size);
// 		}
// 		mcstep_lin(wrt,stop-size,stop,stop);
//   	}

// 	for(int j=0;j<10000;j++)
// 		mcstep(wrt,(int)(gsl_ran_flat(r,1,stop-1)));
// }

// void xtraj::mcstep(ltraj* wrt, int i){
// 	if((MC_chosen==INTX)||(MC_chosen==INTXXYZ)){
// 		mcstep_x(wrt,i);
// 	} else if((MC_chosen==INTV)||(MC_chosen==INTVXYZ)){
// 		mcstep_v(wrt,i);
// 	}

// 	if((MC_chosen==INTXXYZ)||(MC_chosen==INTVXYZ)||(MC_chosen==XYZ)){
// 		mcstep_xyz(wrt,i);
// 	}
// }

// void xtraj::mcstep_x(ltraj* wrt, int i){
// 	double xo=x(i);
// 	double dx=gsl_ran_flat(r,-.2,.2);

// 	xtraj xtt=*this;
// 	xtt.set_x(i,xo+dx);
// 	xtt.calc_vdata_from_xdata(i,xtt.length());

// 	if(accept_move(wrt->LLR(*this,xtt,i-1,i+2)))
// 		set_x(i,xo+dx);
// }

// void xtraj::mcstep_v(ltraj* wrt, int i){
// 	double xo=v(i);
// 	double dx=gsl_ran_flat(r,-.2,.2);

// 	xtraj xtt=*this;
// 	xtt.set_v(i,xo+dx);
// 	xtt.calc_xdata_from_vdata(i,xtt.length());

// 	if(accept_move(wrt->LLR(*this,xtt,i-1,i+2)))
// 		set_v(i,xo+dx);
// }

// void xtraj::mcstep_xyz(ltraj* wrt, int i){
// 	double xo=px(i);
// 	double yo=py(i);
// 	double zo=pz(i);
// 	double dx=gsl_ran_flat(r,-.2,.2);
// 	double dy=gsl_ran_flat(r,-.2,.2);
// 	double dz=gsl_ran_flat(r,-.2,.2);

// 	xtraj xtt=*this;
// 	xtt.set_px(i,xo+dx);
// 	xtt.set_py(i,yo+dx);
// 	xtt.set_pz(i,zo+dx);

// 	if(accept_move(wrt->LLR(*this,xtt,i-1,i+2))){
// 		set_px(i,xo+dx);
// 		set_py(i,yo+dy);
// 		set_pz(i,zo+dz);
// 	}
// }

// void xtraj::mcstep_lin(ltraj* wrt, int i1, int i2, int i3){
// 	if((MC_chosen==INTX)||(MC_chosen==INTXXYZ)){
// 		mcstep_lin_x(wrt, i1, i2, i3);
// 	} else if((MC_chosen==INTV)||(MC_chosen==INTVXYZ)){
// 		mcstep_lin_v(wrt, i1, i2, i3);
// 	}

// 	if((MC_chosen==INTXXYZ)||(MC_chosen==INTVXYZ)||(MC_chosen==XYZ)){
// 		mcstep_lin_xyz(wrt, i1, i2, i3);
// 	}
// }

// void xtraj::mcstep_lin_x(ltraj* wrt, int i1, int i2, int i3){
// 	double range=.3;
// 	double dx=gsl_ran_flat(r,-range,range);
// 	xtraj xtt=*this;

// 	if(i2-i1) for(int i=i1;i<i2;i++){
// 		xtt.set_x(i,x(i)+((double)(i-i1))/((double)(i2-i1))*dx);
// 	}
// 	if(i3-i2) for(int i=i2;i<i3;i++){
// 		xtt.set_x(i,x(i)+((double)(i3-i))/((double)(i3-i2))*dx);
// 	}

// 	if(accept_move(wrt->LLR(*this,xtt,i1,i3)))
// 		gsl_vector_memcpy(xdata,xtt.xdata);
// }

// void xtraj::mcstep_lin_v(ltraj* wrt, int i1, int i2, int i3){
// 	double range=.3;
// 	double dx=gsl_ran_flat(r,-range,range);
// 	xtraj xtt=*this;

// 	if(i2-i1 > 2){
// 		int m1=(i1+i2)/2;
// 			for(int i=i1;i<m1;i++)
// 				xtt.set_v(i,xtt.l(i)+((double)(i-i1))/((double)(m1-i1))*dx);
// 			for(int i=m1;i<i2;i++)
// 				xtt.set_v(i,xtt.l(i)+((double)(i2-i))/((double)(i2-m1))*dx);
// 	}

// 	if(i3-i2 > 2){
// 		int m2=(i2+i3)/2;
// 		for(int i=i2;i<m2;i++)
// 			xtt.set_v(i,xtt.l(i)-((double)(i-i2))/((double)(m2-i2))*dx);
// 		for(int i=m2;i<i3;i++)
// 			xtt.set_v(i,xtt.l(i)-((double)(i3-i))/((double)(i3-m2))*dx);
// 	}

// 	if((i2-i1)!=(i3-i2)) xtt.calc_xdata_from_vdata(i1,xtt.length());
// 	else xtt.calc_xdata_from_vdata(i1,i3);
	
// 	if(accept_move(wrt->LLR(*this,xtt,i1,i3))){
// 		gsl_vector_memcpy(vdata,xtt.vdata);
// 		gsl_vector_memcpy(xdata,xtt.xdata);
// 	}
// }

// void xtraj::mcstep_lin_xyz(ltraj* wrt, int i1, int i2, int i3){
// 	double range=.2;
// 	double dx=gsl_ran_flat(r,-range,range);
// 	double dy=gsl_ran_flat(r,-range,range);
// 	double dz=gsl_ran_flat(r,-range,range);
// 	xtraj xtt=*this;

// 	if(i2-i1) for(int i=i1;i<i2;i++){
// 		xtt.set_px(i,xtt.px(i)+((double)(i-i1))/((double)(i2-i1))*dx);
// 		xtt.set_py(i,xtt.py(i)+((double)(i-i1))/((double)(i2-i1))*dy);
// 		xtt.set_pz(i,xtt.pz(i)+((double)(i-i1))/((double)(i2-i1))*dz);
// 	}
// 	if(i3-i2) for(int i=i2;i<i3;i++){
// 		xtt.set_px(i,xtt.px(i)+((double)(i3-i))/((double)(i3-i2))*dx);
// 		xtt.set_py(i,xtt.py(i)+((double)(i3-i))/((double)(i3-i2))*dy);
// 		xtt.set_pz(i,xtt.pz(i)+((double)(i3-i))/((double)(i3-i2))*dz);
// 	}

// 	if(accept_move(wrt->LLR(*this,xtt,i1,i3))){
// 		gsl_vector_memcpy(pxdata,xtt.pxdata);
// 		gsl_vector_memcpy(pydata,xtt.pydata);
// 		gsl_vector_memcpy(pzdata,xtt.pzdata);
// 	}
// }

double xtraj::LLR(ltraj* wrt, uint i1, uint i2){
	if(wrt!=assoc){
		return 0;
	} else{
		xtraj* xt1=this;
		xtraj* xt2=static_cast<xtraj*>(ltt);
		double llrl=wrt->LLR_Langevin(xt1,xt2,i1,i2);
		double llri=wrt->LLR_Intensity(xt1,xt2,i1,i2);
		return llrl+llri;
	}
}

// double xtraj::LLR_Langevin(ltraj* wrt, uint i1, uint i2){
// 	if(wrt!=assoc) return 0;

// 	xtraj* xt1=this;
// 	xtraj* xt2=static_cast<xtraj*>(ltt);
// 	if(xt1->length()!=xt2->length()) throw dle(dle::LTRAJ_LLP_SIZE);

// 	double lle=0;
// 	double lll=0;
// 	double d1_1=0;
// 	double d1_2=0;
// 	double d2_1=1;
// 	double d2_2=1;
// 	double dt=xt1->get_dt();
// 	if(i2==xt1->length()) i2--;
// 	xtraj* xtrajhold=wrt->xdata;

// 	for(uint i=i1;i<i2;i++){
// 		if(!D1ZERO){
// 			wrt->xdata=xt1;
// 			d1_1=D(1,i,xdata->t(i));
// 			xdata=xt2;
// 			d1_2=D(1,i,xdata->t(i));
// 		}
// 		if(!D2CONST){
// 			xdata=xt1;
// 			d2_1=D(2,i,xdata->t(i));
// 			xdata=xt2;
// 			d2_2=D(2,i,xdata->t(i));
// 			lll-=log(d2_1)-log(d2_2);
// 		}
// 		lle-=gsl_pow_2(xt1->l(i+1)-xt1->l(i)-d1_1*dt)/d2_1
// 			-gsl_pow_2(xt2->l(i+1)-xt2->l(i)-d1_2*dt)/d2_2;
//  	}

// 	if(D2CONST){
// 		lle/=4*D(2,0,0)*dt;
// 		lll=0;//((i2-i1-1)/2)*log(4*M_PI*D(2,0,0)*xt1.get_dt(N));
// 	} else{
// 		lle/=(4*dt);
// 		lll=lll/2-((i2-i1-1)/2)*log(4*M_PI*dt);
// 	}

// 	xdata=xtrajhold;
// 	return lle+lll;
// }

void xtraj::L_select(L_type kind){
	L_chosen=kind;
	switch(L_chosen){
	case POSITION:
		params=xdata;
		break;
	case VELOCITY:
		params=vdata;
		break;
	}
}

xtraj::L_type xtraj::L_get(){
	return L_chosen;
}

// xtraj xtraj::operator+(const xtraj &xa){
// 	xtraj temp=*this;

// 	if((du==xa.du)&&(u1==xa.u1)&&(length()==xa.length())){
// 		gsl_vector_add(temp.xdata,xa.xdata);
// 		return temp;
// 	}

// 	//Find range of this that is intersected by range of xa
// 	double t1=max(u(0),xa.u(0));
// 	double ct;
// 	double t2=min(u(length()-1),xa.u(xa.length())-1);
// 	int n1=n(t1);
// 	int n2=n(t2);
// 	double ta;
// 	double tb;
// 	int ni;
// 	double ix;

// 	for(int i=n1;i<=n2;i++){
// 		ct=t(i);
// 		try{
// 			ni=xa.n(ct);
// 			ta=xa.t(ni);
// 			tb=xa.t(ni+1);
// 			ix=((ct-ta)*xa.x(ni+1)+(tb-ct)*xa.x(ni))/(tb-ta);
// 		} catch(...){
// 			ix=xa.x(ni);
// 		}
// 		temp.set_x(i,temp.x(i)+ix);
// 	}

// 	return temp;
// }

// xtraj xtraj::operator-(const xtraj& xa){
// 	xtraj a=xa;
// 	gsl_vector_scale(a.xdata,-1.0);
// 	return operator+(a);
// }

// xtraj xtraj::operator*(const double fa){
// 	xtraj temp=*this;
// 	gsl_vector_scale(temp.xdata,fa);
// 	return temp;
// }

// xtraj xtraj::operator/(const double fa){
// 	xtraj temp=*this;
// 	gsl_vector_scale(temp.xdata,1.0/fa);
// 	return temp;
// }

xtraj& xtraj::operator=(const xtraj& xa){
	du=xa.du;
	dt1=1/du;
	u1=xa.u1;
	L_select(xa.L_chosen);
	setup(xa.length());
	gsl_vector_memcpy(xdata,xa.xdata);
	gsl_vector_memcpy(vdata,xa.vdata);
	return *this;
}

xtraj* xtraj::clone(){
	return new xtraj(*this);
}

void xtraj::init(){
	dt1=0;
	xdata=NULL;
	vdata=NULL;
	assoc=NULL;
	ltt=NULL;
	L_select(POSITION);
	name="xtraj";
}

void xtraj::free_xdata(){
	if(xdata) gsl_vector_free(xdata);
	xdata=NULL;

	if(vdata) gsl_vector_free(vdata);
	vdata=NULL;
}

void xtraj::alloc_xdata(int size){
	free_xdata();
	xdata=gsl_vector_calloc(size);
	vdata=gsl_vector_calloc(size);
	L_select(L_chosen);
}

void xtraj::setup(int size){
	if(!xdata){
		alloc_xdata(size);
	} else if(xdata->size != uint(size)){
		alloc_xdata(size);
	} else{
		gsl_vector_set_all(xdata,0);
		gsl_vector_set_all(vdata,0);
	}
}
