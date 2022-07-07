#include "langsim.h"

langsim::langsim(){
	rnd=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rnd,time(NULL));
	x0=0;
	v0=0;
	x=0;
	v=0;
	t=0;
	dyn=SLOW_FRICTION;
	dgam=0;
	beta=1;
	dt=1;
	mass=1;
	memlength=1000;
	mem=gsl_vector_calloc(memlength);
	markovian=true;
	Vp=Vp_default;
	gam=gam_default;
}

langsim::~langsim(){
	gsl_vector_free(mem);
}

//Reset the system to time zero
void langsim::Reset(){
	x=x0;
	v=v0;
	t=0;
	gsl_vector_set_zero(mem);
}

//Evolve the system by one time step
void langsim::Evolve(){
	t+=dt;
	if(dyn==FAST_FRICTION){
		Evolve_FF();
	} else if(dyn==SLOW_FRICTION){
		Evolve_SF();
	}
}

//Evolve under fast friction assumption
void langsim::Evolve_FF(){
	//TODO: Check this!!!
	x+=-dt*Vp(x)/gam(x,dgam)+gsl_ran_gaussian(rnd,sqrt(2/gam(x,dgam)/beta*dt));
	v=0;
}

//Evolve under slow friction assumption
void langsim::Evolve_SF(){
	double fric=0;

	if(markovian){
		//no memory
		fric=v;
	} else {
		//integrate the kernel
		for(int i=0;i<memlength;i++)
			fric+=K(dt*i)*gsl_vector_get(mem,i);
	}

	v+=dt*(-Vp(x)/mass-gam(x,dgam)*fric)+gsl_ran_gaussian(rnd,sqrt(2*gam(x,dgam)/beta/mass*dt));
	x+=v*dt;

	//Update memory
	if(!markovian){
		for(int i=memlength-1;i>0;i--)
			gsl_vector_set(mem,i,gsl_vector_get(mem,i-1));
		gsl_vector_set(mem,0,v);
	}
}

//Output current state
void langsim::Print(ostream& f){
	f << scientific << t << '\t' << x << '\t' << v << endl;
}

//Derivative of the potential energy function
double langsim::Vp_default(double x){
	return x;
}

//Derivative of the potential energy function
double langsim::gam_default(double x,double b){
	return b;
}

//Friction kernel
double langsim::K(double t){
	return exp(-t);
}
