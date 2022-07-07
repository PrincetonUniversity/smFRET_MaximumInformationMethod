/*

implementation of trajectory class

*/

#include "traj.h"

trajectory::trajectory(uint32_t N){
	rnd=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rnd,::time(NULL));
	Nt=N;
	traj=gsl_vector_alloc(Nt);
	time=gsl_vector_alloc(Nt);
}

void trajectory::fill(gsl_vector* t, gsl_vector* x, gsl_vector* s){
	gsl_vector_memcpy(time,t);
	dt=gsl_vector_get(time,1)-gsl_vector_get(time,0);
	gsl_vector_memcpy(traj,x);
//	for(int i=0;i<Nt;i++){
//		tx=gsl_vector_get(x,i);
//		ts=gsl_vector_get(s,i);
//		gsl_vector_set(traj,i,gsl_ran_flat(rnd,tx-2*ts,tx+2*ts));
//	}
}

void trajectory::print(std::ostream& dout){
	dout << std::scientific;
	for(uint i=0;i<Nt;i++)
		dout << gsl_vector_get(time,i) << '\t' << gsl_vector_get(traj,i) << std::endl;
}

uint32_t trajectory::findtime(double t){
	uint32_t b=(uint32_t)(t/dt-0.5);
	if(b>=Nt){
		cerr << "traj.cpp: ERROR: time out of range - " << t << endl;
		b=Nt-1;
	}
	return b;
}

double trajectory::distance(double t){
	uint32_t b=findtime(t);
	if(b>=Nt-1){
		return gsl_vector_get(traj,Nt-1);
	}
	double du=gsl_vector_get(time,b+1)-t;
	double dd=t-gsl_vector_get(time,b);
	double xu=gsl_vector_get(traj,b+1);
	double xd=gsl_vector_get(traj,b);

	return (du*xd+dd*xu)/dt;
}

